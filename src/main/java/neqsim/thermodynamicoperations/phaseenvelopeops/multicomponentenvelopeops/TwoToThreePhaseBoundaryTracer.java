package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.IncipientPhaseStabilityAnalyzer.CandidatePhase;

/**
 * Adaptive pressure-parameter continuation of a specified two-to-three-phase boundary.
 *
 * <p>
 * Each accepted point is independently corrected by {@link TwoToThreePhaseBoundaryPointSolver}; pressure is only the
 * continuation parameter and is never used to interpolate a displayed boundary. A secant temperature predictor and
 * composition-continuity gate reject branch jumps, while failed corrections reduce the pressure step. This is the
 * regular-pressure portion of the continuation algorithm; a pseudo-arclength driver is still required at pressure
 * turning points.
 * </p>
 */
public final class TwoToThreePhaseBoundaryTracer {
  private final SystemInterface template;
  private final CandidatePhase retainedPhaseZero;
  private final CandidatePhase retainedPhaseOne;
  private final CandidatePhase incipientPhase;
  private double minimumTemperatureK = 100.0;
  private double maximumTemperatureK = 500.0;
  private double initialPressureStepBara = 5.0;
  private double minimumPressureStepBara = 0.05;
  private double maximumPressureStepBara = 20.0;
  private int maximumAcceptedPoints = 200;
  private double maximumCompositionJump = 0.35;
  private double maximumPredictionErrorK = 25.0;

  /** Creates a regular-pressure continuation tracer for three distinct physical phase families. */
  public TwoToThreePhaseBoundaryTracer(SystemInterface template, CandidatePhase retainedPhaseZero,
      CandidatePhase retainedPhaseOne, CandidatePhase incipientPhase) {
    if (template == null || retainedPhaseZero == null || retainedPhaseOne == null || incipientPhase == null) {
      throw new IllegalArgumentException("template and all three phase families must be specified");
    }
    if (retainedPhaseZero == retainedPhaseOne || retainedPhaseZero == incipientPhase
        || retainedPhaseOne == incipientPhase) {
      throw new IllegalArgumentException("the three physical phase families must be distinct");
    }
    this.template = template.clone();
    this.retainedPhaseZero = retainedPhaseZero;
    this.retainedPhaseOne = retainedPhaseOne;
    this.incipientPhase = incipientPhase;
  }

  /** Sets the temperature domain used by both the incipient and two-to-three-phase correctors. */
  public TwoToThreePhaseBoundaryTracer setTemperatureRange(double minimumTemperatureK, double maximumTemperatureK) {
    if (!Double.isFinite(minimumTemperatureK) || minimumTemperatureK < 50.0 || !Double.isFinite(maximumTemperatureK)
        || maximumTemperatureK <= minimumTemperatureK) {
      throw new IllegalArgumentException("invalid continuation temperature range");
    }
    this.minimumTemperatureK = minimumTemperatureK;
    this.maximumTemperatureK = maximumTemperatureK;
    return this;
  }

  /** Sets adaptive pressure-step controls. */
  public TwoToThreePhaseBoundaryTracer setPressureStepControls(double initialPressureStepBara,
      double minimumPressureStepBara, double maximumPressureStepBara, int maximumAcceptedPoints) {
    if (!Double.isFinite(initialPressureStepBara) || initialPressureStepBara <= 0.0
        || !Double.isFinite(minimumPressureStepBara) || minimumPressureStepBara <= 0.0
        || !Double.isFinite(maximumPressureStepBara) || maximumPressureStepBara < minimumPressureStepBara
        || initialPressureStepBara < minimumPressureStepBara || initialPressureStepBara > maximumPressureStepBara
        || maximumAcceptedPoints < 2) {
      throw new IllegalArgumentException("invalid pressure continuation controls");
    }
    this.initialPressureStepBara = initialPressureStepBara;
    this.minimumPressureStepBara = minimumPressureStepBara;
    this.maximumPressureStepBara = maximumPressureStepBara;
    this.maximumAcceptedPoints = maximumAcceptedPoints;
    return this;
  }

  /** Sets branch-jump rejection tolerances. */
  public TwoToThreePhaseBoundaryTracer setContinuityGates(double maximumCompositionJump,
      double maximumPredictionErrorK) {
    if (!Double.isFinite(maximumCompositionJump) || maximumCompositionJump <= 0.0
        || !Double.isFinite(maximumPredictionErrorK) || maximumPredictionErrorK <= 0.0) {
      throw new IllegalArgumentException("invalid continuation continuity gates");
    }
    this.maximumCompositionJump = maximumCompositionJump;
    this.maximumPredictionErrorK = maximumPredictionErrorK;
    return this;
  }

  /**
   * Traces from start to end pressure, including both ends when correction succeeds.
   *
   * @param startPressureBara first pressure in bara
   * @param endPressureBara requested last pressure in bara
   * @return immutable accepted branch and termination diagnostics
   */
  public Result trace(double startPressureBara, double endPressureBara) {
    if (!Double.isFinite(startPressureBara) || startPressureBara <= 0.0 || !Double.isFinite(endPressureBara)
        || endPressureBara <= 0.0 || startPressureBara == endPressureBara) {
      throw new IllegalArgumentException("continuation pressures must be distinct and positive");
    }
    double direction = Math.signum(endPressureBara - startPressureBara);
    double stepMagnitude = initialPressureStepBara;
    List<Point> points = new ArrayList<Point>();
    int attemptedCorrections = 0;
    int rejectedCorrections = 0;

    Correction start = correct(startPressureBara);
    attemptedCorrections++;
    if (!start.isConverged()) {
      return new Result(points, attemptedCorrections, rejectedCorrections, false, startPressureBara,
          "initial boundary correction failed: " + start.failureMessage);
    }
    points.add(new Point(start.boundary));

    String failureMessage = null;
    double lastAttemptedPressure = startPressureBara;
    while (points.size() < maximumAcceptedPoints) {
      Point last = points.get(points.size() - 1);
      if (hasReached(last.pressureBara, endPressureBara, direction)) {
        break;
      }
      double trialPressure = last.pressureBara + direction * stepMagnitude;
      if (hasPassed(trialPressure, endPressureBara, direction)) {
        trialPressure = endPressureBara;
      }
      lastAttemptedPressure = trialPressure;
      double predictedTemperature = predictTemperature(points, trialPressure);
      Correction correction = correct(trialPressure);
      attemptedCorrections++;
      boolean accepted = correction.isConverged();
      String rejectionReason = correction.failureMessage;
      if (accepted) {
        Point candidate = new Point(correction.boundary);
        double predictionError = Math.abs(candidate.temperatureK - predictedTemperature);
        double compositionJump = compositionJump(last, candidate);
        if (predictionError > maximumPredictionErrorK || compositionJump > maximumCompositionJump) {
          accepted = false;
          rejectionReason = "continuity gate rejected correction: predictionErrorK=" + predictionError
              + ", compositionJump=" + compositionJump;
        } else {
          points.add(candidate);
          if (predictionError < 0.25 * maximumPredictionErrorK && compositionJump < 0.25 * maximumCompositionJump) {
            stepMagnitude = Math.min(maximumPressureStepBara, 1.35 * stepMagnitude);
          }
        }
      }
      if (!accepted) {
        rejectedCorrections++;
        stepMagnitude *= 0.5;
        if (stepMagnitude < minimumPressureStepBara) {
          failureMessage = rejectionReason == null ? "minimum pressure step reached" : rejectionReason;
          break;
        }
      }
    }
    boolean reachedEnd = !points.isEmpty()
        && hasReached(points.get(points.size() - 1).pressureBara, endPressureBara, direction);
    if (!reachedEnd && failureMessage == null && points.size() >= maximumAcceptedPoints) {
      failureMessage = "maximum accepted continuation point count reached";
    }
    return new Result(points, attemptedCorrections, rejectedCorrections, reachedEnd, lastAttemptedPressure,
        failureMessage);
  }

  private Correction correct(double pressureBara) {
    IncipientPhaseBoundaryPointSolver.Result retainedDew = new IncipientPhaseBoundaryPointSolver(template,
        retainedPhaseZero, retainedPhaseOne).setNumericalControls(40, 80, 1.0e-5, 1.0e-8)
        .solve(pressureBara, minimumTemperatureK, maximumTemperatureK);
    if (!retainedDew.isConverged()) {
      return Correction.failure("retained incipient boundary failed: " + retainedDew.getFailureMessage());
    }
    double upperTemperatureK = Math.min(maximumTemperatureK, retainedDew.getTemperatureK() - 0.25);
    if (upperTemperatureK <= minimumTemperatureK) {
      return Correction.failure("retained two-phase temperature interval is empty");
    }
    TwoToThreePhaseBoundaryPointSolver.Result boundary = new TwoToThreePhaseBoundaryPointSolver(template,
        retainedPhaseZero, retainedPhaseOne, incipientPhase).setNumericalControls(36, 70, 2.0e-5, 1.0e-8)
        .solve(pressureBara, minimumTemperatureK, upperTemperatureK, 1.0 - 1.0e-6, overallComposition(),
            retainedDew.getIncipientComposition());
    return boundary.isConverged() ? Correction.success(boundary)
        : Correction.failure("two-to-three-phase correction failed: " + boundary.getFailureMessage());
  }

  private double[] overallComposition() {
    double[] composition = new double[template.getPhase(0).getNumberOfComponents()];
    for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
      composition[componentIndex] = template.getPhase(0).getComponent(componentIndex).getz();
    }
    return composition;
  }

  private static double predictTemperature(List<Point> points, double trialPressureBara) {
    Point last = points.get(points.size() - 1);
    if (points.size() < 2) {
      return last.temperatureK;
    }
    Point previous = points.get(points.size() - 2);
    double pressureDifference = last.pressureBara - previous.pressureBara;
    if (Math.abs(pressureDifference) < 1.0e-12) {
      return last.temperatureK;
    }
    return last.temperatureK
        + (last.temperatureK - previous.temperatureK) / pressureDifference * (trialPressureBara - last.pressureBara);
  }

  private static double compositionJump(Point first, Point second) {
    return Math.max(compositionDistance(first.phaseZeroComposition, second.phaseZeroComposition),
        Math.max(compositionDistance(first.phaseOneComposition, second.phaseOneComposition),
            compositionDistance(first.incipientComposition, second.incipientComposition)));
  }

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int componentIndex = 0; componentIndex < first.length; componentIndex++) {
      distance += Math.abs(first[componentIndex] - second[componentIndex]);
    }
    return distance;
  }

  private static boolean hasReached(double pressure, double endPressure, double direction) {
    return direction > 0.0 ? pressure >= endPressure - 1.0e-12 : pressure <= endPressure + 1.0e-12;
  }

  private static boolean hasPassed(double pressure, double endPressure, double direction) {
    return direction > 0.0 ? pressure > endPressure : pressure < endPressure;
  }

  private static final class Correction {
    private final TwoToThreePhaseBoundaryPointSolver.Result boundary;
    private final String failureMessage;

    private Correction(TwoToThreePhaseBoundaryPointSolver.Result boundary, String failureMessage) {
      this.boundary = boundary;
      this.failureMessage = failureMessage;
    }

    private static Correction success(TwoToThreePhaseBoundaryPointSolver.Result boundary) {
      return new Correction(boundary, null);
    }

    private static Correction failure(String failureMessage) {
      return new Correction(null, failureMessage);
    }

    private boolean isConverged() {
      return boundary != null && boundary.isConverged();
    }
  }

  /** One accepted, fully corrected continuation point. */
  public static final class Point {
    private final double temperatureK;
    private final double pressureBara;
    private final double beta;
    private final double[] phaseZeroComposition;
    private final double[] phaseOneComposition;
    private final double[] incipientComposition;
    private final double retainedFlashResidual;
    private final double tangentPlaneDistance;
    private final double stationarityResidual;

    private Point(TwoToThreePhaseBoundaryPointSolver.Result boundary) {
      this.temperatureK = boundary.getTemperatureK();
      this.pressureBara = boundary.getPressureBara();
      this.beta = boundary.getBeta();
      this.phaseZeroComposition = boundary.getPhaseZeroComposition();
      this.phaseOneComposition = boundary.getPhaseOneComposition();
      this.incipientComposition = boundary.getIncipientComposition();
      this.retainedFlashResidual = boundary.getRetainedFlashResidual();
      this.tangentPlaneDistance = boundary.getTangentPlaneDistance();
      this.stationarityResidual = boundary.getStationarityResidual();
    }

    public double getTemperatureK() {
      return temperatureK;
    }

    public double getPressureBara() {
      return pressureBara;
    }

    public double getBeta() {
      return beta;
    }

    public double[] getPhaseZeroComposition() {
      return phaseZeroComposition.clone();
    }

    public double[] getPhaseOneComposition() {
      return phaseOneComposition.clone();
    }

    public double[] getIncipientComposition() {
      return incipientComposition.clone();
    }

    public double getRetainedFlashResidual() {
      return retainedFlashResidual;
    }

    public double getTangentPlaneDistance() {
      return tangentPlaneDistance;
    }

    public double getStationarityResidual() {
      return stationarityResidual;
    }
  }

  /** Immutable regular-pressure continuation result. */
  public static final class Result {
    private final List<Point> points;
    private final int attemptedCorrections;
    private final int rejectedCorrections;
    private final boolean reachedRequestedEnd;
    private final double lastAttemptedPressureBara;
    private final String failureMessage;

    private Result(List<Point> points, int attemptedCorrections, int rejectedCorrections, boolean reachedRequestedEnd,
        double lastAttemptedPressureBara, String failureMessage) {
      this.points = Collections.unmodifiableList(new ArrayList<Point>(points));
      this.attemptedCorrections = attemptedCorrections;
      this.rejectedCorrections = rejectedCorrections;
      this.reachedRequestedEnd = reachedRequestedEnd;
      this.lastAttemptedPressureBara = lastAttemptedPressureBara;
      this.failureMessage = failureMessage;
    }

    public List<Point> getPoints() {
      return points;
    }

    public int getAttemptedCorrections() {
      return attemptedCorrections;
    }

    public int getRejectedCorrections() {
      return rejectedCorrections;
    }

    public boolean hasReachedRequestedEnd() {
      return reachedRequestedEnd;
    }

    public double getLastAttemptedPressureBara() {
      return lastAttemptedPressureBara;
    }

    public String getFailureMessage() {
      return failureMessage;
    }
  }
}
