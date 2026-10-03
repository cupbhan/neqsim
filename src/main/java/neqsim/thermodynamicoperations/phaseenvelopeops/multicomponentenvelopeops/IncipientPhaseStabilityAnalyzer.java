package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;

/**
 * Non-destructive tangent-plane stability analysis for candidate gas, oil, and aqueous incipient phases.
 *
 * <p>
 * This analyzer is intended for phase-envelope branch monitoring. Unlike {@code TPmultiflash}, it never adds or removes
 * phases from the caller. It minimizes Michelsen's tangent-plane distance from physically distinct gas-, oil-, and
 * water-rich seeds and returns every converged stationary trial. A negative non-trivial stationary point indicates that
 * an additional phase can lower Gibbs energy and therefore brackets a topology transition along a traced boundary.
 * </p>
 */
public final class IncipientPhaseStabilityAnalyzer {
  private static final double MINIMUM_COMPOSITION = 1.0e-100;
  private static final double MINIMUM_LOG_VALUE = -700.0;
  private static final double MAXIMUM_LOG_VALUE = 700.0;
  private static final double DEFAULT_CONVERGENCE_TOLERANCE = 1.0e-9;
  private static final double DEFAULT_INSTABILITY_TOLERANCE = -1.0e-8;
  private static final double EQUIVALENT_STATIONARY_POINT_TOLERANCE = 1.0e-7;

  private final SystemInterface sourceSystem;
  private int maximumIterations = 150;
  private double convergenceTolerance = DEFAULT_CONVERGENCE_TOLERANCE;
  private double instabilityTolerance = DEFAULT_INSTABILITY_TOLERANCE;
  private double dampingFactor = 0.5;

  /** Candidate phase families considered by the analyzer. */
  public enum CandidatePhase {
    GAS(PhaseType.GAS), OIL(PhaseType.OIL), AQUEOUS(PhaseType.AQUEOUS);

    private final PhaseType neqsimPhaseType;

    CandidatePhase(PhaseType neqsimPhaseType) {
      this.neqsimPhaseType = neqsimPhaseType;
    }

    private PhaseType getNeqsimPhaseType() {
      return neqsimPhaseType;
    }
  }

  /**
   * Creates a stability analyzer for an already initialized state.
   *
   * @param sourceSystem thermodynamic state at one PT boundary point
   */
  public IncipientPhaseStabilityAnalyzer(SystemInterface sourceSystem) {
    if (sourceSystem == null) {
      throw new IllegalArgumentException("sourceSystem must not be null");
    }
    this.sourceSystem = sourceSystem;
  }

  /**
   * Sets the fixed-point iteration limit for each physical seed.
   *
   * @param maximumIterations positive iteration limit
   * @return this analyzer
   */
  public IncipientPhaseStabilityAnalyzer setMaximumIterations(int maximumIterations) {
    if (maximumIterations < 5) {
      throw new IllegalArgumentException("maximumIterations must be at least 5");
    }
    this.maximumIterations = maximumIterations;
    return this;
  }

  /**
   * Sets fixed-point convergence and negative-instability tolerances.
   *
   * @param convergenceTolerance maximum logarithmic stationarity residual
   * @param instabilityTolerance negative threshold for {@code 1 - sum(W)}
   * @return this analyzer
   */
  public IncipientPhaseStabilityAnalyzer setTolerances(double convergenceTolerance, double instabilityTolerance) {
    if (!Double.isFinite(convergenceTolerance) || convergenceTolerance <= 0.0 || !Double.isFinite(instabilityTolerance)
        || instabilityTolerance >= 0.0) {
      throw new IllegalArgumentException("convergence tolerance must be positive and instability tolerance negative");
    }
    this.convergenceTolerance = convergenceTolerance;
    this.instabilityTolerance = instabilityTolerance;
    return this;
  }

  /**
   * Sets logarithmic successive-substitution damping.
   *
   * @param dampingFactor value in {@code (0, 1]}
   * @return this analyzer
   */
  public IncipientPhaseStabilityAnalyzer setDampingFactor(double dampingFactor) {
    if (!Double.isFinite(dampingFactor) || dampingFactor <= 0.0 || dampingFactor > 1.0) {
      throw new IllegalArgumentException("dampingFactor must be in (0, 1]");
    }
    this.dampingFactor = dampingFactor;
    return this;
  }

  /**
   * Runs all physically relevant trial phases without mutating the caller.
   *
   * @return immutable analysis result ordered from lowest to highest tangent-plane distance
   */
  public Result analyze() {
    SystemInterface reference = sourceSystem.clone();
    if (reference == null) {
      throw new IllegalStateException("thermodynamic system clone failed");
    }
    reference.init(1);
    ChemicalPotentialReference chemicalPotentialReference = chemicalPotentialReference(reference);
    List<Candidate> trials = new ArrayList<Candidate>();
    trials.add(minimize(reference, chemicalPotentialReference.logFugacity, CandidatePhase.GAS));
    trials.add(minimize(reference, chemicalPotentialReference.logFugacity, CandidatePhase.OIL));
    if (reference.hasComponent("water")) {
      trials.add(minimize(reference, chemicalPotentialReference.logFugacity, CandidatePhase.AQUEOUS));
    }
    List<Candidate> candidates = distinctStationaryPoints(trials);
    Collections.sort(candidates, new Comparator<Candidate>() {
      @Override
      public int compare(Candidate first, Candidate second) {
        return Double.compare(first.getTangentPlaneDistance(), second.getTangentPlaneDistance());
      }
    });
    return new Result(trials, candidates, chemicalPotentialReference.maximumLogFugacitySpread, instabilityTolerance);
  }

  /**
   * Minimizes one requested trial phase from an explicit physical composition seed.
   *
   * <p>
   * This entry point is used by continuation algorithms that already know the incipient composition at the preceding
   * boundary point. It avoids losing a hydrocarbon stationary point when a generic Wilson seed is attracted to a much
   * deeper water-rich minimum.
   * </p>
   *
   * @param candidatePhase requested physical trial family
   * @param composition normalized or unnormalised non-negative mole-fraction seed
   * @return stationary trial result
   */
  public Candidate analyzeCandidate(CandidatePhase candidatePhase, double[] composition) {
    if (candidatePhase == null || composition == null) {
      throw new IllegalArgumentException("candidate phase and composition seed must be specified");
    }
    SystemInterface reference = sourceSystem.clone();
    if (reference == null) {
      throw new IllegalStateException("thermodynamic system clone failed");
    }
    reference.init(1);
    int componentCount = reference.getPhase(0).getNumberOfComponents();
    if (composition.length != componentCount) {
      throw new IllegalArgumentException("composition seed must contain one value per component");
    }
    double[] normalized = composition.clone();
    for (double value : normalized) {
      if (!Double.isFinite(value) || value < 0.0) {
        throw new IllegalArgumentException("composition seed contains an invalid value");
      }
    }
    normalize(normalized);
    double[] initialLogW = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      initialLogW[componentIndex] = Math.log(Math.max(normalized[componentIndex], MINIMUM_COMPOSITION));
    }
    Candidate fixedPoint = minimize(reference, chemicalPotentialReference(reference).logFugacity, candidatePhase,
        initialLogW);
    if (fixedPoint.isConverged()) {
      return fixedPoint;
    }

    double[] newtonSeed = usableComposition(fixedPoint.getComposition()) ? fixedPoint.getComposition() : normalized;
    IncipientPhaseStationaryPointSolver.Result newton = new IncipientPhaseStationaryPointSolver(reference,
        candidatePhase).setNumericalControls(maximumIterations, convergenceTolerance, 2.0e-5).solve(newtonSeed);
    if (newton.isConverged() || newton.getStationarityResidual() < fixedPoint.getStationarityResidual()) {
      double tangentPlaneDistance = newton.getTangentPlaneDistance();
      double stabilityFunction = Double.isFinite(tangentPlaneDistance) ? 1.0 - Math.exp(-tangentPlaneDistance)
          : Double.NaN;
      return new Candidate(candidatePhase, newton.getPhysicalPhase(), newton.getComposition(), stabilityFunction,
          tangentPlaneDistance, newton.getStationarityResidual(), newton.getIterations(), newton.isConverged(),
          newton.isTrivial(), newton.getFailureMessage());
    }
    return fixedPoint;
  }

  /**
   * Evaluates an explicit normalized trial composition without minimizing or changing stationary-root identity.
   *
   * <p>
   * At a zero-TPD boundary root every component value {@code ln(x_i) + ln(phi_i) - ln(f_i^ref)} is zero. At a non-zero
   * stationary point the values are all equal to the normalized TPD. Reporting both their composition-weighted mean and
   * their departure from that mean separates a displaced boundary from an iteration that converged to another
   * stationary root.
   * </p>
   *
   * @param candidatePhase requested physical trial family
   * @param composition normalized or unnormalised non-negative mole-fraction seed
   * @return immutable direct seed evidence
   */
  public SeedEvaluation evaluateCandidateSeed(CandidatePhase candidatePhase, double[] composition) {
    return evaluateCandidateSeed(candidatePhase, composition, 1);
  }

  SeedEvaluation evaluateCandidateSeed(CandidatePhase candidatePhase, double[] composition, int trialPhaseIndex) {
    if (candidatePhase == null || composition == null) {
      throw new IllegalArgumentException("candidate phase and composition seed must be specified");
    }
    if (trialPhaseIndex < 1) {
      throw new IllegalArgumentException("trial phase index must be positive");
    }
    SystemInterface reference = sourceSystem.clone();
    if (reference == null) {
      throw new IllegalStateException("thermodynamic system clone failed");
    }
    reference.init(1);
    int componentCount = reference.getPhase(0).getNumberOfComponents();
    if (composition.length != componentCount) {
      throw new IllegalArgumentException("composition seed must contain one value per component");
    }
    double[] normalized = composition.clone();
    for (double value : normalized) {
      if (!Double.isFinite(value) || value < 0.0) {
        throw new IllegalArgumentException("composition seed contains an invalid value");
      }
    }
    normalize(normalized);
    ChemicalPotentialReference chemicalReference = chemicalPotentialReference(reference);
    SystemInterface trial = reference.clone();
    if (trial.getMaxNumberOfPhases() <= trialPhaseIndex) {
      trial.setMaxNumberOfPhases(trialPhaseIndex + 1);
    }
    if (trial.getNumberOfPhases() <= trialPhaseIndex) {
      trial.setNumberOfPhases(trialPhaseIndex + 1);
    }
    trial.setPhaseType(trialPhaseIndex, candidatePhase.getNeqsimPhaseType());
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      trial.getPhase(trialPhaseIndex).getComponent(componentIndex)
          .setx(Math.max(normalized[componentIndex], MINIMUM_COMPOSITION));
    }
    trial.getPhase(trialPhaseIndex).normalize();
    trial.setPhaseType(trialPhaseIndex, candidatePhase.getNeqsimPhaseType());
    trial.init(1, trialPhaseIndex);

    double[] componentValues = new double[componentCount];
    double tangentPlaneDistance = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      componentValues[componentIndex] = Math.log(Math.max(normalized[componentIndex], MINIMUM_COMPOSITION))
          + trial.getPhase(trialPhaseIndex).getComponent(componentIndex).getLogFugacityCoefficient()
          - chemicalReference.logFugacity[componentIndex];
      tangentPlaneDistance += normalized[componentIndex] * componentValues[componentIndex];
    }
    double maximumBoundaryResidual = 0.0;
    double maximumStationarityDeparture = 0.0;
    for (double value : componentValues) {
      maximumBoundaryResidual = Math.max(maximumBoundaryResidual, Math.abs(value));
      maximumStationarityDeparture = Math.max(maximumStationarityDeparture, Math.abs(value - tangentPlaneDistance));
    }
    return new SeedEvaluation(candidatePhase, normalized, componentValues, tangentPlaneDistance,
        maximumBoundaryResidual, maximumStationarityDeparture, chemicalReference.maximumLogFugacitySpread);
  }

  private Candidate minimize(SystemInterface reference, double[] logReferenceFugacity, CandidatePhase candidatePhase) {
    return minimize(reference, logReferenceFugacity, candidatePhase, initialLogMoles(reference, candidatePhase));
  }

  private Candidate minimize(SystemInterface reference, double[] logReferenceFugacity, CandidatePhase candidatePhase,
      double[] initialLogW) {
    int componentCount = reference.getPhase(0).getNumberOfComponents();
    SystemInterface trialSystem = reference.clone();
    if (trialSystem == null) {
      return Candidate.failure(candidatePhase, componentCount, "thermodynamic system clone failed");
    }
    if (trialSystem.getNumberOfPhases() < 2) {
      trialSystem.setNumberOfPhases(2);
    }
    int trialPhaseIndex = 1;
    trialSystem.setPhaseType(trialPhaseIndex, candidatePhase.getNeqsimPhaseType());
    double[] logW = initialLogW.clone();
    double[] targetLogW = new double[componentCount];
    double residual = Double.POSITIVE_INFINITY;
    int iteration = 0;
    String failureMessage = null;

    while (iteration < maximumIterations) {
      iteration++;
      setNormalizedTrialComposition(trialSystem, trialPhaseIndex, logW);
      trialSystem.setPhaseType(trialPhaseIndex, candidatePhase.getNeqsimPhaseType());
      try {
        trialSystem.init(1, trialPhaseIndex);
      } catch (RuntimeException error) {
        failureMessage = error.getMessage();
        break;
      }

      residual = 0.0;
      for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
        double overall = reference.getPhase(0).getComponent(componentIndex).getz();
        if (overall <= MINIMUM_COMPOSITION) {
          targetLogW[componentIndex] = MINIMUM_LOG_VALUE;
          continue;
        }
        double logPhi = trialSystem.getPhase(trialPhaseIndex).getComponent(componentIndex).getLogFugacityCoefficient();
        if (!Double.isFinite(logPhi) || !Double.isFinite(logReferenceFugacity[componentIndex])) {
          failureMessage = "non-finite fugacity coefficient for component "
              + reference.getPhase(0).getComponent(componentIndex).getComponentName();
          residual = Double.POSITIVE_INFINITY;
          break;
        }
        targetLogW[componentIndex] = clamp(logReferenceFugacity[componentIndex] - logPhi);
        residual = Math.max(residual, Math.abs(targetLogW[componentIndex] - logW[componentIndex]));
      }
      if (failureMessage != null || residual <= convergenceTolerance) {
        break;
      }
      for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
        logW[componentIndex] = clamp(
            logW[componentIndex] + dampingFactor * (targetLogW[componentIndex] - logW[componentIndex]));
      }
    }

    boolean converged = failureMessage == null && residual <= convergenceTolerance;
    if (!converged && failureMessage == null) {
      failureMessage = "maximum fixed-point iteration count reached";
    }
    if (converged) {
      System.arraycopy(targetLogW, 0, logW, 0, componentCount);
      setNormalizedTrialComposition(trialSystem, trialPhaseIndex, logW);
      trialSystem.setPhaseType(trialPhaseIndex, candidatePhase.getNeqsimPhaseType());
      try {
        trialSystem.init(1, trialPhaseIndex);
      } catch (RuntimeException error) {
        converged = false;
        failureMessage = error.getMessage();
      }
    }

    double[] composition = trialComposition(trialSystem, trialPhaseIndex, componentCount);
    double sumW = sumExponentials(reference, logW);
    double stabilityFunction = Double.isFinite(sumW) ? 1.0 - sumW : Double.NEGATIVE_INFINITY;
    double tangentPlaneDistance = sumW > 0.0 && Double.isFinite(sumW) ? -Math.log(sumW) : Double.NEGATIVE_INFINITY;
    boolean trivial = isTrivial(reference, composition);
    CandidatePhase physicalPhase = classifyPhysicalPhase(trialSystem.getPhase(trialPhaseIndex));
    return new Candidate(candidatePhase, physicalPhase, composition, stabilityFunction, tangentPlaneDistance, residual,
        iteration, converged, trivial, failureMessage);
  }

  /**
   * Reads the evaluated EOS root family before interpreting liquid composition.
   *
   * <p>
   * SystemThermo relabels gas roots in nonzero phase slots as oil after initialization. Stability workspaces use slot
   * one for trials, so recover the PhaseEos volume criterion instead of trusting that storage-slot label.
   *
   * @param phase evaluated trial phase
   * @return gas, oil or aqueous identity assigned by the EOS
   * @throws IllegalArgumentException when the EOS returned an unsupported family
   */
  static CandidatePhase classifyPhysicalPhase(neqsim.thermo.phase.PhaseInterface phase) {
    if (phase instanceof neqsim.thermo.phase.PhaseEosInterface) {
      double covolume = ((neqsim.thermo.phase.PhaseEosInterface) phase).getB();
      if (phase.getVolume() / covolume > 1.75) {
        return CandidatePhase.GAS;
      }
      return CandidatePhase.valueOf(neqsim.thermo.phase.LiquidPhaseClassification.classify(phase).name());
    }
    return CandidatePhase.valueOf(phase.getType().name());
  }

  private static List<Candidate> distinctStationaryPoints(List<Candidate> trials) {
    List<Candidate> distinct = new ArrayList<Candidate>();
    for (Candidate trial : trials) {
      int equivalentIndex = equivalentStationaryPointIndex(distinct, trial);
      if (equivalentIndex < 0) {
        distinct.add(trial);
      } else if (preferDuplicate(trial, distinct.get(equivalentIndex))) {
        distinct.set(equivalentIndex, trial);
      }
    }
    return distinct;
  }

  private static int equivalentStationaryPointIndex(List<Candidate> candidates, Candidate trial) {
    if (!trial.isConverged()) {
      return -1;
    }
    for (int candidateIndex = 0; candidateIndex < candidates.size(); candidateIndex++) {
      Candidate candidate = candidates.get(candidateIndex);
      if (!candidate.isConverged() || Math.abs(candidate.getTangentPlaneDistance()
          - trial.getTangentPlaneDistance()) > EQUIVALENT_STATIONARY_POINT_TOLERANCE) {
        continue;
      }
      double compositionDistance = 0.0;
      for (int componentIndex = 0; componentIndex < trial.composition.length; componentIndex++) {
        compositionDistance += Math.abs(candidate.composition[componentIndex] - trial.composition[componentIndex]);
      }
      if (compositionDistance <= EQUIVALENT_STATIONARY_POINT_TOLERANCE) {
        return candidateIndex;
      }
    }
    return -1;
  }

  private static boolean preferDuplicate(Candidate candidate, Candidate current) {
    boolean candidateSeedMatchesPhysicalPhase = candidate.seedPhase == candidate.phase;
    boolean currentSeedMatchesPhysicalPhase = current.seedPhase == current.phase;
    if (candidateSeedMatchesPhysicalPhase != currentSeedMatchesPhysicalPhase) {
      return candidateSeedMatchesPhysicalPhase;
    }
    return candidate.stationarityResidual < current.stationarityResidual;
  }

  private static ChemicalPotentialReference chemicalPotentialReference(SystemInterface reference) {
    int componentCount = reference.getPhase(0).getNumberOfComponents();
    double[] logFugacity = new double[componentCount];
    double maximumSpread = 0.0;
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      double largestComposition = -1.0;
      double selected = Double.NaN;
      double minimum = Double.POSITIVE_INFINITY;
      double maximum = Double.NEGATIVE_INFINITY;
      for (int phaseIndex = 0; phaseIndex < reference.getNumberOfPhases(); phaseIndex++) {
        ComponentInterface component = reference.getPhase(phaseIndex).getComponent(componentIndex);
        double composition = component.getx();
        double value = composition > MINIMUM_COMPOSITION ? Math.log(composition) + component.getLogFugacityCoefficient()
            : Double.NaN;
        if (Double.isFinite(value)) {
          minimum = Math.min(minimum, value);
          maximum = Math.max(maximum, value);
          if (composition > largestComposition) {
            largestComposition = composition;
            selected = value;
          }
        }
      }
      if (!Double.isFinite(selected)) {
        throw new IllegalStateException("no finite reference fugacity for component "
            + reference.getPhase(0).getComponent(componentIndex).getComponentName());
      }
      logFugacity[componentIndex] = selected;
      if (Double.isFinite(minimum) && Double.isFinite(maximum)) {
        maximumSpread = Math.max(maximumSpread, maximum - minimum);
      }
    }
    return new ChemicalPotentialReference(logFugacity, maximumSpread);
  }

  private static double[] initialLogMoles(SystemInterface reference, CandidatePhase candidatePhase) {
    int componentCount = reference.getPhase(0).getNumberOfComponents();
    double[] weights = new double[componentCount];
    double temperatureK = reference.getTemperature();
    double pressureBara = reference.getPressure();
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      ComponentInterface component = reference.getPhase(0).getComponent(componentIndex);
      double overall = Math.max(component.getz(), MINIMUM_COMPOSITION);
      double criticalTemperature = component.getTC();
      double criticalPressure = component.getPC();
      double acentricFactor = component.getAcentricFactor();
      double wilsonK = criticalPressure / pressureBara
          * Math.exp(5.373 * (1.0 + acentricFactor) * (1.0 - criticalTemperature / temperatureK));
      wilsonK = Math.max(1.0e-20, Math.min(1.0e20, wilsonK));
      if (candidatePhase == CandidatePhase.GAS) {
        weights[componentIndex] = overall * wilsonK;
      } else if (candidatePhase == CandidatePhase.OIL) {
        weights[componentIndex] = overall / wilsonK;
      } else {
        String name = component.getComponentName().toLowerCase(Locale.ROOT);
        if (name.equals("water")) {
          weights[componentIndex] = Math.max(overall, 0.99);
        } else if (component.isHydrocarbon()) {
          weights[componentIndex] = overall * 1.0e-6;
        } else {
          weights[componentIndex] = overall;
        }
      }
    }
    normalize(weights);
    double[] logW = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      logW[componentIndex] = Math.log(Math.max(weights[componentIndex], MINIMUM_COMPOSITION));
    }
    return logW;
  }

  private static void setNormalizedTrialComposition(SystemInterface system, int phaseIndex, double[] logW) {
    double maximum = MINIMUM_LOG_VALUE;
    for (double value : logW) {
      maximum = Math.max(maximum, value);
    }
    double total = 0.0;
    double[] scaled = new double[logW.length];
    for (int componentIndex = 0; componentIndex < logW.length; componentIndex++) {
      scaled[componentIndex] = Math.exp(Math.max(MINIMUM_LOG_VALUE, logW[componentIndex] - maximum));
      total += scaled[componentIndex];
    }
    for (int componentIndex = 0; componentIndex < logW.length; componentIndex++) {
      system.getPhase(phaseIndex).getComponent(componentIndex)
          .setx(Math.max(MINIMUM_COMPOSITION, scaled[componentIndex] / total));
    }
    system.getPhase(phaseIndex).normalize();
  }

  private static double[] trialComposition(SystemInterface system, int phaseIndex, int componentCount) {
    double[] composition = new double[componentCount];
    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++) {
      composition[componentIndex] = system.getPhase(phaseIndex).getComponent(componentIndex).getx();
    }
    return composition;
  }

  private static double sumExponentials(SystemInterface reference, double[] logW) {
    double sum = 0.0;
    for (int componentIndex = 0; componentIndex < logW.length; componentIndex++) {
      if (reference.getPhase(0).getComponent(componentIndex).getz() > MINIMUM_COMPOSITION) {
        sum += Math.exp(clamp(logW[componentIndex]));
      }
    }
    return sum;
  }

  private static boolean isTrivial(SystemInterface reference, double[] composition) {
    for (int phaseIndex = 0; phaseIndex < reference.getNumberOfPhases(); phaseIndex++) {
      double distance = 0.0;
      for (int componentIndex = 0; componentIndex < composition.length; componentIndex++) {
        distance += Math
            .abs(composition[componentIndex] - reference.getPhase(phaseIndex).getComponent(componentIndex).getx());
      }
      if (distance < 1.0e-5) {
        return true;
      }
    }
    return false;
  }

  private static void normalize(double[] values) {
    double total = 0.0;
    for (double value : values) {
      total += value;
    }
    if (!(total > 0.0) || !Double.isFinite(total)) {
      throw new IllegalStateException("trial composition cannot be normalized");
    }
    for (int index = 0; index < values.length; index++) {
      values[index] /= total;
    }
  }

  private static boolean usableComposition(double[] values) {
    double total = 0.0;
    for (double value : values) {
      if (!Double.isFinite(value) || value < 0.0) {
        return false;
      }
      total += value;
    }
    return total > 0.0 && Double.isFinite(total);
  }

  private static double clamp(double value) {
    return Math.max(MINIMUM_LOG_VALUE, Math.min(MAXIMUM_LOG_VALUE, value));
  }

  private static final class ChemicalPotentialReference {
    private final double[] logFugacity;
    private final double maximumLogFugacitySpread;

    private ChemicalPotentialReference(double[] logFugacity, double maximumLogFugacitySpread) {
      this.logFugacity = logFugacity;
      this.maximumLogFugacitySpread = maximumLogFugacitySpread;
    }
  }

  /** Direct, non-iterative TPD evidence for one explicit trial composition. */
  public static final class SeedEvaluation {
    private final CandidatePhase phase;
    private final double[] composition;
    private final double[] componentValues;
    private final double tangentPlaneDistance;
    private final double maximumBoundaryResidual;
    private final double maximumStationarityDeparture;
    private final double maximumReferenceLogFugacitySpread;

    private SeedEvaluation(CandidatePhase phase, double[] composition, double[] componentValues,
        double tangentPlaneDistance, double maximumBoundaryResidual, double maximumStationarityDeparture,
        double maximumReferenceLogFugacitySpread) {
      this.phase = phase;
      this.composition = composition.clone();
      this.componentValues = componentValues.clone();
      this.tangentPlaneDistance = tangentPlaneDistance;
      this.maximumBoundaryResidual = maximumBoundaryResidual;
      this.maximumStationarityDeparture = maximumStationarityDeparture;
      this.maximumReferenceLogFugacitySpread = maximumReferenceLogFugacitySpread;
    }

    public CandidatePhase getPhase() {
      return phase;
    }

    public double[] getComposition() {
      return composition.clone();
    }

    public double[] getComponentValues() {
      return componentValues.clone();
    }

    public double getTangentPlaneDistance() {
      return tangentPlaneDistance;
    }

    public double getMaximumBoundaryResidual() {
      return maximumBoundaryResidual;
    }

    public double getMaximumStationarityDeparture() {
      return maximumStationarityDeparture;
    }

    public double getMaximumReferenceLogFugacitySpread() {
      return maximumReferenceLogFugacitySpread;
    }
  }

  /** One stationary trial-phase result. */
  public static final class Candidate {
    private final CandidatePhase seedPhase;
    private final CandidatePhase phase;
    private final double[] composition;
    private final double stabilityFunction;
    private final double tangentPlaneDistance;
    private final double stationarityResidual;
    private final int iterations;
    private final boolean converged;
    private final boolean trivial;
    private final String failureMessage;

    private Candidate(CandidatePhase seedPhase, CandidatePhase phase, double[] composition, double stabilityFunction,
        double tangentPlaneDistance, double stationarityResidual, int iterations, boolean converged, boolean trivial,
        String failureMessage) {
      this.seedPhase = seedPhase;
      this.phase = phase;
      this.composition = composition.clone();
      this.stabilityFunction = stabilityFunction;
      this.tangentPlaneDistance = tangentPlaneDistance;
      this.stationarityResidual = stationarityResidual;
      this.iterations = iterations;
      this.converged = converged;
      this.trivial = trivial;
      this.failureMessage = failureMessage;
    }

    private static Candidate failure(CandidatePhase phase, int componentCount, String failureMessage) {
      return new Candidate(phase, phase, new double[componentCount], Double.NaN, Double.NaN, Double.POSITIVE_INFINITY,
          0, false, false, failureMessage);
    }

    /** @return physical family used to seed the stationary-point iteration */
    public CandidatePhase getSeedPhase() {
      return seedPhase;
    }

    /** @return physical family of the trial phase */
    public CandidatePhase getPhase() {
      return phase;
    }

    /** @return normalized stationary trial composition */
    public double[] getComposition() {
      return composition.clone();
    }

    /** @return Michelsen stability function {@code 1 - sum(W)} */
    public double getStabilityFunction() {
      return stabilityFunction;
    }

    /** @return normalized tangent-plane distance, equal to {@code -ln(sum(W))} at a stationary point */
    public double getTangentPlaneDistance() {
      return tangentPlaneDistance;
    }

    /** @return maximum logarithmic stationarity residual */
    public double getStationarityResidual() {
      return stationarityResidual;
    }

    /** @return iteration count */
    public int getIterations() {
      return iterations;
    }

    /** @return true when the fixed-point equations converged */
    public boolean isConverged() {
      return converged;
    }

    /** @return true when the stationary composition duplicates an existing phase */
    public boolean isTrivial() {
      return trivial;
    }

    /** @return failure diagnostic, or {@code null} */
    public String getFailureMessage() {
      return failureMessage;
    }
  }

  /** Immutable collection of all physical trial results. */
  public static final class Result {
    private final List<Candidate> trials;
    private final List<Candidate> candidates;
    private final double maximumReferenceLogFugacitySpread;
    private final double instabilityTolerance;

    private Result(List<Candidate> trials, List<Candidate> candidates, double maximumReferenceLogFugacitySpread,
        double instabilityTolerance) {
      this.trials = Collections.unmodifiableList(new ArrayList<Candidate>(trials));
      this.candidates = Collections.unmodifiableList(new ArrayList<Candidate>(candidates));
      this.maximumReferenceLogFugacitySpread = maximumReferenceLogFugacitySpread;
      this.instabilityTolerance = instabilityTolerance;
    }

    /** @return candidates sorted from lowest to highest tangent-plane distance */
    public List<Candidate> getCandidates() {
      return candidates;
    }

    /** @return all seed trials before equivalent stationary points are collapsed */
    public List<Candidate> getTrials() {
      return trials;
    }

    /** @return largest log-fugacity inconsistency between existing equilibrium phases */
    public double getMaximumReferenceLogFugacitySpread() {
      return maximumReferenceLogFugacitySpread;
    }

    /**
     * Returns the lowest-TPD distinct stationary point for one physical phase family.
     *
     * @param phase requested gas, oil, or aqueous family
     * @return converged candidate, or {@code null} when no candidate in that family converged
     */
    public Candidate getCandidate(CandidatePhase phase) {
      for (Candidate candidate : candidates) {
        if (candidate.isConverged() && candidate.getPhase() == phase) {
          return candidate;
        }
      }
      return null;
    }

    /** @return lowest-TPD converged, non-trivial stationary point, regardless of its sign */
    public Candidate getLowestTangentPlaneDistanceCandidate() {
      for (Candidate candidate : candidates) {
        if (candidate.isConverged() && !candidate.isTrivial()) {
          return candidate;
        }
      }
      return null;
    }

    /** @return most unstable converged, non-trivial candidate or {@code null} */
    public Candidate getMostUnstableCandidate() {
      Candidate candidate = getLowestTangentPlaneDistanceCandidate();
      return candidate != null && candidate.getStabilityFunction() < instabilityTolerance ? candidate : null;
    }

    /** @return true when an additional thermodynamic phase is indicated */
    public boolean hasNonTrivialInstability() {
      return getMostUnstableCandidate() != null;
    }
  }
}
