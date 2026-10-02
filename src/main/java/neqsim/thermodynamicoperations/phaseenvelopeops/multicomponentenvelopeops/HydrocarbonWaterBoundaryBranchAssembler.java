package neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.AnchorPoint;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.BoundaryFamily;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.Branch;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterBoundaryAnchorDiscoverer.EndpointCandidate;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.HydrocarbonWaterPhaseTopology.SpecialPointType;
import neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.TwoToThreePhaseBoundaryQualityGate.EvidencePoint;

/** Extends clustered strict anchors in both directions while preserving isolated and failed branches. */
public final class HydrocarbonWaterBoundaryBranchAssembler {
  private final SystemInterface template;
  private int maximumCorrectorIterations = 80;
  private double residualTolerance = 1.0e-8;
  private double finiteDifferenceStep = 2.0e-5;
  private double initialArcStep = 0.25;
  private double minimumArcStep = 0.01;
  private double maximumArcStep = 0.5;
  private int maximumRetriesPerPoint = 10;
  private double maximumCompositionJump = 0.35;
  private double maximumInternalLogPressureJump = 0.75;
  private double maximumInternalRelativeTemperatureJump = 0.20;
  private int maximumInternalBridgeDepth = 12;
  private int minimumRegularPressurePrimaryAnchors = 3;
  private double maximumRegularPressurePrimaryLogPressureGap = 0.12;
  private double maximumRegularPressurePrimaryRelativeTemperatureGap = 0.08;
  private boolean singleAnchorContinuationEnabled;
  private boolean regularPressureFallbackEnabled = true;
  private boolean domainBoundsEnabled;
  private double minimumTemperatureK;
  private double maximumTemperatureK;
  private double minimumPressureBara;
  private double maximumPressureBara;

  /** Creates a non-destructive branch assembler. */
  public HydrocarbonWaterBoundaryBranchAssembler(SystemInterface template) {
    if (template == null) {
      throw new IllegalArgumentException("thermodynamic template is required");
    }
    this.template = template.clone();
  }

  /** Sets Newton controls for every pseudo-arclength correction. */
  public HydrocarbonWaterBoundaryBranchAssembler setCorrectorControls(int maximumIterations, double residualTolerance,
      double finiteDifferenceStep) {
    if (maximumIterations < 1 || !Double.isFinite(residualTolerance) || residualTolerance <= 0.0
        || !Double.isFinite(finiteDifferenceStep) || finiteDifferenceStep <= 0.0) {
      throw new IllegalArgumentException("invalid branch-assembler corrector controls");
    }
    this.maximumCorrectorIterations = maximumIterations;
    this.residualTolerance = residualTolerance;
    this.finiteDifferenceStep = finiteDifferenceStep;
    return this;
  }

  /** Sets absolute transformed-space continuation and retry controls. */
  public HydrocarbonWaterBoundaryBranchAssembler setStepControls(double initialArcStep, double minimumArcStep,
      double maximumArcStep, int maximumRetriesPerPoint, double maximumCompositionJump) {
    if (!Double.isFinite(initialArcStep) || !Double.isFinite(minimumArcStep) || !Double.isFinite(maximumArcStep)
        || minimumArcStep <= 0.0 || maximumArcStep < minimumArcStep || initialArcStep < minimumArcStep
        || initialArcStep > maximumArcStep || maximumRetriesPerPoint < 1 || !Double.isFinite(maximumCompositionJump)
        || maximumCompositionJump <= 0.0) {
      throw new IllegalArgumentException("invalid branch-assembler continuation controls");
    }
    this.initialArcStep = initialArcStep;
    this.minimumArcStep = minimumArcStep;
    this.maximumArcStep = maximumArcStep;
    this.maximumRetriesPerPoint = maximumRetriesPerPoint;
    this.maximumCompositionJump = maximumCompositionJump;
    return this;
  }

  /** Sets strict fixed-pressure subdivision controls between sparse discovered anchors. */
  public HydrocarbonWaterBoundaryBranchAssembler setInternalBridgeControls(double maximumLogPressureJump,
      double maximumRelativeTemperatureJump, int maximumBridgeDepth) {
    if (!Double.isFinite(maximumLogPressureJump) || maximumLogPressureJump <= 0.0
        || !Double.isFinite(maximumRelativeTemperatureJump) || maximumRelativeTemperatureJump <= 0.0
        || maximumBridgeDepth < 1) {
      throw new IllegalArgumentException("invalid internal anchor-bridge controls");
    }
    this.maximumInternalLogPressureJump = maximumLogPressureJump;
    this.maximumInternalRelativeTemperatureJump = maximumRelativeTemperatureJump;
    this.maximumInternalBridgeDepth = maximumBridgeDepth;
    return this;
  }

  /**
   * Enables experimental continuation from a single discovered anchor.
   *
   * <p>
   * This is deliberately disabled by default. A singleton may be a genuine branch point, a special endpoint, or a
   * metastable numerical root. Production assembly therefore preserves it as isolated until an independent topology
   * classifier establishes that continuation is physically appropriate.
   */
  public HydrocarbonWaterBoundaryBranchAssembler setSingleAnchorContinuationEnabled(boolean enabled) {
    this.singleAnchorContinuationEnabled = enabled;
    return this;
  }

  /** Enables stable-region fixed-pressure continuation after a purely numerical pseudo-arclength start failure. */
  public HydrocarbonWaterBoundaryBranchAssembler setRegularPressureFallbackEnabled(boolean enabled) {
    this.regularPressureFallbackEnabled = enabled;
    return this;
  }

  /** Sets the evidence density required before pressure is selected as the primary local branch parameter. */
  public HydrocarbonWaterBoundaryBranchAssembler setRegularPressurePrimaryControls(int minimumAnchorCount,
      double maximumLogPressureGap, double maximumRelativeTemperatureGap) {
    if (minimumAnchorCount < 3 || !Double.isFinite(maximumLogPressureGap) || maximumLogPressureGap <= 0.0
        || !Double.isFinite(maximumRelativeTemperatureGap) || maximumRelativeTemperatureGap <= 0.0) {
      throw new IllegalArgumentException("invalid regular-pressure primary selection controls");
    }
    this.minimumRegularPressurePrimaryAnchors = minimumAnchorCount;
    this.maximumRegularPressurePrimaryLogPressureGap = maximumLogPressureGap;
    this.maximumRegularPressurePrimaryRelativeTemperatureGap = maximumRelativeTemperatureGap;
    return this;
  }

  /** Sets the declared PT calculation domain used to stop each continuation direction after its first crossing. */
  public HydrocarbonWaterBoundaryBranchAssembler setDomainBounds(double minimumTemperatureK, double maximumTemperatureK,
      double minimumPressureBara, double maximumPressureBara) {
    if (!Double.isFinite(minimumTemperatureK) || !Double.isFinite(maximumTemperatureK)
        || maximumTemperatureK <= minimumTemperatureK || !Double.isFinite(minimumPressureBara)
        || minimumPressureBara <= 0.0 || !Double.isFinite(maximumPressureBara)
        || maximumPressureBara <= minimumPressureBara) {
      throw new IllegalArgumentException("invalid branch-assembler PT domain bounds");
    }
    this.minimumTemperatureK = minimumTemperatureK;
    this.maximumTemperatureK = maximumTemperatureK;
    this.minimumPressureBara = minimumPressureBara;
    this.maximumPressureBara = maximumPressureBara;
    this.domainBoundsEnabled = true;
    return this;
  }

  /**
   * Extends every continuation-capable discovered branch from both anchor ends.
   *
   * @param discovery strict anchor discovery result
   * @param additionalPointsPerDirection requested corrected points beyond each end anchor
   * @return assembled branches; the output contains exactly one entry for every discovered branch
   */
  public Result assemble(HydrocarbonWaterBoundaryAnchorDiscoverer.Result discovery, int additionalPointsPerDirection) {
    if (discovery == null || additionalPointsPerDirection < 1) {
      throw new IllegalArgumentException("anchor discovery and a positive extension count are required");
    }
    Result branches = assembleBranches(discovery.getBranches(), additionalPointsPerDirection);
    return new Result(branches.branches, discovery.getEndpointCandidates());
  }

  /** Extends a previously clustered strict-anchor inventory while preserving every supplied branch. */
  public Result assembleBranches(List<Branch> discoveredBranches, int additionalPointsPerDirection) {
    if (discoveredBranches == null || additionalPointsPerDirection < 1) {
      throw new IllegalArgumentException("discovered branches and a positive extension count are required");
    }
    List<AssembledBranch> assembled = new ArrayList<AssembledBranch>();
    for (Branch branch : discoveredBranches) {
      if (branch == null) {
        throw new IllegalArgumentException("discovered branch inventory must not contain null entries");
      }
      if (!branch.canSeedContinuation()) {
        AssembledBranch recovered = singleAnchorContinuationEnabled && branch.getPoints().size() == 1
            ? assembleSingleAnchorBranch(branch, additionalPointsPerDirection)
            : null;
        assembled.add(recovered == null ? AssembledBranch.isolated(branch) : recovered);
        continue;
      }
      BoundaryFamily family = branch.getFamily();
      List<AnchorPoint> anchors = branch.getPoints();
      int lastIndex = anchors.size() - 1;
      boolean regularPressurePrimary = isRegularPressurePrimarySequence(anchors);
      TwoToThreePhasePseudoArcLengthTracer.Result backward = regularPressurePrimary ? null
          : tracer(family).trace(anchors.get(1).toContinuationState(), anchors.get(0).toContinuationState(),
              additionalPointsPerDirection);
      TwoToThreePhasePseudoArcLengthTracer.Result forward = regularPressurePrimary ? null
          : tracer(family).trace(anchors.get(lastIndex - 1).toContinuationState(),
              anchors.get(lastIndex).toContinuationState(), additionalPointsPerDirection);
      requireCompleteAcceptedGlobalStabilityEvidenceIfPresent(backward);
      requireCompleteAcceptedGlobalStabilityEvidenceIfPresent(forward);
      HydrocarbonWaterRegularPressureBoundaryTracer.Result regularBackward = regularPressurePrimary
          || shouldUseRegularPressureFallback(backward)
              ? regularPressureTracer(family).trace(anchors.get(1), anchors.get(0), additionalPointsPerDirection)
              : null;
      HydrocarbonWaterRegularPressureBoundaryTracer.Result regularForward = regularPressurePrimary
          || shouldUseRegularPressureFallback(forward)
              ? regularPressureTracer(family).trace(anchors.get(lastIndex - 1), anchors.get(lastIndex),
                  additionalPointsPerDirection)
              : null;
      TwoToThreePhasePseudoArcLengthTracer.Result postRegularBackward = continueAfterPressureParameterLimit(family,
          regularBackward, additionalPointsPerDirection);
      TwoToThreePhasePseudoArcLengthTracer.Result postRegularForward = continueAfterPressureParameterLimit(family,
          regularForward, additionalPointsPerDirection);
      requireCompleteAcceptedGlobalStabilityEvidenceIfPresent(postRegularBackward);
      requireCompleteAcceptedGlobalStabilityEvidenceIfPresent(postRegularForward);
      HydrocarbonWaterBoundaryTerminationClassifier terminationClassifier = new HydrocarbonWaterBoundaryTerminationClassifier(
          template);
      HydrocarbonWaterBoundaryTerminationClassifier.Result backwardTermination = classifyTerminalTrace(
          terminationClassifier, backward, regularBackward, postRegularBackward, anchors.get(0).getBoundaryRoot());
      HydrocarbonWaterBoundaryTerminationClassifier.Result forwardTermination = classifyTerminalTrace(
          terminationClassifier, forward, regularForward, postRegularForward, anchors.get(lastIndex).getBoundaryRoot());
      List<EvidencePoint> evidence = new ArrayList<EvidencePoint>();
      if (backwardTermination != null && backwardTermination.isPhysicalEndpoint()) {
        evidence.add(EvidencePoint.from(backwardTermination.getBoundaryRoot()));
      }
      addPseudoBackwardEvidence(evidence, postRegularBackward);
      addRegularBackwardEvidence(evidence, regularBackward);
      if (backward != null) {
        List<TwoToThreePhaseArcLengthCorrector.Result> backwardCorrections = backward.getAcceptedCorrections();
        for (int correctionIndex = backwardCorrections.size() - 1; correctionIndex >= 0; correctionIndex--) {
          evidence.add(EvidencePoint.from(backwardCorrections.get(correctionIndex)));
        }
      }
      BridgeResult bridge = bridgeAnchors(family, anchors);
      evidence.addAll(bridge.evidence);
      if (forward != null) {
        for (TwoToThreePhaseArcLengthCorrector.Result correction : forward.getAcceptedCorrections()) {
          evidence.add(EvidencePoint.from(correction));
        }
      }
      addRegularForwardEvidence(evidence, regularForward);
      addPseudoForwardEvidence(evidence, postRegularForward);
      if (forwardTermination != null && forwardTermination.isPhysicalEndpoint()) {
        evidence.add(EvidencePoint.from(forwardTermination.getBoundaryRoot()));
      }
      assembled.add(new AssembledBranch(branch.getIdentifier(), family, anchors.size(), 0, bridge.insertedPointCount,
          evidence, backward, forward, regularBackward, regularForward, postRegularBackward, postRegularForward,
          backwardTermination, forwardTermination));
    }
    if (assembled.size() != discoveredBranches.size()) {
      throw new IllegalStateException("discovered branch preservation failed during continuation assembly");
    }
    return new Result(assembled, Collections.<EndpointCandidate>emptyList());
  }

  private static void requireCompleteAcceptedGlobalStabilityEvidence(
      TwoToThreePhasePseudoArcLengthTracer.Result trace) {
    long acceptedGlobal = trace.getGlobalStabilityEvidence().stream()
        .filter(HydrocarbonWaterBoundaryGlobalStabilityGate.Result::isAccepted).count();
    if (acceptedGlobal != trace.getAcceptedCorrections().size()) {
      throw new IllegalStateException("accepted continuation corrections lack complete global-stability evidence");
    }
  }

  private static void requireCompleteAcceptedGlobalStabilityEvidenceIfPresent(
      TwoToThreePhasePseudoArcLengthTracer.Result trace) {
    if (trace != null) {
      requireCompleteAcceptedGlobalStabilityEvidence(trace);
    }
  }

  private boolean shouldUseRegularPressureFallback(TwoToThreePhasePseudoArcLengthTracer.Result trace) {
    if (!regularPressureFallbackEnabled || trace == null || trace.hasCompletedRequestedPoints()) {
      return false;
    }
    // A trace that stops at its minimum arc step stalled for a numerical reason, not because the branch ended, so
    // the fixed-pressure tracer should still get a turn. Requiring an empty correction list meant the most common
    // failure - a few good points and then a stall - never reached the fallback at all, and the branch was reported
    // with an unattached end while three untried tracers sat idle.
    if (trace.getTerminationReason() == TwoToThreePhasePseudoArcLengthTracer.TerminationReason.MINIMUM_ARC_STEP) {
      return true;
    }
    if (!trace.getAcceptedCorrections().isEmpty()) {
      return false;
    }
    TwoToThreePhaseArcLengthCorrector.Result terminal = trace.getTerminalCorrection();
    return terminal == null || !terminal.isConverged();
  }

  private boolean isRegularPressurePrimarySequence(List<AnchorPoint> anchors) {
    if (!regularPressureFallbackEnabled || anchors.size() < minimumRegularPressurePrimaryAnchors) {
      return false;
    }
    for (int index = 1; index < anchors.size(); index++) {
      AnchorPoint previous = anchors.get(index - 1);
      AnchorPoint current = anchors.get(index);
      double logPressureGap = Math.log(current.getPressureBara() / previous.getPressureBara());
      double relativeTemperatureGap = Math.abs(current.getTemperatureK() - previous.getTemperatureK())
          / Math.max(1.0, Math.min(current.getTemperatureK(), previous.getTemperatureK()));
      if (!(logPressureGap > 0.0) || logPressureGap > maximumRegularPressurePrimaryLogPressureGap
          || relativeTemperatureGap > maximumRegularPressurePrimaryRelativeTemperatureGap
          || compositionJump(previous.toContinuationState(), current.toContinuationState()) > maximumCompositionJump) {
        return false;
      }
    }
    return true;
  }

  private static boolean hasNewRegularPoints(HydrocarbonWaterRegularPressureBoundaryTracer.Result trace) {
    return trace != null && trace.getPoints().size() > 2;
  }

  private HydrocarbonWaterRegularPressureBoundaryTracer regularPressureTracer(BoundaryFamily family) {
    HydrocarbonWaterRegularPressureBoundaryTracer tracer = new HydrocarbonWaterRegularPressureBoundaryTracer(template,
        family).setStepControls(0.04, 0.0025, 0.08, maximumRetriesPerPoint)
        .setCorrectionControls(11, 8.0, 0.04, maximumCompositionJump, 30.0);
    return domainBoundsEnabled
        ? tracer.setDomainBounds(minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara)
        : tracer;
  }

  private TwoToThreePhasePseudoArcLengthTracer.Result continueAfterPressureParameterLimit(BoundaryFamily family,
      HydrocarbonWaterRegularPressureBoundaryTracer.Result regular, int requestedPointCount) {
    if (regular == null || regular
        .getTerminationReason() != HydrocarbonWaterRegularPressureBoundaryTracer.TerminationReason.PRESSURE_PARAMETER_LIMIT
        || regular.getPoints().size() < 2) {
      return null;
    }
    int remainingPointCount = requestedPointCount - (regular.getPoints().size() - 2);
    if (remainingPointCount < 1) {
      return null;
    }
    List<AnchorPoint> points = regular.getPoints();
    AnchorPoint previous = points.get(points.size() - 2);
    AnchorPoint current = points.get(points.size() - 1);
    return tracer(family).trace(previous.toContinuationState(), current.toContinuationState(), remainingPointCount);
  }

  private static HydrocarbonWaterBoundaryTerminationClassifier.Result classifyTerminalTrace(
      HydrocarbonWaterBoundaryTerminationClassifier classifier, TwoToThreePhasePseudoArcLengthTracer.Result initial,
      HydrocarbonWaterRegularPressureBoundaryTracer.Result regular,
      TwoToThreePhasePseudoArcLengthTracer.Result postRegular, TwoToThreePhaseBoundaryPointSolver.Result anchorRoot) {
    if (postRegular != null) {
      AnchorPoint lastRegular = regular.getPoints().get(regular.getPoints().size() - 1);
      return classifier.classify(postRegular, lastRegular.getBoundaryRoot());
    }
    if (hasNewRegularPoints(regular)) {
      return null;
    }
    return initial == null ? null : classifier.classify(initial, anchorRoot);
  }

  private static void addRegularBackwardEvidence(List<EvidencePoint> evidence,
      HydrocarbonWaterRegularPressureBoundaryTracer.Result regular) {
    if (regular == null) {
      return;
    }
    for (int index = regular.getPoints().size() - 1; index >= 2; index--) {
      evidence.add(regular.getPoints().get(index).toEvidencePoint());
    }
  }

  private static void addRegularForwardEvidence(List<EvidencePoint> evidence,
      HydrocarbonWaterRegularPressureBoundaryTracer.Result regular) {
    if (regular == null) {
      return;
    }
    for (int index = 2; index < regular.getPoints().size(); index++) {
      evidence.add(regular.getPoints().get(index).toEvidencePoint());
    }
  }

  private static void addPseudoBackwardEvidence(List<EvidencePoint> evidence,
      TwoToThreePhasePseudoArcLengthTracer.Result trace) {
    if (trace == null) {
      return;
    }
    List<TwoToThreePhaseArcLengthCorrector.Result> corrections = trace.getAcceptedCorrections();
    for (int index = corrections.size() - 1; index >= 0; index--) {
      evidence.add(EvidencePoint.from(corrections.get(index)));
    }
  }

  private static void addPseudoForwardEvidence(List<EvidencePoint> evidence,
      TwoToThreePhasePseudoArcLengthTracer.Result trace) {
    if (trace == null) {
      return;
    }
    for (TwoToThreePhaseArcLengthCorrector.Result correction : trace.getAcceptedCorrections()) {
      evidence.add(EvidencePoint.from(correction));
    }
  }

  private AssembledBranch assembleSingleAnchorBranch(Branch branch, int additionalPointsPerDirection) {
    AnchorPoint anchor = branch.getPoints().get(0);
    GeneratedSeed generated = generateSecondSeed(branch.getFamily(), anchor);
    if (generated == null) {
      return null;
    }
    TwoToThreePhaseBoundaryPointSolver.Result anchorRoot = anchor.getBoundaryRoot();
    TwoToThreePhaseBoundaryPointSolver.Result firstRoot = anchorRoot.getPressureBara() <= generated.root
        .getPressureBara() ? anchorRoot : generated.root;
    TwoToThreePhaseBoundaryPointSolver.Result secondRoot = firstRoot == anchorRoot ? generated.root : anchorRoot;
    EvidencePoint firstEvidence = firstRoot == anchorRoot ? anchor.toEvidencePoint() : generated.evidence;
    EvidencePoint secondEvidence = firstRoot == anchorRoot ? generated.evidence : anchor.toEvidencePoint();
    TwoToThreePhasePseudoArcLengthTracer.Result backward = tracer(branch.getFamily()).trace(
        TwoToThreePhaseArcLengthCorrector.State.from(secondRoot),
        TwoToThreePhaseArcLengthCorrector.State.from(firstRoot), additionalPointsPerDirection);
    TwoToThreePhasePseudoArcLengthTracer.Result forward = tracer(branch.getFamily()).trace(
        TwoToThreePhaseArcLengthCorrector.State.from(firstRoot),
        TwoToThreePhaseArcLengthCorrector.State.from(secondRoot), additionalPointsPerDirection);
    requireCompleteAcceptedGlobalStabilityEvidence(backward);
    requireCompleteAcceptedGlobalStabilityEvidence(forward);
    HydrocarbonWaterBoundaryTerminationClassifier classifier = new HydrocarbonWaterBoundaryTerminationClassifier(
        template);
    HydrocarbonWaterBoundaryTerminationClassifier.Result backwardTermination = classifier.classify(backward, firstRoot);
    HydrocarbonWaterBoundaryTerminationClassifier.Result forwardTermination = classifier.classify(forward, secondRoot);
    List<EvidencePoint> evidence = new ArrayList<EvidencePoint>();
    if (backwardTermination.isPhysicalEndpoint()) {
      evidence.add(EvidencePoint.from(backwardTermination.getBoundaryRoot()));
    }
    List<TwoToThreePhaseArcLengthCorrector.Result> backwardCorrections = backward.getAcceptedCorrections();
    for (int index = backwardCorrections.size() - 1; index >= 0; index--) {
      evidence.add(EvidencePoint.from(backwardCorrections.get(index)));
    }
    evidence.add(firstEvidence);
    evidence.add(secondEvidence);
    for (TwoToThreePhaseArcLengthCorrector.Result correction : forward.getAcceptedCorrections()) {
      evidence.add(EvidencePoint.from(correction));
    }
    if (forwardTermination.isPhysicalEndpoint()) {
      evidence.add(EvidencePoint.from(forwardTermination.getBoundaryRoot()));
    }
    return new AssembledBranch(branch.getIdentifier(), branch.getFamily(), 1, 1, 0, evidence, backward, forward, null,
        null, null, null, backwardTermination, forwardTermination);
  }

  private GeneratedSeed generateSecondSeed(BoundaryFamily family, AnchorPoint anchor) {
    double[] logPressureOffsets = new double[] {0.125, -0.125, 0.25, -0.25, 0.0625, -0.0625, 0.5, -0.5};
    TwoToThreePhaseArcLengthCorrector corrector = new TwoToThreePhaseArcLengthCorrector(template,
        family.getRetainedPhaseZero(), family.getRetainedPhaseOne(), family.getIncipientPhase())
        .setNumericalControls(maximumCorrectorIterations, residualTolerance, finiteDifferenceStep);
    HydrocarbonWaterBoundaryGlobalStabilityGate stabilityGate = new HydrocarbonWaterBoundaryGlobalStabilityGate(
        template);
    TwoToThreePhaseArcLengthCorrector.State anchorState = anchor.toContinuationState();
    for (double offset : logPressureOffsets) {
      double pressureBara = anchorState.getPressureBara() * Math.exp(offset);
      if (domainBoundsEnabled && (pressureBara < minimumPressureBara || pressureBara > maximumPressureBara)) {
        continue;
      }
      TwoToThreePhaseArcLengthCorrector.Result correction = corrector.correctAtPressure(anchorState, pressureBara);
      GeneratedSeed generated = acceptedGeneratedSeed(anchorState, correction, stabilityGate);
      if (generated != null) {
        return generated;
      }
    }
    double[] temperatureOffsetsK = new double[] {2.0, -2.0, 5.0, -5.0, 1.0, -1.0, 10.0, -10.0, 20.0, -20.0};
    for (double offsetK : temperatureOffsetsK) {
      double temperatureK = anchorState.getTemperatureK() + offsetK;
      if (domainBoundsEnabled && (temperatureK < minimumTemperatureK || temperatureK > maximumTemperatureK)) {
        continue;
      }
      TwoToThreePhaseArcLengthCorrector.Result correction = corrector.correctAtTemperature(anchorState, temperatureK);
      GeneratedSeed generated = acceptedGeneratedSeed(anchorState, correction, stabilityGate);
      if (generated != null) {
        return generated;
      }
    }
    double[] localArcSteps = new double[] {0.05, 0.10, 0.025, 0.20, 0.01, 0.40};
    for (double arcStep : localArcSteps) {
      for (int orientation : new int[] {1, -1}) {
        TwoToThreePhaseArcLengthCorrector.Result correction = corrector.correctFromLocalTangent(anchorState, arcStep,
            orientation);
        GeneratedSeed generated = acceptedGeneratedSeed(anchorState, correction, stabilityGate);
        if (generated != null) {
          return generated;
        }
      }
    }
    return null;
  }

  private GeneratedSeed acceptedGeneratedSeed(TwoToThreePhaseArcLengthCorrector.State anchorState,
      TwoToThreePhaseArcLengthCorrector.Result correction, HydrocarbonWaterBoundaryGlobalStabilityGate stabilityGate) {
    if (!correction.isConverged() || compositionJump(anchorState, correction.getState()) > maximumCompositionJump
        || physicalDistance(anchorState, correction.getState()) <= 1.0e-8) {
      return null;
    }
    TwoToThreePhaseBoundaryPointSolver.Result root = TwoToThreePhaseBoundaryPointSolver.Result
        .fromContinuationState(correction.getState(), correction.getThermodynamicMaximumResidual());
    HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability = stabilityGate.evaluate(root);
    return stability.isAccepted() ? new GeneratedSeed(root, EvidencePoint.from(correction)) : null;
  }

  private static double physicalDistance(TwoToThreePhaseArcLengthCorrector.State first,
      TwoToThreePhaseArcLengthCorrector.State second) {
    double relativeTemperature = (second.getTemperatureK() - first.getTemperatureK())
        / Math.max(1.0, first.getTemperatureK());
    double logPressure = Math.log(second.getPressureBara() / first.getPressureBara());
    double composition = compositionJump(first, second);
    return Math.sqrt(relativeTemperature * relativeTemperature + logPressure * logPressure + composition * composition);
  }

  private TwoToThreePhasePseudoArcLengthTracer tracer(BoundaryFamily family) {
    TwoToThreePhasePseudoArcLengthTracer tracer = new TwoToThreePhasePseudoArcLengthTracer(template,
        family.getRetainedPhaseZero(), family.getRetainedPhaseOne(), family.getIncipientPhase())
        .setCorrectorControls(maximumCorrectorIterations, residualTolerance, finiteDifferenceStep)
        .setStepControls(initialArcStep, minimumArcStep, maximumArcStep, maximumRetriesPerPoint)
        .setMaximumCompositionJump(maximumCompositionJump);
    return domainBoundsEnabled
        ? tracer.setDomainBounds(minimumTemperatureK, maximumTemperatureK, minimumPressureBara, maximumPressureBara)
        : tracer;
  }

  private BridgeResult bridgeAnchors(BoundaryFamily family, List<AnchorPoint> anchors) {
    List<EvidencePoint> evidence = new ArrayList<EvidencePoint>();
    TwoToThreePhaseArcLengthCorrector corrector = new TwoToThreePhaseArcLengthCorrector(template,
        family.getRetainedPhaseZero(), family.getRetainedPhaseOne(), family.getIncipientPhase())
        .setNumericalControls(maximumCorrectorIterations, residualTolerance, finiteDifferenceStep);
    HydrocarbonWaterBoundaryGlobalStabilityGate stabilityGate = new HydrocarbonWaterBoundaryGlobalStabilityGate(
        template);
    int inserted = 0;
    evidence.add(anchors.get(0).toEvidencePoint());
    for (int index = 1; index < anchors.size(); index++) {
      List<EvidencePoint> subdivision = new ArrayList<EvidencePoint>();
      inserted += bridgeSegment(corrector, stabilityGate, anchors.get(index - 1).toContinuationState(),
          anchors.get(index).toContinuationState(), 0, subdivision);
      evidence.addAll(subdivision);
      evidence.add(anchors.get(index).toEvidencePoint());
    }
    return new BridgeResult(evidence, inserted);
  }

  private int bridgeSegment(TwoToThreePhaseArcLengthCorrector corrector,
      HydrocarbonWaterBoundaryGlobalStabilityGate stabilityGate, TwoToThreePhaseArcLengthCorrector.State first,
      TwoToThreePhaseArcLengthCorrector.State second, int depth, List<EvidencePoint> evidence) {
    if (!requiresInternalBridge(first, second)) {
      return 0;
    }
    if (depth >= maximumInternalBridgeDepth) {
      throw new IllegalStateException("strict internal anchor subdivision exceeded its depth limit");
    }
    double firstPressure = first.getPressureBara();
    double secondPressure = second.getPressureBara();
    if (Math.abs(Math.log(secondPressure / firstPressure)) <= 1.0e-12) {
      throw new IllegalStateException("sparse equal-pressure anchors require a temperature-parameter bridge");
    }
    double midpointPressure = Math.sqrt(firstPressure * secondPressure);
    TwoToThreePhaseArcLengthCorrector.Result midpoint = corrector.correctAtPressure(first, midpointPressure);
    if (!midpoint.isConverged()) {
      midpoint = corrector.correctAtPressure(second, midpointPressure);
    }
    if (!midpoint.isConverged()) {
      throw new IllegalStateException("strict internal anchor subdivision failed: " + midpoint.getFailureMessage());
    }
    TwoToThreePhaseBoundaryPointSolver.Result root = TwoToThreePhaseBoundaryPointSolver.Result
        .fromContinuationState(midpoint.getState(), midpoint.getThermodynamicMaximumResidual());
    HydrocarbonWaterBoundaryGlobalStabilityGate.Result stability = stabilityGate.evaluate(root);
    if (!stability.isAccepted()) {
      throw new IllegalStateException(
          "strict internal anchor subdivision failed global stability: " + stability.getFailureMessage());
    }
    int inserted = bridgeSegment(corrector, stabilityGate, first, midpoint.getState(), depth + 1, evidence);
    evidence.add(EvidencePoint.from(midpoint));
    inserted++;
    inserted += bridgeSegment(corrector, stabilityGate, midpoint.getState(), second, depth + 1, evidence);
    return inserted;
  }

  private boolean requiresInternalBridge(TwoToThreePhaseArcLengthCorrector.State first,
      TwoToThreePhaseArcLengthCorrector.State second) {
    double logPressureJump = Math.abs(Math.log(second.getPressureBara() / first.getPressureBara()));
    double relativeTemperatureJump = Math.abs(second.getTemperatureK() - first.getTemperatureK())
        / Math.max(1.0, Math.min(first.getTemperatureK(), second.getTemperatureK()));
    return logPressureJump > maximumInternalLogPressureJump
        || relativeTemperatureJump > maximumInternalRelativeTemperatureJump
        || compositionJump(first, second) > maximumCompositionJump;
  }

  private static double compositionJump(TwoToThreePhaseArcLengthCorrector.State first,
      TwoToThreePhaseArcLengthCorrector.State second) {
    return Math.max(compositionDistance(first.getPhaseZeroComposition(), second.getPhaseZeroComposition()),
        Math.max(compositionDistance(first.getPhaseOneComposition(), second.getPhaseOneComposition()),
            compositionDistance(first.getIncipientComposition(), second.getIncipientComposition())));
  }

  private static double compositionDistance(double[] first, double[] second) {
    double distance = 0.0;
    for (int index = 0; index < first.length; index++) {
      distance += Math.abs(first[index] - second[index]);
    }
    return distance;
  }

  private static final class BridgeResult {
    private final List<EvidencePoint> evidence;
    private final int insertedPointCount;

    private BridgeResult(List<EvidencePoint> evidence, int insertedPointCount) {
      this.evidence = evidence;
      this.insertedPointCount = insertedPointCount;
    }
  }

  private static final class GeneratedSeed {
    private final TwoToThreePhaseBoundaryPointSolver.Result root;
    private final EvidencePoint evidence;

    private GeneratedSeed(TwoToThreePhaseBoundaryPointSolver.Result root, EvidencePoint evidence) {
      this.root = root;
      this.evidence = evidence;
    }
  }

  /** One preserved branch with strict anchor and continuation evidence. */
  public static final class AssembledBranch {
    private final String identifier;
    private final BoundaryFamily family;
    private final int anchorCount;
    private final int generatedSeedPointCount;
    private final int internalBridgePointCount;
    private final List<EvidencePoint> evidencePoints;
    private final TwoToThreePhasePseudoArcLengthTracer.Result backwardTrace;
    private final TwoToThreePhasePseudoArcLengthTracer.Result forwardTrace;
    private final HydrocarbonWaterRegularPressureBoundaryTracer.Result regularBackwardTrace;
    private final HydrocarbonWaterRegularPressureBoundaryTracer.Result regularForwardTrace;
    private final TwoToThreePhasePseudoArcLengthTracer.Result postRegularBackwardTrace;
    private final TwoToThreePhasePseudoArcLengthTracer.Result postRegularForwardTrace;
    private final HydrocarbonWaterBoundaryTerminationClassifier.Result backwardTermination;
    private final HydrocarbonWaterBoundaryTerminationClassifier.Result forwardTermination;

    private AssembledBranch(String identifier, BoundaryFamily family, int anchorCount, int generatedSeedPointCount,
        int internalBridgePointCount, List<EvidencePoint> evidencePoints,
        TwoToThreePhasePseudoArcLengthTracer.Result backwardTrace,
        TwoToThreePhasePseudoArcLengthTracer.Result forwardTrace,
        HydrocarbonWaterRegularPressureBoundaryTracer.Result regularBackwardTrace,
        HydrocarbonWaterRegularPressureBoundaryTracer.Result regularForwardTrace,
        TwoToThreePhasePseudoArcLengthTracer.Result postRegularBackwardTrace,
        TwoToThreePhasePseudoArcLengthTracer.Result postRegularForwardTrace,
        HydrocarbonWaterBoundaryTerminationClassifier.Result backwardTermination,
        HydrocarbonWaterBoundaryTerminationClassifier.Result forwardTermination) {
      this.identifier = identifier;
      this.family = family;
      this.anchorCount = anchorCount;
      this.generatedSeedPointCount = generatedSeedPointCount;
      this.internalBridgePointCount = internalBridgePointCount;
      this.evidencePoints = Collections.unmodifiableList(new ArrayList<EvidencePoint>(evidencePoints));
      this.backwardTrace = backwardTrace;
      this.forwardTrace = forwardTrace;
      this.regularBackwardTrace = regularBackwardTrace;
      this.regularForwardTrace = regularForwardTrace;
      this.postRegularBackwardTrace = postRegularBackwardTrace;
      this.postRegularForwardTrace = postRegularForwardTrace;
      this.backwardTermination = backwardTermination;
      this.forwardTermination = forwardTermination;
    }

    private static AssembledBranch isolated(Branch branch) {
      List<EvidencePoint> evidence = new ArrayList<EvidencePoint>();
      for (AnchorPoint anchor : branch.getPoints()) {
        evidence.add(anchor.toEvidencePoint());
      }
      return new AssembledBranch(branch.getIdentifier(), branch.getFamily(), branch.getPoints().size(), 0, 0, evidence,
          null, null, null, null, null, null, null, null);
    }

    public String getIdentifier() {
      return identifier;
    }

    public BoundaryFamily getFamily() {
      return family;
    }

    public int getAnchorCount() {
      return anchorCount;
    }

    /** @return strict fixed-pressure seed points generated from a single discovered anchor */
    public int getGeneratedSeedPointCount() {
      return generatedSeedPointCount;
    }

    /** @return strict globally stable fixed-pressure points inserted between sparse discovery anchors */
    public int getInternalBridgePointCount() {
      return internalBridgePointCount;
    }

    public List<EvidencePoint> getEvidencePoints() {
      return evidencePoints;
    }

    public boolean isIsolated() {
      return anchorCount + generatedSeedPointCount < 2;
    }

    public TwoToThreePhasePseudoArcLengthTracer.Result getBackwardTrace() {
      return backwardTrace;
    }

    public TwoToThreePhasePseudoArcLengthTracer.Result getForwardTrace() {
      return forwardTrace;
    }

    /** @return fixed-pressure backward fallback evidence, or {@code null} when pseudo-arclength started normally */
    public HydrocarbonWaterRegularPressureBoundaryTracer.Result getRegularBackwardTrace() {
      return regularBackwardTrace;
    }

    /** @return fixed-pressure forward fallback evidence, or {@code null} when pseudo-arclength started normally */
    public HydrocarbonWaterRegularPressureBoundaryTracer.Result getRegularForwardTrace() {
      return regularForwardTrace;
    }

    /** @return pseudo-arclength trace started only after backward pressure-parameter collapse */
    public TwoToThreePhasePseudoArcLengthTracer.Result getPostRegularBackwardTrace() {
      return postRegularBackwardTrace;
    }

    /** @return pseudo-arclength trace started only after forward pressure-parameter collapse */
    public TwoToThreePhasePseudoArcLengthTracer.Result getPostRegularForwardTrace() {
      return postRegularForwardTrace;
    }

    /** @return classified backward endpoint, or {@code null} for an isolated branch */
    public HydrocarbonWaterBoundaryTerminationClassifier.Result getBackwardTermination() {
      return backwardTermination;
    }

    /** @return classified forward endpoint, or {@code null} for an isolated branch */
    public HydrocarbonWaterBoundaryTerminationClassifier.Result getForwardTermination() {
      return forwardTermination;
    }

    /** Builds a gate input only after the caller supplies independently classified physical endpoints. */
    public TwoToThreePhaseBoundaryQualityGate.Branch toQualityGateBranch(SpecialPointType startPointType,
        SpecialPointType endPointType) {
      if (isIsolated() || evidencePoints.size() < 2) {
        throw new IllegalStateException("an isolated anchor cannot be promoted to a boundary curve");
      }
      return new TwoToThreePhaseBoundaryQualityGate.Branch(identifier, family.getDefinition(),
          family.getRetainedPhaseZero(), family.getRetainedPhaseOne(), family.getIncipientPhase(), evidencePoints,
          startPointType, endPointType, backwardTermination, forwardTermination);
    }
  }

  /** Immutable branch assembly result preserving one entry per discovered branch. */
  public static final class Result {
    private final List<AssembledBranch> branches;
    private final List<EndpointCandidate> endpointCandidates;

    private Result(List<AssembledBranch> branches, List<EndpointCandidate> endpointCandidates) {
      this.branches = Collections.unmodifiableList(new ArrayList<AssembledBranch>(branches));
      this.endpointCandidates = Collections.unmodifiableList(new ArrayList<EndpointCandidate>(endpointCandidates));
    }

    public List<AssembledBranch> getBranches() {
      return branches;
    }

    /** @return corrected special-point candidates preserved outside ordinary continuation branches */
    public List<EndpointCandidate> getEndpointCandidates() {
      return endpointCandidates;
    }

    public int getEndpointCandidateCount() {
      return endpointCandidates.size();
    }

    public int getIsolatedBranchCount() {
      int isolated = 0;
      for (AssembledBranch branch : branches) {
        if (branch.isIsolated()) {
          isolated++;
        }
      }
      return isolated;
    }
  }
}
