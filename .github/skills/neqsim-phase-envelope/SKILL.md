---
name: neqsim-phase-envelope
description: "Generate, plot, interpret, validate, and troubleshoot NeqSim PT phase envelopes. USE WHEN: calculating phase envelopes, dew and bubble curves, cricondenbar, cricondentherm, critical points, retrograde regions, envelope segments, or fixing Michelsen continuation and singular-Jacobian failures caused by zero or trace components."
last_verified: "2026-09-12"
---

# NeqSim PT Phase Envelopes

## Overview

Use this skill to generate physically interpretable PT phase envelopes and to modify the
Michelsen continuation solver safely. Prefer the structured segment API, verify branch
identity from physics, and preserve the caller's thermodynamic system.

## Route the Task

- Use `@thermo-fluid` for envelope generation, plotting, properties, and solver defects.
- Chain to `neqsim-eos-regression` or `@pvt-simulation` when matching lab data or tuning EOS parameters.
- Chain to `neqsim-flow-assurance` for operating-path, hydrate, wax, or pipeline assessments.
- Chain to `neqsim-ccs-hydrogen` for CO2/H2 impurity envelopes.
- Chain to `neqsim-troubleshooting` when continuation does not converge or output is incomplete.

## Workflow

### 1. Verify the API and inputs

Read these sources before generating or changing code:

- `ThermodynamicOperations.calcPTphaseEnvelope(...)` and `getEnvelopeSegments()`
- `PTPhaseEnvelopeMichelsen` result getters and continuation settings
- `EnvelopeSegment` getters and `PhaseType`
- Neighboring tests in `PTPhaseEnvelopeMichelsenTest`, `PTPhaseEnvelopeRobustnessTest`, and
  `PTPhaseEnvelopeSegmentsTest`

Validate that temperature is in K, pressure is in bara, component amounts are finite and
non-negative, and the mixing rule is set. Choose the EOS from the fluid chemistry, not from
the desired plot shape.

### 2. Handle zero and trace components correctly

The Michelsen Jacobian contains terms divided by overall mole fraction $z_i$. Components with
zero or numerically negligible $z_i$ can therefore make the Jacobian singular.

- Current `PTPhaseEnvelopeMichelsen` excludes components with $z_i < 10^{-12}$ from a private
  clone after composition initialization.
- Do not remove components from the caller's `SystemInterface`.
- Do not discard physically meaningful trace impurities at or above the threshold merely to
  obtain convergence. Their impact can be important for CO2 quality, water behavior, and H2S.
- When changing the threshold or filtering logic, add a regression proving both finite envelope
  output and unchanged caller composition.

### 3. Generate the envelope

```java
SystemInterface fluid = new SystemSrkEos(273.15 + 25.0, 50.0);
fluid.addComponent("methane", 0.85);
fluid.addComponent("ethane", 0.10);
fluid.addComponent("propane", 0.05);
fluid.setMixingRule("classic");

ThermodynamicOperations ops = new ThermodynamicOperations(fluid);
ops.calcPTphaseEnvelope(true, 1.0);

List<EnvelopeSegment> segments = ops.getEnvelopeSegments();
double[] cricondenbar = ops.get("cricondenbar");
double[] cricondentherm = ops.get("cricondentherm");
```

Use `calcPTphaseEnvelope(true, lowPressureBara)` when the start side and lower pressure
matter. Use the no-argument overload for the standard calculation.

### 4. Consume structured segments by default

Use `getEnvelopeSegments()` for plotting, JSON export, and AI-generated analysis. Each segment
is contiguous and contains no `NaN` branch separators.

```java
for (EnvelopeSegment segment : ops.getEnvelopeSegments()) {
  double[] temperaturesK = segment.getTemperatures();
  double[] pressuresBara = segment.getPressures();
  EnvelopeSegment.PhaseType storedType = segment.getPhaseType();
  // Plot each segment independently; do not connect separate segments.
}
```

The legacy `ops.get("dewT")`, `ops.get("dewP")`, `ops.get("bubT")`, and
`ops.get("bubP")` arrays remain supported. They may contain intentional `NaN` values that mark
branch breaks. Never reject an envelope solely because these flat arrays contain `NaN`; reject
infinities and require finite physical points.

### 5. Classify physical branches

For `calcPTphaseEnvelope(true, 1.0)`, stored dew/bubble labels can be swapped by the historical
Michelsen tracing order. This applies to legacy array names and can affect `EnvelopeSegment`
stored `PhaseType` labels. Do not infer physical identity from a getter or enum name alone.

Classify the physical dew side as the branch containing the cricondentherm, normally the branch
with the highest finite temperature. The other side is the physical bubble branch. Preserve
segment boundaries while grouping segments by physical side.

```python
finite_a = branch_a_t[np.isfinite(branch_a_t)]
finite_b = branch_b_t[np.isfinite(branch_b_t)]
if finite_a.max() > finite_b.max():
    dew_t, dew_p = branch_a_t, branch_a_p
    bubble_t, bubble_p = branch_b_t, branch_b_p
else:
    dew_t, dew_p = branch_b_t, branch_b_p
    bubble_t, bubble_p = branch_a_t, branch_a_p
```

If the branches have overlapping maxima or an unusual topology, use the reported
cricondentherm point and critical-point continuity instead of relying only on maximum
temperature.

### 6. Validate the result

Require all applicable checks:

1. At least one non-empty structured segment.
2. At least one finite point on each expected physical branch.
3. No infinite temperatures or pressures.
4. Positive pressures and physically plausible temperatures.
5. Cricondentherm equals the maximum finite envelope temperature within numerical tolerance.
6. Cricondenbar equals the maximum finite envelope pressure within numerical tolerance.
7. Bubble and dew branches meet continuously near the critical point.
8. The original fluid still has the same component count and composition.
9. Benchmark key points against lab, literature, or another trusted EOS implementation when
   results drive an engineering decision.

Do not require every flat-array value to be finite because `NaN` is the branch-break sentinel.

## Solver Change Protocol

When modifying `PTPhaseEnvelopeMichelsen` or `SysNewtonRhapsonPhaseEnvelope`:

1. Add a focused regression that fails before the code change.
2. Assert physical outputs, not only `assertDoesNotThrow`.
3. Cover zero fraction, near-zero fraction, ordinary gas, and a heavy/TBP mixture as relevant.
4. Verify the operation may use a filtered clone while the caller's system is unchanged.
5. Preserve two-pass continuation, critical-point handling, and segment construction.
6. Run the focused regression immediately after the first edit.
7. Run the full Michelsen and robustness test classes.
8. Run Java 8 compilation, Spotless, Checkstyle, and JavaDoc gates.

Focused Windows commands:

```powershell
.\mvnw.cmd test "-Dtest=PTPhaseEnvelopeMichelsenTest#testMethod"
.\mvnw.cmd test "-Dtest=PTPhaseEnvelopeMichelsenTest,PTPhaseEnvelopeRobustnessTest,PTPhaseEnvelopeSegmentsTest"
.\mvnw.cmd --file pomJava8.xml test "-Dtest=PTPhaseEnvelopeMichelsenTest#testMethod" "-Djacoco.skip=true"
.\mvnw.cmd spotless:apply
.\mvnw.cmd spotless:check
```

Use explicit Java 8 types. Do not use `var`, `List.of`, `String.repeat`, records, or other Java
9+ APIs. Use Log4j2 rather than `System.out` or `System.err`.

## Troubleshooting

| Symptom | Likely cause | Action |
|---|---|---|
| Singular Jacobian or division by $z_i$ | Zero/negligible component reached continuation equations | Confirm the current private-clone filter runs after `init(0)`; never mutate the caller |
| Envelope arrays contain `NaN` | Intentional break between disjoint traced segments | Use `getEnvelopeSegments()` or split flat arrays on `NaN` |
| Bubble and dew labels look reversed | Historical stored-label behavior with bubble-first tracing | Classify from cricondentherm and physical topology |
| Empty or very short branch | Poor start point, extreme composition, unsuitable EOS, or continuation failure | Check logs, lower start pressure, simplify a clone for diagnosis, and compare with neighboring robustness tests |
| Unrealistic critical point | Wrong EOS, uncharacterized heavy end, bad composition, or unit error | Validate inputs, characterize plus fractions, and benchmark independently |
| Trace impurity disappears | Fraction is below the numerical filter threshold | Decide whether it is physically relevant; if so, use a defensible non-negligible composition and document sensitivity |
| Envelope reports a cricondentherm far colder than a direct flash finds liquid at | Michelsen continuation truncated the dew branch | Verify with a flash scan (below) before trusting any envelope on a lean gas with a small heavy tail |
| `dewPointTemperatureFlash` returns the initial temperature guess unchanged | Degenerate incipient-liquid seed | Fixed for zero-fraction water (see below). Otherwise reseed the flash near the expected root |
| Point dew point sits above the cricondentherm | Flash converged on the low-temperature retrograde root | Seed the flash at the cricondentherm temperature and assert `T_dew <= T_cricondentherm` |

### Always Cross-Check a Lean Gas With a Flash Scan

On a lean gas with a small heavy tail the Michelsen continuation can **truncate the dew
branch and still return a plausible-looking envelope**. Observed on a 90.5 mol% methane
gas-cap gas with 0.25 mol% C7+: `calcPTphaseEnvelope(true, 1.0)` reported a
cricondentherm of −66 °C, while a direct `TPflash` at 6 °C and 45 bara found 0.9 mol%
liquid — the true cricondentherm was +41 °C. Nothing in the envelope output flagged the
truncation, and the error would have hidden the whole liquid-handling and slugging issue
for a subsea tie-back.

Cheap guard: bisect on the phase count at one or two pressures and compare.

```python
def dew_point_temperature(fluid, pressure_bara, t_high=90.0, t_low=-60.0):
    """Upper (retrograde) dew temperature by bisection on the phase count."""
    lo, hi = t_low, t_high
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        w = fluid.clone()
        w.setTemperature(273.15 + mid)
        w.setPressure(pressure_bara)
        ThermodynamicOperations(w).TPflash()
        if int(w.getNumberOfPhases()) > 1:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
```

If the scan and the envelope disagree by more than a few kelvin, trust the scan and
report the envelope as unreliable for that fluid.

### Point Dew-Point Flashes on Wet Gas

`ThermodynamicOperations.dewPointTemperatureFlash()` seeds an aqueous incipient
liquid whenever water carries moles, so on a wet gas it returns the **water** dew
point. For a hydrocarbon dew point, clone the fluid and
`removeComponent("water")` first.

Seed the flash at the cricondentherm temperature so it descends onto the upper
(physical) dew branch rather than a low-temperature retrograde root:

```java
double[] cct = envOps.get("cricondentherm"); // [T (K), P (bara)]
hcFluid.setPressure(pBara, "bara");
hcFluid.setTemperature(cct[0]);
hcFluid.init(0);
new ThermodynamicOperations(hcFluid).dewPointTemperatureFlash();
```

Zero-fraction water no longer changes the result: the aqueous seed is gated on
`ConstantDutyTemperatureFlash.hasSignificantWater`, so a component with
$z_i = 0$ behaves like an absent one.

## Output Convention

Report:

- EOS and mixing rule
- normalized composition and any filtered numerical-zero components
- critical point, cricondenbar, and cricondentherm with K/bara units
- physically classified dew and bubble segments
- operating points or paths overlaid on the envelope when relevant
- convergence or filtering warnings
- benchmark source and deviations for decision-critical work

## Silent Truncation of the Cricondenbar (READ THIS)

`calcPTphaseEnvelope` marches along the saturation curve and reports the highest pressure it
visits. For mixtures carrying several trace heavy components the march can **terminate early**,
and the endpoint of the partial trace is then reported as the cricondenbar — with no exception, no
warning, and a plausible-looking number.

Observed case: a rich natural gas whose true cricondenbar is 123.4 bara at −6 °C was reported at
62.0 bara at −65.8 °C — **49.8 % low**. Used for a dense-phase pipeline design that would have set
the minimum operating pressure at 77 bara instead of 138 bara, and the line would have run
two-phase while believed single-phase.

**Always cross-check any cricondenbar that feeds a design decision.** Use
`neqsim.thermodynamicoperations.phaseenvelopeops.multicomponentenvelopeops.RobustPhaseEnvelope`,
which finds the boundary by a flash grid plus bisection and is independent of continuation:

```java
RobustPhaseEnvelope env = new RobustPhaseEnvelope(fluid);
env.setTemperatureRange(150.0, 350.0);
env.setPressureRange(1.0, 300.0);
env.calculate();
double pcb = env.getCricondenbarPressure();
boolean suspect = env.continuationLooksTruncated(continuationValue, 0.10);
```

Symptoms that a continuation result is truncated:
- The cricondenbar temperature sits at an implausible extreme (well below −50 °C for a
  hydrocarbon gas).
- The cricondentherm is far colder than the heaviest component's boiling point suggests.
- Adding or removing a trace component changes the cricondenbar by more than a few bar.

## Explicit specified-phase research APIs

- `SpecifiedMultiphaseFlashSolver` and `RetainedPhaseBifurcationCorrector` separate numerical roots from physical
  acceptance. `isConverged()` is insufficient: call `validateEquilibrium` or `toValidatedThermodynamicSystem`
  before consuming a finite-phase equilibrium.
- The shared validator replays balance and fugacity, rejects coincident/vanishing phases and identity mismatches,
  and requires completed gas/oil/water stability searches. Finite searches do not prove global stability or
  laboratory accuracy. Critical and incipient-phase roots require boundary-specific checks.
- `HydrocarbonWaterBoundaryModelProfile` is clone-only and opt-in. Per-boundary water-Kij multipliers remain
  research settings, never a default production calibration inferred from software benchmark agreement.
- Preserve archive originals and migrate assertions, including rejected high-pressure roots. See
  `docs/development/shared-thermo-research-integration.md` and the shared-engine regression policy.
- Liquid slots must use `LiquidPhaseClassification`, matching this branch's mass-based EOS convention;
  never infer aqueous identity solely from a water mole-fraction threshold. Read gas identity from the evaluated
  EOS root, including water vapor. `SystemThermo.init(1, phase)` relabels gas in nonzero slots as oil;
  stability trials in slot one therefore recover the `PhaseEos` volume/covolume criterion before classifying
  a liquid. Phase labels are separate from equilibrium acceptance.
- Upstream merges require official flash regressions as well as personal envelope regressions, including files
  changed only upstream. Retire superseded recovery patches when differential tests show interference. Keep
  historical light-pseudo/water parameters explicit through the model profile; CPA and heavy-fraction data are
  separate models. See `docs/development/shared-thermo-compatibility-audit.md`.

## Independent PVTsim comparisons

- Run the installed licensed Open Structure engine afresh. Record the actual DLL file version/hash and the
  selected NeqSim release JAR hash; an adapter's display name is not version evidence. Keep task fluids and
  raw point results in the task folder, outside code commits.
- Audit input units independently of field names. In the installed FluidLw SDK, component critical pressure
  uses atmospheres, while NeqSim uses bara (1 atm = 1.01325 bar); flash-state pressure uses Pa. Inspect the
  loaded component's unit metadata before conversion. Legacy water-component exports can mix critical-volume
  units. Do not alter frozen delivery resources to make a benchmark agree.
- Preserve the requested polar model and explicit interactions. Verify the half-table orientation and complete
  pair count after import; an ignored interaction must not count as a successful same-parameter comparison.
  SDK input-range rejection is missing comparison coverage, not numerical agreement. In particular, legacy
  N2/water tables below -1 cannot pass the installed SDK's greater-than-minus-one BIP requirement unchanged.
- Keep same-EOS solver verification separate from the accuracy of a model fitted to a different native PVTsim
  reference. Transferring fitted interactions to both solvers does not independently validate that calibration.
- Refine the union of both solvers' phase-transition brackets and keep all topology islands. Count unique
  states separately from repeated coarse/refined evaluations. Gate each point on phase-fraction sum,
  component balance and fugacity, not only returned phase labels or lack of an exception.
- Reject invalid final multiphase statuses. Exclude failed points from accepted boundary estimates and report
  any resulting wider brackets. A high coarse-grid topology agreement does not close boundary failures.
- Promote newly reproduced conservation failures into mandatory merge/build regressions immediately; a release
  gate is not a reason to defer a known solver defect. Check cold starts and pressure round trips. Internal phase
  labels can be provisional before density ordering, so verify compositions and balances rather than assuming
  that absence of an AQUEOUS or GAS label excludes a water-bearing gas/oil state.
- Replay failed states on the old personal release and an isolated official release before attributing them
  to a new upstream merge. Compare numerical residuals as well as pass/fail counts; a retained personal patch
  can leave much larger conservation errors than the official result. See
  `docs/development/shared-thermo-compatibility-audit.md` for the bounded rc.3 findings.

## Known Limitations

- Stored branch labels can differ from physical branch identity for bubble-first tracing.
- Cubic-EOS envelope accuracy depends on heavy-end characterization and binary interaction
  parameters; numerical convergence does not prove physical accuracy.
- PT phase envelopes do not replace hydrate, wax, solid, electrolyte, or reactive-equilibrium
  boundaries. Add those analyses through their dedicated skills.
