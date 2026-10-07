---
title: Shared explicit thermodynamic engine
description: Standalone Java and JSON API for the explicit-component flash, enthalpy and transport methods migrated from fluid-platform.
---

# Shared explicit thermodynamic engine

The maintained `cupbhan/neqsim` source owns the calculation algorithms previously embedded in
`fluid-platform/engine`. Consumers pin a tested shared distribution; the platform supplies its UI,
parameter mapping, task management, provenance and JSON transport. No PVTsim installation, SDK,
platform checkout or sealed reference data is needed to calculate with this API.

## Entry points

Package: `neqsim.thermo.util.explicit`.

| Entry point | Purpose |
| --- | --- |
| `ThermodynamicEngine.compositionCharacterize(String)` / `compositionCharacterize(JsonNode)` | Independent raw SCN/GC characterization, effective mole-basis input and bulk-closure diagnostics |
| `ThermodynamicEngine.compositionFlash(String)` / `compositionFlash(JsonNode)` | Independent PR/SRK database fluid from composition and optional user TBP characterization |
| `ThermodynamicEngine.compositionEnvelope(String)` / `compositionEnvelope(JsonNode)` | Physical envelope of the independently built database fluid, with explicit projection audit |
| `ThermodynamicEngine.calculate(String)` | Explicit-component grid flash with validated enthalpy, isobaric heat capacity and viscosity; JSON in and JSON out |
| `ThermodynamicEngine.envelope(String)` / `envelope(JsonNode)` | Real continuation envelope with physical branch labels, projection audit and incomplete-coverage status |
| `ThermodynamicEngine.calculate(JsonNode, boolean)` | Same grid; `false` selects flash-only evaluation |
| `ThermodynamicEngine.diagnose(JsonNode)` | Fixed-composition phase diagnostics, without another equilibrium calculation |
| `ThermodynamicEngine.standardFlash(String)` | Database-component PR/SRK flash, constant BIPs and PFCT properties |
| `PureWaterViscosity.densityKgM3(double, double)` | IF97 region-1 density from temperature in K and pressure in MPa |
| `PureWaterViscosity.viscosityPaS(double, double)` | Industrial IAPWS viscosity from temperature in K and density in kg/m3 |
| `CspViscosity.calculate(JsonNode, JsonNode, JsonNode)` | Explicit CSP mixture viscosity at a specified phase composition |

The JSON string entry points can also be called through a JVM bridge such as JPype. This does not
add a new MCP tool or replace NeqSim's default system getters: the explicit compatibility profile is
selected deliberately. In particular, the enthalpy convention and the 95 K CSP transition hypothesis
must not silently become the default for unrelated calculations.

The complete public synthetic request is
[`example-request.json`](../../src/test/resources/neqsim/thermo/explicit/example-request.json).
It is executed by `ThermodynamicEngineTest.documentationExampleRunsWithoutPlatform`:

```java
String response = neqsim.thermo.util.explicit.ThermodynamicEngine.calculate(request);
```

Use the core JAR and dependency POM from the shared distribution, or its standalone runner JAR as a
classpath dependency. The core distribution requires Java 17; the standalone MCP runner requires
Java 21. Source code remains compatible with the repository's Java 8 syntax policy.

## Explicit request contract

`components` gives ordered, normalized mole fractions `z`, names and explicit EOS data:
`tcK`, `pcBar`, `mwKgMol`, `tbK`, `acentricFactor`, `omegaA`, `omegaB`, `alpha`,
`penelouxLmol` and `penelouxTLmolK`. Thermal evaluation additionally requires `vcCm3Mol`,
four finite `cpJMolK` polynomial coefficients and `mwWeightKgMol`. The heat-capacity polynomial
uses absolute temperature in kelvin. Component database templates allocate component objects;
the supplied parameters replace the thermodynamic values used by the explicit calculation.

`eos` accepts `SrkPeneloux`, `SrkPenelouxTemperatureDependent`, `PengRobinsonPeneloux` or
`PengRobinson78Peneloux`. `mixingRule` is `Classic` (default) or `HV`. `interactions` contains
zero-based component indices `i`, `j`, `kij` and `kijT`; the latter uses a 288.15 K reference.
Huron-Vidal pairs also supply `mixingRule: "HuronVidal"`, dimensionless `alpha`, `gij`, `gji`,
`gijT`, `gjiT`. The `gij`/`gji` inputs are energy divided by the gas constant, in K; their
temperature coefficients `gijT`/`gjiT` are dimensionless. Unlisted pairs are zero.

`grid` is a nonempty array of `{temperatureC, pressureMPa}`; pressure is absolute.
Set `nativeProperties.enabled` to `true`, `nativeProperties.profile` to `native-properties-v1`
and `propertyAuditVersion` to `4` to select the migrated thermal profile. These historical field
names are retained for request compatibility. Supply `transport.normalFactors` and `heavyFactors`
as two positive finite numbers each, `nativeCandidate` as
`native-csp-mccarty1974-sdktransition95-v3`, and `waterCandidate` as
`iapws2008-if97-region1-pure-water-v1`. The example contains all required fields.

Results contain `points`; each point has `status`, phase compositions and equilibrium diagnostics.
Phase values are molar fraction, density in kg/m3, molar mass in kg/mol, dimensionless Z,
enthalpy in J/mol, `cpJMolK` in J/(mol K) and dynamic viscosity in Pa s. Failed states remain separate from successful
neighbors. Unsupported properties are `null` with explicit availability reasons, never zero.
Requests are not modified; each call constructs its own thermodynamic state.

## Explicit envelope contract

The envelope entry point consumes the same explicit components and complete unordered interaction table.
It never reads a simulator reference or a sealed curve. Optional `envelope` settings are
`temperatureMinC` (default -120), `temperatureMaxC` (250), `pressureMinMPa` (0.01),
`pressureMaxMPa` (40), `maximumIterations` (1200; range 100–5000) and
`hydrocarbonOnly` (false). Temperature bounds set the initial tracing temperature; the pressure bounds
are continuation limits. The actual envelope can extend beyond the requested temperature interval.
Continuation chooses its own steps; a regular-grid `pointCount` is not an accepted setting.

Without explicit projection, positive water, inhibitor or salt fractions produce an unavailable result.
With `hydrocarbonOnly: true`, these components are removed from a private model and the remaining
composition is normalized. `sourceFeed`, `feed` and `omittedComponents` preserve source indices,
names, original fractions and omission reasons. Numerical-zero fractions below $10^{-12}$ are also
reported; meaningful trace components are retained. A PVTsim comparison must apply the **same**
projection to its independently calculated reference, rather than comparing dry and wet envelopes.

Results contain contiguous `segments` with physical `dew`/`bubble` types and arrays of
`{temperatureC, pressureMPa}` points. Characteristic points, termination reasons, iteration limits,
engine/model provenance and unavailable water/three-phase boundaries are explicit. `status: partial`
retains calculated points while identifying incomplete continuation. A complete result requires both
physical sides, a critical point, normal termination, and independent flashes above the reported
cricondentherm without a phase split. This check can detect truncation; it is not a proof of global
stability. No schematic or interpolated substitute curve is returned when calculation fails.

A projected HV model whose remaining pairs are all Classic uses the exactly equivalent Classic
submodel. Classic linear temperature-dependent interactions are bound to the continuation evaluator.
HV interaction terms with nonzero Classic `kijT` remain explicitly unavailable for envelopes; they
are never silently frozen at the initial temperature. Single-component vapor-pressure curves and
aqueous/three-phase boundaries are outside this entry point's current supported coverage.

## Algorithms and qualification limits

- `ExplicitFlash`: PR/SRK, explicit OmegaA/B, Classic/HV interactions, temperature-dependent BIPs,
  Peneloux translation, material/fugacity balance checks and the qualified NaCl path. NaCl feed and
  reported phases use formula-unit moles; its internal EOS uses ionic moles. Supersaturation remains
  a failed state, not a precipitated-solid prediction.
- `ExplicitProperties` and `ExplicitPropertyCalculator`: explicit ideal heat capacities,
  temperature-dependent EOS residual derivatives and volume-translation enthalpy correction.
  `native-enthalpy-sdk-convention-v1` retains the 273.15 K reference and compatibility gas constant;
  derivative step agreement and density/Cp closure remain mandatory.
- `ExplicitHeatCapacity`: the derivative of the qualified phase enthalpy with temperature at fixed
  composition and pressure. Independent 0.05 K and 0.025 K steps must agree; root changes,
  negative/nonfinite Cp and failed enthalpy checks produce unavailable metadata. One-sided
  second-order differences are used at the thermal domain edges. This is phase Cp, not the
  derivative of total equilibrium enthalpy across changing phase fractions.
- `CspViscosity`: McCarty (1974), *Cryogenics* 14, 276–280, table 2; Hanley et al. (1977),
  *JPCRD* 6, 597–610, tables 1, 2 and 4. Explicit mixture/heavy-oil equations follow the recorded
  PVTsim Nova 2026 method-manual section 11.1. The 95 K transition is a separately identified SDK
  parity hypothesis, not a physical freezing point or a value claimed to be specified by that manual.
- `PureWaterViscosity`: IAPWS R7-97(2012), equations 7 and 30, and R12-08, equations 11 and 12.
  Critical enhancement is omitted. Aqueous mixture use requires at least 99 mol% water and excludes
  appreciable salt or inhibitor; output explicitly identifies a pure-water approximation.

The migrated thermal profile remains limited to 25–150 °C and 1–20 MPa. Salt enthalpy/viscosity and
inhibitor-containing aqueous viscosity remain unavailable. Moving code does not expand the validated
domain or constitute a new full simulator comparison.

## Validation and release

`ThermodynamicEngineTest` runs published water/methane verification values, original synthetic
PR/SRK × Classic/HV migration results, external-consumer JSON calls, repeatability and failure
isolation. The synthetic snapshots are regression evidence, not independent accuracy evidence.
Private platform datasets remain outside this repository and are used separately for consumer parity.

`ExplicitRealtimePropertiesTest` compares the new Cp derivative to independent analytic constant-kij
EOS derivatives, covers HV and temperature-dependent interactions and thermal domain edges, and
checks real envelope branches, extrema, unchanged inputs, explicit dry projection and missing coverage.

The enhancement policy requires this suite for every shared candidate. Follow
[shared runtime management](shared-thermo-repository.md) to build a clean immutable candidate,
validate consumers and activate it. Do not overwrite an existing release or silently retarget a
consumer's fixed version. The original platform `v1.0.0` tag and sealed comparisons remain historical
baselines; a new consumer run records the new JAR and thin-adapter hashes.

## Independent composition modeling

New application calculations and realtime cross-engine comparisons use `compositionFlash` and
`compositionEnvelope`. A sample library supplies component identities and mole fractions only.
NeqSim constructs a new PR or SRK fluid with its own component database, ordinary property methods,
Classic mixing rule (2) and database binary interactions. A reference simulator independently constructs
its own fluid from the same composition. Their parameter sets can differ; numerical differences are
model differences as well as solver differences. The older explicit-parameter entry points remain
available for deliberately requested parameter-controlled research and historical replay.

Example composition-only request:

```json
{
  "eos": "PR",
  "volumeCorrection": false,
  "temperatureC": 30,
  "pressureMPa": 5,
  "components": [
    {"id": "gas", "name": "C1", "z": 0.7},
    {"id": "middle", "name": "nC4", "z": 0.2},
    {"id": "heavy", "name": "C10", "z": 0.1}
  ]
}
```

Allowed root fields are exactly `eos`, `volumeCorrection`, `temperatureC`, `pressureMPa`,
`components`, `interactions` and `envelope`. Components accept exactly `id`, `name`, optional
`formula`, optional `kind`, `z` and optional `tbp`. Aliases include N2, CO2, H2O, MeOH, C1–C3, iC4/nC4,
iC5/nC5 and normal C6–C20; other unambiguous NeqSim catalog names are also accepted.
No source critical constants, molar masses, heat capacities, viscosity parameters or mixing-rule
profiles are accepted. Contradictory identities and duplicate database bindings fail explicitly.

Unknown fractions such as C10-C13 or C7+ require **user-supplied**
`tbp: {"molarMassKgMol": 0.17, "densityRelative": 0.79}`. NeqSim's TBP characterization then
generates its own parameters. These two characterization inputs must not be obtained implicitly
from the reference model. Supplying TBP parameters for a recognized pure database component fails.
Ions require a separately qualified electrolyte model and are rejected by this cubic composition API.
`kind` accepts `pure`, `pseudo`, `plus` or `custom`. Declared `pseudo`/`plus` fractions always require
user TBP data for a positive mole fraction, even when their name is a pure-compound alias such as C7.
Exactly zero declared fractions are omitted without inventing physical parameters; their original
indices remain in `omittedComponents` and their zero entries remain in each phase composition.
Omitted rows have no entry in the actual modeled `components` bindings. Missing kind, `pure` and
`custom` retain ordinary catalog lookup semantics.
Optional `interactions` entries are explicit local constant overrides `{i, j, value}` with zero-based
input indices. All other pairs retain NeqSim database defaults; the complete applied table and its
source are returned. Returned pair `i`/`j` index the actual modeled bindings, while `sourceI`/`sourceJ`
index the caller's original rows. Requests and component order are preserved.

Flash output has `ok`, `phases`, component `bindings` in `components`, `diagnostics` and
`interactions`. Phase values use density kg/m³, viscosity cP, enthalpy kJ/kg, heat capacity
J/(kg K), molar mass kg/mol and dimensionless compressibility factor. Formation enthalpy is
disabled, retaining NeqSim's sensible-enthalpy reference at 273.15 K. `propertyAvailability`
describes each independently evaluated field; missing fields are null, while remaining valid fields
survive. `unavailableProperties` retains compatibility entries for density, viscosity and enthalpy.
Property method classes, independent component provenance, EOS, volume correction, material balance
and fugacity diagnostics accompany results. This is normal NeqSim property evaluation; the old
explicit thermal profile's 25–150 °C qualification does not establish accuracy bounds for this API.

Composition envelopes use the same envelope settings and physical segment representation described
above. Positive water and common inhibitors require `hydrocarbonOnly: true`, which removes them
from a private normalized copy and records all original indices, including zero-feed entries.
Pure-fluid saturation, water and three-phase boundaries remain explicitly unavailable. Any reference
scan must apply the returned `feed` projection exactly before independent fluid construction.

`CompositionThermodynamicEngineTest` is required by the shared release gate. It compares PR/SRK
results to separately assembled NeqSim database systems, audits aliases and zero-row indices,
requires user TBP inputs, rejects imported parameters, preserves default BIPs, checks aqueous
property availability, and traces physical dry envelope branches without modifying caller input.

### rc.7 consumer acceptance

`3.23.0-cupbhan.1-rc.7` was built from clean commit
`d2869c34a0269b40ff50cdd9f20374e0e01a2fc8` and selected in the shared runtime lock after
consumer validation. The release gate checked 485 core tests (484 passed, one pre-existing upstream
disabled pure-component envelope test), four MCP tests and a real STDIO runtime smoke test.
Spotless, Checkstyle, Javadoc and isolated Java 8 source compilation also passed.

The platform replayed 11 sealed production records / 99 states using the released JAR. All 19,836
numeric comparisons across 199 phases were identical to the original implementation: 98 successful
states and one expected NaCl supersaturation failure. This is migration parity, not a new SDK run or
an expansion of the qualification domain. Private models and phase data remain in the consumer's
ignored local-data directory. The consumer report SHA-256 is
`ec8a47fcd153fcf518248ac157260307f746904eb303a642785c4fb5f8d1ad0b`.

## Independent raw SCN and GC characterization

The composition entry points also accept a declared raw heavy-end assay. No source database,
PVTsim installation or generated PVTsim fraction properties are required. `compositionCharacterize`
returns `effectiveInput` and a `characterization` audit; `compositionFlash` and
`compositionEnvelope` independently resolve the same raw input before calculating and return the audit.

| Field | Meaning |
| --- | --- |
| `eos` | `PR` or `SRK`; pure-component and EOS data come from NeqSim |
| `basis` | `mole` or `mass`; ordered component `z` values sum to one on this raw basis |
| `components` | Identities, raw fractions, `pure`/`pseudo` kinds and identifiable carbon numbers; up to 100 components |
| `characterization.method` | `pedersen-scn-bulk-closure` |
| `characterization.startCarbonNumber` / `terminalCarbonNumber` | Declared SCN heavy region, C7 through C79; the terminal entry is aggregate Cn+ |
| `characterization.molarMassKgMol` | Measured bulk heavy-end molar mass, kg/mol |
| `characterization.densityRelative` | Measured bulk standard liquid density, numerically g/cm3 |
| `characterization.includeDefinedHeavy` | Include NeqSim's independently resolved defined heavy components in GC bulk constraints |
| `characterization.molarMassClosure` | Optional `scale-unmeasured` for the declared GC case; otherwise close with the terminal aggregate |
| `components[].tbp` | Measured pre-plus fraction `molarMassKgMol` and `densityRelative`; these inputs are retained |

NeqSim generates missing SCN molar masses and densities from its own correlations, applies the bulk
constraints and obtains thermodynamic parameters from its TBP method. Raw fractions are preserved.
Mass-basis amounts are converted with NeqSim's own resolved molar masses. Foreign critical properties,
interaction coefficients, transport parameters and undeclared fields are rejected.

The terminal Cn+ is currently one aggregate TBP fraction; there is no terminal-distribution
extrapolation. A reference engine may independently extrapolate and lump, so equal raw inputs do not
imply equal characterized component lists or mole fractions. Each engine must calculate its envelope
with its own characterized feed. The audit records the terminal convention, generated cuts, bulk
molar-mass/density closure residuals and `sourceParametersUsed: false`.

Synthetic standalone verification is in `CompositionCharacterizationTest`: both input bases,
raw and resolved PR/SRK flashes, retained measured cuts, defined GC components, bulk closure,
physical lower bounds, rejected foreign parameters and input immutability. These five tests are
required by the shared-distribution regression policy.
