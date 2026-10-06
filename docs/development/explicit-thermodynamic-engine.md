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
| `ThermodynamicEngine.calculate(String)` | Explicit-component grid flash with validated enthalpy and viscosity; JSON in and JSON out |
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
enthalpy in J/mol and dynamic viscosity in Pa s. Failed states remain separate from successful
neighbors. Unsupported properties are `null` with explicit availability reasons, never zero.
Requests are not modified; each call constructs its own thermodynamic state.

## Algorithms and qualification limits

- `ExplicitFlash`: PR/SRK, explicit OmegaA/B, Classic/HV interactions, temperature-dependent BIPs,
  Peneloux translation, material/fugacity balance checks and the qualified NaCl path. NaCl feed and
  reported phases use formula-unit moles; its internal EOS uses ionic moles. Supersaturation remains
  a failed state, not a precipitated-solid prediction.
- `ExplicitProperties` and `ExplicitPropertyCalculator`: explicit ideal heat capacities,
  temperature-dependent EOS residual derivatives and volume-translation enthalpy correction.
  `native-enthalpy-sdk-convention-v1` retains the 273.15 K reference and compatibility gas constant;
  derivative step agreement and density/Cp closure remain mandatory.
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

The enhancement policy requires this suite for every shared candidate. Follow
[shared runtime management](shared-thermo-repository.md) to build a clean immutable candidate,
validate consumers and activate it. Do not overwrite an existing release or silently retarget a
consumer's fixed version. The original platform `v1.0.0` tag and sealed comparisons remain historical
baselines; a new consumer run records the new JAR and thin-adapter hashes.

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
