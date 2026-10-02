---
name: neqsim-java8-rules
description: "Java compatibility boundaries for NeqSim. USE WHEN: writing or reviewing Java code or tests. Core src/ keeps Java 8 source compatibility; the independent MCP module uses Java 21. Covers forbidden core APIs, replacement patterns, API verification, and JavaDoc requirements."
last_verified: "2026-10-02"
---

# Java 8 Compatibility Rules for NeqSim

Core Java code under `src/` — including test classes in `src/test/java/` — **MUST** remain Java 8 source-compatible. The default core artifact targets Java 17; `pomJava8.xml` defines the separate Java 8 artifact. The independent `neqsim-mcp-server/` module requires Java 21 and uses its own POM and Spotless configuration. Do not apply the core API ban to that module.

Run formatting checks separately for the core and MCP module when both are involved. A root Spotless pass does not validate the MCP sources. See the authoritative boundaries in [AGENTS.md](../../../AGENTS.md#critical-constraint-java-compatibility-boundaries). The personal distribution workflow is documented in [shared thermodynamics repository management](../../../docs/development/shared-thermo-repository.md).

## Forbidden Java 9+ Features

| Forbidden | Java 8 Alternative |
|-----------|-------------------|
| `var x = ...` | Explicit type: `String x = ...`, `Map<String, Object> map = ...` |
| `List.of(a, b)` | `Arrays.asList(a, b)` or `Collections.singletonList(a)` |
| `Set.of(a, b)` | `new HashSet<>(Arrays.asList(a, b))` |
| `Map.of(k, v)` | `Collections.singletonMap(k, v)` or `new HashMap<>()` |
| `"str".repeat(n)` | `StringUtils.repeat("str", n)` (Apache Commons) |
| `str.isBlank()` | `str.trim().isEmpty()` |
| `str.strip()` | `str.trim()` |
| `str.lines()` | `str.split("\\R")` or BufferedReader |
| `Optional.isEmpty()` | `!optional.isPresent()` |
| Text blocks `"""..."""` | Regular strings with `\n` |
| Records | Regular class with fields, constructor, getters |
| Pattern matching `instanceof` | Traditional `instanceof` + cast |
| `Stream.toList()` | `.collect(Collectors.toList())` |

## Common `var` Replacements

```java
// WRONG (Java 10+):
var map = someMethod.toMap();
var list = getItems();
var result = calculate();

// CORRECT (Java 8):
Map<String, Object> map = someMethod.toMap();
List<String> list = getItems();
CalculationResult result = calculate();
```

## Required Import for String Repeat

```java
import org.apache.commons.lang3.StringUtils;
// Usage: StringUtils.repeat("=", 70)
```

## Code Formatting (Spotless) — MANDATORY

AI-generated Java is NOT auto-formatted. After creating or editing ANY `.java`
file, reformat it before committing — do not rely on local pre-commit hooks
being installed:

```bash
./mvnw spotless:apply    # reformats Java to the project style (Eclipse profile)
./mvnw spotless:check    # verifies formatting — this is what CI runs
```

- Formatter profile: `.config/neqsim_formatter.xml` (configured in `pom.xml`),
  applied to `src/main/java` and `src/test/java`.
- CI runs `./mvnw spotless:check` and FAILS the build on any unformatted file.
- Run `spotless:apply`, then `git add` the reformatted files, then commit.
- NEVER bypass the gate with `git commit --no-verify`.

## API Verification (MANDATORY)

Before using any NeqSim class in code or examples:

1. **Search** for the class: `file_search("**/ClassName.java")`
2. **Read** constructor and method signatures from the actual source
3. **Use only methods that actually exist** with correct parameter types
4. Do NOT assume convenience overloads — check first

Common API mistakes:
- Assuming `getXxx95()` when actual is `getXxx(int percentile)`
- Assuming 1-arg constructors when 2+ args are required
- Calling methods on wrong class
- Assuming `calculate()` when actual is `calculateRisk()` or `run()`

## JavaDoc Requirements

All classes and methods (public, protected, AND private) require complete JavaDoc:
- Class-level: description, `@author`, `@version`
- Method-level: description, `@param` for every parameter, `@return` for non-void, `@throws` for each exception
- HTML5 compatible: use `<caption>` in tables (no `summary` attribute)
- No `@see` with plain text — only valid Java references
- No lambda arrows (`->`) in JavaDoc code examples

## Build Commands

```bash
./mvnw install                            # full build
./mvnw test -Dtest=ClassName              # single test class
./mvnw test -Dtest=ClassName#methodName   # single method
./mvnw checkstyle:check spotbugs:check pmd:check  # static analysis
./mvnw javadoc:javadoc                    # verify JavaDoc
```

## Serialization — SE_BAD_FIELD Rule (MANDATORY)

SpotBugs enforces that all instance fields in `Serializable` classes are either
serializable themselves or marked `transient`. This applies to any class extending
`ProcessEquipmentBaseClass`, `MeasurementDeviceBaseClass`, `MechanicalDesign`,
thermo phase classes, or any other `Serializable` class.

### When to use `transient`

Mark a field `transient` when its type does NOT implement `Serializable`:
- Functional interfaces: `Function`, `BiConsumer`, `Consumer`, `Supplier`
- JDBC: `Connection`, `Statement`, `ResultSet`
- Threads: `Thread`, `ExecutorService`
- Apache Commons Math: `BicubicInterpolator`, `BicubicInterpolatingFunction`, `LinearInterpolator`
- Inner classes that don't implement `Serializable` (e.g., `NetworkNode`, `GibbsComponent`)
- External library types not designed for serialization

### Correct modifier order

```java
// private fields
private transient MyType field;
private final transient List<NonSerializableInner> items = new ArrayList<>();

// package-private fields
transient SomeType field;
```

### Verify with SpotBugs

```bash
./mvnw spotbugs:check 2>&1 | Select-String "SE_BAD_FIELD"  # should return empty
```
