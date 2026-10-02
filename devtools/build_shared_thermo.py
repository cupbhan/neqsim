"""Build a traceable personal NeqSim candidate from a clean commit.

Uses the existing Maven wrapper and caches, but a distinct version for both modules.
No upload, tag creation, consumer replacement, or public Maven deployment occurs.
"""

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from shared_thermo_smoke import run_smoke

ROOT = Path(__file__).resolve().parents[1]
CORE_TESTS = ("ComponentQueryTest,FieldFluidRunnerTest,WaterIF97RunnerTest,"
              "HeavyOilMultimediaFluidTest,TPmultiflashSolveStatusTest,HydrocarbonWater*Test")
VERSION_PATTERN = re.compile(r"\d+\.\d+\.\d+-cupbhan\.\d+(?:-rc\.\d+)?\Z")


def digest(path):
    """Hash an artifact without reading the complete file into memory."""
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def git(*args, root=ROOT):
    return subprocess.check_output(["git", *args], cwd=root, text=True, encoding="utf-8").strip()


def require_clean(root=ROOT):
    status = git("status", "--porcelain", "--untracked-files=normal", root=root)
    if status:
        raise RuntimeError("Candidate builds require a clean committed tree:\n" + status)
    return git("rev-parse", "HEAD", root=root)


def test_summary(directory):
    """Reject missing, skipped, failing or empty selected test reports."""
    files = sorted(Path(directory).glob("TEST-*.xml"))
    summary = {"suites": [], "tests": 0, "failures": 0, "errors": 0, "skipped": 0}
    for path in files:
        suite = ET.parse(path).getroot()
        summary["suites"].append(suite.attrib["name"])
        for key in ("tests", "failures", "errors", "skipped"):
            summary[key] += int(suite.attrib.get(key, 0))
    if not summary["tests"] or any(summary[k] for k in ("failures", "errors", "skipped")):
        raise RuntimeError("Selected regression gate did not fully pass: " + json.dumps(summary))
    return summary


def verify_jar_version(path, artifact, expected):
    entry = f"META-INF/maven/com.equinor.neqsim/{artifact}/pom.properties"
    with zipfile.ZipFile(path) as jar:
        properties = dict(line.split("=", 1) for line in jar.read(entry).decode().splitlines()
                          if "=" in line and not line.startswith("#"))
    if properties.get("version") != expected:
        raise RuntimeError(f"{artifact} version mismatch: {properties.get('version')} != {expected}")


def build(version):
    if not VERSION_PATTERN.fullmatch(version):
        raise ValueError("Use a personal version such as 3.17.0-cupbhan.1-rc.1")
    commit = require_clean()
    base = ROOT / "build/shared-thermo"
    destination = base / "releases" / version
    if destination.exists():
        raise FileExistsError("Refusing to replace an existing distribution: " + str(destination))
    logs = base / "validation" / version
    logs.mkdir(parents=True, exist_ok=True)
    wrapper = str(ROOT / ("mvnw.cmd" if os.name == "nt" else "mvnw"))
    java_home = os.environ.get("JAVA_HOME")
    java = str(Path(java_home) / "bin" / ("java.exe" if os.name == "nt" else "java")) if java_home else "java"
    common = [wrapper, "-B", "-ntp", f"-Drevision={version}", f"-Dmcp.revision={version}",
              "-Dmaven.javadoc.skip=true", "-Djacoco.skip=true"]
    commands = []

    def run(label, command):
        print(label + " ...", flush=True)
        with (logs / (label + ".log")).open("w", encoding="utf-8") as stream:
            result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                                    creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        commands.append({"step": label, "command": command, "exitCode": result.returncode})
        if result.returncode:
            raise RuntimeError(f"{label} failed; inspect {logs / (label + '.log')}")

    run("core", common + ["clean", "spotless:check", "install", "-Dneqsim.shade.skip=true",
                          "-Dtest=" + CORE_TESTS])
    core_tests = test_summary(ROOT / "target/surefire-reports")
    run("mcp", common + ["-f", "neqsim-mcp-server/pom.xml", "clean", "spotless:check", "package"])
    mcp_tests = test_summary(ROOT / "neqsim-mcp-server/target/surefire-reports")
    core = ROOT / "target" / f"neqsim-{version}.jar"
    runner = ROOT / "neqsim-mcp-server/target" / f"neqsim-mcp-server-{version}-runner.jar"
    verify_jar_version(core, "neqsim", version)
    verify_jar_version(runner, "neqsim-mcp-server", version)
    verify_jar_version(runner, "neqsim", version)
    print("stdio-runtime ...", flush=True)
    smoke = run_smoke(runner, version, logs, ROOT, java)
    if require_clean() != commit:
        raise RuntimeError("Source commit changed during the build")
    baseline = json.loads((ROOT / "distribution/cupbhan/source-baseline.json").read_text(encoding="utf-8"))
    report = {"schemaVersion": 1, "version": version, "sourceCommit": commit,
              "scope": "Selected shared-engine regressions and real STDIO runtime; not the full repository suite",
              "core": core_tests, "mcp": mcp_tests, "runtime": smoke, "commands": commands,
              "pvtsimRerun": False, "productRuntimeSwitched": False}
    (logs / "validation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    destination.mkdir(parents=True)
    for source in (core, runner, ROOT / "LICENSE", ROOT / "distribution/cupbhan/source-baseline.json",
                   logs / "validation.json"):
        shutil.copy2(source, destination / source.name)
    artifacts = [{"file": p.name, "sizeBytes": p.stat().st_size, "sha256": digest(p)}
                 for p in sorted(destination.iterdir()) if p.is_file()]
    manifest = {"schemaVersion": 1, "version": version, "status": "tested-candidate",
                "repository": baseline["canonicalRepository"], "commit": commit, "dirty": False,
                "branch": git("branch", "--show-current"),
                "upstreamCommit": baseline["upstreamBaseline"],
                "builtAt": datetime.now(timezone.utc).isoformat(), "artifacts": artifacts,
                "java": subprocess.check_output([java, "-version"], stderr=subprocess.STDOUT, text=True).strip(),
                "validation": "validation.json", "publicMavenPublished": False,
                "knownLimitations": ["See source-baseline.json for model qualification limits",
                                     "Dedicated heavy-oil multimedia MCP endpoints are not yet exposed",
                                     "Full PVTsim and consumer integration gates remain separate"]}
    (destination / "SOURCE.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    checksum_files = sorted(p for p in destination.iterdir() if p.is_file())
    (destination / "SHA256SUMS").write_text(
        "".join(digest(p) + "  " + p.name + "\n" for p in checksum_files), encoding="utf-8")
    print(str(destination), flush=True)
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()
    try:
        build(args.version)
    except (RuntimeError, ValueError, OSError, subprocess.CalledProcessError) as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
