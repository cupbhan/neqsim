"""Prevent upstream synchronization from dropping personal contracts or model data."""

import hashlib
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]


def verify(root=ROOT):
    root = Path(root)
    baseline = json.loads((root / "distribution/cupbhan/upstream.json").read_text(encoding="utf-8"))
    policy = json.loads((root / "distribution/cupbhan/personal-enhancements.json").read_text(encoding="utf-8"))
    revision = ET.parse(root / "pom.xml").getroot().findtext("{*}properties/{*}revision")
    if revision != baseline["version"]:
        raise RuntimeError("Core POM does not match the recorded official version")
    subprocess.run(["git", "merge-base", "--is-ancestor", baseline["commit"], "HEAD"],
                   cwd=root, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for name in policy["requiredFiles"]:
        if not (root / name).is_file():
            raise RuntimeError("Missing personal enhancement: " + name)
    for name, expected in policy["frozenResources"].items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != expected:
            raise RuntimeError("Frozen personal model data changed: " + name)
    source = (root / "neqsim-mcp-server/src/main/java/neqsim/mcp/server/NeqSimTools.java").read_text(encoding="utf-8")
    methods = set(re.findall(r"public String (\w+)\s*\(", source))
    missing = set(policy["requiredTools"]) - methods
    if missing:
        raise RuntimeError("Missing personal MCP contracts: " + ", ".join(sorted(missing)))
    return baseline, policy


if __name__ == "__main__":
    baseline, policy = verify()
    print(f"Verified {baseline['tag']}: {len(policy['requiredTools'])} required tools, "
          f"{len(policy['frozenResources'])} frozen model resources")
