"""Resolve the pinned shared engine without importing or starting Java."""

import hashlib
import json
import os
from pathlib import Path


def resolve_runtime():
    configured = os.environ.get("NEQSIM_HOME")
    if configured:
        root = Path(configured).expanduser().resolve()
    else:
        root = next((parent for parent in Path(__file__).resolve().parents
                     if (parent / "distribution/cupbhan/runtime.json").is_file()), None)
        if root is None:
            raise RuntimeError("Set NEQSIM_HOME to the shared NeqSim runtime directory")
    lock = json.loads((root / "distribution/cupbhan/runtime.json").read_text(encoding="utf-8"))
    jar = (root / lock["jar"]).resolve()
    if lock.get("schemaVersion") != 1 or not jar.is_relative_to(root):
        raise RuntimeError("Invalid shared NeqSim runtime lock")
    checksum = hashlib.sha256()
    with jar.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(chunk)
    if checksum.hexdigest() != lock["sha256"]:
        raise RuntimeError(f"NeqSim SHA-256 mismatch: {jar}")
    return jar, lock
