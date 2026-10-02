"""Install a verified candidate and atomically select it for local consumers.

Usage: python devtools/activate_shared_thermo.py build/shared-thermo/releases/VERSION
The runtime lock is tracked; binaries and historical directories remain local.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def activate(source, root):
    source = source.resolve()
    metadata = json.loads((source / "SOURCE.json").read_text(encoding="utf-8"))
    version = metadata["version"]
    if not re.fullmatch(r"\d+\.\d+\.\d+-cupbhan\.\d+(?:-rc\.\d+)?", version):
        raise ValueError("Invalid release version")
    if metadata.get("dirty") is not False or metadata.get("status") != "tested-candidate":
        raise ValueError("Only clean, tested candidates can be activated")
    for artifact in metadata["artifacts"]:
        name = artifact["file"]
        if Path(name).name != name or "/" in name or "\\" in name:
            raise ValueError("Invalid artifact name")
        if digest(source / name) != artifact["sha256"]:
            raise ValueError(f"Checksum mismatch: {name}")
    runner = f"neqsim-mcp-server-{version}-runner.jar"
    entry = next(a for a in metadata["artifacts"] if a["file"] == runner)
    target = root / "runtime" / "releases" / version
    if target.exists():
        for file in source.iterdir():
            if file.is_file() and (not (target / file.name).is_file()
                                   or digest(file) != digest(target / file.name)):
                raise ValueError(f"Refusing to overwrite existing release: {target}")
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, target)
    if digest(target / runner) != entry["sha256"]:
        raise ValueError("Installed runner checksum mismatch")
    lock = {
        "schemaVersion": 1, "version": version, "status": metadata["status"],
        "commit": metadata["commit"],
        "jar": (target / runner).relative_to(root).as_posix(),
        "sha256": entry["sha256"],
        "source": (target / "SOURCE.json").relative_to(root).as_posix(),
    }
    lock_path = root / "distribution" / "cupbhan" / "runtime.json"
    temporary = lock_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")
    temporary.replace(lock_path)
    return lock


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()
    print(json.dumps(activate(args.candidate, Path(__file__).resolve().parents[1]), indent=2))
