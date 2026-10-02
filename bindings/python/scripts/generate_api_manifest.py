"""Generate the offline NeqSim API manifest consumed by ``neqsim.discovery``.

Runs once at release time. It reflects the whole Java API (via ``neqsim.discovery``
live reflection) and writes ``src/neqsim/data/neqsim_api.json`` containing every
class with its constructors and public method names. Shipping that file makes
``discovery.list_classes`` / ``find_classes`` / ``describe`` work instantly and
JVM-free for end users, and it can drive a docs site.

Usage
-----
    python scripts/generate_api_manifest.py [--out PATH] [--packages neqsim ...]
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir():
    sys.path.insert(0, str(_SRC))

from neqsim import discovery  # noqa: E402


def _class_entry(fqcn: str) -> Dict[str, List[str]]:
    """Build the manifest entry (constructors + methods) for one class.

    Args:
        fqcn: Fully-qualified Java class name.

    Returns:
        A dict with ``constructors`` and ``methods`` lists (possibly empty).
    """
    constructors: List[str] = []
    methods: List[str] = []
    try:
        java_class = discovery.get_class(fqcn).class_
        constructors = sorted(
            discovery._simplify_signature(str(c)) for c in java_class.getConstructors()
        )
        method_names = set()
        for m in java_class.getMethods():
            name = str(m.getName())
            if name.isidentifier():
                method_names.add(name)
        methods = sorted(method_names)
    except Exception:
        pass
    return {"constructors": constructors, "methods": methods}


def _detect_version() -> str:
    """Return the NeqSim jar version if derivable, else ``"unknown"``."""
    try:
        for entry in discovery._classpath_entries():
            name = Path(entry).name
            if name.startswith("neqsim-") and name.endswith(".jar"):
                return name[len("neqsim-") : -len(".jar")]
    except Exception:
        pass
    return "unknown"


def generate(out_path: Path, packages: List[str]) -> int:
    """Generate the manifest for the requested packages.

    Args:
        out_path: Destination JSON file.
        packages: Package prefixes to include.

    Returns:
        The number of classes written.
    """
    class_names: List[str] = []
    for prefix in packages:
        class_names.extend(discovery.list_classes(prefix, recursive=True))
    class_names = sorted(set(class_names))

    classes: Dict[str, Dict[str, List[str]]] = {}
    for fqcn in class_names:
        classes[fqcn] = _class_entry(fqcn)

    manifest = {
        "schema": "neqsim-api-manifest/1.0",
        "version": _detect_version(),
        "generated": datetime.now(timezone.utc).isoformat(),
        "class_count": len(classes),
        "classes": classes,
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(manifest, indent=1, sort_keys=True)
    if out_path.suffix == ".gz":
        with gzip.open(out_path, "wt", encoding="utf-8") as handle:
            handle.write(payload)
    else:
        with open(out_path, "w", encoding="utf-8") as handle:
            handle.write(payload)
    return len(classes)


def main(argv: List[str] | None = None) -> int:
    """CLI entry point.

    Args:
        argv: Optional argument list (defaults to ``sys.argv``).

    Returns:
        Process exit code.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        default=str(_SRC / "neqsim" / "data" / "neqsim_api.json.gz"),
        help="Output path for the manifest JSON (use .gz to gzip).",
    )
    parser.add_argument(
        "--packages",
        nargs="*",
        default=["neqsim"],
        help="Package prefixes to include (default: the whole neqsim tree).",
    )
    args = parser.parse_args(argv)

    count = generate(Path(args.out), args.packages)
    print(f"Wrote manifest with {count} classes to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
