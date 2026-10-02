"""Generate a fully-typed Python client for the NeqSim REST API (JVM-free path).

For users who consume NeqSim through the cloud REST service (NeqSimAPI) rather
than the in-process JPype bridge, a generated client gives real Pydantic models
with autocomplete and validation and needs no local JVM.

This wraps the ``openapi-python-client`` generator. It does not vendor an API
spec: pass the URL or file of the OpenAPI/Swagger document for the deployment you
target.

Prerequisites
-------------
    pip install openapi-python-client

Usage
-----
    python scripts/generate_openapi_client.py --spec https://<host>/openapi.json
    python scripts/generate_openapi_client.py --spec openapi.json --out clients/

The generated package is importable and typed, e.g.::

    from neqsim_api_client import Client
    from neqsim_api_client.api.default import some_endpoint
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List

_REPO_ROOT = Path(__file__).resolve().parent.parent


def _has_generator() -> bool:
    """Return whether ``openapi-python-client`` is available on PATH.

    Returns:
        True if the generator CLI is installed.
    """
    return shutil.which("openapi-python-client") is not None


def generate(spec: str, out_dir: Path, overwrite: bool) -> int:
    """Run ``openapi-python-client`` for the given spec.

    Args:
        spec: A URL or filesystem path to the OpenAPI document.
        out_dir: Directory to generate the client into.
        overwrite: Whether to overwrite an existing generated client.

    Returns:
        The subprocess exit code.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    is_url = spec.startswith(("http://", "https://"))
    if is_url:
        source_flag, source_value = "--url", spec
    else:
        # Resolve to an absolute path so it still works with cwd=out_dir.
        source_flag, source_value = "--path", str(Path(spec).resolve())
    cmd: List[str] = [
        "openapi-python-client",
        "generate",
        source_flag,
        source_value,
    ]
    if overwrite:
        cmd.append("--overwrite")
    print("Running:", " ".join(cmd), f"(cwd={out_dir})")
    return subprocess.call(cmd, cwd=str(out_dir))


def main(argv: List[str] | None = None) -> int:
    """CLI entry point.

    Args:
        argv: Optional argument list (defaults to ``sys.argv``).

    Returns:
        Process exit code.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spec",
        required=True,
        help="URL or path to the OpenAPI/Swagger JSON for the NeqSim REST API.",
    )
    parser.add_argument(
        "--out",
        default=str(_REPO_ROOT / "clients"),
        help="Directory to generate the client package into (default: clients/).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing generated client.",
    )
    args = parser.parse_args(argv)

    if not _has_generator():
        print(
            "openapi-python-client is not installed.\n"
            "Install it first:  pip install openapi-python-client",
            file=sys.stderr,
        )
        return 2

    return generate(args.spec, Path(args.out), args.overwrite)


if __name__ == "__main__":
    raise SystemExit(main())
