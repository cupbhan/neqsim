"""Shared runtime integrity tests; no JVM or optional Python dependencies needed."""

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from activate_shared_thermo import activate

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "shared_resolver", ROOT / "bindings/python/src/neqsim/_shared_runtime.py")
resolver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(resolver)


class SharedRuntimeTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "distribution/cupbhan").mkdir(parents=True)
        self.source = self.root / "candidate"
        self.source.mkdir()
        self.version = "3.17.0-cupbhan.1-rc.1"
        self.name = f"neqsim-mcp-server-{self.version}-runner.jar"
        self.bytes = b"fixture, not an executable JAR"
        (self.source / self.name).write_bytes(self.bytes)
        self.metadata = {"version": self.version, "dirty": False, "status": "tested-candidate",
                         "commit": "a" * 40, "artifacts": [{"file": self.name,
                         "sha256": hashlib.sha256(self.bytes).hexdigest()}]}
        self.write_metadata()
        self.env = patch.dict(os.environ, {"NEQSIM_HOME": str(self.root)})
        self.env.start()
        self.addCleanup(self.env.stop)

    def write_metadata(self):
        (self.source / "SOURCE.json").write_text(json.dumps(self.metadata), encoding="utf-8")

    def test_install_and_resolve_repeatably(self):
        expected = activate(self.source, self.root)
        self.assertEqual(activate(self.source, self.root), expected)
        jar, lock = resolver.resolve_runtime()
        self.assertEqual(jar.read_bytes(), self.bytes)
        self.assertEqual(lock["commit"], self.metadata["commit"])

    def test_corrupt_candidate_cannot_replace_lock(self):
        lock = activate(self.source, self.root)
        (self.source / self.name).write_bytes(b"corrupt")
        with self.assertRaisesRegex(ValueError, "Checksum"):
            activate(self.source, self.root)
        self.assertEqual(resolver.resolve_runtime()[1], lock)

    def test_existing_release_is_immutable(self):
        activate(self.source, self.root)
        (self.source / self.name).write_bytes(b"different valid artifact")
        self.metadata["artifacts"][0]["sha256"] = hashlib.sha256(b"different valid artifact").hexdigest()
        self.write_metadata()
        with self.assertRaisesRegex(ValueError, "overwrite"):
            activate(self.source, self.root)

    def test_corrupt_runtime_is_rejected(self):
        lock = activate(self.source, self.root)
        (self.root / lock["jar"]).write_bytes(b"corrupt")
        with self.assertRaisesRegex(RuntimeError, "SHA-256"):
            resolver.resolve_runtime()

    def test_path_escape_is_rejected(self):
        lock = activate(self.source, self.root)
        lock["jar"] = "../outside.jar"
        (self.root / "distribution/cupbhan/runtime.json").write_text(json.dumps(lock))
        with self.assertRaisesRegex(RuntimeError, "Invalid"):
            resolver.resolve_runtime()

    def test_dirty_source_and_invalid_version_rejected(self):
        self.metadata["dirty"] = True
        self.write_metadata()
        with self.assertRaisesRegex(ValueError, "clean"):
            activate(self.source, self.root)
        self.metadata["dirty"] = False
        self.metadata["version"] = ".."
        self.write_metadata()
        with self.assertRaisesRegex(ValueError, "version"):
            activate(self.source, self.root)


if __name__ == "__main__":
    unittest.main()
