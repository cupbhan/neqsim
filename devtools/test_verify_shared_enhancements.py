"""The personal compatibility inventory must fail closed when data or contracts disappear."""

import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from verify_shared_enhancements import verify


class EnhancementInventoryTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.write("pom.xml", '<project><properties><revision>3.23.0</revision></properties></project>')
        self.write("distribution/cupbhan/upstream.json", json.dumps({"version": "3.23.0", "commit": "a" * 40}))
        self.write("model.json", '{"calibration":"frozen"}')
        self.write("test.java", "fixture")
        policy = {"requiredFiles": ["test.java"], "requiredTools": ["runFieldFluid"],
                  "frozenResources": {"model.json": hashlib.sha256((self.root / "model.json").read_bytes()).hexdigest()}}
        self.write("distribution/cupbhan/personal-enhancements.json", json.dumps(policy))
        self.tool_path = "neqsim-mcp-server/src/main/java/neqsim/mcp/server/NeqSimTools.java"
        self.write(self.tool_path, 'public String runFieldFluid(String input) { return input; }')
        self.git_check = patch("verify_shared_enhancements.subprocess.run")
        self.git_check.start()
        self.addCleanup(self.git_check.stop)

    def write(self, name, value):
        p = self.root / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(value, encoding="utf-8")

    def test_retained_inventory_passes(self):
        baseline, policy = verify(self.root)
        self.assertEqual("3.23.0", baseline["version"])
        self.assertEqual(["runFieldFluid"], policy["requiredTools"])

    def test_lost_personal_test_is_rejected(self):
        (self.root / "test.java").unlink()
        with self.assertRaisesRegex(RuntimeError, "Missing personal enhancement"):
            verify(self.root)

    def test_model_parameter_changes_are_rejected(self):
        self.write("model.json", '{"calibration":"silently replaced"}')
        with self.assertRaisesRegex(RuntimeError, "Frozen personal model data"):
            verify(self.root)

    def test_dropped_interface_is_rejected(self):
        self.write(self.tool_path, 'public String runOther(String input) { return input; }')
        with self.assertRaisesRegex(RuntimeError, "Missing personal MCP contracts"):
            verify(self.root)

    def test_stale_provenance_is_rejected(self):
        self.write("pom.xml", '<project><properties><revision>3.17.0</revision></properties></project>')
        with self.assertRaisesRegex(RuntimeError, "official version"):
            verify(self.root)


if __name__ == "__main__":
    unittest.main()
