"""Tests for distribution provenance failures that must prevent a candidate release."""

import subprocess
import tempfile
import unittest
import zipfile
from pathlib import Path

from build_shared_thermo import VERSION_PATTERN, require_clean, test_summary, verify_jar_version


class DistributionGateTest(unittest.TestCase):
    def test_rejects_official_snapshot_and_path_versions(self):
        for value in ("3.17.0", "3.17.0-SNAPSHOT", "../3.17.0-cupbhan.1", "3.17.0-cupbhan.1\n"):
            self.assertIsNone(VERSION_PATTERN.fullmatch(value), value)
        self.assertIsNotNone(VERSION_PATTERN.fullmatch("3.17.0-cupbhan.1-rc.1"))

    def test_rejects_old_core_hidden_in_new_runner(self):
        with tempfile.TemporaryDirectory() as temporary:
            jar = Path(temporary) / "candidate.jar"
            with zipfile.ZipFile(jar, "w") as archive:
                archive.writestr("META-INF/maven/com.equinor.neqsim/neqsim/pom.properties", "version=3.16.0\n")
            with self.assertRaisesRegex(RuntimeError, "version mismatch"):
                verify_jar_version(jar, "neqsim", "3.17.0-cupbhan.1-rc.1")

    def test_accepts_matching_embedded_version(self):
        with tempfile.TemporaryDirectory() as temporary:
            jar = Path(temporary) / "candidate.jar"
            with zipfile.ZipFile(jar, "w") as archive:
                archive.writestr("META-INF/maven/com.equinor.neqsim/neqsim/pom.properties",
                                 "# generated\nversion=3.17.0-cupbhan.1-rc.1\n")
            verify_jar_version(jar, "neqsim", "3.17.0-cupbhan.1-rc.1")

    def test_missing_reports_cannot_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(RuntimeError):
                test_summary(temporary)

    def test_failed_or_skipped_reports_cannot_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            report = Path(temporary) / "TEST-regression.xml"
            for field in ("failures", "errors", "skipped"):
                report.write_text(f'<testsuite name="regression" tests="1" {field}="1"/>')
                with self.assertRaises(RuntimeError):
                    test_summary(temporary)

    def test_reports_preserve_selected_suite_names_and_counts(self):
        with tempfile.TemporaryDirectory() as temporary:
            report = Path(temporary) / "TEST-regression.xml"
            report.write_text('<testsuite name="regression" tests="7" failures="0" errors="0" skipped="0"/>')
            summary = test_summary(temporary)
            self.assertEqual(7, summary["tests"])
            self.assertEqual(["regression"], summary["suites"])

    def test_untracked_sources_and_modified_sources_prevent_release(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            def run(*args):
                subprocess.run(["git", *args], cwd=root, check=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
            run("init")
            # This identity belongs only to an ephemeral test fixture, never the user's repository.
            run("config", "user.name", "Distribution test fixture")
            run("config", "user.email", "fixture@example.invalid")
            source = root / "source.txt"
            source.write_text("committed source\n")
            run("add", "source.txt")
            run("commit", "-m", "test fixture")
            self.assertEqual(40, len(require_clean(root)))
            extra = root / "untracked.java"
            extra.write_text("untracked source\n")
            with self.assertRaisesRegex(RuntimeError, "clean committed tree"):
                require_clean(root)
            extra.unlink()
            source.write_text("changed source\n")
            with self.assertRaisesRegex(RuntimeError, "clean committed tree"):
                require_clean(root)


if __name__ == "__main__":
    unittest.main()
