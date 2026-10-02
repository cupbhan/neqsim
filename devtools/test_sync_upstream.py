"""Exercise upstream synchronization with real disposable Git histories, offline."""

import json
from pathlib import Path
import tempfile
import unittest

from sync_upstream import apply, git, inspect, value


class UpstreamSyncTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        git(self.root, "init", "-b", "master")
        for key, val in [("user.name", "Sync test fixture"), ("user.email", "fixture@example.invalid"),
                         ("core.autocrlf", "false"), ("commit.gpgsign", "false")]:
            git(self.root, "config", key, val)
        self.write(".gitignore", "build/\n")
        self.write("pom.xml", self.pom("3.17.0"))
        self.write("shared.txt", "baseline\n")
        self.write("distribution/cupbhan/runtime.json", '{"version":"old-tested-runtime"}\n')
        self.commit("baseline")
        git(self.root, "switch", "-c", "official")
        self.write("pom.xml", self.pom("3.23.0"))
        self.write("shared.txt", "official implementation\n")
        self.commit("official release")
        self.official = value(self.root, "rev-parse", "HEAD")
        git(self.root, "tag", "v3.23.0")
        git(self.root, "switch", "master")
        self.write("personal.txt", "personal improvement\n")
        self.commit("personal improvement")

    @staticmethod
    def pom(version):
        return f'<project><properties><revision>{version}</revision></properties></project>\n'

    def write(self, name, content):
        p = self.root / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content, encoding="utf-8", newline="\n")

    def commit(self, message):
        git(self.root, "add", ".")
        git(self.root, "commit", "-m", message)

    def plan(self):
        return inspect(self.root, "v3.23.0", "v3.23.0")

    def test_inspection_does_not_modify_branch_or_files(self):
        before = value(self.root, "rev-parse", "HEAD")
        report = self.plan()
        self.assertEqual(1, report["upstreamCommitsToMerge"])
        self.assertEqual(before, value(self.root, "rev-parse", "HEAD"))
        self.assertEqual("master", value(self.root, "branch", "--show-current"))
        self.assertEqual("", value(self.root, "status", "--porcelain"))

    def test_merge_retains_personal_work_history_and_runtime_pin(self):
        report = self.plan()
        runtime = (self.root / "distribution/cupbhan/runtime.json").read_bytes()
        apply(self.root, report)
        self.assertEqual("merged-awaiting-validation", report["status"])
        self.assertEqual("personal improvement\n", (self.root / "personal.txt").read_text())
        self.assertEqual(runtime, (self.root / "distribution/cupbhan/runtime.json").read_bytes())
        for ancestor in (report["before"], self.official):
            self.assertEqual(0, git(self.root, "merge-base", "--is-ancestor", ancestor, "HEAD").returncode)
        lock = json.loads((self.root / "distribution/cupbhan/upstream.json").read_text())
        self.assertEqual(self.official, lock["commit"])
        self.assertEqual("", value(self.root, "status", "--porcelain"))

    def test_dirty_work_is_refused_and_preserved(self):
        report = self.plan()
        self.write("personal.txt", "unsaved work\n")
        with self.assertRaisesRegex(RuntimeError, "local edits"):
            apply(self.root, report)
        self.assertEqual("unsaved work\n", (self.root / "personal.txt").read_text())

    def test_existing_sync_branch_is_never_overwritten(self):
        report = self.plan()
        git(self.root, "branch", "sync/upstream-v3.23.0")
        with self.assertRaisesRegex(RuntimeError, "not be overwritten"):
            apply(self.root, report)
        self.assertEqual(report["before"], value(self.root, "rev-parse", "HEAD"))

    def test_conflict_is_left_visible_without_selecting_either_side(self):
        self.write("shared.txt", "personal implementation\n")
        self.commit("personal conflict")
        report = self.plan()
        with self.assertRaisesRegex(RuntimeError, "requires review"):
            apply(self.root, report)
        self.assertEqual(["shared.txt"], report["conflicts"])
        self.assertIn("<<<<<<<", (self.root / "shared.txt").read_text())
        self.assertEqual(report["before"], value(self.root, "rev-parse", "HEAD"))
        self.assertEqual(report["before"], value(self.root, "rev-parse", report["backupBranch"]))

    def test_repeated_sync_is_noop_after_integration(self):
        apply(self.root, self.plan())
        report = self.plan()
        before = value(self.root, "rev-parse", "HEAD")
        apply(self.root, report)
        self.assertEqual("already-integrated", report["status"])
        self.assertEqual(before, value(self.root, "rev-parse", "HEAD"))

    def test_wrong_release_tag_or_version_is_rejected(self):
        with self.assertRaises(ValueError):
            inspect(self.root, "v3.23.0", "v3.24.0")
        with self.assertRaises(ValueError):
            inspect(self.root, "v3.23.0", "../../master")

    def test_changed_head_invalidates_plan(self):
        report = self.plan()
        self.write("new.txt", "later work\n")
        self.commit("later work")
        with self.assertRaisesRegex(RuntimeError, "HEAD changed"):
            apply(self.root, report)


if __name__ == "__main__":
    unittest.main()
