"""Inspect an official release, or merge it into a new review branch.

Default: fetch and report the latest stable GitHub release; no checkout or merge.
--apply: require a clean tree, preserve the current commit, and merge on sync/*.
Conflicts are left for review. Never reset, rebase, force-push, or change runtime.json.
"""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET

OFFICIAL = "https://github.com/equinor/neqsim.git"
TAG_PATTERN = re.compile(r"v\d+\.\d+\.\d+\Z")
ROOT = Path(__file__).resolve().parents[1]


def git(root, *args, check=True):
    return subprocess.run(["git", *args], cwd=root, check=check, text=True,
                          encoding="utf-8", stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def value(root, *args):
    return git(root, *args).stdout.strip()


def official_release(tag=None):
    if tag and not TAG_PATTERN.fullmatch(tag):
        raise ValueError("Use a stable official release tag such as v3.23.0")
    suffix = f"tags/{tag}" if tag else "latest"
    headers = {"User-Agent": "cupbhan-neqsim-upstream-sync", "Accept": "application/vnd.github+json"}
    token = os.environ.get("GH_TOKEN")
    if token:
        headers["Authorization"] = "Bearer " + token
    request = urllib.request.Request(f"https://api.github.com/repos/equinor/neqsim/releases/{suffix}",
                                     headers=headers)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                release = json.load(response)
            break
        except urllib.error.HTTPError as error:
            if attempt == 2 or error.code not in (429, 500, 502, 503, 504):
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError):
            if attempt == 2:
                raise
        time.sleep(0.5 * (2 ** attempt))
    if release.get("draft") or release.get("prerelease") or not TAG_PATTERN.fullmatch(release["tag_name"]):
        raise ValueError("Only published stable releases are accepted")
    return {key: release[key] for key in ("tag_name", "published_at", "html_url")}


def fetch_release(root, tag):
    if "upstream" not in value(root, "remote").splitlines():
        git(root, "remote", "add", "upstream", OFFICIAL)
    url = value(root, "remote", "get-url", "upstream")
    if url.rstrip("/") not in (OFFICIAL, OFFICIAL.removesuffix(".git"), "git@github.com:equinor/neqsim.git"):
        raise ValueError("The upstream remote must point to equinor/neqsim")
    # Separate namespace prevents an official tag from overwriting any personal/local tag.
    ref = f"refs/remotes/upstream-tags/{tag}"
    result = git(root, "fetch", "--no-tags", "upstream", f"refs/tags/{tag}:{ref}")
    if result.stderr:
        print(result.stderr.strip(), file=sys.stderr)
    return ref


def inspect(root, upstream_ref, tag):
    if not TAG_PATTERN.fullmatch(tag):
        raise ValueError("Invalid release tag")
    before = value(root, "rev-parse", "HEAD")
    upstream = value(root, "rev-parse", upstream_ref + "^{commit}")
    base = value(root, "merge-base", before, upstream)
    official_version = ET.fromstring(value(root, "show", upstream + ":pom.xml")).findtext(
        "{*}properties/{*}revision")
    if official_version != tag[1:]:
        raise ValueError(f"Release tag/POM version mismatch: {tag} vs {official_version}")
    ours = set(value(root, "diff", "--name-only", base, before).splitlines())
    theirs = set(value(root, "diff", "--name-only", base, upstream).splitlines())
    count = int(value(root, "rev-list", "--count", before + ".." + upstream))
    return {"schemaVersion": 1, "repository": OFFICIAL, "tag": tag, "version": tag[1:],
            "before": before, "commit": upstream, "mergeBase": base,
            "upstreamCommitsToMerge": count, "alreadyIntegrated": count == 0,
            "personalChangedFileCount": len(ours), "upstreamChangedFileCount": len(theirs),
            "overlappingFiles": sorted(ours & theirs),
            "workingTreeDirty": bool(value(root, "status", "--porcelain"))}


def apply(root, report):
    if value(root, "status", "--porcelain"):
        raise RuntimeError("Commit or preserve local edits before applying an upstream merge")
    if value(root, "rev-parse", "HEAD") != report["before"]:
        raise RuntimeError("HEAD changed after inspection; generate a fresh report")
    if report["alreadyIntegrated"]:
        report["status"] = "already-integrated"
        return
    branch = "sync/upstream-" + report["tag"]
    if git(root, "show-ref", "--verify", "--quiet", "refs/heads/" + branch, check=False).returncode == 0:
        raise RuntimeError(f"Review existing branch {branch}; it will not be overwritten")
    backup = f"archive/pre-upstream-{report['tag']}-{report['before'][:10]}"
    git(root, "branch", backup, report["before"])
    git(root, "switch", "-c", branch)
    report.update(branch=branch, backupBranch=backup)
    result = git(root, "merge", "--no-ff", "--no-edit", report["commit"], check=False)
    if result.returncode:
        report.update(status="conflicts", conflicts=value(root, "diff", "--name-only", "--diff-filter=U").splitlines())
        raise RuntimeError("Merge requires review; inspect the conflict list. No conflict resolution was guessed.")
    lock = {"schemaVersion": 1, "repository": OFFICIAL, "tag": report["tag"],
            "version": report["version"], "commit": report["commit"], "previousHead": report["before"]}
    path = root / "distribution/cupbhan/upstream.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8", newline="\n")
    git(root, "add", "distribution/cupbhan/upstream.json")
    git(root, "commit", "-m", "Record official upstream baseline " + report["tag"])
    report.update(status="merged-awaiting-validation", mergedCommit=value(root, "rev-parse", "HEAD"))


def write_report(root, report):
    folder = root / "build/upstream-sync" / report["tag"]
    folder.mkdir(parents=True, exist_ok=True)
    if report.get("conflicts"):
        # CI runners disappear after the job; retain the exact index stages for local resolution.
        (folder / "conflicts.diff").write_text(value(root, "diff", "--cc"), encoding="utf-8")
        for index, filename in enumerate(report["conflicts"]):
            item = folder / "conflicts" / str(index)
            item.mkdir(parents=True, exist_ok=True)
            (item / "path.txt").write_text(filename + "\n", encoding="utf-8")
            for stage, name in ((1, "base"), (2, "personal"), (3, "official")):
                result = subprocess.run(["git", "show", f":{stage}:{filename}"], cwd=root,
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if result.returncode == 0:
                    (item / (name + ".bin")).write_bytes(result.stdout)
    (folder / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    lines = [f"# Official NeqSim {report['tag']} synchronization", "",
             f"Status: {report.get('status', 'inspection-only')}", "",
             f"Official commit: `{report['commit']}`", f"Personal starting commit: `{report['before']}`",
             f"Official commits to merge: {report['upstreamCommitsToMerge']}",
             f"Files changed on both sides: {len(report['overlappingFiles'])}", "",
             "## Files requiring focused review", ""]
    lines += [f"- `{p}`" for p in report.get("conflicts", report["overlappingFiles"])]
    lines += ["", "Runtime selection is unchanged. Build and test a candidate before promoting it.", ""]
    (folder / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({key: report.get(key) for key in ("tag", "commit", "status", "alreadyIntegrated",
          "upstreamCommitsToMerge", "branch")}, indent=2))
    print(f"Report: {folder}")
    if os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as stream:
            for key, val in {"tag": report["tag"], "version": report["version"],
                             "branch": report.get("branch", ""),
                             "changed": str(report.get("status") == "merged-awaiting-validation").lower()}.items():
                stream.write(f"{key}={val}\n")
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as stream:
            stream.write("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", help="Published stable release; default is GitHub latest")
    parser.add_argument("--apply", action="store_true", help="Merge on a new sync branch")
    args = parser.parse_args()
    release = official_release(args.tag)
    ref = fetch_release(ROOT, release["tag_name"])
    report = inspect(ROOT, ref, release["tag_name"])
    report["release"] = release
    try:
        if args.apply:
            apply(ROOT, report)
    finally:
        write_report(ROOT, report)


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
