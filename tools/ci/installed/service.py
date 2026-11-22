"""Read-only GitHub identities/downloads; tokens never enter scientific children.

The terminal receiver runs this trusted checkout code. Downloaded source/tools
are data only and are never imported, invoked or checked out for execution.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import zipfile
import time
import urllib.error
import urllib.parse
import urllib.request

# Only this trusted checkout's physical helper directory enters sys.path.
CI = Path(__file__).resolve().parent
if CI.is_symlink() or not CI.is_dir():
    raise ValueError("Trusted receiver source directory required")
sys.path.insert(0, str(CI))
from evidence import directory_bounds, member_names, EVIDENCE_LIMIT, MIB, document, json_value, require
from verify import digest, write_json

API = "https://api.github.com"
API_VERSION = "2026-03-10"
WORKFLOW_PATH = ".github/workflows/installed-package.yaml"
JOB_NAMES = {"linux-gold": "linux-gold", "windows": "installed-windows"}


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def repository(value):
    require(isinstance(value, str) and re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", value),
            "Explicit GitHub repository identity required")
    return value


def numeric(value):
    require(re.fullmatch("[1-9][0-9]*", str(value)) is not None, "Positive GitHub service identity required")
    return str(value)


def now():
    return datetime.now(timezone.utc).isoformat()


def get(path, *, raw=False):
    require(path.startswith("/repos/") and not any(char in path for char in "\r\n#"),
            "Read-only repository API path required")
    request = urllib.request.Request(API + path, headers={
        "Accept": "application/vnd.github+json", "X-GitHub-Api-Version": API_VERSION,
        "User-Agent": "MungeSumstats-installed-CI/1"})
    token = os.environ.get("GITHUB_TOKEN")
    if token:
        request.add_header("Authorization", "Bearer " + token)
    with urllib.request.build_opener(NoRedirect).open(request, timeout=30) as response:
        require(response.status == 200 and response.url == API + path, "Unexpected GitHub API response identity")
        data = response.read(16 * MIB + 1)
    require(len(data) <= 16 * MIB, "Bounded service metadata required")
    value = json_value(data)
    return (data, value) if raw else value


def event_identity(env, payload):
    repo = repository(env["GITHUB_REPOSITORY"])
    require(payload["repository"]["full_name"] == repo, "GitHub event repository differs")
    kind, checkout = env["GITHUB_EVENT_NAME"], env["GITHUB_SHA"]
    require(kind in {"push", "pull_request", "workflow_dispatch"} and re.fullmatch("[0-9a-f]{40}", checkout),
            "Unsupported triggering source event")
    selected = dict(event=kind, repository=repo, checkout_commit=checkout, ref=env["GITHUB_REF"])
    if kind == "pull_request":
        pr = payload["pull_request"]
        selected.update(pull_request=int(payload["number"]), head_commit=pr["head"]["sha"],
            base_commit=pr["base"]["sha"], head_repository=pr["head"]["repo"]["full_name"])
        require(selected["ref"] == "refs/pull/" + str(selected["pull_request"]) + "/merge"
                and all(re.fullmatch("[0-9a-f]{40}", selected[name]) for name in ("head_commit", "base_commit")),
                "Expected exact tested pull-request merge identity")
    else:
        selected["head_commit"] = checkout
        require(selected["ref"].startswith(("refs/heads/", "refs/tags/")), "Expected explicit push/dispatch ref")
    return selected


def observe(path, lane):
    env = os.environ
    repo = repository(env["GITHUB_REPOSITORY"])
    run_id, attempt = numeric(env["GITHUB_RUN_ID"]), numeric(env["GITHUB_RUN_ATTEMPT"])
    payload = document(Path(env["GITHUB_EVENT_PATH"]))
    event = event_identity(env, payload)
    endpoint = "/repos/" + repo + "/actions/runs/" + run_id + "/attempts/" + attempt
    run = get(endpoint)
    require(str(run["id"]) == run_id and str(run["run_attempt"]) == attempt
            and run["repository"]["full_name"] == repo and run["path"] == WORKFLOW_PATH,
            "Actual workflow run/attempt/source differs")
    jobs_path = endpoint + "/jobs?per_page=100"
    jobs = get(jobs_path)
    require(jobs["total_count"] <= 100 and len(jobs["jobs"]) == jobs["total_count"], "Bounded complete producer Jobs census")
    matches = [job for job in jobs["jobs"] if job["name"] == JOB_NAMES[lane] and job["status"] == "in_progress"]
    require(len(matches) == 1, "Exact current lane Jobs API record required")
    job = matches[0]
    require(job["run_id"] == int(run_id) and job["head_sha"] == run["head_sha"]
            and run["head_sha"] in {event["head_commit"], event["checkout_commit"]} and env["GITHUB_JOB"] == JOB_NAMES[lane]
            and job["runner_name"] == env["RUNNER_NAME"], "Jobs/event/current-runner identity differs")
    expected_label = "windows-2022" if lane == "windows" else "ubuntu-22.04"
    require(sorted(job["labels"]) == [expected_label], "Standard hosted lane label required")
    if lane == "windows":
        import windows_job
        host = windows_job.windows_host()
        require(env["RUNNER_OS"] == "Windows" and env["RUNNER_ARCH"] == "X64" and env["ImageOS"] == "win22",
                "Native standard-hosted Windows2022x64 required")
    else:
        require(platform.system() == "Linux" and platform.machine().lower() == "x86_64"
                and env["RUNNER_OS"] == "Linux" and env["RUNNER_ARCH"] == "X64", "Native Linux x64 required")
        host = dict(name=platform.node(), os="Linux", build=platform.release(), os_arch="X64", process_arch="X64")
    observed_tick = time.perf_counter()
    observed = now()
    record = dict(schema="mss-current-job-start-v1", lane=lane, repository=repo, run_id=run_id,
        run_attempt=attempt, job_id=numeric(job["id"]), started_utc=job["started_at"],
        observed_utc=observed, observed_monotonic_seconds=observed_tick, host=host,
        image=dict(os=env.get("ImageOS"), version=env.get("ImageVersion")),
        event=event, jobs=jobs, run=run, request_url=API + jobs_path, API_version=API_VERSION)
    write_json(path, record)
    return record


def bind_running(start, env, checkout):
    require(start["schema"] == "mss-current-job-start-v1", "Actual current first observation required")
    repo, run_id, attempt = start["repository"], start["run_id"], start["run_attempt"]
    require(env["GITHUB_REPOSITORY"] == repo and env["GITHUB_RUN_ID"] == run_id
            and env["GITHUB_RUN_ATTEMPT"] == attempt and env["GITHUB_SHA"] == start["event"]["checkout_commit"]
            and env["GITHUB_EVENT_NAME"] == start["event"]["event"], "Current GitHub environment differs")
    actual = subprocess.check_output(["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True, timeout=10).strip()
    parents = subprocess.check_output(["git", "-C", str(checkout), "show", "-s", "--format=%P", "HEAD"],
                                     text=True, timeout=10).strip().split()
    require(actual == start["event"]["checkout_commit"], "Actual triggering checkout differs")
    if start["event"]["event"] == "pull_request":
        require(parents == [start["event"]["base_commit"], start["event"]["head_commit"]],
                "Tested merge parents do not bind event base/head")
    else:
        require(actual == start["run"]["head_sha"], "Push/dispatch checkout differs from service head")
    job, = [row for row in start["jobs"]["jobs"] if str(row["id"]) == start["job_id"]]
    require(job["name"] == JOB_NAMES[start["lane"]] and job["runner_name"] == env["RUNNER_NAME"]
            and job["status"] == "in_progress", "Actual initial job identity differs")
    return dict(schema="mss-current-runner-v1", repository=repo, run_id=run_id, run_attempt=attempt,
        job_id=start["job_id"], job_name=job["name"], runner_name=job["runner_name"], labels=job["labels"],
        event=start["event"], checkout_commit=actual, checkout_parents=parents, host=start["host"], image=start["image"])


def validate_event(run, event):
    require(event["repository"] == run["repository"]["full_name"] and event["event"] == run["event"]
            and re.fullmatch("[0-9a-f]{40}", event["checkout_commit"]), "Event/service source identity differs")
    if event["event"] == "pull_request":
        rows = [row for row in run["pull_requests"] if row["number"] == event["pull_request"]]
        require(len(rows) == 1, "Exact service pull-request provenance required")
        pr = rows[0]
        require(pr["head"]["sha"] == event["head_commit"] and pr["base"]["sha"] == event["base_commit"]
                and run["head_sha"] in {event["head_commit"], event["checkout_commit"]}
                and event["ref"] == "refs/pull/" + str(event["pull_request"]) + "/merge",
                "Original service PR head/base/tested merge differs")
    else:
        require(event["event"] in {"push", "workflow_dispatch"} and
                event["checkout_commit"] == event["head_commit"] == run["head_sha"],
                "Push/dispatch service checkout differs")


def validate_commit(event, commit):
    require(commit["sha"] == event["checkout_commit"], "Tested Git commit differs")
    if event["event"] == "pull_request":
        require([row["sha"] for row in commit["parents"]] == [event["base_commit"], event["head_commit"]],
                "Tested merge parents do not bind original PR base/head")


def download_artifact(repo, item, target):
    repo, artifact_id = repository(repo), numeric(item["id"])
    if os.name == "posix":
        import resource
        resource.setrlimit(resource.RLIMIT_CORE, (0,0))
    require(not target.exists() and not target.is_symlink(), "Fresh exact original ZIP destination required")
    endpoint = API + "/repos/" + repo + "/actions/artifacts/" + artifact_id + "/zip"
    require(item["archive_download_url"] == endpoint, "Artifact endpoint identity differs")
    headers = {"Accept": "application/vnd.github+json", "X-GitHub-Api-Version": API_VERSION,
               "User-Agent": "MungeSumstats-installed-CI/1"}
    if os.environ.get("GITHUB_TOKEN"):
        headers["Authorization"] = "Bearer " + os.environ["GITHUB_TOKEN"]
    request = urllib.request.Request(endpoint, headers=headers)
    try:
        response = urllib.request.build_opener(NoRedirect).open(request, timeout=30)
    except urllib.error.HTTPError as error:
        require(error.code == 302, "Expected artifact's signed download redirect")
        location = error.headers["Location"]
        parsed = urllib.parse.urlsplit(location)
        require(parsed.scheme == "https" and parsed.hostname and not parsed.username and not parsed.password,
                "Artifact download HTTPS required")
        # Deliberately omit the GitHub token from the signed public-byte request.
        error.close()
        # No further redirects: neither a bearer token nor a signed URL may be
        # forwarded to another authority by this trusted acquisition client.
        try:
            response = urllib.request.build_opener(NoRedirect).open(
                urllib.request.Request(location, headers={"User-Agent": headers["User-Agent"]}), timeout=60)
        except Exception:
            raise ValueError("Artifact signed-byte acquisition failed") from None
    require(response.status == 200, "Artifact bytes download failed")
    size = item["size_in_bytes"]
    require(type(size) is int and 0 < size <= EVIDENCE_LIMIT + 96 * MIB, "Bounded exact service artifact size required")
    count, h = 0, hashlib.sha256()
    deadline = time.monotonic() + 600
    with response, target.open("xb") as stream:
        while True:
            require(time.monotonic() < deadline, "Bounded artifact download deadline")
            block = response.read(min(MIB, size + 1 - count))
            if not block:
                break
            count += len(block)
            require(count <= size, "Artifact exceeds exact service byte bound")
            h.update(block)
            stream.write(block)
    require(count == size and item["digest"] == "sha256:" + h.hexdigest(), "Exact ZIP bytes/service digest differ")
    return dict(bytes=count, sha256=h.hexdigest(), artifact_id=artifact_id, request_url=endpoint)


def terminal(repo, run_id, attempt, output):
    repo, run_id, attempt = repository(repo), numeric(run_id), numeric(attempt)
    require(not output.exists(), "Fresh terminal-service custody root")
    output.mkdir()
    base = "/repos/" + repo
    endpoint = base + "/actions/runs/" + run_id + "/attempts/" + attempt
    documents = {}
    for name, path in {"run.json": endpoint, "jobs.json": endpoint + "/jobs?per_page=100",
                       "artifacts.json": base + "/actions/runs/" + run_id + "/artifacts?per_page=100"}.items():
        raw, value = get(path, raw=True)
        (output / name).write_bytes(raw)
        documents[name] = value
    run, jobs, artifacts = (documents[name] for name in ("run.json", "jobs.json", "artifacts.json"))
    require(run["id"] == int(run_id) and run["run_attempt"] == int(attempt) and run["status"] == "completed"
            and run["repository"]["full_name"] == repo and run["path"] == WORKFLOW_PATH, "Exact terminal producer run required")
    require(jobs["total_count"] == len(jobs["jobs"]) == 2
            and {row["name"] for row in jobs["jobs"]} == set(JOB_NAMES.values())
            and all(row["status"] == "completed" and row["completed_at"] and row["run_id"] == int(run_id)
                    and row["head_sha"] == run["head_sha"] for row in jobs["jobs"]), "Complete exact terminal producer Jobs required")
    require(artifacts["total_count"] == len(artifacts["artifacts"]) <= 100, "Complete bounded service artifact list required")
    retained, missing = {}, {}
    for lane in JOB_NAMES:
        name = "mss-installed-" + lane + "-" + run_id + "-" + attempt
        matches = [row for row in artifacts["artifacts"] if row["name"] == name]
        require(len(matches) <= 1, "Ambiguous exact lane artifact")
        if not matches:
            missing[lane] = "absent_exact_lane_artifact"
            continue
        item = matches[0]
        require(item["expired"] is False and item["workflow_run"]["id"] == int(run_id)
                and item["workflow_run"]["head_sha"] == run["head_sha"], "Artifact service/source/expiry differs")
        retained[lane] = dict(name=name, service=item,
            original_zip=download_artifact(repo, item, output / (lane + ".zip")))
    sources = {}
    for lane in retained:
        with (output / (lane + ".zip")).open("rb") as stream:
            directory_bounds(stream, retained[lane]["original_zip"]["bytes"])
            with zipfile.ZipFile(stream) as archive:
                files, _, _ = member_names(archive, lane=lane)
                require("job-start.json" in files and files["job-start.json"].file_size <= MIB,
                        "Bounded exact original first observation required")
                start = json_value(archive.read(files["job-start.json"]))
        event = start["event"]
        validate_event(run, event)
        commit = event["checkout_commit"]
        if commit not in sources:
            commit_raw, commit_doc = get(base + "/git/commits/" + commit, raw=True)
            require(commit_doc["sha"] == commit, "Git commit service identity differs")
            tree = commit_doc["tree"]["sha"]
            tree_raw, tree_doc = get(base + "/git/trees/" + tree + "?recursive=1", raw=True)
            require(tree_doc["sha"] == tree and tree_doc["truncated"] is False,
                    "Complete exact Git tree required")
            validate_commit(event, commit_doc)
            for name, raw in (("git-commit-" + commit + ".json", commit_raw),
                              ("git-tree-" + tree + ".json", tree_raw)):
                (output / name).write_bytes(raw)
            sources[commit] = {"tree": tree, "commit_document": "git-commit-" + commit + ".json",
                               "tree_document": "git-tree-" + tree + ".json"}
    receipt = dict(schema="mss-terminal-service-v1", repository=repo, run_id=run_id, run_attempt=attempt,
        fetched_utc=now(), run_conclusion=run["conclusion"], qualified=False, sources=sources,
        missing_artifacts=missing, acquisition_complete=not missing,
        terminal_metadata={path.name:digest(path) for path in sorted(output.glob("*.json"))}, artifacts=retained,
        artifact_code_execution=False, receiver_verification_still_required=True)
    write_json(output / "service-receipt.json", receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    first = sub.add_parser("observe")
    first.add_argument("--lane", choices=JOB_NAMES, required=True)
    first.add_argument("--out", type=Path, required=True)
    fetch = sub.add_parser("terminal")
    fetch.add_argument("--repository", required=True)
    fetch.add_argument("--run-id", required=True)
    fetch.add_argument("--attempt", required=True)
    fetch.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "observe":
        value = observe(args.out, args.lane)
        print(json.dumps({key:value[key] for key in ("lane", "run_id", "run_attempt", "job_id", "started_utc")}))
    else:
        value = terminal(args.repository, args.run_id, args.attempt, args.out)
        print(json.dumps({key:value[key] for key in ("repository", "run_id", "run_attempt", "run_conclusion")}))


if __name__ == "__main__":
    main()
