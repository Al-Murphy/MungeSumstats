"""Build an explicit pair of local commits with the accepted Darwin build producer."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import os
from pathlib import Path
import re
import subprocess
import sys


RUNTIME = Path("/private/tmp/mss-fullfile-20261121/runtime")
CALLER = RUNTIME / "matched-builds.py"
PRODUCER = RUNTIME / "namespace-builds.py"
PRODUCER_SHA = "615535da8374c81d3d7f23f0300c313d26d06ebea638e73239492243f851a65f"


def commit_id(value):
    if not re.fullmatch(r"[0-9a-f]{40}", value) or len(set(value)) == 1:
        raise argparse.ArgumentTypeError("expected a nonplaceholder 40-character lowercase Git commit ID")
    return value


def load_producer():
    if PRODUCER.resolve(strict=True) != PRODUCER or not PRODUCER.is_file():
        raise ValueError(f"Expected the physical accepted producer: {PRODUCER}")
    with PRODUCER.open("rb") as stream:
        data = stream.read(1024**2 + 1)
    if len(data) > 1024**2 or hashlib.sha256(data).hexdigest() != PRODUCER_SHA:
        raise ValueError(f"Accepted producer digest changed: {PRODUCER}")
    spec = importlib.util.spec_from_file_location("matched_frozen_producer", PRODUCER)
    producer = importlib.util.module_from_spec(spec)
    exec(compile(data, str(PRODUCER), "exec"), producer.__dict__)
    return producer


def run(root, commits):
    started = datetime.now(timezone.utc)
    if sys.platform != "darwin" or not (sys.flags.isolated and sys.flags.dont_write_bytecode and sys.flags.no_site):
        raise ValueError("Use the pinned Darwin Python with -I -B -S")
    producer = load_producer()
    require, physical, write_json = producer.require, producer.physical, producer.write_json
    require(Path(sys.executable) == physical(producer.PYTHON), "Use the pinned physical Python executable")
    require(physical(Path(__file__)) == CALLER, f"Use the physical caller source: {CALLER}")
    require(root.is_absolute() and root == physical(root.parent) / root.name and not root.is_relative_to(producer.REPO),
            "Expected an absolute physical output path outside the source repository")
    os.umask(0o077)
    root.mkdir(mode=0o700)
    info = root.lstat()
    binding, warm = producer.load_binding()
    require(warm.SOURCES[warm.GUARD] == producer.GUARD_SHA, "Unexpected authenticated guard")
    guard_pin = warm.pin(warm.GUARD, producer.GUARD_SHA)
    native = warm.d.module(guard_pin, "matched_build_guard")
    floor, cap = producer.FLOOR, producer.CAP
    guard = native.Guard(native.GuardConfig(root, started, producer.WALL, producer.RSS, cap, floor))
    require(guard._root_id == (info.st_dev, info.st_ino), "Output root changed since creation")
    identity = {"producer": warm.pin(CALLER, hashlib.sha256(producer.source_bytes(CALLER)).hexdigest()),
                "reused_producer": warm.pin(PRODUCER, PRODUCER_SHA),
                "binding": warm.pin(producer.BINDING, producer.BINDING_SHA), "guard": guard_pin,
                "argv": sys.orig_argv, "repository": str(producer.REPO), "commits": commits}
    with guard.shutdown_signals():
        available = guard._space(walk=False)["free_bytes"]
        admitted = available >= floor + cap
        if not admitted:
            print(f"Available bytes {available} cannot fund the continuing floor and build output cap", file=sys.stderr)
        admission = write_json(root / "admission.json", {
            **identity, "scope": "matched normal installs and runtime metadata only", "admitted": admitted,
            "issued_utc": started.isoformat(), "observed_utc": datetime.now(timezone.utc).isoformat(),
            "root": {"path": str(root), "dev": info.st_dev, "ino": info.st_ino}, "supervisor_pid": os.getpid(),
            "available_bytes": available, "minimum_free_bytes": floor, "owned_output_cap_bytes": cap,
            "sampled_group_rss_limit_bytes": producer.RSS, "wall_seconds": producer.WALL,
        }, guard)
        if not admitted:
            return 1
        result = {**identity, "status": "failed", "admission": admission, "phases": [], "roles": {}, "validated_commits": {}}
        result_path = root / "operator-result.json"
        try:
            repo = physical(producer.REPO)
            for role, commit in commits.items():
                guard._stop_check()
                observed = binding.git(repo, "rev-parse", "--verify", f"{commit}^{{commit}}")
                require(observed == commit, f"{role} must identify the actual Git commit object: {commit}")
                result["validated_commits"][role] = observed
            guard._stop_check()
            require(not os.path.lexists(physical(binding.COMMON) / "MungeSumstats"), "Common library contains a selected package")
            for name in ("home", "tmp"):
                (root / name).mkdir(mode=0o700)
            operator = producer.BuildRun(binding, warm, guard, result)
            toolchain, runtime, transcripts = operator.captures()
            baseline, candidate = (operator.role(role, commit, toolchain) for role, commit in commits.items())
            require(baseline["source_tree"] != candidate["source_tree"]
                    and baseline["compile_and_link_argv"] == candidate["compile_and_link_argv"],
                    "Expected distinct source trees with identical native commands")
            result["provenance"] = write_json(root / "PROVENANCE.json", {
                "repository": str(repo), **result["roles"], "runtime_receipt": runtime, "runtime_transcripts": transcripts,
                "minimum_free_bytes": floor, "estimated_output_bytes": producer.ESTIMATED_OUTPUT,
                "contention_note": (f"Independent local build and metadata admission observed {available} free bytes. "
                                    "The build phases ran sequentially under one guard. Benchmark admission, host "
                                    "contention and cache observations, and timings are pending."),
            }, guard)
            write_json(result_path, {**result, "status": "completed",
                                     "finished_utc": datetime.now(timezone.utc).isoformat()}, guard)
            return 0
        except (OSError, ValueError, RuntimeError, KeyError, TypeError, csv.Error, subprocess.SubprocessError,
                KeyboardInterrupt) as error:
            result.update(error=f"{type(error).__name__}: {error}", finished_utc=datetime.now(timezone.utc).isoformat())
            print(result["error"], file=sys.stderr)
            try:
                write_json(result_path, result, guard, check_stop=False)
            except (OSError, ValueError, RuntimeError) as receipt_error:
                print(f"Could not retain operator result: {receipt_error}", file=sys.stderr)
            return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_commit", metavar="BASELINE_COMMIT", type=commit_id)
    parser.add_argument("candidate_commit", metavar="CANDIDATE_COMMIT", type=commit_id)
    parser.add_argument("output", metavar="FRESH_OUTPUT", type=Path, help="absent physical directory for private build evidence")
    args = parser.parse_args()
    if args.baseline_commit == args.candidate_commit:
        parser.error("baseline and candidate commits must be distinct")
    try:
        if os.path.lexists(args.output):
            raise FileExistsError(f"Output already exists: {args.output}")
        return run(args.output, {"baseline": args.baseline_commit, "candidate": args.candidate_commit})
    except (OSError, ValueError, RuntimeError) as error:
        print(error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
