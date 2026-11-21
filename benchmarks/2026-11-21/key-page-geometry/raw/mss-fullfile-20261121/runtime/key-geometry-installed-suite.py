"""Qualify the fixed Darwin key-page build with one unfiltered installed suite."""
from collections import Counter
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import stat
import sys
import tarfile
from types import ModuleType


RUNTIME = Path("/private/tmp/mss-fullfile-20261121/runtime")
SOURCE = RUNTIME / "key-geometry-installed-suite.py"
ROOT = Path("/private/tmp/mss-key-geometry-installed-suite-20261121")
PYTHON = Path("/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15")
HISTORY = Path("/private/tmp/mss-performance-20261119")
RHOME = Path("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources")
PRIOR = Path("/private/tmp/mss-namespace-installed-suite-20261121")
BUILD = Path("/private/tmp/mss-key-geometry-builds-20261121/candidate")
LIBRARY = BUILD / "library"
PACKAGE = LIBRARY / "MungeSumstats"
REFERENCE = Path("/private/tmp/mss-reference-closure-20261121/attempt-02")
UNION = REFERENCE / "union"
CONTROL = HISTORY / "tools/darwin-functional-stage178/control.R"
QUALIFY = HISTORY / "repo/benchmarks/2026-11-19/performance-upgrade/operations/qualify.R"
NATIVE = HISTORY / "tools/darwin-full-stage185/guard.py"
PRODUCER = RUNTIME / "namespace-builds.py"
COMMIT = "8b85aa5ffc45ad7632cfc5319e85dd5eb3089be6"
ADAPTED_TESTS = {"test-check_impute_se_beta.R", "test-read_vcf.R"}
LABEL = "candidate-unfiltered"
MIB = 1024**2
LIMITS = dict(wall_seconds=3600, rss_bytes=8 * 1024**3,
              disk_bytes=2 * 1024**3, floor_bytes=186 * 1024**3)
PINS = {
    PYTHON: "6cc907e4a022fefc843bcd93954eaf5a2cbfc77a095640bfadeeed96154190c5",
    NATIVE: "58e0818aa20801c30194bdce26bd117112b4b4d5d10f74ea556360c0867b56a2",
    PRODUCER: "615535da8374c81d3d7f23f0300c313d26d06ebea638e73239492243f851a65f",
    CONTROL: "1cc32ff00c4d8455a6d1ec7d05316f11ef0a4d6822607432effcd8c338de3c39",
    QUALIFY: "087e55847f6a5333781f7526625ea5e59776cfd81f18f83a3c7492d55cd712b9",
    BUILD / "build.json": "4f4cc66733417e431281dc1892d71f5acc3fe155b2000d2404814f2c5793b74a",
    REFERENCE / "union-inventory.json": "d30550f4c01897764c858eefe03a1d3f037d3cdf0ee90be1224a419c7a2a3213",
    PRIOR / "admission.json": "cbab726393fca23c7a53265bd07383d263798774848d0791b85da0034efe3ed2",
    PRIOR / "RESULT.json": "383561ec242ddbdb66d619b9ce3b826f46c783fbfb40c31de226fcf108284090",
    PRIOR / "deferred-unfiltered.jsonl": "91c9e317497311468927b73f354b113006edb44d9a9673bb9e4dba63b8275def",
    PRIOR / "receipts/deferred/tests.tsv": "f1f878ae3c82dd9a43a538269bb76204ede85f6c70b04e76f4cdecdc8a24c333",
}


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read_file(path, expected=None, *, limit=MIB):
    require(path.resolve(strict=True) == path, f"Expected physical input {path}")
    flags = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC
    with os.fdopen(os.open(path, flags), "rb") as stream:
        before = os.fstat(stream.fileno())
        require(stat.S_ISREG(before.st_mode) and before.st_size <= limit, f"Invalid input {path}")
        data = stream.read(limit + 1)
        after = os.fstat(stream.fileno())
    require(len(data) == before.st_size == after.st_size
            and before.st_mtime_ns == after.st_mtime_ns and before.st_ctime_ns == after.st_ctime_ns,
            f"Input changed while reading {path}")
    require(expected is None or sha(data) == expected, f"Pin differs for {path}")
    return data


def load_pinned(path, name):
    data = read_file(path, PINS[path])
    module = ModuleType(name)
    module.__file__ = str(path)
    sys.modules[name] = module
    exec(compile(data, str(path), "exec"), module.__dict__)
    return module


def source_inventory(build, tests):
    path = BUILD / "MungeSumstats.tar"
    require(build["archive"]["path"] == str(path), "Build expects a different staged archive")
    contents = read_file(path, build["archive"]["sha256"], limit=32 * MIB)
    product, archive_tests = {}, {}
    with tarfile.open(fileobj=io.BytesIO(contents), mode="r:") as archive:
        require(archive.pax_headers.get("comment") == COMMIT, "Archived source revision differs")
        for member in archive:
            name = PurePosixPath(member.name).relative_to("MungeSumstats")
            require(".." not in name.parts and (member.isdir() or member.isfile()),
                    f"Unexpected staged archive member: {member.name}")
            if member.isdir():
                continue
            if name.parts[0] in {"R", "src", "DESCRIPTION", "NAMESPACE"}:
                table, key = product, name.as_posix()
            elif name.is_relative_to("tests/testthat"):
                table, key = archive_tests, name.relative_to("tests/testthat").as_posix()
            else:
                continue
            require(key not in table, f"Duplicate staged archive member: {member.name}")
            table[key] = sha(archive.extractfile(member).read())
    require(len(product) == 136 and archive_tests.keys() == tests.keys(), "Staged source census differs")
    require({key for key in tests if tests[key] != archive_tests[key]} == ADAPTED_TESTS,
            "Staged test digests differ outside the reviewed ALS adaptations")
    return {"commit": COMMIT, "tree": build["source_tree"], "archive": build["archive"],
            "product_files": product, "archive_test_sources": archive_tests}


def installed_inventory(build, check_stop):
    check_stop()
    require(LIBRARY.resolve(strict=True) == LIBRARY and list(LIBRARY.iterdir()) == [PACKAGE]
            and PACKAGE.resolve(strict=True) == PACKAGE, "Expected a one-package candidate library")
    files, identity, total = {}, hashlib.sha256(), 0
    paths = []
    for path in PACKAGE.rglob("*"):
        check_stop()
        paths.append(path)
    for path in sorted(paths):
        check_stop()
        info = path.lstat()
        if stat.S_ISDIR(info.st_mode):
            continue
        total += info.st_size
        require(len(files) < 1000 and total <= 32 * MIB, "Candidate inventory exceeds input bounds")
        data = read_file(path, limit=32 * MIB)
        name, digest = path.relative_to(PACKAGE).as_posix(), sha(data)
        files[name] = {"bytes": len(data), "sha256": digest}
        identity.update(name.encode("utf-8") + b"\0" + digest.encode("ascii") + b"\n")
    check_stop()
    require(identity.hexdigest() == build["package_tree_sha256"], "Installed candidate tree differs from build")
    dll = "libs/MungeSumstats.so"
    require(build["dll"]["path"] == str(PACKAGE / dll)
            and files[dll]["sha256"] == build["dll"]["sha256"] == build["dll_sha256"],
            "Candidate DLL differs from build")
    return {"path": str(PACKAGE), "tree_sha256": identity.hexdigest(), "files": files}


def test_sources(directory, pins, check_stop):
    check_stop()
    require(directory.resolve(strict=True) == directory
            and {path.name for path in directory.iterdir()} == pins.keys(), f"Test census differs in {directory}")
    sources = {}
    for name, expected in sorted(pins.items()):
        check_stop()
        sources[name] = read_file(directory / name, expected)
    check_stop()
    return sources


def verify_union(manifest, check_stop):
    check_stop()
    targets = manifest["physical_targets"]
    require(manifest["overlay"] == str(UNION) and UNION.resolve(strict=True) == UNION
            and len(targets) == manifest["entry_count"] == 102
            and "MungeSumstats" not in targets and {path.name for path in UNION.iterdir()} == targets.keys(),
            "Reference union census differs")
    for name, target in targets.items():
        check_stop()
        entry = UNION / name
        require(entry.is_symlink() and str(entry.resolve(strict=True)) == target
                and (entry / "DESCRIPTION").is_file(), f"Reference union target differs for {name}")


def environment(prior_plan):
    data = {
        "PATH": "/usr/bin:/bin:/usr/sbin:/sbin", "HOME": str(ROOT / "home"),
        **dict.fromkeys(("TMPDIR", "TMP", "TEMP"), str(ROOT / "tmp")),
        "LANG": "en_US.UTF-8", "LC_ALL": "en_US.UTF-8", "TZ": "UTC",
        "R_LIBS": f"{LIBRARY}:{UNION}",
        **dict.fromkeys(("R_LIBS_SITE", "R_LIBS_USER", "R_ENVIRON_USER", "R_PROFILE_USER"), "/dev/null"),
        "OMP_DYNAMIC": "FALSE", "MKL_DYNAMIC": "FALSE",
        **dict.fromkeys(("BLIS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "OMP_NUM_THREADS",
                         "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS", "R_DATATABLE_NUM_THREADS",
                         "VECLIB_MAXIMUM_THREADS"), "1"),
    }
    require(sorted(data) == prior_plan["environment_keys"], "Environment keys differ from reviewed suite")
    return data


def table(data, names, check_stop):
    check_stop()
    rows = csv.DictReader(io.StringIO(data.decode()), delimiter="\t", strict=True)
    require(rows.fieldnames == names, f"Unexpected receipt header {rows.fieldnames}")
    for row in rows:
        check_stop()
        require(None not in row and None not in row.values(), "Malformed receipt row")
        yield row


def suite_counts(data, check_stop):
    columns = "file test nb passed skipped failed error warning".split()
    cases, summary = Counter(), dict(tests=0, passed=0, skipped=0, failed=0, errors=0, test_warnings=0)
    for row in table(data, columns, check_stop):
        nb, passed, failed, warning = (int(row[key]) for key in ("nb", "passed", "failed", "warning"))
        require(min(nb, passed, failed, warning) >= 0 and nb >= passed + failed + warning
                and row["skipped"] in ("TRUE", "FALSE") and row["error"] in ("TRUE", "FALSE"),
                "Invalid test counts")
        cases[row["file"], row["test"], passed] += 1
        summary["tests"] += 1
        summary["passed"] += passed
        summary["skipped"] += row["skipped"] == "TRUE"
        summary["failed"] += failed
        summary["errors"] += row["error"] == "TRUE"
        summary["test_warnings"] += warning
    summary["active_files"] = sorted({file for file, _, _ in cases})
    return cases, summary


def collect_suite(directory, check_stop):
    paths = {name: directory / f"{name}.tsv" for name in ("tests", "conditions", "qualification")}
    contents = {}
    for name, path in paths.items():
        check_stop()
        contents[name] = read_file(path)
    cases, summary = suite_counts(contents["tests"], check_stop)
    events = table(contents["conditions"], ["file", "test", "type", "message"], check_stop)
    condition_counts = Counter(row["type"] for row in events)
    summary["suite_warnings"] = condition_counts["suite_warning"]
    summary["errors"] += condition_counts["suite_error"]
    columns = "qualified tests passed skipped failed errors test_warnings suite_warnings".split()
    qualified, = table(contents["qualification"], columns, check_stop)
    require(all(qualified[key] == str(summary[key]) for key in columns[1:])
            and condition_counts["expectation_warning"] == summary["test_warnings"], "Qualification receipts disagree")
    pins = {name: {"path": str(path), "sha256": sha(contents[name])} for name, path in paths.items()}
    return cases, summary, qualified["qualified"] == "TRUE", pins


def guard_receipt(receipt, argv, env, paths, limits):
    planned, running, terminal = [json.loads(line) for line in receipt.splitlines()]
    require(all(row["label"] == LABEL for row in (planned, running, terminal))
            and planned["status"] == "planned" and planned["argv"] == argv and planned["cwd"] == str(ROOT)
            and planned["limits"] == limits and planned["allowed_exits"] == [0]
            and all(planned[key] == str(value) for key, value in paths.items())
            and running["status"] == "running" and running["environment_keys"] == sorted(env)
            and running["pid"] == running["expected_pgid"] == running["expected_sid"] > 1
            and running["uid"] == os.geteuid() and planned["utc"] <= running["utc"] <= terminal["utc"]
            and terminal["status"] == "completed" and terminal["closure"] == "confirmed"
            and terminal["raw_exit"] == 0 and terminal["leader_reaped"]
            and terminal["primary"] is None and not terminal["cleanup_errors"]
            and terminal.get("receipt_error") is None and terminal["missed_samples"] == 0, "Guard receipt is incomplete")
    return terminal


def run():
    start = datetime.now(timezone.utc)
    require(len(sys.argv) == 1 and sys.platform == "darwin" and Path(sys.executable) == PYTHON
            and sys.flags.isolated and sys.flags.dont_write_bytecode and sys.flags.no_site
            and Path(__file__) == SOURCE, "Use the fixed physical launch command")
    native = load_pinned(NATIVE, "key_geometry_suite_guard")
    write_json = load_pinned(PRODUCER, "key_geometry_suite_producer").write_json
    os.umask(0o077)
    ROOT.mkdir(mode=0o700)
    root_info = ROOT.lstat()
    guard = native.Guard(native.GuardConfig(ROOT, start, **LIMITS))
    check_stop = guard._stop_check
    require(guard._root_id == (root_info.st_dev, root_info.st_ino), "Output root identity changed")
    result = {"scope": "one unfiltered installed Darwin testthat suite", "status": "failed",
              "scientific_completion": False, "performance_claim": False}
    with guard.shutdown_signals():
        try:
            free = guard._space(walk=False)["free_bytes"]
            require(free >= LIMITS["floor_bytes"] + LIMITS["disk_bytes"],
                    f"Only {free} free bytes at admission")
            inputs = {path: read_file(path, digest) for path, digest in PINS.items()}
            prior_admission = json.loads(inputs[PRIOR / "admission.json"])
            prior_result = json.loads(inputs[PRIOR / "RESULT.json"])
            prior_records = [json.loads(line) for line in inputs[PRIOR / "deferred-unfiltered.jsonl"].splitlines()]
            prior_counts, previous = suite_counts(inputs[PRIOR / "receipts/deferred/tests.tsv"], check_stop)
            pins, census = prior_admission["test_sources"], prior_admission["test_census"]
            require(len(pins) == 92 and len(census["test_files"]) == 90 and len(previous["active_files"]) == 89
                    and previous["tests"] == 269 and previous["passed"] == 3077
                    and sorted(set(census["test_files"]) - set(census["comment_only"])) == previous["active_files"]
                    and prior_result["phases"]["deferred-unfiltered"]["raw_exit"] == prior_records[-1]["raw_exit"] == 0,
                    "Reviewed suite census differs")
            build = json.loads(inputs[BUILD / "build.json"])
            require(build["library"] == str(LIBRARY) and build["source_commit"] == COMMIT
                    and build["terminal"]["status"] == "completed" and build["terminal"]["closure"] == "confirmed"
                    and build["terminal"]["raw_exit"] == 0 and build["terminal"]["leader_reaped"],
                    "Build identity or installation outcome differs")
            manifest = json.loads(inputs[REFERENCE / "union-inventory.json"])
            verify_union(manifest, check_stop)
            runner = RHOME / "bin/Rscript"
            r_pins = {path: manifest["input_pins"][str(path)] for path in
                      (runner, RHOME / "bin/exec/R", RHOME / "lib/libR.dylib")}
            for path, expected in r_pins.items():
                read_file(path, expected, limit=16 * MIB)
            staged = test_sources(PRIOR / "tests", pins, check_stop)
            for name in ("home", "tmp", "tests"):
                (ROOT / name).mkdir(mode=0o700)
            guard._check_space(guard._space(walk=True), evidence_bytes=sum(map(len, staged.values())),
                               metadata_blocks=len(staged), new_entries=len(staged))
            for name, data in staged.items():
                guard._stop_check()
                with os.fdopen(guard._new_file(ROOT / "tests" / name), "wb") as stream:
                    stream.write(data)

            def integrity():
                tested = test_sources(ROOT / "tests", pins, check_stop)
                return {"installed": installed_inventory(build, check_stop),
                        "tests": {name: sha(data) for name, data in tested.items()}}

            before = integrity()
            result["before"] = write_json(ROOT / "before.json", before, guard)
            paths = {"receipt": ROOT / f"{LABEL}.jsonl", "stdout": ROOT / f"{LABEL}.stdout",
                     "stderr": ROOT / f"{LABEL}.stderr"}
            suite = ROOT / "suite"
            argv = list(map(str, (runner, "--vanilla", CONTROL, "suite", LIBRARY, UNION, QUALIFY, ROOT / "tests", suite)))
            env = environment(prior_records[1])
            limits = {key: val for key, val in asdict(guard.config).items() if key not in {"root", "admitted_utc"}}
            result["admission"] = write_json(ROOT / "admission.json", {
                "issued_utc": start.isoformat(), "supervisor_pid": os.getpid(), "root": str(ROOT),
                "root_dev_ino": guard._root_id, "available_bytes": free, "limits": limits,
                "monotonic_deadline": guard.deadline, "operator": {"path": str(SOURCE), "sha256": sha(read_file(SOURCE))},
                "controls_sha256": {str(path): digest for path, digest in (PINS | r_pins).items()},
                "source": source_inventory(build, pins), "build": build, "prior_census": census,
                "prior_case_totals": previous, "reference_union": str(UNION), "reference_packages": manifest["entry_count"],
                "cwd": str(ROOT), "environment": env, "argv": argv, "allowed_exits": [0], "before": result["before"],
                "phase_paths": {key: str(value) for key, value in paths.items()},
            }, guard)
            outcome = guard.run_phase(LABEL, argv, cwd=ROOT, env=env, allowed_exits=(0,), **paths)
            result["outcome"] = asdict(outcome)
            result["post_verification"] = "incomplete" if outcome.closure == "confirmed" else "not_started"
            check_stop()
            if paths["receipt"].is_file():
                receipt = read_file(paths["receipt"])
                result["guard_receipt"] = {"path": str(paths["receipt"]), "sha256": sha(receipt)}
            # Preserve the post-run fence and any usable test reports after a closed failure.
            if outcome.closure == "confirmed":
                after = integrity()
                check_stop()
                result["after"] = write_json(ROOT / "after.json", after, guard)
                require(before == after, "Candidate or staged tests changed during the suite")
                require(test_sources(PRIOR / "tests", pins, check_stop) == staged, "Prior staged tests changed")
                read_file(REFERENCE / "union-inventory.json", PINS[REFERENCE / "union-inventory.json"])
                verify_union(manifest, check_stop)
                for path in (CONTROL, QUALIFY, BUILD / "build.json", *r_pins):
                    check_stop()
                    read_file(path, (PINS | r_pins)[path], limit=16 * MIB)
                if suite.is_dir():
                    cases, summary, qualified, result["suite_receipts"] = collect_suite(suite, check_stop)
                    result["observed"] = summary
                check_stop()
                result["post_verification"] = "complete"
            require(outcome.status == "completed" and outcome.closure == "confirmed" and outcome.raw_exit == 0
                    and outcome.primary is None and not outcome.cleanup_errors and outcome.receipt_error is None,
                    f"Guard did not complete successfully: {outcome}")
            require("guard_receipt" in result, "Successful phase has no guard receipt")
            result["terminal"] = guard_receipt(receipt, argv, env, paths, limits)
            require("observed" in result and qualified and cases == prior_counts
                    and summary["active_files"] == previous["active_files"]
                    and not any(summary[key] for key in ("failed", "skipped", "errors")),
                    "Installed suite differs from the full reviewed case census or reported failures")
            write_json(ROOT / "RESULT.json", {
                **result, "status": "completed", "scientific_completion": True,
                "finished_utc": datetime.now(timezone.utc).isoformat(),
            }, guard)
            return 0
        except (OSError, ValueError, RuntimeError, KeyError, TypeError, csv.Error, tarfile.TarError, KeyboardInterrupt) as error:
            result.update(error=f"{type(error).__name__}: {error}", finished_utc=datetime.now(timezone.utc).isoformat())
            print(result["error"], file=sys.stderr)
            try:
                write_json(ROOT / "RESULT.json", result, guard, check_stop=False)
            except (OSError, ValueError, RuntimeError) as receipt_error:
                print(f"Could not retain suite result: {receipt_error}", file=sys.stderr)
            return 1


if __name__ == "__main__":
    try:
        sys.exit(run())
    except (OSError, ValueError, RuntimeError) as error:
        print(error, file=sys.stderr)
        sys.exit(1)
