#!/usr/bin/env python3
"""Run a pinned baseline/candidate comparison in serial fresh R processes."""

import argparse
import csv
from datetime import datetime, timezone
from decimal import InvalidOperation
import hashlib
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import sys
import time

from process_metrics import r_command, run_process
from verify import (
    compare_bytes, digest, read_json, tsv_pairs, validate_hard_case_oracle, verify_editions,
    write_json,
)


HERE = Path(__file__).resolve().parent
BASELINE = "a86e8c94d753630eff8c6a695ed2f240d1123d0d"
ROLES = ("baseline", "candidate")
REQUIRED_REFERENCE_FILES = {
    "manifest.rds", "core-version", "fasta.fai", "keys", "offsets", "records",
    "rejected.tsv", "source.tsv",
}
REQUIRED_SUPPLEMENTS = {"distinct_shuffled", "coordinate_only", "id_only",
                        "negative", "dense", "mixed", "hard_cases"}
SETTINGS = {
    "batch_size": 100000, "batch_bytes": 8388608, "max_record_bytes": 67108864,
    "candidate_limit": 10000, "frequency_type": "unknown", "effect_type": "beta",
    "delimiter": "tab", "verify_acquisition_receipt": False,
}
RUNTIME_LIBRARIES = {"r_shared": "R_shared", "blas": "BLAS_library", "lapack": "LAPACK_library"}
RUNTIME_FIELDS = (
    "R_version", "R_platform", "R_home", "R_bin", "R_engine", "R_shared", "R_shared_resolution",
    "locale", "BLAS", "BLAS_library", "LAPACK", "LAPACK_library", "threads",
)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def valid_hash(text, size=64):
    return isinstance(text, str) and re.fullmatch(f"[a-f0-9]{{{size}}}", text) is not None


def integer(value, minimum, maximum, name):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f"Invalid {name}")


def directory(value):
    path = Path(value).resolve(strict=True)
    if not path.is_dir():
        raise ValueError(f"Expected directory: {path}")
    return path


def plain_name(value):
    return isinstance(value, str) and re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9._-]*", value)


def package_tree(path):
    if not (path / "DESCRIPTION").is_file():
        raise ValueError(f"Not an installed package directory: {path}")
    files = {}
    identity = hashlib.sha256()
    for item in sorted(path.rglob("*")):
        if item.is_symlink():
            raise ValueError(f"Installed package contains an unpinned symlink: {item}")
        if item.is_file():
            name = item.relative_to(path).as_posix()
            entry = digest(item)
            files[name] = entry
            identity.update(name.encode("utf-8") + b"\0" + entry["sha256"].encode("ascii") + b"\n")
    return {"tree_sha256": identity.hexdigest(), "files": files}


def resolve_configuration(config, suite):
    if config["schema"] != 1 or not config["contention_note"].strip():
        raise ValueError("Expected schema 1 and a recorded contention/admission condition")
    hashes = {}

    def register(pin):
        path = Path(pin["path"]).resolve(strict=True)
        sha256 = pin["sha256"]
        if not path.is_file() or not valid_hash(sha256) or hashes.get(path, sha256) != sha256:
            raise ValueError(f"Invalid or conflicting asset pin: {path}")
        hashes[path] = sha256
        pin["path"] = str(path)
        return path

    register(config["r_executable"])
    register(config["r_engine"])
    if not isinstance(config["runtime_files"], list):
        raise ValueError("runtime_files must be an array of role-labelled runtime pins")
    runtime_roles = set()
    for pin in config["runtime_files"]:
        role = pin["role"]
        if role not in RUNTIME_LIBRARIES or role in runtime_roles:
            raise ValueError(f"Unknown or repeated runtime file role: {role}")
        runtime_roles.add(role)
        register(pin)
    config["dependency_libraries"] = [str(directory(path)) for path in config["dependency_libraries"]]
    integer(config["minimum_free_bytes"], 0, 2**63 - 1, "minimum_free_bytes")
    integer(config["estimated_output_bytes"], 1, 2**63 - 1, "estimated_output_bytes")
    packages = {}
    receipts = {}
    for role in ROLES:
        item = config[role]
        if not valid_hash(item["commit"], 40) or not valid_hash(item["package_tree_sha256"]):
            raise ValueError(f"Invalid source or package tree pin: {role}")
        library = directory(item["library"])
        item["library"] = str(library)
        package = directory(library / "MungeSumstats")
        dll = register(item["dll"])
        if package not in dll.parents:
            raise ValueError(f"DLL is outside selected package: {role}")
        receipts[role] = register(item["build_receipt"])
        packages[role] = package
    if config["baseline"]["commit"] != BASELINE or config["candidate"]["commit"] == BASELINE:
        raise ValueError("Expected the exact task baseline and a distinct candidate source commit")
    if packages["baseline"] == packages["candidate"]:
        raise ValueError("Provide independent installed package trees")

    references = config["references"]
    for reference in references.values():
        refdir = directory(reference["directory"])
        reference["directory"] = str(refdir)
        if reference["format"] not in (2, 3):
            raise ValueError("Unsupported reference format")
        if not isinstance(reference["bind_fasta"], bool):
            raise ValueError("bind_fasta must be boolean")
        register(reference["fasta"])
        names = set(reference["files"])
        if not REQUIRED_REFERENCE_FILES.issubset(names):
            raise ValueError("Incomplete reference inventory")
        for name, sha256 in reference["files"].items():
            if not plain_name(name):
                raise ValueError(f"Invalid reference inventory entry: {name}")
            if register({"path": str(refdir / name), "sha256": sha256}) != refdir / name:
                raise ValueError(f"Reference artifact is an indirect path: {name}")

    ids = set()
    scales = {}
    by_format = {2: set(), 3: set()}
    expanded_cases = []
    for case in config["cases"]:
        if "reference_formats" in case:
            if suite == "legacy" or "reference" in case or case.get("target"):
                raise ValueError("reference_formats is for supplemental preservation cases")
            for refname in case["reference_formats"]:
                suffix = references[refname]["format"]
                replica = {key: value for key, value in case.items() if key != "reference_formats"}
                expanded_cases.append({**replica, "id": f'{case["id"]}-f{suffix}', "reference": refname})
        else:
            expanded_cases.append(case)
    config["cases"] = expanded_cases
    for case in config["cases"]:
        name = case["id"]
        if not plain_name(name) or name.casefold() in ids:
            raise ValueError(f"Invalid or repeated case identifier: {name}")
        ids.add(name.casefold())
        integer(case["rows"], 1, 2**53 - 1, f"{name} rows")
        integer(case["pairs"], 1, 3, f"{name} pairs")
        register(case["input"])
        register(case["cardinality_receipt"])
        if "oracle" in case:
            register(case["oracle"])
        if suite == "legacy":
            register(case["legacy_panel"])
            continue
        reference = references[case["reference"]]
        settings = case.get("settings", {})
        if set(settings) - set(SETTINGS) - {"build", "columns", "source_strand"}:
            raise ValueError(f"Unknown public call setting: {name}")
        case["settings"] = options = {**SETTINGS, **settings}
        for option in ("batch_size", "batch_bytes", "max_record_bytes", "candidate_limit"):
            integer(options[option], 1, 2**31 - 1, option)
        if (options["delimiter"] not in ("tab", "comma") or
                options["frequency_type"] not in ("unknown", "EAF", "MAF") or
                options["effect_type"] not in ("beta", "log_odds", "none") or
                not isinstance(options["verify_acquisition_receipt"], bool)):
            raise ValueError(f"Invalid public call setting: {name}")
        columns = options.get("columns", {})
        if not isinstance(columns, dict) or any(not plain_name(key) or not isinstance(val, str)
                                                 for key, val in columns.items()):
            raise ValueError(f"Invalid source column mappings: {name}")
        if ((options.get("source_strand") is not None or
             {"source_ref", "source_alt", "source_strand"} & columns.keys()) and
                not case.get("source_evidence", "").strip()):
            raise ValueError(f"Explicit source declarations require source evidence: {name}")
        if case.get("target"):
            if case["rows"] in scales or reference["format"] != 2:
                raise ValueError("Historical scale targets must be unique and use reference format 2")
            scales[case["rows"]] = case["pairs"]
        tags = set(case.get("tags", []))
        if not case.get("target"):
            by_format[reference["format"]].update(tags)
        if "hard_cases" in tags:
            if "oracle" not in case:
                raise ValueError(f"Hard-case workload requires an independent complete-row oracle: {name}")
            validate_hard_case_oracle(Path(case["oracle"]["path"]))
    if not ids:
        raise ValueError("At least one case is required")
    if suite == "matrix":
        if scales != {1_000_000: 3, 10_000_000: 3, 40_000_000: 1}:
            raise ValueError("Expected historical 1M/10M/40M cases with three/three/one pairs")
        for format_number, tags in by_format.items():
            if REQUIRED_SUPPLEMENTS - tags:
                raise ValueError(f"Missing supplemental workloads in format {format_number}")
    return hashes, packages, receipts


def asset_inventory(hashes):
    assets = {}
    for path, expected in hashes.items():
        actual = digest(path)
        if actual["sha256"] != expected:
            raise ValueError(f"Asset SHA256 mismatch: {path}")
        assets[str(path)] = actual
    return assets


def read_provenance(config, packages, receipts, assets):
    installed = {}
    builds = {}
    for role in ROLES:
        expected = config[role]
        installed[role] = package_tree(packages[role])
        if installed[role]["tree_sha256"] != expected["package_tree_sha256"]:
            raise ValueError(f"Installed package tree SHA256 mismatch: {role}")
        build = read_json(receipts[role])
        if (build["source_commit"] != expected["commit"] or
                build["dll_sha256"] != assets[expected["dll"]["path"]]["sha256"] or
                build["package_tree_sha256"] != installed[role]["tree_sha256"]):
            raise ValueError(f"Build receipt does not bind source, DLL and installation: {role}")
        if not {"compiler", "compiler_flags", "R_version", "R_platform"}.issubset(build["toolchain"]):
            raise ValueError(f"Incomplete declared build toolchain: {role}")
        builds[role] = build
    if builds["baseline"]["toolchain"] != builds["candidate"]["toolchain"]:
        raise ValueError("Build toolchain declarations differ")
    return installed, builds["baseline"]["toolchain"]


def verify_runtime(config, metrics, toolchain):
    observed = {name: metrics[name] for name in RUNTIME_FIELDS}
    if (observed["threads"] != "1" or observed["R_version"] != toolchain["R_version"] or
            observed["R_platform"] != toolchain["R_platform"]):
        raise ValueError("Observed runtime disagrees with build provenance")
    engine = Path(observed["R_engine"]).resolve(strict=True)
    r_home = directory(observed["R_home"])
    r_bin = directory(observed["R_bin"])
    if engine != Path(config["r_engine"]["path"]) or r_bin not in engine.parents or r_home not in r_bin.parents:
        raise ValueError("Observed R engine disagrees with the engine pin or R installation")
    if sys.platform == "win32" and engine != Path(config["r_executable"]["path"]):
        raise ValueError("The observed Windows R engine must be the pinned Rterm launcher")
    libraries = {pin["role"]: Path(pin["path"]) for pin in config["runtime_files"]}
    for role, field in RUNTIME_LIBRARIES.items():
        path = observed[field]
        actual = None if path == "unavailable" else Path(path).resolve(strict=True)
        if actual != libraries.get(role):
            raise ValueError(f"Observed R runtime library disagrees with its role pin: {role}")
    return observed


def cardinality(case):
    receipt = read_json(Path(case["cardinality_receipt"]["path"]))
    if receipt["input_sha256"] != case["input"]["sha256"] or receipt["input_rows"] != case["rows"]:
        raise ValueError(f"Cardinality receipt does not bind the input: {case['id']}")
    integer(receipt["distinct_keys"], 1, case["rows"], "distinct_keys")
    if not receipt["key_definition"].strip() or not receipt["method"].strip():
        raise ValueError("Cardinality receipt must define the counted key and provenance")
    if "distinct_shuffled" in case.get("tags", ()):
        if receipt["distinct_keys"] != case["rows"]:
            raise ValueError(f"Expected one distinct key per source row: {case['id']}")
        derivation = receipt.get("shuffle_derivation")
        if not isinstance(derivation, str) or not derivation.strip():
            raise ValueError(f"Distinct workload requires shuffle derivation or seed: {case['id']}")
    return receipt


def write_settings(path, settings):
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", quoting=csv.QUOTE_ALL, lineterminator="\n")
        writer.writerow(("key", "value"))
        writer.writerows(settings.items())


def run_one(config, case, run_dir, role, *, suite, profile, threshold, assets):
    settings = {
        "mode": suite, "profile": profile, "allocation_threshold": threshold,
        "input": case["input"]["path"], "rows": case["rows"], "output": str(run_dir),
        "dll": config[role]["dll"]["path"], "library.0": config[role]["library"],
    }
    settings["dll_md5"] = assets[settings["dll"]]["md5"]
    for index, path in enumerate(config["dependency_libraries"], start=1):
        settings[f"library.{index}"] = path
    if suite == "legacy":
        settings["legacy_panel"] = case["legacy_panel"]["path"]
    else:
        reference = config["references"][case["reference"]]
        options = case["settings"]
        settings.update(reference=reference["directory"], core_version=f'mss-preservation-{reference["format"]}')
        settings["fasta_pin_path"] = reference["fasta"]["path"]
        if reference["bind_fasta"]:
            settings["fasta"] = reference["fasta"]["path"]
        settings.update({key: value for key, value in options.items()
                         if value is not None and key not in {"columns", "verify_acquisition_receipt"}})
        settings.update({f"column.{key}": value for key, value in options.get("columns", {}).items()})
        if options["verify_acquisition_receipt"]:
            raw = assets[case["input"]["path"]]
            settings.update(acquisition_md5=raw["md5"], acquisition_bytes=raw["bytes"])
    run_dir.mkdir()
    ini = run_dir / "settings.tsv"
    write_settings(ini, settings)
    command = r_command(config["r_executable"]["path"], HERE / "worker.R", "run", str(ini))
    start = utc_now()
    process = run_process(command, run_dir)
    process.update(started_utc=start, finished_utc=utc_now(), role=role, case=case["id"])
    write_json(run_dir / "process.json", process)
    if process["exit_code"] != 0:
        raise RuntimeError(f"R worker failed with exit {process['exit_code']}: {run_dir}")
    metrics = tsv_pairs(run_dir / "metrics.tsv", ("key", "value"))
    return process, metrics


def dependency_state(run_dir, config, role):
    with (run_dir / "namespaces.tsv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if reader.fieldnames != ["package", "version", "path"]:
            raise ValueError("Unexpected namespace inventory")
        entries = {}
        for row in reader:
            package = row["package"]
            if package in entries:
                raise ValueError("Duplicate package in namespace inventory")
            location = str(Path(row["path"]).resolve(strict=True))
            if package == "MungeSumstats":
                expected = str((Path(config[role]["library"]) / package).resolve(strict=True))
                if location != expected:
                    raise ValueError("Unexpected MungeSumstats namespace")
                location = "<selected installation>"
            entries[package] = [row["version"], location]
    if "MungeSumstats" not in entries:
        raise ValueError("MungeSumstats namespace was not observed")
    return entries


def verify_manifests(config, work, runs):
    output = work / "stable-manifests.txt"
    manifests = [str(run / "edition" / "manifest.rds") for run in runs]
    command = r_command(config["r_executable"]["path"], HERE / "worker.R",
                        "manifests", str(output), *manifests)
    process = run_process(command, work, name="manifest-verification")
    write_json(work / "manifest-verification.json", process)
    if process["exit_code"] != 0 or output.read_text(encoding="ascii").strip() != "pass":
        raise ValueError("Stable RDS manifests disagree")


def summarize_case(case, observations, *, suite):
    def phase_values(entries, name):
        return {phase: statistics.median(float(item["metrics"][f"{name}.{phase}"]) for item in entries)
                for phase in ("elapsed", "user.self", "sys.self", "user.child", "sys.child")}

    paired = [{role: next(item for item in observations
                         if item["pair"] == pair and item["process"]["role"] == role)
               for role in ROLES} for pair in range(1, case["pairs"] + 1)]
    measures = {}
    for role in ROLES:
        entries = [pair[role] for pair in paired]
        processes = [entry["process"] for entry in entries]
        rss = [entry["peak_rss_bytes"] for entry in processes if entry["peak_rss_bytes"]]
        measures[role] = {
            "package_setup": phase_values(entries, "package_setup"),
            "reference_open": phase_values(entries, "reference_open") if suite != "legacy" else None,
            "processing": phase_values(entries, "processing"),
            "process_wall_median_seconds": statistics.median(entry["wall_seconds"] for entry in processes),
            "process_user_median_seconds": statistics.median(entry["user_seconds"] for entry in processes),
            "process_system_median_seconds": statistics.median(entry["system_seconds"] for entry in processes),
            "median_peak_rss_bytes": statistics.median(rss) if len(rss) == len(entries) else None,
            "max_peak_rss_bytes": max(rss, default=None),
            "rss_complete": bool(rss) and all(entry["rss_complete"] for entry in processes),
        }
    base, candidate = (measures[role] for role in ROLES)
    a, b = (measure["processing"]["elapsed"] for measure in (base, candidate))
    pair_times = [(float(pair["baseline"]["metrics"]["processing.elapsed"]),
                   float(pair["candidate"]["metrics"]["processing.elapsed"])) for pair in paired]
    # Millisecond resolution from proc.time() can make a tiny hard-case call zero.
    ratio = a / b if a > 0 and b > 0 else None
    rss_ok = (candidate["median_peak_rss_bytes"] <= base["median_peak_rss_bytes"]
              if base["rss_complete"] and candidate["rss_complete"] else None)
    return {
        "medians": measures,
        "ratio_of_median_processing_seconds": ratio,
        "baseline_over_candidate_by_pair": [x / y if x > 0 and y > 0 else None for x, y in pair_times],
        "median_peak_rss_nonregression": rss_ok,
        "target": bool(case.get("target") and case["rows"] != 40_000_000 and suite == "matrix"),
    }


def benchmark(config, output, *, suite, profile, threshold, configuration_pin):
    hashes, packages, receipts = resolve_configuration(config, suite)
    output.mkdir()
    if shutil.disk_usage(output).free < config["minimum_free_bytes"] + config["estimated_output_bytes"]:
        raise ValueError("Free space does not cover estimated retained outputs and caller-supplied reserve")
    started = utc_now()
    before = asset_inventory(hashes)
    installed, toolchain = read_provenance(config, packages, receipts, before)
    cardinalities = {case["id"]: cardinality(case) for case in config["cases"]}
    write_json(output / "pins.json", {
        "configuration": config, "configuration_file": configuration_pin,
        "assets": before, "installed_packages": installed,
        "r_toolchain": toolchain, "cardinalities": cardinalities,
        "harness": {name: digest(HERE / name)
                    for name in ("run.py", "worker.R", "verify.py", "process_metrics.py")},
        "host": {"platform": platform.platform(), "architecture": platform.machine(),
                 "logical_cpus": os.cpu_count(), "python": platform.python_version()},
        "suite": suite, "profile": profile, "started_utc": started,
        "cache": ("Fresh R process and verified reference handle per run in preservation. "
                  "OS page cache uncontrolled; preread by asset pins, never flushed."),
    })
    runs = {}
    runtime = {}
    dependencies = {}
    for case in config["cases"]:
        work = output / case["id"]
        work.mkdir()
        observations = runs[case["id"]] = []
        for repetition in range(1, case["pairs"] + 1):
            order = ROLES if repetition % 2 else tuple(reversed(ROLES))
            for role in order:
                if shutil.disk_usage(output).free < config["minimum_free_bytes"]:
                    raise RuntimeError("Caller-supplied free-space reserve was reached")
                run_dir = work / f"pair-{repetition:02d}-{role}"
                process, metrics = run_one(config, case, run_dir, role, suite=suite,
                                           profile=profile, threshold=threshold, assets=before)
                relevant = verify_runtime(config, metrics, toolchain)
                state = dependency_state(run_dir, config, role)
                if not runtime:
                    runtime = relevant
                if runtime != relevant:
                    raise ValueError("R runtime, locale or thread settings vary across runs")
                if case["id"] in dependencies and state != dependencies[case["id"]]:
                    raise ValueError(f"Namespace versions or locations vary within {case['id']}")
                dependencies[case["id"]] = state
                observations.append({"directory": str(run_dir), "pair": repetition,
                                     "process": process, "metrics": metrics})

    # Data and installed package rechecks and all artifact work follow the final timed run.
    verified_at = time.monotonic()
    after = asset_inventory(hashes)
    after_installed, _ = read_provenance(config, packages, receipts, after)
    if before != after or installed != after_installed:
        raise ValueError("Post-run asset or installation pins changed")
    write_json(output / "post-pins.json", {"assets": after, "installed_packages": after_installed})
    summaries = {}
    for case in config["cases"]:
        work = output / case["id"]
        observations = runs[case["id"]]
        directories = [Path(item["directory"]) for item in observations]
        if suite == "legacy":
            entries, offsets = compare_bytes([path / "legacy.tsv" for path in directories])
            equality = {"legacy_tsv_digests": dict(zip(map(str, directories), entries)),
                        "differing_blocks": offsets}
            write_json(work / "verification.json", equality)
            if offsets or len({item["metrics"]["legacy_returned_rows"] for item in observations}) != 1:
                raise ValueError("Legacy adapter outputs or counts disagree")
        else:
            reference = config["references"][case["reference"]]
            equality = verify_editions(case, reference, directories, work, before,
                                       config.get("edition_encoding", "utf-8"))
            verify_manifests(config, work, directories)
            equality["stable_rds_manifest_equality"] = True
            write_json(work / "verification.json", equality)
        summaries[case["id"]] = {
            **summarize_case(case, observations, suite=suite),
            "rows": case["rows"], "cardinality": cardinalities[case["id"]],
            "observations": observations, "verification": equality,
        }
    targets = [summary for summary in summaries.values() if summary["target"]]
    outcome = None
    target_status = "not_applicable"
    exit_code = 0
    if targets:
        checks = []
        for item in targets:
            speedup = item["ratio_of_median_processing_seconds"]
            checks.extend((None if speedup is None else speedup >= 2,
                           item["median_peak_rss_nonregression"]))
        if False in checks:
            target_status, outcome, exit_code = "miss", False, 2
        elif None in checks:
            target_status, exit_code = "indeterminate", 3
        else:
            target_status, outcome = "pass", True
    write_json(output / "summary.json", {
        "scope": ("Adapted legacy format_sumstats with historical matched-SNP reference loader"
                  if suite == "legacy" else "Complete installed public harmonise_file calls"),
        "suite": suite, "profile": profile, "target_status": target_status,
        "processing_speed_and_rss_target_met": outcome,
        "fresh_processes": True, "verification_seconds": time.monotonic() - verified_at,
        "started_utc": started, "finished_utc": utc_now(), "cases": summaries,
        "target_rule": "Ratio of processing medians >= 2 and no increase in median OS peak RSS.",
    })
    return exit_code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    experiment = commands.add_parser("run")
    experiment.add_argument("configuration", type=Path)
    experiment.add_argument("output", type=Path)
    experiment.add_argument("--suite", choices=("matrix", "diagnostic", "legacy"), default="matrix")
    experiment.add_argument("--profile", choices=("none", "cpu", "alloc"), default="none")
    experiment.add_argument("--allocation-threshold", type=int, default=65536)
    pin = commands.add_parser("pin", help="Calculate hashes for admitted explicit files or installed trees")
    pin.add_argument("files", type=Path, nargs="+")
    args = parser.parse_args()
    if args.action == "pin":
        for path in args.files:
            path = path.resolve(strict=True)
            identity = package_tree(path)["tree_sha256"] if path.is_dir() else digest(path)["sha256"]
            print(f"{identity}  {path}")
        return 0
    if args.suite == "matrix" and args.profile != "none":
        parser.error("Use a separate diagnostic suite for profiled runs")
    integer(args.allocation_threshold, 0, 2**31 - 1, "allocation threshold")
    output = args.output.resolve()
    if output.exists():
        parser.error("Output must be a new directory under an existing parent")
    try:
        config_path = args.configuration.resolve(strict=True)
        configuration_pin = {"path": str(config_path), **digest(config_path)}
        config = read_json(config_path)
        return benchmark(config, output, suite=args.suite, profile=args.profile,
                         threshold=args.allocation_threshold, configuration_pin=configuration_pin)
    except (OSError, RuntimeError, ValueError, KeyError, TypeError, csv.Error, InvalidOperation) as error:
        if output.is_dir():
            write_json(output / "failure.json", {"error": str(error), "utc": utc_now()})
        print(error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
