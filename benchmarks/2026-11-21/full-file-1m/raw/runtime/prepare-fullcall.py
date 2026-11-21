#!/usr/bin/env python3
"""Prepare schema-1 smoke metadata from successful, pinned Darwin builds."""
import argparse
import copy
import csv
import hashlib
import json
import os
import shlex
import stat
import subprocess
import sys
import tarfile
from pathlib import Path

TOP = Path(__file__).resolve().parents[1]
BENCH = Path("/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade")
STAGE176 = Path("/private/tmp/mss-parser-qualification-stage176/build/build.json")
STAGE178 = Path("/private/tmp/mss-functional-stage178")
COMMON = TOP / "restored-01/common"
INSTALL = COMMON.parent / "receipts"
RHOME = Path("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources")
CLT = Path("/Library/Developer/CommandLineTools")
OTOOL = CLT / "usr/bin/otool"
ROLES = {
    "baseline": ("a86e8c94d753630eff8c6a695ed2f240d1123d0d", "e728d244e98190a4df9f410d7ad651276fc473e4",
                 "c4917d6e4276afaf8863d91754e77c00640315d61bef2528fca68bcb1b17698a"),
    "candidate": ("bfa192c83bb1fe15043c781d3b330f0f71bafb70", "28a48dda81b644bc0999cfc3922f1a4a26b4b75d",
                  "8b47c8164251f5380f282e92a21a1ce0686af468e83583c344bb780e95a9cd5d"),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def checked_source(path, sha):
    require(hashlib.sha256(path.read_bytes()).hexdigest() == sha, f"Control pin changed: {path}")
    return path


def guarded_phase(root, prefix, label, command):
    files = tuple(prefix.with_suffix(f".{suffix}") for suffix in ("jsonl", "stdout", "stderr"))
    planned, running, terminal = map(json.loads, files[0].read_text().splitlines())
    require(all(row["label"] == label for row in (planned, running, terminal))
            and planned["status"] == "planned" and running["status"] == "running"
            and planned["argv"] == list(map(str, command)) and planned["cwd"] == str(root)
            and [planned[key] for key in ("receipt", "stdout", "stderr")] == list(map(str, files))
            and running["pid"] == running["expected_pgid"] == running["expected_sid"],
            f"{label}: unexpected persisted phase binding")
    require(terminal["status"] == "completed" and terminal["closure"] == "confirmed"
            and type(terminal["raw_exit"]) is int and terminal["raw_exit"] == 0
            and terminal["leader_reaped"] is True and terminal["primary"] is None
            and not terminal["cleanup_errors"] and not terminal.get("receipt_error"),
            f"{label}: no successful closed terminal record")
    return files, terminal


def guarded_role(root, role, sdk):
    work = root / role
    archive, library = work / "MungeSumstats.tar", work / "library"
    command = ["/usr/bin/time", "-l", "/usr/bin/env", f"SDKROOT={sdk}", RHOME / "bin/R",
               "CMD", "INSTALL", "--use-vanilla", f"--library={library}", archive]
    return archive, library, *guarded_phase(root, work / "install", f"normal-install-{role}", command)


def native_inventory(paths, pin):
    images, resolved = {}, {}
    for path in sorted(set(paths)):
        text = subprocess.run([OTOOL, "-arch", "arm64", "-L", path], check=True, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10,
                              env={"PATH": "/usr/bin:/bin", "LC_ALL": "C"}).stdout
        lines = text.splitlines()
        require(len(lines) >= 2 and lines[0] == f"{path}:", f"No arm64 link record: {path}")
        images[str(path)] = {"file": pin(path), "arm64_link_entries": lines[1:]}
        for line in lines[1:]:
            link = line.strip().split(" (compatibility version", 1)[0]
            target = (path.parent / link.removeprefix("@loader_path/")) if link.startswith("@loader_path/") else Path(link)
            if target.is_absolute() and target.is_file() and str(target) not in resolved:
                resolved[str(target)] = pin(target)
    return {"images": images, "resolved_link_files": resolved}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--admission-note", type=Path, required=True)
    parser.add_argument("--check", action="store_true", help="Read and validate the inputs without writing metadata")
    args = parser.parse_args()
    root, output = args.build_root.resolve(strict=True), args.output.resolve(strict=True)
    require(root == TOP / "builds-01" and output.parent == TOP and output != root
            and output.stat().st_uid == os.getuid() and stat.S_IMODE(output.stat().st_mode) == 0o700,
            "Expected builds-01 and a distinct private output directory")
    names = ("baseline-build.json", "candidate-build.json", "dependency-inventory.json", "installations.json", "smoke.json")
    require(not any(os.path.lexists(output / name) for name in names), "Metadata output already exists")
    note = args.admission_note.read_text().strip()
    require(note, "Record actual local smoke admission conditions")
    historical = json.loads(checked_source(STAGE176, "2a8c81d9ba5f7cf1f9055fdbe04c2f386355adae01dca967e53c5444f4e0337e").read_text())
    old = json.loads(checked_source(STAGE178 / "smoke-baseline-build.json", "0c76263d1b4fb9a649b5ebb92d80298c195cca1f4d7838901f638644153c9dc2").read_text())
    smoke = json.loads(checked_source(STAGE178 / "smoke.json", "b395193366af7a2d3c27457b22b486d05529e574b09744ed37deaebba14ff685").read_text())
    toolchain, sdk = old["toolchain"], historical["toolchain"]["sdk"]
    builds = {role: guarded_role(root, role, sdk) for role in ROLES}
    for key, sha in {
        "run": "c0c5f2f72a334b2b77721eccca59a00d1e47fb1663f513ddcbe51700a13ec043",
        "verify": "cb0f9774aa9523514ce9fc75868e438272606de3eb3564699addca103de1930b",
        "process_metrics": "f402957540febe6e0f9e8c56422b7acc5f68bcc32afe4377705347de199a05d7",
        "operations/build": "a17b8a9f9c998a0ae27c35df9da31e08f7a6cba56a8b870d6af022350fb3a3f8",
    }.items():
        checked_source(BENCH / f"{key}.py", sha)
    worker = checked_source(BENCH / "worker.R", "f8477b5b0047b0361d636e4af63b5ed9bf2be87229d7cf9c682331a51147e7df")
    sys.path.insert(0, str(BENCH / "operations"))
    from build import inventory_library, pin, package_tree
    from run import asset_inventory, read_provenance, resolve_configuration
    from verify import tsv_pairs, write_json

    upstream = old["upstream"]["full_toolchain"]
    physical = [*historical["toolchain"]["pins"], upstream["Makeconf"], upstream["R_wrapper"]]
    runtime = historical["runtime"]
    require(all(pin(Path(p["path"]))["sha256"] == p["sha256"]
                for p in [*physical, runtime["r_executable"], runtime["r_engine"], *runtime["runtime_files"]])
            and not Path(historical["toolchain"]["absent_makevars_site"]).exists(),
            "R, compiler, SDK or runtime binary pins changed")
    runtime_tsv = output / "runtime.tsv"
    probe_files, _ = guarded_phase(output, output / "runtime", "runtime",
                                   [RHOME / "bin/Rscript", "--vanilla", worker, "runtime", runtime_tsv])
    observed = tsv_pairs(runtime_tsv, ("key", "value"))
    fields = {"r_shared": "R_shared", "blas": "BLAS_library", "lapack": "LAPACK_library"}
    require(observed["R_version"] == toolchain["R_version"] and observed["R_platform"] == toolchain["R_platform"]
            and Path(observed["R_home"]).resolve(strict=True) == RHOME.resolve(strict=True)
            and observed["R_engine"] == runtime["r_engine"]["path"]
            and all(observed[fields[p["role"]]] == p["path"] for p in runtime["runtime_files"]),
            "Observed runtime differs from verified build toolchain")
    installations = {**runtime, "dependency_libraries": [str(COMMON)]}
    tree_details, role_receipts = {}, {}
    expected_commands = [shlex.split(line.replace("/private/tmp/harvest-Rlib/Rcpp/include", str(COMMON / "Rcpp/include")))
                         for line in upstream["actual_compiler_and_link_commands"]]
    for role, (archive, library, files, terminal) in builds.items():
        commit, tree, sha = ROLES[role]
        source_pin = pin(archive)
        with tarfile.open(archive, "r:") as package:
            require(package.pax_headers.get("comment") == commit and source_pin["sha256"] == sha,
                    f"{role}: source archive pin or commit changed")
        log = files[2].read_text()
        required = [toolchain["compiler"][0], f"using SDK: ‘{Path(sdk).name}’", "using C++17",
                    "** testing if installed package can be loaded from temporary location",
                    "** testing if installed package can be loaded from final location", "* DONE (MungeSumstats)"]
        actual = [shlex.split(line) for line in files[1].read_text().splitlines() if line.startswith("clang++ ")]
        require(all(item in log for item in required) and actual == expected_commands
                and set(library.iterdir()) == {library / "MungeSumstats"}, f"{role}: incomplete or unexpected install")
        package = library / "MungeSumstats"
        dll = pin(package / "libs/MungeSumstats.so")
        results = package_tree(package)
        tree_details[role] = results
        role_receipts[role] = {"source_commit": commit, "source_tree": tree, "archive": source_pin,
                               "dll_sha256": dll["sha256"], "package_tree_sha256": results["tree_sha256"],
                               "toolchain": toolchain, "compile_and_link_argv": actual,
                               "install_transcripts": [pin(p) for p in files], "terminal": terminal}
        installations[role] = {"commit": commit, "library": str(library), "dll": dll,
                               "package_tree_sha256": results["tree_sha256"]}

    packages_tsv = INSTALL / "installation/packages.tsv"
    with packages_tsv.open() as stream:
        package_rows = list(csv.DictReader(stream, delimiter="\t"))
    lib = (RHOME / "library").resolve(strict=True)
    dependency_trees = {str(p): inventory_library(p) for p in (COMMON, lib)}
    require((INSTALL / "install.status").read_text().strip() == "0" and len(package_rows) == 94
            and set(dependency_trees[str(COMMON)]) == {p["package"] for p in package_rows}
            and all(Path(p["path"]) == COMMON / p["package"]
                    and pin(Path(p["description"]))["sha256"] == pin(Path(p["path"]) / "DESCRIPTION")["sha256"]
                    for p in package_rows), "Common library differs from the successful restoration inventory")
    paths = [Path(item["path"]) for item in runtime["runtime_files"]]
    for directory in (COMMON, lib, RHOME / "lib", *(Path(installations[role]["library"]) for role in ROLES)):
        paths.extend(p for p in directory.rglob("*") if p.suffix in {".so", ".dylib"}
                     and p.is_file() and not any(part.endswith(".dSYM") for part in p.parts))
    native = native_inventory(paths, pin)
    audit = {"runtime_receipt": pin(runtime_tsv), "runtime_observed": observed,
             "runtime_transcripts": [pin(p) for p in probe_files],
             "build_toolchain_pins": physical, "benchmark_worker": pin(worker),
             "restore_install_receipt": pin(packages_tsv), "installed_dependency_tree_sha256": dependency_trees,
             "role_trees": tree_details, "native_linkage": native, "otool": pin(OTOOL)}
    print(f"Verified both build trees, {len(package_rows)} common packages, "
          f"{len(dependency_trees[str(lib)])} framework packages and {len(native['images'])} native images")
    if args.check:
        hashes, _, _ = resolve_configuration(smoke, "diagnostic")
        earlier_roles = {Path(smoke[role][key]["path"]) for role in ROLES for key in ("dll", "build_receipt")}
        asset_inventory({path: sha for path, sha in hashes.items() if path not in earlier_roles})
        return
    write_json(output / "dependency-inventory.json", audit)
    installations["dependency_inventory"] = pin(output / "dependency-inventory.json")
    for role, receipt in role_receipts.items():
        path = output / f"{role}-build.json"
        write_json(path, receipt)
        installations[role]["build_receipt"] = pin(path)
    smoke.update(installations, contention_note=note, minimum_free_bytes=186 * 1024**3)
    resolved = copy.deepcopy(smoke)
    hashes, pkgs, receipts = resolve_configuration(resolved, "diagnostic")
    read_provenance(resolved, pkgs, receipts, asset_inventory(hashes))
    write_json(output / "installations.json", installations)
    write_json(output / "smoke.json", smoke)


if __name__ == "__main__":
    os.umask(0o077)
    main()
