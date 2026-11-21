"""Guarded current-path handoff to the frozen C9 same-handle companion."""
import argparse
import hashlib
import importlib.util
import io
import json
import math
import os
import stat
import sys
import time
from contextlib import redirect_stdout
from decimal import Decimal
from pathlib import Path

BASE = Path("/private/tmp/mss-fullmatrix-20261121")
PREP, COLD, WARM = (BASE / name for name in ("prep-01", "cold-01", "warm-01"))
MATRIX = COLD / "matrix"
PACKET = WARM / "cold-handoff.json"
OBS = BASE / "cold-observation.json"
CONFIG, PROFILE = PREP / "matrix.json", PREP / "snp-oracle.json"
HIST = Path("/private/tmp/mss-performance-20261119")
STAGE = HIST / "tools/darwin-matrix-stage185"
BENCH = HIST / "repo/benchmarks/2026-11-19/performance-upgrade"
WORKER, PRIME = BENCH / "worker.R", STAGE / "prime.R"
ORACLE = HIST / "repo/inst/verification/fai-safety-stage8/whole-file/harness/verify.py"
GUARD = HIST / "tools/darwin-full-stage185/guard.py"
BUILD = BENCH / "operations/build.py"
PYTHON = Path("/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15")
COLD_PY = "/opt/homebrew/opt/python@3.15/bin/python3.15"
SOURCES = {
    STAGE / "driver.py": "16e82e1e581dd92e911f8947f2465f4092759124bc157d1a1b0f8133c10761b4",
    STAGE / "warm.py": "46dac35e0452add60a41a0b5923b2096d0cce14309421bf0ca4585702f4c6382",
    PRIME: "28d96fb69242596bd5c9adb73013f5de70f8cb6f326f3900ed780c8ba5118661",
    BENCH / "run.py": "c0c5f2f72a334b2b77721eccca59a00d1e47fb1663f513ddcbe51700a13ec043",
    BENCH / "verify.py": "cb0f9774aa9523514ce9fc75868e438272606de3eb3564699addca103de1930b",
    BENCH / "process_metrics.py": "f402957540febe6e0f9e8c56422b7acc5f68bcc32afe4377705347de199a05d7",
    WORKER: "f8477b5b0047b0361d636e4af63b5ed9bf2be87229d7cf9c682331a51147e7df",
    ORACLE: "d4bbe46a9037597203e779b9e0f830a293b15a476421973fc6a0b69b2dc45c95",
    GUARD: "58e0818aa20801c30194bdce26bd117112b4b4d5d10f74ea556360c0867b56a2",
    BUILD: "a17b8a9f9c998a0ae27c35df9da31e08f7a6cba56a8b870d6af022350fb3a3f8",
    PYTHON: "6cc907e4a022fefc843bcd93954eaf5a2cbfc77a095640bfadeeed96154190c5",
}
ACCEPTED = {
    CONFIG: "156143e169744ead89eb7a01721ac433da0d0983be68924b417bcc72f00426c7",
    PROFILE: "b387920dbab4d53af316a8b2b320d577be463abdad7ce80dbf8b9de411e303d0",
    PREP / "preparation.json": "e913abf35404a3b9851fb74a845e11921b0011f016ab28b611a4b0b367c29922",
    PREP / "COMPLETE.json": "ef3ab62f840dde88841e91583d8f12ca99741ae1f4cb89068d7dd443f2d110a4",
    COLD / "admission.json": "d668e0ed1b02641e6e7a125dd6452948d8ad4cadd45982fbf0b014a615d70bab",
}
OPEN = ["adapted legacy", "CPU and allocation profiles", "current lookup pressure",
        "full installed functional qualification", "Windows qualification", "reviewed formal Stage185 successor"]
CONTROLS = ("settings.tsv", "metrics.tsv", "namespaces.tsv", "process.json")
ROLE_LOGS = ("process.stdout.log", "process.stderr.log", "r-heap.rds")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pin(path, sha):
    return {"path": str(path), "sha256": sha}


def modules():
    global d, w, r, build, v
    require(Path(sys.executable).resolve(strict=True) == PYTHON and sys.flags.isolated
            and sys.dont_write_bytecode and sys.flags.no_site, "Require pinned isolated Python -I -B -S")
    # Authenticate the whole source list before executing even the bootstrap driver.
    for path, sha in SOURCES.items():
        require(path.resolve(strict=True) == path and path.is_file() and path.stat().st_size < 1024**2,
                f"Unexpected frozen source {path}")
        with path.open("rb") as stream:
            data = stream.read(1024**2 + 1)
        require(len(data) <= 1024**2 and hashlib.sha256(data).hexdigest() == sha, f"Changed frozen source {path}")
        if path == STAGE / "driver.py":
            driver_bytes = data
    path = STAGE / "driver.py"
    spec = importlib.util.spec_from_file_location("driver", path)
    d = importlib.util.module_from_spec(spec)
    sys.modules["driver"] = d
    exec(compile(driver_bytes, str(path), "exec"), d.__dict__)
    for path, name in [(BENCH / (n + ".py"), n) for n in ("process_metrics", "verify", "run")] + [
            (STAGE / "warm.py", "matrix_warm"), (BUILD, "matrix_build")]:
        require(name not in sys.modules, f"Unexpected loaded module {name}")
        d.module(pin(path, SOURCES[path]), name)
    r, v, w, build = (sys.modules[n] for n in ("run", "verify", "matrix_warm", "matrix_build"))
    sys.path.remove(str(BENCH))  # build.py's legacy convenience import directory


def root_identity(root):
    info = root.stat()
    require(root.resolve(strict=True) == root and stat.S_ISDIR(info.st_mode)
            and stat.S_IMODE(info.st_mode) == 0o700 and info.st_uid == os.getuid(), f"Expected private {root}")
    return {"path": str(root), "dev": info.st_dev, "ino": info.st_ino}


class Controls:
    """Record the digest of each small retained control as its actual bytes are consumed."""
    def __init__(self):
        self.files = {}

    def data(self, path):
        path = Path(path)
        require(path.resolve(strict=True) == path, f"Indirect control {path}")
        data = d.control_bytes(path)
        key, sha = str(path), d.sha(data)
        require(self.files.setdefault(key, sha) == sha, f"Changed control {path}")
        return data

    def pin(self, path):
        self.data(path)
        return pin(path, self.files[str(path)])

    def obj(self, path):
        return json.loads(self.data(path))


def check_warm_environment(guard_keys, actual):
    expected = d.environment(WARM)
    # Darwin may add this key at interpreter startup, after the guard records its environment.
    startup_keys = {"__CF_USER_TEXT_ENCODING"} if sys.platform == "darwin" else set()
    require(set(guard_keys) == expected.keys() and (actual.keys() - expected.keys()) <= startup_keys
            and all(actual.get(key) == value for key, value in expected.items()),
            "Warm worker environment differs from frozen single-thread settings")


def dispatch(label):
    root = root_identity(WARM)
    admission = json.loads(d.private_grant_bytes(WARM / "admission.json"))
    require(admission["root"] == root and admission["guard_sha256"] == SOURCES[GUARD]
            and admission["adapter_sha256"] == d.pin_control(Path(__file__))["sha256"]
            and Path.cwd() == WARM, "Warm admission/root differs")
    path = WARM / (label + ".jsonl")
    # The native guard appends running after launching the process.
    deadline = time.monotonic() + 5
    while True:
        try:
            data = d.control_bytes(path)
        except ValueError as exc:
            if not str(exc).startswith("Control changed during read:"):
                raise
            data = b""
        if data.endswith(b"\n") and data.count(b"\n") >= 2:
            break
        require(time.monotonic() < deadline, "Missing live dispatch record")
        time.sleep(0.01)
    rows = [json.loads(line) for line in data.splitlines()]
    request = {"label": label, "allowed": [0], "root": str(WARM), "parent_pid": admission["supervisor_pid"]}
    d.validate_dispatch(request, rows, sys.orig_argv,
                        (os.getppid(), os.getpid(), os.getpgrp(), os.getsid(0), os.geteuid()), time.monotonic())
    planned = rows[0]
    check_warm_environment(rows[1]["environment_keys"], dict(os.environ))
    require(planned["receipt"] == str(path) and all(planned[k] == str(WARM / f"{label}.{k}")
            for k in ("stdout", "stderr")), "Unexpected warm receipt/log paths")
    limits = planned["limits"]
    require(0 < limits["rss_bytes"] <= 24_000_000_000 and limits["floor_bytes"] >= 186 * 1024**3
            and 0 < limits["disk_bytes"] <= 50_000_000_000, "Unexpected warm guard envelope")
    return root


def cold_phases(c):
    observer = c.obj(OBS)
    admission = c.obj(COLD / "admission.json")
    require(observer["schema"] == "current-cold-observation-v1" and observer["root"] == root_identity(COLD)
            and observer["guard"] == pin(GUARD, SOURCES[GUARD])
            and observer["admission"] == c.pin(COLD / "admission.json"), "Unbound root observation")
    phases = observer["phases"]
    cold_argv = [COLD_PY, "-E", "-s", "-B", str(BENCH / "run.py"), "run", str(CONFIG), str(MATRIX),
                 "--suite", "matrix", "--profile", "none"]
    require(admission["argv"] == cold_argv and Path(COLD_PY).resolve(strict=True) == PYTHON
            and admission["guard_sha256"] == SOURCES[GUARD] and admission["config_sha256"] == ACCEPTED[CONFIG],
            "Cold launch source/config differs")
    commands = {"full-cold-matrix": cold_argv}
    for fixture in d.SUPPLEMENTS:
        for role in d.ROLES:
            commands[f"cross-{fixture}-{role}"] = ["/usr/bin/cmp", *(str(
                MATRIX / f"{fixture}-f{fmt}" / f"pair-01-{role}/edition/normalized.tsv") for fmt in (2, 3))]
    for size in d.SCALES:
        edition = MATRIX / f"matched-{size}/pair-01-candidate/edition"
        commands[f"snp-oracle-{size}"] = [str(PYTHON), "-I", "-B", "-S", str(ORACLE), "--profile", str(PROFILE),
            "--stage", size, "--edition", f"mss7={edition}", "--receipt", str(COLD / f"snp-oracle-{size}.json")]
    require(set(phases) == commands.keys(), "Require matrix, sixteen comparisons and three SNP phases")
    post_path = COLD / "post-cold-admission.json"
    post = c.obj(post_path)
    require(post["guard_sha256"] == SOURCES[GUARD] and post["config_sha256"] == ACCEPTED[CONFIG]
            and post["root"] == root_identity(COLD), "Unbound post-cold admission")
    matrix_end = None
    for label, argv in commands.items():
        stem = "matrix" if label == "full-cold-matrix" else label
        allowed = [0, 2, 3] if stem == "matrix" else [0]
        path = COLD / (stem + ".jsonl")
        seen = phases[label]
        require(seen["receipt"] == c.pin(path), f"Wrong observed receipt {label}")
        planned, running, terminal = [json.loads(line) for line in c.data(path).splitlines()]
        raw = d.phase_status(terminal, allowed)
        grant = admission if stem == "matrix" else post
        require(seen["admission"] == c.pin(COLD / "admission.json" if stem == "matrix" else post_path)
                and planned["limits"]["rss_bytes"] == grant["rss_group_limit_bytes"]
                and planned["limits"]["floor_bytes"] == grant["guard_floor_bytes"], f"Wrong admission {label}")
        if stem == "matrix":
            matrix_end = terminal["utc"]
        else:
            require(matrix_end <= planned["utc"], f"Phase predates matrix closure {label}")
        require(all(x["label"] == label for x in (planned, running, terminal))
                and planned["status"] == "planned" and running["status"] == "running"
                and planned["receipt"] == str(path) and planned["argv"] == argv and planned["cwd"] == str(COLD)
                and planned["allowed_exits"] == allowed and running["uid"] == os.geteuid()
                and running["pid"] == running["expected_pgid"] == running["expected_sid"] > 1
                and seen["pgid"] == seen["sid"] == running["pid"]
                and seen["group_absent"] is True and seen["session_absent"] is True
                and type(seen["observed_raw_exit"]) is int and seen["observed_raw_exit"] == raw
                and type(seen["observed_supervisor_exit"]) is int and seen["observed_supervisor_exit"] == 0,
                f"Parent observation or phase differs {label}")
        for key in ("stdout", "stderr"):
            log = COLD / f"{stem}.{key}"
            require(planned[key] == str(log), f"Incorrect {label} {key}")
            c.pin(log)
        if stem == "matrix":
            require(set(running["environment_keys"]) == set(d.environment(COLD)), "Cold environment keys")
            for field, bound in (("wall_seconds", "wall_seconds"), ("rss_bytes", "rss_group_limit_bytes"),
                                 ("disk_bytes", "cold_pool_bytes"), ("floor_bytes", "guard_floor_bytes")):
                require(planned["limits"][field] == admission[bound], f"Cold {field} differs")
    return phases["full-cold-matrix"]["observed_raw_exit"]


def dependencies(config):
    deps = json.loads(d.check_control(config["dependency_inventory"]))
    trees = deps["installed_dependency_tree_sha256"]
    require(len(trees) == 2 and sorted(map(len, trees.values())) == [30, 94]
            and all(build.inventory_library(Path(path)) == packages for path, packages in trees.items()),
            "Common/framework dependency inventory changed")
    return deps


def expected_settings(cfg, case, role, dest, assets):
    ref, options = cfg["references"][case["reference"]], case["settings"]
    dll = cfg[role]["dll"]["path"]
    values = dict(mode="matrix", profile="none", allocation_threshold=65536, input=case["input"]["path"],
                  rows=case["rows"], output=str(dest), dll=dll, dll_md5=assets[dll]["md5"],
                  reference=ref["directory"], core_version=f'mss-preservation-{ref["format"]}',
                  fasta_pin_path=ref["fasta"]["path"])
    for i, lib in enumerate([cfg[role]["library"], *cfg["dependency_libraries"]]):
        values[f"library.{i}"] = lib
    if ref["bind_fasta"]:
        values["fasta"] = ref["fasta"]["path"]
    values.update({k: val for k, val in options.items() if val is not None
                   and k not in {"columns", "verify_acquisition_receipt"}})
    values.update({f"column.{k}": val for k, val in options.get("columns", {}).items()})
    if options["verify_acquisition_receipt"]:
        src = assets[values["input"]]
        values.update(acquisition_md5=src["md5"], acquisition_bytes=src["bytes"])
    return {k: str(val) for k, val in values.items()}


def checked_process(proc, argv, *, complete_rss=True):
    require(proc["command"] == argv and type(proc["exit_code"]) is int and proc["exit_code"] == 0
            and type(proc["pid"]) is int and proc["pid"] > 1 and proc["rss_source"] == "wait4 ru_maxrss"
            and type(proc["peak_rss_bytes"]) is int and proc["peak_rss_bytes"] >= 0
            and proc["rss_complete"] is (proc["peak_rss_bytes"] > 0) and (proc["rss_complete"] or not complete_rss)
            and proc["wall_seconds"] > 0 and all(math.isfinite(proc[k]) and proc[k] >= 0
                for k in ("wall_seconds", "user_seconds", "system_seconds")), f"Failed process or accounting {argv}")


def checked_run(c, cfg, case, dest, role, assets, deps, toolchain, worker):
    controls = {name: c.pin(dest / name) for name in CONTROLS}
    exp = expected_settings(cfg, case, role, dest, assets)
    require(w.settings(controls["settings.tsv"]) == exp, f"Settings differ {dest}")
    process = c.obj(dest / "process.json")
    checked_process(process, r.r_command(cfg["r_executable"]["path"], worker, "run", str(dest / "settings.tsv")),
                    complete_rss=worker == WORKER)
    require(process["case"] == case["id"] and process["role"] == role
            and process["started_utc"] <= process["finished_utc"], f"Process identity differs {dest}")
    metrics = w.settings(controls["metrics.tsv"])
    require(all(Decimal(metrics[k]) == case["rows"] for k in ("source_rows", "normalized_rows"))
            and all(exp[k] == metrics[k] for k in ("dll", "dll_md5")), f"Worker rows or DLL differ {dest}")
    require(r.verify_runtime(cfg, metrics, toolchain) == {**deps["runtime_observed"], "threads": "1"},
            f"Runtime differs {dest}")
    state = r.dependency_state(dest, cfg, role)
    trees = deps["installed_dependency_tree_sha256"]
    for pkg, (version, location) in state.items():
        if pkg != "MungeSumstats":
            require(version and any(pkg in names and Path(location) == Path(lib) / pkg for lib, names in trees.items()),
                    f"Unknown namespace {location}")
    return controls, process, metrics, state


def check_retained_digest(path, expected):
    require(r.digest(path) == expected, f"Retained edition digest differs {path}")


def checked_artifacts(c, cfg, case, dirs, work, assets, *, rehash):
    result = c.obj(work / "verification.json")
    inventory = c.obj(work / "artifact-digests.json")
    cols = list(v.NORMALIZED_COLUMNS)
    rows = case["rows"]
    require(result["edition_count"] == len(dirs) and result["normalized_columns"] == cols
            and result["normalized_rows_checked"] == sum(result["status_counts"].values()) == rows
            and result["stable_byte_equality"] is True and result["stable_rds_manifest_equality"] is True
            and result["source_artifact"] == assets[case["input"]["path"]]
            and set(inventory) == set(map(str, dirs)), f"Incomplete edition check {work}")
    gold = "oracle" in case
    require(result["oracle_rows_checked"] == (rows if gold else 0) and
            (set(result["oracle_fields"]) == set(cols) if gold else result["oracle_fields"] is None),
            f"Literal gold coverage differs {work}")
    decoded_name = "source.tsv" if case["settings"]["delimiter"] == "tab" else "source.csv"
    required = {"source-artifact", decoded_name, "normalized.tsv", "manifest.rds", "audit.tsv",
                "metadata.tsv", "columns.tsv", "checksums.tsv", "reference_checksums.tsv"}
    first = inventory[str(dirs[0])]
    require(required <= first.keys() <= required | {"reference_contigs.tsv"}
            and first["source-artifact"] == result["source_artifact"] and first[decoded_name] == result["decoded_source"],
            f"Wrong edition inventory {work}")
    for dest in dirs:
        members = inventory[str(dest)]
        edition = dest / "edition"
        require(set(members) == set(first) == {p.name for p in edition.iterdir()}, f"Edition membership {dest}")
        for name, dig in members.items():
            f = edition / name
            require(f.resolve(strict=True) == f and f.is_file() and f.stat().st_size == dig["bytes"]
                    and r.valid_hash(dig["sha256"]) and r.valid_hash(dig["md5"], 32)
                    and (name == "manifest.rds" or dig == first[name]), f"Edition bytes/receipt {f}")
            if rehash:
                check_retained_digest(f, dig)
    stable = work / "stable-manifests.txt"
    require(c.data(stable).strip() == b"pass", f"Failed RDS equality {work}")
    checked_process(c.obj(work / "manifest-verification.json"), r.r_command(cfg["r_executable"]["path"], WORKER,
                    "manifests", str(stable), *(str(p / "edition/manifest.rds") for p in dirs)))
    for stream in ("stdout", "stderr"):
        c.pin(work / f"manifest-verification.{stream}.log")
    return result, inventory


def check_oracle_counts(report, cold_statuses, rows):
    statuses = report["statuses"]
    resolved_labels = {"match", "swap", "complement", "swapped_complement", "sequence_only"}
    require(statuses == cold_statuses and statuses.keys() <= resolved_labels | {
                "ambiguous", "alt_alt_contrast", "allele_mismatch"}
            and all(type(n) is int and n > 0 for census in (statuses, cold_statuses) for n in census.values())
            and sum(statuses.values()) == rows, "Oracle status census differs")
    resolved = sum(n for status, n in statuses.items() if status in resolved_labels)
    expected = {
        "source_rows": rows, "normalized_rows": rows, "all_19_fields_correct_rows": rows,
        "expected_resolved_rows": resolved, "observed_resolved_rows": resolved,
        "effect_nonmissing_denominator": resolved, "effect_missing_correct_rows": rows - resolved,
        "eaf_nonmissing_denominator": resolved, "eaf_missing_correct_rows": rows - resolved,
    }
    counts = report["counts"]
    require(counts == expected and all(type(n) is int for n in counts.values()), "Oracle count scope differs")


def oracle(c, size, case, ref, assets, result, inventory):
    profile = c.obj(PROFILE)
    refdir = Path(ref["directory"])
    require(profile["stages"][size] == {"rows": case["rows"], "input": case["input"]}
            and ref["format"] == 2 and profile["reference"] == {
                "build": "GRCh37", "source": str(refdir / "source.tsv"), "fai": str(refdir / "fasta.fai"),
                "fasta": ref["fasta"]["path"]} and all(assets[str(refdir / f)]["sha256"] == ref["files"][f]
                for f in ("source.tsv", "fasta.fai")) and assets[ref["fasta"]["path"]]["sha256"] == ref["fasta"]["sha256"],
            f"Oracle profile differs {size}")
    label = f"snp-oracle-{size}"
    raw = COLD / f"{label}.json"
    obj = c.obj(raw)
    rows = case["rows"]
    require(obj == c.obj(COLD / (label + ".stdout")) and obj["schema"] == 1 and obj["stage"] == size
            and obj["passed"] is True and obj["all_rows_examined"] is True and obj["rows_expected"] == rows
            and obj["source_bytes_sha256"] == case["input"]["sha256"] and obj["source_copies_exact_each"] == 2
            and obj["scope"] == ("Repeated real chr22 unique biallelic non-palindromic SNP panel; every source and "
                "normalized row/all19fields verified; indels/multiallelic contexts are validated in separate "
                "scientific suites, not timed here") and obj["deltas"] == {} and obj["unaccounted_deltas"] == 0 and set(obj["versions"]) == {"mss7"},
            f"Oracle receipt/guard output {size}")
    panel = obj["reference_inventory"]
    require(all(panel[k] == 437458 for k in ("source_rows", "loci", "eligible_panel_records"))
            and panel["ids"] >= 437458 and all(panel[k] == 0 for k in (
                "multiallelic_loci", "quadallelic_loci", "ids_with_multiple_events")), f"Oracle panel {size}")
    check_oracle_counts(obj["versions"]["mss7"], result["status_counts"], rows)
    dest = MATRIX / case["id"] / "pair-01-candidate"
    edition = dest / "edition"
    binding = c.obj(COLD / f"oracle-binding-{size}.json")
    snapshots = {case["input"]["path"]: assets[case["input"]["path"]], **{
        str(edition / name): inventory[str(dest)][name] for name in ("source.tsv", "source-artifact", "normalized.tsv")}}
    require(binding["schema"] == "current-oracle-binding-v1" and binding["stage"] == size
            and binding["phase"] == c.pin(COLD / f"{label}.jsonl")
            and binding["oracle_source"] == pin(ORACLE, SOURCES[ORACLE]) and binding["profile"] == c.pin(PROFILE)
            and binding["edition"] == str(edition) and binding["exclusive_retained_window"] is True
            and binding["before"] == binding["after"] == snapshots, f"Oracle byte window {size}")
    return c.pin(raw)


def cold_matrix(c, raw):
    config = c.obj(CONFIG)
    prep = c.obj(PREP / "preparation.json")
    require(c.obj(PREP / "COMPLETE.json") == {"status": "prepared", "receipt": c.pin(PREP / "preparation.json")},
            "Wrong preparation seal")
    opening = c.obj(MATRIX / "pins.json")
    closing = c.obj(MATRIX / "post-pins.json")
    summary = c.obj(MATRIX / "summary.json")
    require(opening["suite"] == summary["suite"] == "matrix" and opening["profile"] == summary["profile"] == "none"
            and summary["fresh_processes"] is True and opening["configuration_file"] == {
                "path": str(CONFIG), **r.digest(CONFIG)}, "Wrong cold configuration or scope")
    hashes, packages, builds = r.resolve_configuration(config, "matrix")
    assets = r.asset_inventory(hashes)
    require(config == opening["configuration"] and config["references"] == prep["reference_formats"]
            and assets == opening["assets"] == closing["assets"] == prep["assets"]
            and opening["cardinalities"] == prep["cardinalities"], "Cold opening/closure/preparation differ")
    installed, toolchain = r.read_provenance(config, packages, builds, assets)
    deps = dependencies(config)
    require(installed == opening["installed_packages"] == closing["installed_packages"] == prep["installed_packages"]
            == deps["role_trees"] and toolchain == opening["r_toolchain"] == prep["toolchain"], "Role trees/toolchain")
    require({role: config[role]["commit"] for role in d.ROLES} == {
        "baseline": "a86e8c94d753630eff8c6a695ed2f240d1123d0d",
        "candidate": "bfa192c83bb1fe15043c781d3b330f0f71bafb70"}, "Role commits")
    prov = prep["provenance"]
    for field in ("installations", "dependency_inventory"):
        obj = json.loads(d.check_control(prov[field]))
        if field == "installations":
            require(all(config[k] == val for k, val in obj.items()), "Prepared installations")
        else:
            require(prov[field] == config[field], "Prepared dependency inventory")
        c.pin(Path(prov[field]["path"]))
    for role, path in builds.items():
        b = c.obj(path)
        expected = prov["builds"][role]
        require(c.pin(path) == expected["build_receipt"] and b["source_commit"] == expected["source_commit"]
                and b["source_tree"] == expected["source_tree"] and b["archive"] == expected["archive"],
                f"Prepared source build {role}")
    harness = ("run.py", "worker.R", "verify.py", "process_metrics.py")
    require(opening["harness"] == {name: r.digest(BENCH / name) for name in harness}, "Cold runner harness differs")
    roster = d.expected_roster()
    cases = {case["id"]: case for case in config["cases"]}
    require(set(roster) == set(cases) == set(summary["cases"]) and len(config["cases"]) == len(roster) == 19
            and prep["roster"] == {"expanded_cases": 19, "pairs": 23, "editions": 46}, "Cold matrix roster")
    primary = {"matrix_config": c.pin(CONFIG), "matrix_pins": c.pin(MATRIX / "pins.json"),
               "matrix_post_pins": c.pin(MATRIX / "post-pins.json"), "scale_controls": {}, "snp_oracles": {}}
    rss = {role: [] for role in d.ROLES}
    checks = []
    for name, (rows, pairs, fmt) in roster.items():
        case, row, work = cases[name], summary["cases"][name], MATRIX / name
        ref = config["references"][case["reference"]]
        require(row["rows"] == case["rows"] == rows and case["pairs"] == pairs and ref["format"] == fmt
                and row["cardinality"] == opening["cardinalities"][name] == r.cardinality(case), f"Case scope {name}")
        order = [(pair, role) for pair in range(1, pairs + 1) for role in (d.ROLES if pair % 2 else d.ROLES[::-1])]
        dirs = []
        comparator = {}
        namespaces = None
        for (pair, role), observed in zip(order, row["observations"], strict=True):
            dest = work / f"pair-{pair:02d}-{role}"
            dirs.append(dest)
            controls, process, metrics, state = checked_run(c, config, case, dest, role, assets, deps, toolchain, WORKER)
            require(observed == {"directory": str(dest), "pair": pair, "process": process, "metrics": metrics}
                    and (namespaces is None or state == namespaces), f"Unbound observation {dest}")
            namespaces = state
            rss[role].append(process)
            if pair == 1:
                comparator[role] = controls
        result, inventory = checked_artifacts(c, config, case, dirs, work, assets, rehash=True)
        measures = r.summarize_case(case, row["observations"], suite="matrix")
        require(row["verification"] == result and all(row[k] == val for k, val in measures.items())
                and measures["target"] == (name in {"matched-1M", "matched-10M"}), f"Cold case summary {name}")
        if measures["target"]:
            ratio = measures["ratio_of_median_processing_seconds"]
            checks.extend((None if ratio is None else ratio >= 2, measures["median_peak_rss_nonregression"]))
        if name.startswith("matched-"):
            size = name.removeprefix("matched-")
            primary["snp_oracles"][size] = oracle(c, size, case, ref, assets, result, inventory)
            primary["scale_controls"][size] = {**comparator, "artifact_digests": c.pin(work / "artifact-digests.json")}
    target = "miss" if False in checks else "indeterminate" if None in checks else "pass"
    require(summary["target_status"] == target == {0: "pass", 2: "miss", 3: "indeterminate"}[raw]
            and summary["processing_speed_and_rss_target_met"] is {"pass": True, "miss": False, "indeterminate": None}[target],
            "Median targets and raw matrix exit disagree")
    return primary, target, {role: d.summarize_rss(obs) for role, obs in rss.items()}


def seal(root):
    c = Controls()
    for path, sha in ACCEPTED.items():
        require(c.pin(path) == pin(path, sha), f"Accepted bytes changed {path}")
    c.pin(WARM / "admission.json")
    raw = cold_phases(c)
    primary, target, rss = cold_matrix(c, raw)
    rss["baseline"].pop("limit_status")
    receipt = {"schema": "current-cold-scientific-handoff-v1", "scientific_completion": True,
               "observed_cold_supervisor_exit": c.obj(OBS)["phases"]["full-cold-matrix"]["observed_supervisor_exit"],
               "root_observation": c.pin(OBS), "runner_raw_exit": raw, "runner_target_status": target, "candidate_8GiB_status": rss["candidate"]["limit_status"],
               "rss": rss, "primary": primary, "cold_root": root_identity(COLD), "warm_root": root, "pending": OPEN,
               "adapter": d.pin_control(Path(__file__)), "sources": {str(p): h for p, h in SOURCES.items()},
               "controls": c.files}
    for p, h in c.files.items():
        d.check_control(pin(p, h))
    r.write_json(PACKET, receipt)
    print(d.encode(d.pin_control(PACKET)).decode(), end="")


def prior_warm_phase(controls, label, args):
    receipt = controls.pin(WARM / (label + ".jsonl"))
    planned, _, _ = d.closed_phase(receipt, WARM, label, [0])
    expected = [str(PYTHON), "-I", "-B", "-S", str(Path(__file__)), *args]
    require(planned["argv"] == expected and planned["allowed_exits"] == [0], f"Wrong warm predecessor {label}")
    for suffix in ("stdout", "stderr"):
        path = WARM / (label + "." + suffix)
        require(planned[suffix] == str(path), f"Wrong predecessor log {label}")
        controls.pin(path)


def scale_action(args, root):
    handoff = pin(PACKET, args.handoff_sha256)
    receipt = json.loads(d.check_control(handoff))
    require(receipt["schema"] == "current-cold-scientific-handoff-v1" and receipt["scientific_completion"] is True
            and receipt["adapter"] == d.pin_control(Path(__file__)) and receipt["cold_root"] == root_identity(COLD)
            and receipt["warm_root"] == root and receipt["pending"] == OPEN
            and receipt["sources"] == {str(p): h for p, h in SOURCES.items()}, "Handoff provenance differs")
    for path, sha in receipt["controls"].items():
        d.check_control(pin(path, sha))
    c = Controls()
    prior_warm_phase(c, "seal-cold", ["seal"])
    require(c.obj(WARM / "seal-cold.stdout") == handoff, "Guarded eligibility output differs from the pinned handoff")
    if args.action == "verify":
        prior_warm_phase(c, f"warm-run-{args.size}", ["measure", args.size, "--handoff-sha256", args.handoff_sha256])
    primary = receipt["primary"]
    cfg = json.loads(d.check_control(primary["matrix_config"]))
    deps = dependencies(cfg)
    scale = WARM / f"matched-{args.size}"
    call = (r, primary, WARM, COLD, args.size, WORKER, PRIME)
    if args.action == "measure":
        status = w.measure(*call)
        for path, sha in SOURCES.items():
            d.check_control(pin(path, sha))
        return status
    require(scale.resolve(strict=True) == scale and not (scale / "verified-outputs.json").exists(),
            "Expected physical unverified scale")
    # Keep the exact verifier's JSON, including its stdout bytes. Only publish the
    # supplementary pin list after every process, path and evidence check succeeds.
    stdout = io.StringIO()
    with redirect_stdout(stdout):
        status = w.verify(*call)
    text = stdout.getvalue()
    print(text, end="", flush=True)
    require(status == 0, f"Frozen verification failed {status}")
    result = json.loads(text)
    with (scale / "verifier.stdout.json").open("xb") as stream:
        stream.write(text.encode())
    for name in ("inputs-before.json", "inputs-after.json", "warm-worker.R", "verifier.stdout.json"):
        c.pin(scale / name)
    identity = c.obj(scale / "inputs-after.json")
    cfg["cases"] = [case for case in cfg["cases"] if case["id"] == scale.name]
    r.resolve_configuration(cfg, "diagnostic")
    case, = cfg["cases"]
    assets = identity["assets"]
    dirs = w.edition_roots(COLD, WARM, args.size)
    verified, _ = checked_artifacts(c, cfg, case, dirs, scale, assets, rehash=False)
    require(result["verification"] == verified and result["size"] == args.size and result["pairs"] == 1
            and result["independent_oracle_transferred_by_byte_parity"] == primary["snp_oracles"][args.size]
            and result["inputs_before"] == c.pin(scale / "inputs-before.json")
            and result["inputs_after"] == c.pin(scale / "inputs-after.json"), "Frozen verifier output differs")
    for role in d.ROLES:
        dest = scale / role
        _, process, metrics, _ = checked_run(c, cfg, case, dest, role, assets, deps, identity["toolchain"], scale / "warm-worker.R")
        require(result["roles"][role]["process"] == process and result["roles"][role]["warm_metrics"] == metrics
                and process["worker"] == c.pin(scale / "warm-worker.R")
                and process["application_state"] == "same_verified_handle_after_complete_input_prime"
                and process["process_peak_rss_includes_prime"] is True and process["os_page_cache"] == "uncontrolled",
                f"Warm role process {dest}")
        require(result["roles"][role]["prime_metrics"] == w.settings(c.pin(dest / "prime/metrics.tsv"))
                and result["roles"][role]["cold_comparator"] == primary["scale_controls"][args.size][role],
                f"Warm prime/comparator differs {dest}")
        if role == "candidate":
            require(result["candidate_8GiB_status"] == d.summarize_rss([process])["limit_status"],
                    "Candidate RSS evaluation differs")
        for filename in (*ROLE_LOGS, "prime/metrics.tsv"):
            c.pin(dest / filename)
    dependencies(cfg)
    for path, sha in SOURCES.items():
        d.check_control(pin(path, sha))
    r.write_json(scale / "verified-outputs.json", {"scope": "verified scale outputs pending root phase closure",
                 "size": args.size, "handoff": handoff, "files": c.files, "adapter": receipt["adapter"],
                 "sources": receipt["sources"], "candidate_8GiB_status": result["candidate_8GiB_status"]})
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("seal")
    for action in ("measure", "verify"):
        p = sub.add_parser(action)
        p.add_argument("size", choices=("1M", "10M", "40M"))
        p.add_argument("--handoff-sha256", required=True)
    args = parser.parse_args()
    os.umask(0o077)
    modules()
    if args.action == "seal":
        seal(dispatch("seal-cold"))
        return 0
    label = "warm-run-" if args.action == "measure" else "warm-verify-"
    return scale_action(args, dispatch(label + args.size))


if __name__ == "__main__":
    sys.exit(main())
