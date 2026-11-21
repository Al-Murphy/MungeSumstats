"""Produce guarded normal staged Darwin builds and namespace benchmark provenance."""
import argparse
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shlex
import stat
import sys
import time


RUNTIME = Path("/private/tmp/mss-fullfile-20261121/runtime")
SOURCE = RUNTIME / "namespace-builds.py"
BINDING = RUNTIME / "namespace-bench.py"
BINDING_SHA = "a61cbab2a7af12e0befd8d4fb6c0e9f81d786160db1ad391c95389717a002e94"
GUARD_SHA = "58e0818aa20801c30194bdce26bd117112b4b4d5d10f74ea556360c0867b56a2"
PYTHON = Path("/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15")
REPO = Path("/private/tmp/mss-lazy-variantannotation-20261121")
COMMITS = {"baseline": "14ea7af51a19272c7792fbe582a3d42b291b4a1a",
           "candidate": "f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b"}
CLT = Path("/Library/Developer/CommandLineTools/usr/bin")
FLOOR, CAP, RSS, WALL = 186 * 1024**3, 512 * 1024**2, 4 * 1024**3, 900
ESTIMATED_OUTPUT = 3 * 1024**3
CONTROL_LIMIT = 1024**2
REQUIRED_FLAGS = ("CC", "CFLAGS", "CPPFLAGS", "CXX17", "CXX17FLAGS", "LDFLAGS")
OPTIONAL_FLAGS = ("CXX17STD", "SHLIB_CXX17LD", "SHLIB_CXX17LDFLAGS")
GIT_SETTINGS = {"GIT_CONFIG_NOSYSTEM": "1", "GIT_CONFIG_GLOBAL": "/dev/null",
                "GIT_NO_REPLACE_OBJECTS": "1", "GIT_NO_LAZY_FETCH": "1"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def physical(path):
    require(path.is_absolute() and path.resolve(strict=True) == path, f"Expected a physical path: {path}")
    return path


def source_bytes(path):
    physical(path)
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC), "rb") as stream:
        info = os.fstat(stream.fileno())
        require(stat.S_ISREG(info.st_mode) and 0 < info.st_size <= CONTROL_LIMIT, f"Invalid Python source: {path}")
        data = stream.read(CONTROL_LIMIT + 1)
    require(len(data) == info.st_size, f"Python source size changed: {path}")
    return data


def load_binding():
    data = source_bytes(BINDING)
    require(hashlib.sha256(data).hexdigest() == BINDING_SHA, f"Reviewed binding changed: {BINDING}")
    spec = importlib.util.spec_from_file_location("namespace_build_binding", BINDING)
    binding = importlib.util.module_from_spec(spec)
    exec(compile(data, str(BINDING), "exec"), binding.__dict__)
    return binding, binding.load_frozen()


def write_json(path, value, guard, *, check_stop=True):
    data = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    require(len(data) <= CONTROL_LIMIT, f"Oversized operator control: {path}")
    pin = {"path": str(path), "sha256": hashlib.sha256(data).hexdigest()}
    guard._check_space(guard._space(walk=True), evidence_bytes=len(data), metadata_blocks=1, new_entries=1)
    if check_stop:
        guard._stop_check()
    with guard._directory(path.parent) as parent, os.fdopen(guard._new_file(path), "wb") as stream:
        opened = os.fstat(stream.fileno())
        try:
            stream.write(data)
            stream.close()
            if check_stop:
                guard._stop_check()
        except BaseException:
            current = os.stat(path.name, dir_fd=parent, follow_symlinks=False)
            require((current.st_dev, current.st_ino) == (opened.st_dev, opened.st_ino),
                    f"Operator control changed before cleanup: {path}")
            os.unlink(path.name, dir_fd=parent)
            raise
    return pin


def environment(root, common):
    return {
        "PATH": f"{CLT}:/usr/bin:/bin:/usr/sbin:/sbin",
        "HOME": str(root / "home"),
        **{key: str(root / "tmp") for key in ("TMPDIR", "TMP", "TEMP")},
        "LANG": "en_US.UTF-8", "LC_ALL": "en_US.UTF-8", "TZ": "UTC",
        "R_LIBS": str(common), "R_LIBS_USER": "", "R_LIBS_SITE": "",
        **dict.fromkeys(("R_ENVIRON_USER", "R_PROFILE_USER", "R_MAKEVARS_USER"), "/dev/null"),
        "MAKEFLAGS": "-j1", "OMP_DYNAMIC": "FALSE", "MKL_DYNAMIC": "FALSE",
        **dict.fromkeys(("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS",
                         "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
                         "NUMEXPR_NUM_THREADS", "R_DATATABLE_NUM_THREADS"), "1"),
    }


def parse_runtime(data, expected):
    rows = iter(csv.reader(io.StringIO(data.decode()), delimiter="\t", strict=True))
    require(next(rows, None) == ["key", "value"], "Unexpected runtime TSV header")
    observed = {}
    for row in rows:
        require(len(row) == 2 and row[0] not in observed, "Invalid or repeated runtime TSV key")
        observed[row[0]] = row[1]
    require(observed == expected, "Fresh runtime differs from the pinned matrix runtime")
    return observed


def compiler_flags(r_config):
    captured = {}
    for line in r_config.splitlines():
        key, separator, value = line.partition("=")
        if separator:
            captured[key.strip()] = value.strip()
    require(all(name in captured for name in REQUIRED_FLAGS), "Incomplete R CMD config --all")
    return {name: captured[name] for name in (*REQUIRED_FLAGS, *OPTIONAL_FLAGS) if name in captured}


class BuildRun:
    def __init__(self, binding, warm, guard, result):
        self.binding, self.warm, self.guard, self.result = binding, warm, guard, result
        self.root = guard.root
        self.env = environment(self.root, binding.COMMON)
        self.sdk_env = {"SDKROOT": str(binding.PINNED_SDK)}
        self.r = str(binding.RHOME / "bin/R")

    def phase(self, label, argv, *, stem=None, archive=None, overrides=None):
        self.guard._stop_check()
        stem = self.root / label if stem is None else stem
        paths = [stem.with_suffix(".jsonl"), archive or stem.with_suffix(".stdout"), stem.with_suffix(".stderr")]
        effective = {**self.env, **(overrides or {})}
        dispatched = ["/usr/bin/env", *(f"{k}={v}" for k, v in overrides.items()), *argv] if overrides else argv
        record = {"label": label, "argv": dispatched, "command": argv, "cwd": str(self.root),
                  "environment": self.env, "effective_environment": effective, "allowed_exits": [0],
                  "paths": dict(zip(("receipt", "stdout", "stderr"), map(str, paths)))}
        self.result["phases"].append(record)
        outcome = asdict(self.guard.run_phase(label, dispatched, cwd=self.root, env=self.env, allowed_exits=(0,),
                                              receipt=paths[0], stdout=paths[1], stderr=paths[2]))
        record["outcome"] = outcome
        rows = []
        if paths[0].exists():
            evidence = self.warm.d.control_bytes(physical(paths[0]))
            record["receipt"] = {"path": str(paths[0]), "sha256": hashlib.sha256(evidence).hexdigest()}
            rows = [json.loads(line) for line in evidence.splitlines()]
            if rows and rows[-1]["status"] not in {"planned", "running"}:
                record["terminal"] = rows[-1]
        require(outcome["status"] == "completed" and outcome["receipt_error"] is None,
                f"{label}: guard outcome {outcome}")
        planned, running, terminal = rows
        require(all(row["label"] == label for row in rows)
                and planned["status"] == "planned" and running["status"] == "running"
                and planned["argv"] == dispatched and planned["cwd"] == str(self.root)
                and planned["allowed_exits"] == [0] and all(planned[key] == value for key, value in record["paths"].items())
                and running["environment_keys"] == sorted(self.env)
                and running["pid"] == running["expected_pgid"] == running["expected_sid"] > 1
                and running["uid"] == os.geteuid() and planned["utc"] <= running["utc"] <= terminal["utc"]
                and outcome["raw_exit"] == self.warm.d.phase_status(terminal, [0]),
                f"{label}: retained guard records differ from the completed phase")
        self.guard._stop_check()
        stdout_pin = self.warm.build.pin(physical(paths[1])) if archive else self.warm.d.pin_control(physical(paths[1]))
        record["transcripts"] = [record["receipt"], stdout_pin, self.warm.d.pin_control(physical(paths[2]))]
        return record

    def captures(self):
        warm, binding = self.warm, self.binding
        config_pin = warm.pin(warm.CONFIG, warm.ACCEPTED[warm.CONFIG])
        config = json.loads(warm.d.check_control(config_pin))
        previous = json.loads(warm.d.check_control(config["dependency_inventory"]))
        toolchain = {"sdk": str(physical(binding.PINNED_SDK)),
                     "build_toolchain_pins": previous["build_toolchain_pins"]}
        toolchain_pins = binding.toolchain_pins(toolchain, previous["build_toolchain_pins"])
        require(not os.path.lexists(binding.RHOME / "etc/Makevars.site"), "Unexpected R site Makevars")
        runtime_pins = [{key: pin[key] for key in ("path", "sha256")}
                        for pin in [config["r_executable"], config["r_engine"], *config["runtime_files"]]]
        for pin in [*toolchain_pins, *runtime_pins]:
            self.guard._stop_check()
            binding.checked_file(warm, pin)
        self.result.update(matrix_configuration=config_pin, historical_dependency_inventory=config["dependency_inventory"],
                           verified_toolchain_files=toolchain["build_toolchain_pins"], verified_runtime_files=runtime_pins)
        compiler = self.phase("compiler", [str(CLT / "clang"), "--version"])["transcripts"][1]
        r_config = self.phase("r-config", [self.r, "CMD", "config", "--all"],
                              overrides=self.sdk_env)["transcripts"][1]
        runtime_path = self.root / "runtime.tsv"
        command = warm.r.r_command(config["r_executable"]["path"], warm.WORKER, "runtime", str(runtime_path))
        runtime = self.phase("runtime", command)["transcripts"]
        runtime_receipt = warm.d.pin_control(physical(runtime_path))
        observed = parse_runtime(warm.d.check_control(runtime_receipt), previous["runtime_observed"])
        self.result["runtime_receipt"] = runtime_receipt
        compiler_lines = warm.d.check_control(compiler).decode().strip().splitlines()
        require(compiler_lines and compiler_lines[0].startswith("Apple clang version "), "Unexpected compiler version capture")
        config_lines = warm.d.check_control(r_config).decode()
        toolchain.update(R_version=observed["R_version"], R_platform=observed["R_platform"],
                         compiler=compiler_lines, compiler_flags=compiler_flags(config_lines),
                         compiler_receipt=compiler, r_config_receipt=r_config)
        binding.verify_compiler_capture(toolchain, compiler_lines, config_lines.splitlines())
        warm.r.verify_runtime(config, {**observed, "threads": "1"}, toolchain)
        return toolchain, runtime_receipt, runtime

    def role(self, role, commit, toolchain):
        warm = self.warm
        self.guard._stop_check()
        directory = self.root / role
        directory.mkdir(mode=0o700)
        library = directory / "library"
        library.mkdir(mode=0o700)
        git = ["/usr/bin/git", "--no-replace-objects", "-C", str(REPO)]
        revision = self.phase(f"source-tree-{role}", [*git, "rev-parse", "--verify", f"{commit}^{{tree}}"],
                              stem=directory / "source-tree", overrides=GIT_SETTINGS)
        tree = warm.d.check_control(revision["transcripts"][1]).decode().strip()
        require(warm.r.valid_hash(tree, 40), f"{role}: invalid Git tree ID")
        archive_path = directory / "MungeSumstats.tar"
        archive_phase = self.phase(f"git-archive-{role}", [*git, "archive", "--format=tar", "--prefix=MungeSumstats/", commit],
                                   stem=directory / "archive", archive=archive_path, overrides=GIT_SETTINGS)
        archive = archive_phase["transcripts"][1]
        install = [self.r, "CMD", "INSTALL", "--use-vanilla", f"--library={library}", str(archive_path)]
        phase = self.phase(f"normal-install-{role}", install, stem=directory / "install", overrides=self.sdk_env)
        stdout, stderr = (warm.d.check_control(pin).decode() for pin in phase["transcripts"][1:])
        commands = [shlex.split(line) for line in stdout.splitlines() if line.startswith("clang++ ")]
        require(len(commands) == 3 and all(f"-I{self.binding.COMMON}/Rcpp/include" in c and "-c" in c for c in commands[:2])
                and "-dynamiclib" in commands[2] and "MungeSumstats.so" in commands[2],
                f"{role}: unexpected compiler or linker commands")
        markers = ("** using staged installation", toolchain["compiler"][0], "using C++17",
                   f"using SDK: ‘{self.binding.PINNED_SDK.name}’",
                   "** testing if installed package can be loaded from temporary location",
                   "** testing if installed package can be loaded from final location", "* DONE (MungeSumstats)")
        offsets = [stderr.find(marker) for marker in markers]
        require(offsets[0] >= 0 and offsets == sorted(set(offsets)), f"{role}: incomplete normal staged install evidence")
        package = library / "MungeSumstats"
        require(set(library.iterdir()) == {package}, f"{role}: library must contain only the newly installed package")
        dll = warm.build.pin(physical(package / "libs/MungeSumstats.so"))
        tree_sha = warm.build.package_tree(physical(package))["tree_sha256"]
        self.guard._stop_check()
        receipt = {"source_commit": commit, "source_tree": tree, "archive": archive, "library": str(physical(library)),
                   "dll": dll, "dll_sha256": dll["sha256"], "package_tree_sha256": tree_sha,
                   "install_argv": phase["argv"], "install_transcripts": phase["transcripts"], "terminal": phase["terminal"],
                   "compile_and_link_argv": commands, "toolchain": toolchain}
        self.result["roles"][role] = write_json(directory / "build.json", receipt, self.guard)
        return receipt


def run(root):
    started = datetime.now(timezone.utc)
    require(sys.platform == "darwin" and physical(Path(sys.executable)) == PYTHON and sys.flags.isolated
            and sys.flags.dont_write_bytecode and sys.flags.no_site, "Use the pinned Darwin Python with -I -B -S")
    require(root.is_absolute() and root == physical(root.parent) / root.name
            and not root.is_relative_to(REPO), "Expected an absolute physical output path outside the source repository")
    require(physical(Path(__file__)) == SOURCE, f"Use the physical producer source: {SOURCE}")
    os.umask(0o077)
    root.mkdir(mode=0o700)
    info = root.lstat()
    require(stat.S_IMODE(info.st_mode) == 0o700 and info.st_uid == os.geteuid(), "Expected a private owner output root")
    source = {"path": str(SOURCE), "sha256": hashlib.sha256(source_bytes(SOURCE)).hexdigest()}
    binding, warm = load_binding()
    require(warm.SOURCES[warm.GUARD] == GUARD_SHA and binding.BASELINE == COMMITS["baseline"],
            "Unexpected authenticated guard or baseline")
    guard_pin = warm.pin(warm.GUARD, GUARD_SHA)
    native = warm.d.module(guard_pin, "namespace_build_guard")
    guard = native.Guard(native.GuardConfig(root, started, WALL, RSS, CAP, FLOOR))
    require(guard._root_id == (info.st_dev, info.st_ino), "Output root changed since creation")
    with guard.shutdown_signals():
        available = guard._space(walk=False)["free_bytes"]
        admitted = available >= FLOOR + CAP
        if not admitted:
            print(f"Available bytes {available} cannot fund the continuing floor and build output cap", file=sys.stderr)
        admission = write_json(root / "admission.json", {
            "scope": "matched normal installs and runtime metadata only", "admitted": admitted,
            "issued_utc": started.isoformat(), "observed_utc": datetime.now(timezone.utc).isoformat(),
            "root": {"path": str(root), "dev": info.st_dev, "ino": info.st_ino}, "supervisor_pid": os.getpid(),
            "producer": source, "binding": {"path": str(BINDING), "sha256": BINDING_SHA}, "guard_sha256": GUARD_SHA,
            "available_bytes": available, "minimum_free_bytes": FLOOR, "owned_output_cap_bytes": CAP,
            "sampled_group_rss_limit_bytes": RSS, "wall_seconds": WALL,
        }, guard)
        if not admitted:
            return 1
        result = {"status": "failed", "producer": source, "admission": admission, "repository": str(REPO),
                  "guard": guard_pin, "commits": COMMITS, "phases": [], "roles": {}}
        try:
            physical(REPO)
            common = physical(binding.COMMON)
            require(not os.path.lexists(common / "MungeSumstats"), "Common library contains a selected package")
            for name in ("home", "tmp"):
                (root / name).mkdir(mode=0o700)
            operator = BuildRun(binding, warm, guard, result)
            toolchain, runtime, transcripts = operator.captures()
            baseline, candidate = (operator.role(role, commit, toolchain) for role, commit in COMMITS.items())
            require(baseline["source_tree"] != candidate["source_tree"]
                    and baseline["compile_and_link_argv"] == candidate["compile_and_link_argv"],
                    "Expected distinct source trees with identical native commands")
            guard._stop_check()
            result["provenance"] = write_json(root / "PROVENANCE.json", {
                "repository": str(REPO), **result["roles"], "runtime_receipt": runtime, "runtime_transcripts": transcripts,
                "minimum_free_bytes": FLOOR, "estimated_output_bytes": ESTIMATED_OUTPUT,
                "contention_note": (f"Independent build/metadata admission observed {available} free bytes. "
                                    "Build phases ran sequentially under the shared guard. "
                                    "Benchmark host admission, contention and cache observations, and timings are pending."),
            }, guard)
            completed = {**result, "status": "completed", "finished_utc": datetime.now(timezone.utc).isoformat()}
            write_json(root / "operator-result.json", completed, guard)
            return 0
        except (OSError, ValueError, RuntimeError, KeyError, TypeError, csv.Error, KeyboardInterrupt) as error:
            result.update(error=f"{type(error).__name__}: {error}", finished_utc=datetime.now(timezone.utc).isoformat())
            print(result["error"], file=sys.stderr)
            try:
                write_json(root / "operator-result.json", result, guard, check_stop=False)
            except (OSError, ValueError, RuntimeError) as receipt_error:
                print(f"Could not retain operator result: {receipt_error}", file=sys.stderr)
            return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="absent physical directory for private build evidence")
    try:
        return run(parser.parse_args().output)
    except (OSError, ValueError, RuntimeError) as error:
        print(error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
