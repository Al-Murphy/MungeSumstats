"""Bind the matched Darwin namespace experiment to the frozen diagnostic runner."""
import argparse
import copy
import csv
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys


BOOTSTRAP = Path(__file__).parent / "warm-current.py"
BOOTSTRAP_SHA = "c58b14cda2ad650cefa078e4a73cd2a2a4d5ab8ffc4f7477e5375e0aaa794220"
BASELINE = "14ea7af51a19272c7792fbe582a3d42b291b4a1a"
COMMON = Path("/private/tmp/mss-parallel-support-20261121/common")
COMMON_PIN = {"path": str(COMMON.parent / "inventory.json"),
              "sha256": "181de043a869c040e467039cd533d936efe3a69734d7c9c27c6a2a80d4e84472"}
RHOME = Path("/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources")
FRAMEWORK = RHOME / "library"
PINNED_SDK = Path("/Library/Developer/CommandLineTools/SDKs/MacOSX27.2.sdk")
ROLES = ("baseline", "candidate")
CASE_NAMES = ("hard-cases", "declared", "matched-1M")
CASE_IDS = ("hard-cases-f2", "hard-cases-f3", "declared-f2", "declared-f3", "matched-1M")
HISTORICAL_BASELINE = "a86e8c94d753630eff8c6a695ed2f240d1123d0d"
EAGER_ONLY = frozenset({"AnnotationDbi", "DBI", "GenomicFeatures", "KEGGREST", "RSQLite",
                       "VariantAnnotation", "bit", "bit64", "blob", "cachem", "fastmap", "memoise", "png"})
SELECTED = "<selected installation>"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def physical(path):
    path = Path(path)
    require(path.is_absolute() and path.resolve(strict=True) == path, f"Expected physical path: {path}")
    return path


def load_frozen():
    with physical(BOOTSTRAP).open("rb") as stream:
        data = stream.read(1024**2 + 1)
    require(len(data) <= 1024**2 and hashlib.sha256(data).hexdigest() == BOOTSTRAP_SHA,
            "Frozen warm bootstrap digest changed")
    spec = importlib.util.spec_from_file_location("namespace_frozen_bootstrap", BOOTSTRAP)
    warm = importlib.util.module_from_spec(spec)
    exec(compile(data, str(BOOTSTRAP), "exec"), warm.__dict__)
    require(Path(sys.executable) == warm.PYTHON, "Use the physical pinned Python executable")
    warm.modules()
    require(warm.r.BASELINE == HISTORICAL_BASELINE, "Frozen baseline differs")
    return warm


class NamespacePolicy:
    def __init__(self, installed):
        self.installed = installed
        self.observed = {}

    def observe(self, case, role, state):
        for package, identity in state.items():
            require(identity == self.installed[role].get(package),
                    f"{case}/{role}: wrong version or origin for {package}: {identity}")
        require("MungeSumstats" in state, f"{case}/{role}: missing selected installation")
        previous = self.observed.get((case, role))
        require(previous is None or previous == state, f"{case}/{role}: changed across repetitions")
        if role == "baseline":
            require(EAGER_ONLY <= state.keys(), f"{case}: baseline missing {sorted(EAGER_ONLY - state.keys())}")
        else:
            baseline = self.observed.get((case, "baseline"))
            require(baseline is not None, f"{case}: baseline must be observed before candidate")
            require(not (state.keys() & EAGER_ONLY),
                    f"{case}: candidate loaded eager-only {sorted(state.keys() & EAGER_ONLY)}")
            require(state.keys() <= baseline.keys(),
                    f"{case}: candidate-only namespaces {sorted(state.keys() - baseline.keys())}")
            require(baseline.keys() - state.keys() == EAGER_ONLY,
                    f"{case}: additional baseline-only namespaces {sorted(baseline.keys() - state.keys() - EAGER_ONLY)}")
            require(all(identity == baseline[name] for name, identity in state.items()),
                    f"{case}: shared namespace metadata differs")
        self.observed[case, role] = {package: list(identity) for package, identity in state.items()}
        return {package: identity for package, identity in state.items() if package not in EAGER_ONLY}


def role_commit(warm, commit):
    require(warm.r.valid_hash(commit, 40) and len(set(commit)) > 1 and commit != BASELINE,
            "The actual distinct candidate Git commit (40 lowercase hex characters) is required")
    return commit


def checked_file(warm, pin):
    path = physical(pin["path"])
    require(warm.build.pin(path) == pin, f"File pin changed: {path}")
    return path


def guarded_phase(warm, pins, argv, label):
    require(len(pins) == 3, f"{label}: require actual guard JSONL, stdout and stderr")
    guard, stdout, stderr = [warm.d.check_control(pin) for pin in pins]
    planned, running, terminal = [json.loads(line) for line in guard.splitlines()]
    require(all(row["label"] == label for row in (planned, running, terminal))
            and planned["status"] == "planned" and running["status"] == "running"
            and planned["argv"] == argv and planned["allowed_exits"] == [0]
            and [planned[key] for key in ("receipt", "stdout", "stderr")] == [pin["path"] for pin in pins]
            and running["pid"] == running["expected_pgid"] == running["expected_sid"] > 1
            and running["uid"] == os.geteuid() and planned["utc"] <= running["utc"] <= terminal["utc"],
            f"{label}: guard phase does not bind the actual command and records")
    warm.d.phase_status(terminal, [0])
    return planned, running, terminal, stdout.decode(), stderr.decode()


def git(repo, *args, stdin=None):
    return subprocess.run(["/usr/bin/git", "--no-replace-objects", "-C", str(repo), *args],
                          stdin=stdin, stdout=subprocess.PIPE, check=True, text=True, timeout=15,
                          env={"PATH": "/usr/bin:/bin", "LC_ALL": "C", "GIT_NO_LAZY_FETCH": "1",
                               "GIT_CONFIG_NOSYSTEM": "1", "GIT_CONFIG_GLOBAL": "/dev/null"}).stdout.strip()


def verified_role(warm, role, commit, receipt_pin, repo):
    receipt = json.loads(warm.d.check_control(receipt_pin))
    require(receipt["source_commit"] == commit and git(repo, "rev-parse", "--verify", f"{commit}^{{commit}}") == commit
            and git(repo, "rev-parse", "--verify", f"{commit}^{{tree}}") == receipt["source_tree"],
            f"{role}: build source does not match the actual Git commit and tree")
    archive = checked_file(warm, receipt["archive"])
    with archive.open("rb") as stream:
        require(git(repo, "get-tar-commit-id", stdin=stream) == commit, f"{role}: wrong Git archive commit")
    library = physical(receipt["library"])
    package = library / "MungeSumstats"
    require(set(library.iterdir()) == {package} and physical(package) == package,
            f"{role}: expected a separate, clean normal install")
    dll = checked_file(warm, receipt["dll"])
    require(dll == package / "libs/MungeSumstats.so" and receipt["dll_sha256"] == receipt["dll"]["sha256"],
            f"{role}: wrong installed DLL")
    require(warm.r.valid_hash(receipt["package_tree_sha256"]), f"{role}: missing installation digest")
    toolchain = receipt["toolchain"]
    sdk = physical(toolchain["sdk"])
    argv = receipt["install_argv"]
    install = ["/usr/bin/env", f"SDKROOT={sdk}", str(RHOME / "bin/R"), "CMD", "INSTALL",
               "--use-vanilla", f"--library={library}", str(archive)]
    require(argv in (install, ["/usr/bin/time", "-l", *install]), f"{role}: expected normal staged R installation")
    _, running, terminal, stdout, stderr = guarded_phase(warm, receipt["install_transcripts"], argv,
                                                         f"normal-install-{role}")
    markers = ("** using staged installation", toolchain["compiler"][0], "using C++17", f"using SDK: ‘{sdk.name}’",
               "** testing if installed package can be loaded from temporary location",
               "** testing if installed package can be loaded from final location", "* DONE (MungeSumstats)")
    offsets = [stderr.find(marker) for marker in markers]
    require(offsets == sorted(set(offsets)) and offsets[0] >= 0 and terminal == receipt["terminal"],
            f"{role}: incomplete guarded build, temporary or final load check")
    commands = [shlex.split(line) for line in stdout.splitlines() if line.startswith("clang++ ")]
    require(len(commands) == 3 and commands == receipt["compile_and_link_argv"]
            and all(f"-I{COMMON}/Rcpp/include" in command and "-c" in command for command in commands[:2])
            and "-dynamiclib" in commands[2] and "MungeSumstats.so" in commands[2],
            f"{role}: missing or inconsistent compile/link argv against the common library")
    config = {"commit": commit, "library": str(library), "dll": receipt["dll"],
              "package_tree_sha256": receipt["package_tree_sha256"], "build_receipt": receipt_pin}
    return config, receipt, set(running["environment_keys"])


def select_config(original):
    cases = {case["id"]: case for case in original["cases"]}
    return {**original,
            "references": {name: original["references"][name] for name in ("matched_2", "scientific_2", "scientific_3")},
            "cases": [{**cases[name], "pairs": 3, "target": False} for name in CASE_NAMES]}


def toolchain_pins(toolchain, historical_pins):
    pins = {pin["path"]: pin for pin in toolchain["build_toolchain_pins"]}
    historical = {pin["path"]: pin for pin in historical_pins}
    settings_path = str(PINNED_SDK / "SDKSettings.json")
    settings_pin = historical.pop(settings_path, None)
    require(toolchain["sdk"] == str(PINNED_SDK) and settings_pin is not None
            and pins.get(settings_path) == settings_pin
            and all(pins.get(path) == pin for path, pin in historical.items()),
            "Expected the pinned compiler, SDK and R Makeconf")
    return pins.values()


def verify_compiler_capture(toolchain, compiler, r_config):
    flags = {part[0].strip(): part[2].strip() for line in r_config if (part := line.partition("="))[1]}
    required = {"CC", "CFLAGS", "CPPFLAGS", "CXX17", "CXX17FLAGS", "LDFLAGS"}
    declared = toolchain["compiler_flags"]
    require(compiler == toolchain["compiler"] and required <= declared.keys() <= flags.keys()
            and all(isinstance(value, str) and flags[key] == value for key, value in declared.items()),
            "Compiler version or effective R flags disagree with actual captures")


def inputs(warm, provenance_pin, candidate):
    root = json.loads(warm.d.check_control(provenance_pin))
    config = select_config(json.loads(warm.d.check_control(warm.pin(warm.CONFIG, warm.ACCEPTED[warm.CONFIG]))))
    previous = json.loads(warm.d.check_control(config["dependency_inventory"]))
    common = json.loads(warm.d.check_control(COMMON_PIN))
    inherited, framework = (previous["installed_dependency_tree_sha256"][str(path)]
                            for path in (Path("/private/tmp/mss-fullfile-20261121/restored-01/common"), FRAMEWORK))
    require(common["common"] == str(COMMON) and len(common["packages"]) == 96 and len(inherited) == 94
            and {name: digest for name, digest in common["packages"].items() if name in inherited} == inherited
            and common["packages"].keys() - inherited.keys() == {"BiocBaseUtils", "GenomicFiles"}
            and len(framework) == 30, "Common/framework archive inventory differs from the accepted controls")
    expected = {str(COMMON): common["packages"], str(FRAMEWORK): framework}
    repo = physical(root["repository"])
    builds, environments = {}, {}
    for role, commit in zip(ROLES, (BASELINE, candidate)):
        config[role], builds[role], environments[role] = verified_role(warm, role, commit, root[role], repo)
    a, b = (builds[role] for role in ROLES)
    toolchain = a["toolchain"]
    require(a["archive"]["path"] != b["archive"]["path"] and a["source_tree"] != b["source_tree"]
            and a["install_argv"][:-2] == b["install_argv"][:-2]
            and a["compile_and_link_argv"] == b["compile_and_link_argv"] and toolchain == b["toolchain"]
            and environments["baseline"] == environments["candidate"]
            and {"R_LIBS", "R_LIBS_SITE", "R_LIBS_USER", "R_MAKEVARS_USER", "MAKEFLAGS"} <= environments["baseline"],
            "Role archives, guard environments, commands or toolchain are not comparable")
    pins = toolchain_pins(toolchain, previous["build_toolchain_pins"])
    require(not (RHOME / "etc/Makevars.site").exists(), "Unexpected R site Makevars")
    for pin in [*pins, config["r_executable"], config["r_engine"],
                *({key: p[key] for key in ("path", "sha256")} for p in config["runtime_files"])]:
        checked_file(warm, pin)
    compiler = warm.d.check_control(toolchain["compiler_receipt"]).decode().strip().splitlines()
    r_config = warm.d.check_control(toolchain["r_config_receipt"]).decode().splitlines()
    verify_compiler_capture(toolchain, compiler, r_config)
    runtime_path = root["runtime_receipt"]["path"]
    runtime = warm.d.check_control(root["runtime_receipt"]).decode()
    guarded_phase(warm, root["runtime_transcripts"], warm.r.r_command(config["r_executable"]["path"], warm.WORKER,
                  "runtime", runtime_path), "runtime")
    header, *records = csv.reader(io.StringIO(runtime), delimiter="\t", strict=True)
    observed = dict(records)
    require(header == ["key", "value"] and len(observed) == len(records) and observed == previous["runtime_observed"],
            "Fresh frozen-worker runtime differs from the pinned matrix runtime")
    warm.r.verify_runtime(config, {**observed, "threads": "1"}, toolchain)
    config.update(dependency_libraries=[str(COMMON)], dependency_inventory=COMMON_PIN,
                  **{name: root[name] for name in ("contention_note", "minimum_free_bytes", "estimated_output_bytes")})
    resolved = copy.deepcopy(config)
    warm.r.resolve_configuration(resolved, "diagnostic")
    require(tuple(case["id"] for case in resolved["cases"]) == CASE_IDS
            and all(case["pairs"] == 3 and case["target"] is False for case in resolved["cases"]),
            "Unexpected diagnostic case or pair roster")
    return config, expected


def snapshot(warm, config, expected):
    libraries = {path: warm.build.inventory_library(physical(path)) for path in expected}
    require(libraries == expected, "Current common or framework package inventory differs")
    roles = {role: warm.build.package_tree(Path(config[role]["library"]) / "MungeSumstats") for role in ROLES}
    require(all(tree["tree_sha256"] == config[role]["package_tree_sha256"] for role, tree in roles.items()),
            "Current installed role package tree differs")
    return {"installed_dependency_tree_sha256": libraries, "role_trees": roles}


def namespace_identities(warm, config, libraries):
    def identity(path):
        desc = warm.d.control_bytes(path / "DESCRIPTION").decode("latin-1")
        versions = re.findall(r"^Version:[ \t]+([0-9]+(?:[.-][0-9]+)+)[ \t]*$", desc, re.MULTILINE)
        require(len(versions) == 1, f"Missing or invalid installed DESCRIPTION Version: {path}")
        return [".".join(str(int(part)) for part in re.split(r"[.-]", versions[0])), str(physical(path))]

    shared = {}
    for root, names in libraries.items():
        require(shared.keys().isdisjoint(names), f"Package shadowing in library: {root}")
        shared.update({name: identity(Path(root) / name) for name in names})
    require("MungeSumstats" not in shared, "MungeSumstats shadows a role installation")
    return {role: {**shared, "MungeSumstats": identity(Path(config[role]["library"]) / "MungeSumstats")}
            for role in ROLES}


def source_identity(warm):
    return {"binding": warm.d.pin_control(physical(__file__)), "bootstrap": warm.pin(BOOTSTRAP, BOOTSTRAP_SHA),
            "frozen_sources": {str(path): digest for path, digest in warm.SOURCES.items()},
            "matrix_config": warm.pin(warm.CONFIG, warm.ACCEPTED[warm.CONFIG]),
            "baseline_override": {"frozen": HISTORICAL_BASELINE, "actual": BASELINE}}


def prepare(warm, args, candidate):
    output = args.output
    if not os.path.lexists(output):
        physical(output.parent)
        output.mkdir(mode=0o700)
    warm.root_identity(output)
    require(not any(output.iterdir()), "Expected a fresh metadata directory")
    provenance = warm.d.pin_control(physical(args.provenance))
    config, expected = inputs(warm, provenance, candidate)
    inventory_path = output / "preparation-inventory.json"
    warm.r.write_json(inventory_path, snapshot(warm, config, expected))
    config_path = output / "diagnostic.json"
    warm.r.write_json(config_path, config)
    warm.r.write_json(output / "binding.json", {**source_identity(warm), "provenance": provenance,
                      "candidate": candidate, "configuration": warm.d.pin_control(config_path),
                      "preparation_inventory": warm.d.pin_control(inventory_path)})
    print(config_path)
    return 0


def run(warm, args):
    config_path = physical(args.configuration)
    warm.root_identity(config_path.parent)
    output = args.output
    warm.root_identity(output.parent)
    require(not os.path.lexists(output) and config_path.parent not in output.parents,
            "Use a new measurement output under a separate private root")
    binding_path = physical(config_path.parent / "binding.json")
    binding_pin = warm.d.pin_control(binding_path)
    binding = json.loads(warm.d.check_control(binding_pin))
    identity = source_identity(warm)
    require(binding["configuration"]["path"] == str(config_path)
            and all(binding[key] == value for key, value in identity.items()),
            "Measurement config or binding sources differ from preparation")
    config, expected = inputs(warm, binding["provenance"], role_commit(warm, binding["candidate"]))
    require(config == json.loads(warm.d.check_control(binding["configuration"])), "Prepared configuration differs")
    before = snapshot(warm, config, expected)
    require(before == json.loads(warm.d.check_control(binding["preparation_inventory"])),
            "Preflight inventory differs from the prepared role installs")
    installed = namespace_identities(warm, config, expected)
    policy = NamespacePolicy({role: {**names, "MungeSumstats": [names["MungeSumstats"][0], SELECTED]}
                              for role, names in installed.items()})
    original_state, original_argv = warm.r.dependency_state, sys.argv
    observations = []

    def dependency_state(run_dir, cfg, role):
        raw = warm.d.pin_control(run_dir / "namespaces.tsv")
        state = original_state(run_dir, cfg, role)
        records = csv.DictReader(io.StringIO(warm.d.check_control(raw).decode()), delimiter="\t", strict=True)
        raw_names = set()
        for row in records:
            name = row["package"]
            require(name not in raw_names and None not in row
                    and installed[role].get(name) == [row["version"], row["path"]]
                    and state.get(name) == [row["version"], SELECTED if name == "MungeSumstats" else row["path"]],
                    f"Nonphysical or unexpected raw namespace metadata: {run_dir}/{name}")
            raw_names.add(name)
        require(raw_names == state.keys(), f"Raw and frozen-parser namespace names differ: {run_dir}")
        case = run_dir.parent.name
        comparison = policy.observe(case, role, state)
        observations.append({"case": case, "role": role, "worker": run_dir.name, "raw_inventory": raw})
        return comparison

    try:
        warm.r.dependency_state = dependency_state
        sys.argv = [warm.r.__file__, "run", str(config_path), str(output), "--suite", "diagnostic", "--profile", "none"]
        status = warm.r.main()
    finally:
        warm.r.dependency_state, sys.argv = original_state, original_argv
    if status != 0:
        return status
    summary = warm.r.read_json(output / "summary.json")
    expected_workers = [(case, role, f"pair-{pair:02d}-{role}") for case in CASE_IDS for pair in (1, 2, 3)
                        for role in (ROLES if pair % 2 else ROLES[::-1])]
    require(summary["target_status"] == "not_applicable" and summary["suite"] == "diagnostic"
            and summary["profile"] == "none" and list(summary["cases"]) == list(CASE_IDS)
            and [(r["case"], r["role"], r["worker"]) for r in observations] == expected_workers,
            "Incomplete or unexpected diagnostic observations")
    for name, result in summary["cases"].items():
        verification = result["verification"]
        require(verification["stable_byte_equality"] is True and verification["stable_rds_manifest_equality"] is True
                and verification["edition_count"] == 6 and verification["normalized_rows_checked"] == result["rows"]
                and (name == "matched-1M" or (verification["oracle_rows_checked"] == result["rows"]
                     and set(verification["oracle_fields"]) == set(warm.v.NORMALIZED_COLUMNS))),
                f"Incomplete edition or literal gold parity: {name}")
        require(len(result["observations"]) == 6 and all(item["process"]["rss_complete"] is True
                    and item["process"]["rss_source"] == "wait4 ru_maxrss"
                    and item["process"]["peak_rss_bytes"] > 0 for item in result["observations"]),
                f"Incomplete Darwin RSS observations: {name}")
    for worker in observations:
        warm.d.check_control(worker["raw_inventory"])
    after = snapshot(warm, config, expected)
    require(before == after and installed == namespace_identities(warm, config, expected),
            "Namespace installation inventories changed during measurement")
    warm.r.write_json(output / "namespace-inputs-before.json", before)
    warm.r.write_json(output / "namespace-inputs-after.json", after)
    cases = {}
    for case in CASE_IDS:
        loaded = {role: policy.observed[case, role] for role in ROLES}
        cases[case] = {"baseline_only": sorted(loaded["baseline"].keys() - loaded["candidate"].keys()),
                       "candidate_only": sorted(loaded["candidate"].keys() - loaded["baseline"].keys()),
                       "common_comparison_state": loaded["candidate"],
                       "actual_role_states": {role: {**loaded[role], "MungeSumstats": installed[role]["MungeSumstats"]}
                                              for role in ROLES}}
    for pin in [binding_pin, binding["configuration"], binding["provenance"], identity["bootstrap"],
                *(warm.pin(path, digest) for path, digest in warm.SOURCES.items())]:
        warm.d.check_control(pin)
    require(source_identity(warm) == identity, "Binding source changed during measurement")
    warm.r.write_json(output / "namespace-policy.json", {**identity,
                      "preparation": binding_pin, "configuration": binding["configuration"],
                      "source_commits": {role: config[role]["commit"] for role in ROLES},
                      "required_baseline_only": sorted(EAGER_ONLY), "cases": cases, "workers": observations})
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    p = commands.add_parser("prepare")
    p.add_argument("provenance", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--candidate", required=True)
    r = commands.add_parser("run")
    r.add_argument("configuration", type=Path)
    r.add_argument("output", type=Path)
    r.add_argument("--suite", choices=("diagnostic",), default="diagnostic")
    r.add_argument("--profile", choices=("none",), default="none")
    args = parser.parse_args()
    os.umask(0o077)
    try:
        warm = load_frozen()
        warm.r.BASELINE = BASELINE
        if args.action == "prepare":
            return prepare(warm, args, role_commit(warm, args.candidate))
        return run(warm, args)
    except (OSError, ValueError, KeyError, TypeError, IndexError, csv.Error, subprocess.SubprocessError) as error:
        print(error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
