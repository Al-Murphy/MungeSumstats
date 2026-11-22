"""Offline trusted terminal receiver. Artifact contents are always data only.

The workflow runs this file from its trusted workflow checkout. Only that
checkout supplies validators/contracts; downloaded Python/R is never executed.
A failed or incomplete producer returns an explicit unqualified result.
"""
import argparse
from datetime import datetime, timezone
from email.parser import Parser
import hashlib
import json
import math
from pathlib import Path, PurePosixPath, PureWindowsPath
import re
import sys
import zipfile

CI = Path(__file__).resolve().parent
sys.path.insert(0, str(CI))
import evidence as E
import service as S
from legacy import compare, operation, OPERATIONS, version_maps
from package import included_source
from verify import digest, write_json, tsv_pairs


def sha1_blob(raw):
    return hashlib.sha1(b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw).hexdigest()


def npath(value, lane):
    if lane == "windows":
        return E.canonical_work_path(value)
    E.require(isinstance(value, str) and value.startswith("/") and not any(
        part in {"", ".", ".."} or "\0" in part for part in value.split("/")[1:]), "Canonical Linux path required")
    return PurePosixPath(value)


def local_native(value, native, extracted):
    relative = value.relative_to(native).as_posix()
    E.relative_name(relative)
    result = extracted / relative
    E.plain_path(result, False)
    return result


def pin_native(pin, native, extracted, lane):
    file = local_native(npath(pin["path"], lane), native, extracted)
    actual = digest(file)
    E.require(all(pin[key] == actual[key] for key in ("bytes", "sha256", "md5") if key in pin),
              "Native retained file pin differs")
    return file


def same_pin(pin, path):
    actual = digest(path)
    E.require(all(pin[key] == actual[key] for key in ("bytes", "sha256", "md5") if key in pin), "Retained digest differs")


def read_zip_json(path, name, lane):
    with path.open("rb") as stream:
        E.directory_bounds(stream, path.stat().st_size)
        with zipfile.ZipFile(stream) as archive:
            files, _, _ = E.member_names(archive, lane=lane)
            E.require(name in files and files[name].file_size <= E.MIB, "Missing bounded outer receipt")
            return E.json_value(archive.read(files[name]))


def source_identity(service, run, start, checkout, root):
    event = start["event"]
    S.validate_event(run, event)
    E.require(checkout["schema"] == "mss-checkout-v1" and checkout["normal_current_checkout"] is True
              and checkout["commit"] == event["checkout_commit"], "Actual triggering checkout receipt differs")
    entry = service["sources"][checkout["commit"]]
    E.require(entry["commit_document"] == "git-commit-" + checkout["commit"] + ".json"
              and entry["tree_document"] == "git-tree-" + entry["tree"] + ".json", "Exact source service document leaves required")
    commit = E.document(root / entry["commit_document"])
    tree = E.document(root / entry["tree_document"])
    S.validate_commit(event, commit)
    E.require(commit["tree"]["sha"] == tree["sha"] == entry["tree"] == checkout["tree"]
              and tree["truncated"] is False and len(tree["tree"]) <= 50000, "Exact complete Git tree differs")
    expected = {}
    for row in tree["tree"]:
        E.relative_name(row["path"])
        if row["type"] == "tree":
            E.require(row["mode"] == "040000", "Invalid Git directory mode")
            continue
        E.require(row["type"] == "blob" and row["mode"] in {"100644", "100755"}
                  and re.fullmatch("[0-9a-f]{40}", row["sha"]), "Unsupported tracked source kind")
        E.require(row["path"] not in expected, "Duplicate Git tree path")
        expected[row["path"]] = row
    E.require(set(expected) == set(checkout["files"]) and 0 < len(expected) <= 50000, "Full Git source membership differs")
    return expected


def verify_checkout(root, service_root, service, run, start):
    receipt = E.document(root / "checkout/checkout.json")
    rows = source_identity(service, run, start, receipt, service_root)
    source = root / "checkout/source"
    actual = set()
    total = 0
    for path in source.rglob("*"):
        E.plain_path(path, path.is_dir())
        if path.is_dir():
            continue
        name = path.relative_to(source).as_posix()
        E.require(name in rows, "Untracked producer source file")
        E.require(path.stat().st_size <= 64 * E.MIB, "Oversized source blob")
        raw = path.read_bytes()
        row, pin = rows[name], receipt["files"][name]
        E.require(len(raw) <= 64 * E.MIB and sha1_blob(raw) == row["sha"] == pin["git_blob"]
                  and len(raw) == row["size"] == pin["bytes"] and row["mode"] == pin["git_mode"]
                  and hashlib.sha256(raw).hexdigest() == pin["sha256"], "Actual source bytes/Git blob differ")
        total += len(raw)
        E.require(total <= 512 * E.MIB, "Source total exceeds contract")
        actual.add(name)
    E.require(actual == set(rows) and receipt["bytes"] == total, "Incomplete source freeze")
    # These are independently reviewed prospective contracts, not an artifact's
    # choice of the test counts, closure or fixture assertions it wishes to pass.
    for name in ("suite-contract.json", "closure.json", "legacy-fixtures.json"):
        E.require((source / "tools/ci/installed" / name).read_bytes() == (CI / name).read_bytes(),
                  "Test/dependency contract differs from trusted receiver revision")
    return receipt, source


def command_records(evidence, lane, checkout):
    records = {}
    logs = evidence / "logs"
    E.plain_path(logs, True)
    for file in logs.glob("*.json"):
        row = E.document(file)
        if "argv" not in row:
            continue
        E.require(row["phase"] not in records and row["exit_code"] == 0
                  and isinstance(row["argv"], list) and row["argv"] and row["elapsed_seconds"] >= 0,
                  "Failed/duplicate/malformed command receipt")
        E.require(row["source_pins"]["checkout_commit"] == checkout["commit"]
                  and row["source_pins"]["checkout_tree"] == checkout["tree"], "Command used another source checkout")
        same_pin(row["source_pins"]["checkout_receipt"], evidence / "checkout/checkout.json")
        if lane == "windows":
            E.require(row["stopped"] is None, "Stopped Windows command")
        else:
            E.require(row["cleanup_signaled"] is False and row["error"] is None
                      and row["rss_limit_mode"] == "sampled_not_hard_Job"
                      and row["maximum_sampled_group_rss"] <= 7 * 1024**3,
                      "Linux command lacks native unsignaled closure")
        for suffix in (".stdout", ".stderr"):
            E.plain_path(file.with_suffix(suffix), False)
        records[row["phase"]] = (file, row)
    E.require(records, "Missing complete command receipts")
    return records


def command_spelling(value):
    return str(value).replace("\\", "/")


def check_argv(row, expected):
    E.match([command_spelling(x) for x in row["argv"]],
            [command_spelling(x) for x in expected], "actual command argv")


def command_recipes(checkout, installation, tools, native, lane):
    """Trusted single-checkout recipes, independent of artifact-selected argv."""
    source = native / "evidence/checkout/source"
    ci_source = source / "tools/ci/installed"
    benchmark = source / "benchmarks/2026-11-19/performance-upgrade"
    ops = benchmark / "operations"
    dependencies = native / "dependencies"
    out = native / "evidence"
    library = installation["library"]
    r, rscript = tools["R_launcher"]["path"], tools["Rscript"]["path"]
    worker = benchmark / "worker.R"
    observation = 'stopifnot(getRversion()=="4.5.3",.Machine$sizeof.pointer==8);writeLines(normalizePath(c(R.home(),R.home("bin")),winslash="/",mustWork=TRUE))'
    recipes = {
        "R-observation": [rscript, "--vanilla", "-e", observation],
        "compiler-version": [tools["compiler"]["path"], "--version"],
        "R-config": [r, "CMD", "config", "--all"],
        "current-R-CMD-build": [r, "CMD", "build", "--no-build-vignettes", "--no-manual", native / "package-build/source"],
        "current-R-CMD-INSTALL": [r, "CMD", "INSTALL", "--install-tests", "--preclean", "--library=" + library,
            out / "build" / ("MungeSumstats_" + tools["source_package_version"] + ".tar.gz")],
        "current-installed-identity": [rscript, "--vanilla", ci_source / "install-receipt.R", library,
            dependencies, out / "checkout/checkout.json", ci_source / "closure.json", lane, out / "build/installation.json"],
    }
    if lane == "windows":
        engine = tools["native_runtime"]["R_engine"]["path"]
        recipes["native-R-runtime"] = [engine, "--vanilla", "--slave", "--file=" + str(worker),
                                       "--args", "runtime", out / "runtime.tsv"]
    else:
        recipes["native-R-runtime"] = [rscript, "--vanilla", worker, "runtime", out / "runtime.tsv"]
        recipes["native-Linux-ownership-controls"] = [tools["Python"]["command"], "-I", "-S", "-B", ci_source / "linux-control.py"]
    for name, version in installation["dependencies"].items():
        archive = name + "_" + version + ".tar.gz"
        recipes["dependency-" + name] = [r, "CMD", "INSTALL", "--preclean", "--library=" + str(dependencies),
                                       native / "dependency-archives" / archive]
    for name in OPERATIONS:
        recipes["legacy-" + name] = [rscript, "--vanilla", ci_source / "legacy.R", name, lane, library,
            dependencies, out / "legacy-inputs", out / "legacy" / name, out / "build/installation-with-payloads.json"]
    if lane == "windows":
        small, probes, scratch = out / "scientific-inputs", out / "native", native / "suite-scratch/native"
        recipes["complete-current-installed-suite"] = [rscript, "--vanilla", ci_source / "installed-suite.R",
            library, dependencies, npath(installation["package_path"], lane) / "tests/testthat",
            out / "installed-suite", native / "home"]
        recipes["scientific-current-references"] = [rscript, "--vanilla", ops / "references.R",
            library, dependencies, small, small / "unused-matched", "scientific"]
        recipes["current-VCF-write-index-query"] = [rscript, "--vanilla", ci_source / "vcf.R",
            library, dependencies, checkout["commit"], out / "vcf"]
        for name in ("hard-cases", "declared"):
            for number in (2, 3):
                run = out / "science" / (name + "-f" + str(number))
                recipes["science-" + name + "-f" + str(number)] = [engine, "--vanilla", "--slave",
                    "--file=" + str(worker), "--args", "run", run / "settings.tsv"]
        for name in ("preservation", "registered"):
            recipes["checkout-" + name] = [rscript, "--vanilla", source / "tests/standalone" / (name + ".R"), source]
        for filename, arg in (("test_core.cpp", "fixtures"), ("test_blocks.cpp", "blocks"),
                              ("projection.cpp", "projection.chain"), ("test_format3.cpp", "format3"), ("test_fai.cpp", "fai")):
            stem = filename.removesuffix(".cpp")
            executable = probes / (stem + ".exe")
            recipes["checkout-compile-" + stem] = [tools["compiler"]["path"], "-std=c++17", "-Wall", "-Wextra",
                "-pedantic", "-O2", "-UNDEBUG", source / "tests/native" / filename, "-o", executable]
            recipes["checkout-native-" + stem] = [executable, scratch / arg]
        options = [tools["compiler"]["path"], "-std=c++17", "-O1", "-UNDEBUG", "-I", source / "src"]
        recipes["checkout-compile-lookup"] = [*options, "-DMSS_LOOKUP_COUNTERS",
            source / "tests/standalone/lookup_parity.cpp", "-o", probes / "lookup.exe"]
        for name in ("parity", "pressure"):
            recipes["checkout-lookup-" + name] = [probes / "lookup.exe", name, scratch / ("lookup-" + name)]
        recipes["checkout-compile-normalize"] = [*options, source / "tests/native/test_normalize_reuse.cpp",
                                               "-o", probes / "normalize-reuse.exe"]
        recipes["checkout-normalize-reuse"] = [probes / "normalize-reuse.exe", scratch / "normalize-reuse"]
        for mode in ("failure", "width", "precedence"):
            driver = "failure_parity.R" if mode == "failure" else "parser_error_order.R"
            recipes["checkout-parser-" + mode] = [rscript, "--vanilla", ops / driver, library, dependencies,
                small, probes / mode, *([] if mode == "failure" else [mode])]
        for mode in ("startup", "bp", "parallel-vcf"):
            recipes["checkout-" + mode] = [rscript, "--vanilla", source / "tests/standalone/test_variantannotation_deferral.R",
                mode, library, dependencies, probes / (mode + ".rds")]
    return recipes


def verify_recipe_argv(logs, recipes):
    E.require(set(logs) == set(recipes), "Complete exact current-checkout command recipe roster differs")
    for label, recipe in recipes.items():
        check_argv(logs[label][1], recipe)


def verify_call_settings(actual, expected, path_fields, lane):
    E.require(set(actual) == set(expected), "Scientific worker settings membership differs")
    for key,value in expected.items():
        if key in path_fields:
            E.require(npath(actual[key],lane)==npath(value,lane), "Scientific call path differs: "+key)
        else:
            E.match(actual[key],value,"scientific call setting "+key)


def verify_command_recipes(logs, checkout, installation, tools, native, lane, evidence, extracted):
    E.require(npath(checkout["source"], lane) == native / "evidence/checkout/source",
              "Frozen source native location differs from retained checkout")
    python = tools["Python"]
    E.require(python["version"] == [3, 12, 10] and python["implementation"] == "CPython"
              and npath(python["command"], lane) == npath(python["executable"]["path"], lane),
              "Actual Python identity incomplete")
    identity = python["executable"]
    E.require(type(identity["bytes"]) is int and identity["bytes"] > 0
              and re.fullmatch("[0-9a-f]{64}", identity["sha256"])
              and re.fullmatch("[0-9a-f]{32}", identity["md5"]), "Actual Python executable pin missing")
    recipes = command_recipes(checkout, installation, tools, native, lane)
    verify_recipe_argv(logs, recipes)
    if lane == "windows":
        phases = E.document(evidence / "native/probes.json")
        for name, phase in phases.items():
            if name.startswith("compile-"):
                expected = recipes["checkout-" + name][-1]
                E.require(npath(phase["compiled_executable"]["path"], lane) == npath(str(expected), lane),
                          "Compiled PE is not the binary named by its exact recipe")
                E.check_probe_executable(pin_native(phase["compiled_executable"], native, extracted, lane))
    return recipes

def dependencies(evidence, lane, native, tools, logs):
    import windows_job as w
    from dependencies import plan, url, SIX_REFERENCE_DISK_FLOOR
    rows = plan(E.document(CI / "closure.json"), lane)
    area = evidence / "dependency-plan"
    E.require({p.name for p in area.iterdir()} == {"verified-get.jsonl", "source-attempts.jsonl", "acquisition-capacity.json"},
              "Exact successful acquisition metadata roster required")
    def lines(name):
        path = area / name
        E.plain_path(path, False)
        E.require(0 < path.stat().st_size <= E.MAX_DIRECTORY, "Bounded acquisition ledger required")
        raw = path.read_bytes()
        E.require(raw.endswith(b"\n"), "Incomplete acquisition ledger")
        return [E.json_value(line) for line in raw.splitlines()]
    gets, attempts = lines("verified-get.jsonl"), lines("source-attempts.jsonl")
    E.require(len(gets) == len(rows) and len(rows) <= len(attempts) <= 2*len(rows), "Complete exact download ledger required")
    position = 0
    for row, got in zip(rows,gets):
        E.match({key:got[key] for key in row},row,"frozen dependency source/version/bytes")
        E.require(got["verified"] is True and got["observed"]["bytes"] == int(row["bytes"])
                  and got["observed"]["sha256"] == row["sha256"], "Pinned source acquisition failed")
        first = attempts[position]; position += 1
        E.require(first["package"] == row["package"] and first["source"] == row["source"]
                  and first["attempted_url"] == row["request_url"], "Original exact source request differs")
        accepted = first
        if first["outcome"] == "HTTP_error":
            source = S.urllib.parse.urlsplit(row["source"])
            E.require(first["http_status"] in {404,410} and row["kind"] == "code"
                      and source.hostname in {"cran.r-project.org","cloud.r-project.org"}
                      and source.path == "/src/contrib/" + row["archive"], "Only unavailable exact CRAN source may fall back")
            accepted = attempts[position]; position += 1
            expected = "https://" + source.hostname + "/src/contrib/Archive/" + row["package"] + "/" + row["archive"]
            E.require(accepted["package"] == row["package"] and accepted["source"] == row["source"]
                      and accepted["attempted_url"] == expected, "Fallback changed source/version or repository")
        E.require(accepted["outcome"] == "opened" and accepted["http_status"] == 200
                  and got["actual_request_url"] == accepted["attempted_url"]
                  and got["final_url"] == accepted["final_url"], "Actual accepted pinned source differs")
        url(got["final_url"])
        check_argv(logs["dependency-" + row["package"]][1], [tools["R_launcher"]["path"], "CMD", "INSTALL", "--preclean",
                    "--library=" + str(native / "dependencies"), native / "dependency-archives" / row["archive"]])
        E.match(logs["dependency-" + row["package"]][1]["provenance"],
                {"package":row["package"], "version":row["version"], "archive":got["observed"]}, "ordinary dependency install source")
    E.require(position == len(attempts), "Unconsumed acquisition source attempts")
    capacity = E.document(area / "acquisition-capacity.json")
    archive_bytes = sum(int(row["bytes"]) for row in rows)
    expansion = sum(row["unpacked_bytes"] for row in rows)
    growth = archive_bytes + 2*expansion + max(row["unpacked_bytes"] for row in rows) + 2_000_000_000
    if lane == "windows": growth = max(growth,SIX_REFERENCE_DISK_FLOOR-w.RESERVE)
    E.require(capacity["lane"]==lane and capacity["archive_bytes"]==archive_bytes
              and capacity["total_expansion_bytes"]==expansion and capacity["additional_growth_bytes"]==growth,
              "Acquisition capacity omitted a pinned archive/expansion")
    admission = capacity["admission"]
    replay = w.capacity_record(admission["total_ram_bytes"], admission["available_ram_bytes"], admission["free_bytes"], growth,
                               "six-reference-minimum" if lane=="windows" else "two-reference-minimum")
    E.match({key:admission[key]for key in replay},replay,"actual capacity arithmetic")
    if lane == "windows": E.verify_reference_admission(evidence/"logs/six-reference-admission.jsonl",admission)
    return {"packages":len(rows),"same_version_sources":True,"exact_cran_fallback_only":True}


def ordinary_install(evidence, source, checkout, lane, native, extracted, logs):
    build = E.document(evidence / "build/build.json")
    install = E.document(evidence / "build/installation.json")
    current = E.document(evidence / "build/installation-with-payloads.json")
    E.require(build["schema"] == "mss-normal-build-v1" and build["historical_archive_used"] is False
              and build["tests_installed_by_R_CMD_INSTALL"] is True, "Normal current R build/install required")
    for value in (build, install, current):
        E.require(value["checkout_commit"] == checkout["commit"] and value["checkout_tree"] == checkout["tree"],
                  "Built/installed source differs")
    E.require(current["lane"] == lane and current["source_files"] == checkout["files"]
              and current["normal_current_checkout"] is True and current["development"] is False
              and current["reference_payload_inventory_still_required"] is False
              and current["R_version"] == "4.5.3", "Current installation closure differs")
    E.match({key:value for key,value in current.items() if key not in {"reference_payloads", "reference_payload_inventory_still_required"}},
            {key:value for key,value in install.items() if key != "reference_payload_inventory_still_required"}, "installation extension")
    E.require(all(npath(current[key],lane) == native/"home" for key in ("startup_home","startup_R_USER","startup_HOME")),
              "R startup HOME/R_USER must bind private current process home")
    E.require(npath(current["library"], lane).parent == native and
              npath(current["package_path"], lane) == npath(current["library"], lane) / "MungeSumstats"
              and npath(current["package_path"], lane) in npath(current["dll"], lane).parents,
              "Installed package/DLL outside selected private library")
    version = Parser().parsestr((source / "DESCRIPTION").read_text())["Version"]
    E.require(current["package_version"] == version.replace("-", "."), "Installed current package version differs")
    archive = pin_native(build["generated_current_archive"], native, extracted, lane)
    binding = included_source(archive, checkout, "MungeSumstats")
    E.match(binding, build["source_binding"], "normal build included-source binding")
    for field in ("build_stdout", "build_stderr", "build_execution", "install_stdout", "install_stderr", "install_execution", "installed_identity"):
        pin_native(build[field], native, extracted, lane)
    tools = E.document(evidence / "toolchain.json")
    E.require(tools["schema"] == "mss-current-toolchain-v1" and re.match(r"^R version 4\.5\.3(?: |$)", tools["runtime"]["R_version"]) is not None
              and set(tools["native_runtime"]) == {"R_engine", "R_shared", "BLAS_library", "LAPACK_library"}, "Native toolchain identity incomplete")
    E.require(npath(tools["runtime"]["R_home"],lane) == npath(current["R_home"],lane), "R native process home differs")
    for name,pin in {**tools["native_runtime"], "compiler":tools["compiler"], "Rscript":tools["Rscript"], "R_launcher":tools["R_launcher"]}.items():
        npath(pin["path"],lane)
        E.require(type(pin["bytes"]) is int and pin["bytes"] > 0 and re.fullmatch("[0-9a-f]{64}",pin["sha256"])
                  and re.fullmatch("[0-9a-f]{32}",pin["md5"]), "Incomplete native R/compiler binary identity")
        if name in tools["native_runtime"]: E.require(npath(pin["path"],lane)==npath(tools["runtime"][name],lane), "Native runtime selected another binary")
    E.require((lane == "windows" and tools["Rtools"] == "45") or (lane == "linux-gold" and tools["Rtools"] is None), "Toolchain lane differs")
    build_cmd = logs["current-R-CMD-build"][1]
    check_argv(build_cmd, [tools["R_launcher"]["path"], "CMD", "build", "--no-build-vignettes", "--no-manual", native / "package-build/source"])
    check_argv(logs["current-R-CMD-INSTALL"][1], [tools["R_launcher"]["path"], "CMD", "INSTALL", "--install-tests", "--preclean",
               "--library=" + current["library"], build["generated_current_archive"]["path"]])
    check_argv(logs["current-R-CMD-INSTALL"][1], build["install_argv"])
    before = E.document(evidence / "installed-before.json")
    E.match(before, E.document(evidence / "installed-after.json"), "installed package before/after")
    h = hashlib.sha256()
    for name, pin in sorted(before["files"].items()):
        E.relative_name(name)
        h.update(name.encode() + b"\0" + pin["sha256"].encode() + b"\n")
    E.require(h.hexdigest() == before["tree_sha256"], "Installed package inventory hash differs")
    dll_relative = npath(current["dll"], lane).relative_to(npath(current["package_path"], lane)).as_posix()
    E.require(before["files"][dll_relative]["md5"] == current["dll_md5"], "Installed DLL identity differs")
    for prefix in ("dependency-library", "base-library"):
        E.match(E.document(evidence / (prefix + "-before.json")), E.document(evidence / (prefix + "-after.json")), prefix + " immutable inventory")
    closure = E.document(CI / "closure.json")
    version_maps(current, closure)
    E.require(set(E.document(evidence / "dependency-library-before.json")) == set(current["dependencies"]), "Dependency library membership differs")
    dependencies(evidence,lane,native,tools,logs)
    return current, dict(tools, source_package_version=version)


def legacy_contract(receipt, name):
    fixed = {
        "allele-flip": {"historical_full_output_line_equality", "all_decoded_columns_and_values_equal", "nonempty_unique_variant_membership"},
        "indels": {"known_indel_retained", "two_expected_rows", "indel_alleles_and_anchor_retained",
                   "indel_effect_and_frequency_preserved", "one_unidentified_indel_removed_and_logged",
                   "logged_indel_effect_and_frequency_preserved"},
        "multiple-rsids": {"first_identifier_and_conversion_flag", "historical_full_decoded_equality", "five_historical_log_categories", "historical_log_row_counts"},
        "tabix": {"bgzf_and_index_exist", "historical_nine_column_header", "historical_full_read_equality", "sort-GenomicRanges", "sort-data.table"},
        "headers": {"frozen_vcf_528_lines"},
    }[name]
    wanted = {"frozen_education_fixture_93_records"} | fixed
    if name == "tabix":
        wanted |= {label + "-range-count" for label in ("ordered", "unordered")}
        wanted |= {label + "-range-" + str(i) + "-decoded-equality" for label in ("ordered", "unordered") for i in (1,2)}
    if name == "headers":
        wanted |= {prefix + compression + "-" + skip for prefix in ("tsv-first-", "tsv-all-") for compression in ("plain", "gzip", "bgzf") for skip in ("FALSE", "TRUE")}
        wanted |= {prefix + compression for prefix in ("vcf-metadata-", "vcf-rows-") for compression in ("plain", "gzip", "bgzf")}
    initial = len(wanted)
    names = ("ALSvcf.vcf", "ALSvcf.vcf.gz", "eduAttainOkbay.txt", "eduAttainOkbay.txt", "ALSvcf.vcf", "ALSvcf.vcf.gz")
    wanted |= {"unchanged-" + value + "-" + str(initial+i) for i,value in enumerate(names)}
    wanted.add("loaded_dll_unchanged")
    E.require(set(receipt["checks"]) == wanted and all(value is True for value in receipt["checks"].values()),
              "Complete prospective legacy scientific assertion roster differs")


def verify_qpc_context(start, job, control, qualification, sealer, admissions, finished_utc):
    """Replay the accepted original shared-QPC chronology using generic steps.

    Only named terminal Jobs service timestamps and the independently bound
    first observation select the interval. No producer-supplied phase limits.
    """
    import windows_job as w
    begin = w.receipt_epoch(job["started_at"])
    observed = w.receipt_epoch(start["observed_utc"])
    first_tick = start["observed_monotonic_seconds"]
    E.require(type(first_tick) in (int,float) and math.isfinite(first_tick) and first_tick >= 0,
              "Invalid first shared-QPC observation")
    steps = {step["name"]:step for step in job["steps"]}
    E.require(len(steps) == len(job["steps"]), "Duplicate terminal service step name")
    def bounds(name):
        step = steps[name]
        E.require(step["status"] == "completed" and step["conclusion"] == "success", "Incomplete service qualification step")
        return w.receipt_epoch(step["started_at"]), w.receipt_epoch(step["completed_at"]) + 1
    pre_begin, pre_end = bounds("Observe actual job start")
    phase_begin, phase_end = bounds("Run complete current installed package")
    E.require(pre_begin <= observed < pre_end and begin <= observed <= begin + 180*60,
              "First observation outside original service step/job clock")
    deadlines = {name:begin+minutes*60 for name,minutes in (("compute",300),("drain",305),("seal",316),("final",317))}
    def check_anchor(point, phase):
        E.require(set(point) == {"local_monotonic","utc_epoch","qpc","qpc_end"}
                  and all(type(v) in (int,float) and math.isfinite(v) for v in point.values()), "Invalid native clock anchor")
        wall,lo,hi = point["utc_epoch"],point["qpc"],point["qpc_end"]
        limit = min(phase_end,deadlines[phase])
        E.require(phase_begin <= wall < limit and first_tick+phase_begin-observed <= lo <= hi < first_tick+limit-observed,
                  "Native anchor outside selected service step or original deadline")
    def precedes(a,b):
        E.require(a["utc_epoch"] <= b["utc_epoch"] and a["qpc_end"] <= b["qpc"], "Native shared-QPC phases out of order")
    anchors = control["anchors"]
    E.require(set(anchors) == {"launch","cutoff","drain","complete"}, "Incomplete Job-control anchors")
    origin = anchors["launch"]
    E.require(origin["local_monotonic"] == control["controller_started_monotonic"], "Job-control start clock differs")
    prior = None
    for label,elapsed in (("launch",0),("cutoff",control["deadline_trigger"]["elapsed_seconds"]),
                          ("drain",control["drain_elapsed_seconds"]),("complete",control["elapsed_seconds"])):
        point = anchors[label]
        check_anchor(point,"compute" if label in {"launch","cutoff"} else "drain")
        E.require(point["local_monotonic"]-origin["local_monotonic"] == elapsed, "Control anchor local elapsed differs")
        if prior is not None: precedes(prior,point)
        # Unchanged32ms allowance is ONLY Windows monotonic tick quantization;
        # it is not scientific/numeric tolerance or an extension of deadlines.
        E.require(point["qpc"]-origin["qpc_end"]-.032 <= elapsed <= point["qpc_end"]-origin["qpc"]+.032,
                  "Job-control local elapsed disagrees with shared QPC")
        prior = point
    completion = {}
    launch = qualification["launch_anchor"]
    for mode,stop in (("qualify",qualification),("seal",sealer)):
        wall,tick = stop["drained_utc_epoch"],stop["drained_qpc"]
        E.require(all(type(v)in(int,float)and math.isfinite(v)for v in(wall,tick)), "Invalid drain clock")
        limit = deadlines["drain"if mode=="qualify"else"final"]
        E.require(phase_begin <= wall < min(limit,phase_end) and first_tick <= tick < first_tick+limit-observed,
                  "Late or inconsistent owned worker drain")
        completion[mode] = wall,tick
    check_anchor(launch,"compute")
    precedes(anchors["complete"],launch)
    E.require(launch["utc_epoch"] <= completion["qualify"][0] and launch["qpc_end"] <= completion["qualify"][1],
              "Qualification launched after drain")
    finish = w.receipt_epoch(finished_utc)
    E.require(completion["qualify"][0] <= completion["seal"][0] and completion["qualify"][1] <= completion["seal"][1]
              and launch["utc_epoch"] <= finish <= completion["qualify"][0] and observed <= finish < deadlines["compute"]
              and finish < phase_end and finish <= w.receipt_epoch(job["completed_at"]), "Qualification finish/drain chronology differs")
    E.require(len(admissions) == 2, "Both original commitment gates required")
    for row in admissions: check_anchor(row["anchor"],"compute")
    precedes(launch,admissions[0]["anchor"])
    precedes(admissions[0]["anchor"],admissions[1]["anchor"])
    E.require(admissions[1]["anchor"]["utc_epoch"] <= finish, "Commitment release follows qualification finish")
    return {"shared_QPC_anchors_verified":True,"original_wall_and_monotonic_deadlines_verified":True}


def verify_controls(extracted, evidence, start, job, lane, logs):
    import windows_job as w
    beginning = w.receipt_epoch(start["started_utc"])
    E.require(job["labels"] == ["windows-2022" if lane == "windows" else "ubuntu-22.04"]
              and start["host"]["os_arch"] == start["host"]["process_arch"] == "X64", "Standard native hosted lane differs")
    if lane == "windows": E.require(start["image"]["os"] == "win22", "Windows image differs")
    E.require(start["schema"] == "mss-current-job-start-v1" and start["lane"] == lane
              and start["started_utc"] == job["started_at"] and job["runner_name"] == start["jobs"]["jobs"][[str(r["id"]) for r in start["jobs"]["jobs"]].index(start["job_id"])]["runner_name"], "First actual Job observation differs")
    deadlines = E.document(extracted / "q/deadlines.json")
    E.match(deadlines["absolute_utc_epoch"], {name:beginning + minutes*60 for name,minutes in (("compute",300),("drain",305),("seal",316),("final",317))}, "original phase clocks")
    E.require(beginning <= w.receipt_epoch(start["observed_utc"]) <= w.receipt_epoch(job["completed_at"])
              <= beginning + 360*60, "Actual hosted Job exceeded its finite upload-inclusive policy")
    if lane == "windows":
        for mode in ("job-control", "qualify", "seal"):
            row = E.document(extracted / ("q/" + mode + "-stop.json"))
            E.require(row["mode"] == mode and row["assigned"] is True and row["drained"] is True
                      and row["active_after"] == 0, "Native Job ownership/drain incomplete")
            if mode in {"job-control", "qualify"}:
                E.match(row["attempt"], {key:start[key]for key in ("run_id","run_attempt","job_id")}, "owned Job original attempt")
            if mode == "job-control":
                E.check_job_control(row, evidence / "job-control")
            else:
                E.require(row["exit_code"] == 0 and row["reason"] is None and row["active_before"] == 0,
                          "Worker or sealer failed/required cleanup")
                E.require(row["drained_utc_epoch"] <= beginning + (305 if mode == "qualify" else 317)*60,
                          "Worker/sealer closing deadline exceeded")
        admissions = E.verify_job_commit_admission(evidence)
        control = E.document(extracted / "q/job-control-stop.json")
        qualification = E.document(extracted / "q/qualify-stop.json")
        sealer = E.document(extracted / "q/seal-stop.json")
        E.match(control,E.document(evidence / "job-control-stop.json"),"Job-control duplicate receipt")
        E.match(qualification,E.document(evidence / "qualify-stop.json"),"qualifier duplicate receipt")
        E.match(sealer,E.document(extracted / "q/artifact-status.json")["sealer"],"sealer duplicate receipt")
        verify_qpc_context(start,job,control,qualification,sealer,admissions,
                           E.document(evidence / "qualification-candidate.json")["finished_utc"])
    else:
        stop = E.document(extracted / "q/linux-stop.json")
        E.require(stop["schema"] == "mss-linux-gold-stop-v1" and stop["owned_tree_quiescent"] is True
                  and stop["cleanup_signaled"] is False and stop["error"] is None
                  and w.receipt_epoch(stop["finished_utc"]) <= beginning + 300*60, "Native Linux final drain/compute clock failed")
        file, row = logs["native-Linux-ownership-controls"]
        text = file.with_suffix(".stderr").read_text()
        E.require("test_transitive_workers_and_adoptees_only" in text and "test_new_session_adoptee_gets_killed_and_reaped" in text
                  and re.search(r"Ran 2 tests\s+\n\s*OK\s*$", text) and "skipped" not in text, "Both actual Linux process controls must pass")


def verify_native(evidence, source, checkout, install, tools, logs, native, extracted):
    phases = E.document(evidence / "native/probes.json")
    labels = {"preservation", "registered", "compile-lookup", "lookup-parity", "lookup-pressure", "compile-normalize", "normalize-reuse",
              "parser-failure", "parser-width", "parser-precedence", "startup", "bp", "parallel-vcf"}
    labels |= {prefix + name for name in ("test_core", "test_blocks", "projection", "test_format3", "test_fai") for prefix in ("compile-", "native-")}
    E.require(set(phases) == labels, "Complete single-checkout native/probe roster differs")
    for label, phase in phases.items():
        driver = phase["driver"]
        relative = npath(driver["path"], "windows").relative_to(npath(checkout["source"], "windows")).as_posix()
        E.require(phase["role"] == phase["driver_role"] == "checkout" and phase["commit"] == phase["driver_commit"] == checkout["commit"]
                  and checkout["files"][relative]["sha256"] == driver["sha256"], "Probe driver/checkout identity differs")
        tag = "checkout-" + label
        log, command = logs[tag]
        E.match(command["provenance"], {key:phase[key] for key in ("role", "commit", "driver_role", "driver_commit", "driver")}, "probe command provenance")
        for field, suffix in (("execution", ".json"), ("stdout", ".stdout"), ("stderr", ".stderr")):
            file = pin_native(phase[field], native, extracted, "windows")
            E.require(file == log.with_suffix(suffix), "Probe phase log binding differs")
        if label.startswith("compile-"):
            executable = pin_native(phase["compiled_executable"], native, extracted, "windows")
            E.check_probe_executable(executable)
            E.require(command["argv"][0] == tools["compiler"]["path"] and command["argv"][-2:] == ["-o", phase["compiled_executable"]["path"]]
                      and driver["path"] in command["argv"] and "-UNDEBUG" in command["argv"], "Native compiler recipe differs")
        elif label in {"lookup-parity", "lookup-pressure", "normalize-reuse"}:
            path = evidence / "native" / (label + ".stdout")
            E.match(digest(path), digest(log.with_suffix(".stdout")), "Copied native transcript")
            if label == "lookup-parity": E.match(phase["typed_boundary_sha256"], E.verify_typed_boundaries(path), "typed boundary gold")
            elif label == "lookup-pressure": E.match(phase["formats"], E.verify_lookup_pressure(path), "pressure arithmetic")
            else: E.normalization_phases(path)
        elif label.startswith("parser-"):
            mode = label.removeprefix("parser-")
            E.match(phase["gold"], E.verify_parser_gold(evidence / "native" / mode / "errors.tsv", mode), "parser full-field gold")
        elif label in {"startup", "bp", "parallel-vcf"}:
            pin_native(phase["receipt"], native, extracted, "windows")
            E.deferral_result(evidence / "native" / (label + ".rds"), label)
            E.require(log.with_suffix(".stdout").read_text().splitlines()[-1] == label + " passed", "Incomplete R deferral result")
    return sorted(labels)


def verify_science(evidence, checkout, install, native, extracted, tools):
    from run import SETTINGS
    # Read literal column maps from the trusted receiver checkout, never from
    # downloaded scientific.py or an artifact's chosen options.
    import importlib.util
    path = CI.parents[2] / "benchmarks/2026-11-19/performance-upgrade/operations/scientific.py"
    spec = importlib.util.spec_from_file_location("mss_trusted_literal_receiver", path)
    literal = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(literal)
    fixture = evidence / "scientific-inputs"
    provenance = E.document(fixture / "fixture-provenance.json")
    E.require(provenance["schema"] == "mss-current-scientific-v1" and provenance["normal_current_install"] is True
              and provenance["checkout_commit"] == checkout["commit"], "Current scientific fixture provenance differs")
    native_fixture = native / "evidence/scientific-inputs"
    assets = {str(native_fixture / path.relative_to(fixture).as_posix()): digest(path) for path in fixture.rglob("*") if path.is_file()}
    def local(value): return local_native(npath(value, "windows"), native, extracted)
    row = E.document(evidence / "science/scientific-editions.json")
    E.require(row["single_current_install"] is True and row["source_commit"] == checkout["commit"]
              and row["complete_independent_oracle_fields"] == 19 and row["immutable_fixtures"] is True
              and row["performance_comparison_claimed"] is False and len(row["editions"]) == 4, "Four full scientific editions required")
    expected = {(name, number) for name in ("hard-cases", "declared") for number in (2,3)}
    E.require({(r["case"], r["format"]) for r in row["editions"]} == expected, "Scientific edition membership differs")
    base_namespaces = None
    for item in row["editions"]:
        name, number = item["case"], item["format"]
        count, gold = {"hard-cases":(32,"ca4f4d2f59d70f6052a5eb149e4e72011fefc175bef074897f8f15578824a565"),
                       "declared":(9,"36248164266feb550b628be89212864e80b3d9c0c3b07c71cf61aedbfcf3319e")}[name]
        oracle = fixture / (name + ".gold.tsv")
        E.require(digest(oracle)["sha256"] == item["oracle_sha256"] == gold and item["rows"] == count, "Independent19-field oracle differs")
        ref = provenance["references"][str(number)]
        case = {"rows":count,"input":{"path":str(native_fixture / (name + (".tsv.gz" if name=="hard-cases" else ".tsv")))},
                "oracle":{"path":str(native_fixture / (name+".gold.tsv"))},"settings":dict(SETTINGS, frequency_type="EAF", verify_acquisition_receipt=True,
                    columns=literal.DECLARED_MAPPING if name=="declared" else literal.MAPPING)}
        edition = evidence / "science" / (name + "-f" + str(number)) / "edition"
        run = native / "evidence/science" / (name + "-f" + str(number))
        opt = case["settings"]
        settings = dict(mode="diagnostic", profile="none", allocation_threshold="10000",
            input=case["input"]["path"], rows=str(count), output=str(run), dll=install["dll"],
            dll_md5=install["dll_md5"], reference=str(native_fixture/("scientific_"+str(number))),
            core_version="mss-preservation-"+str(number), fasta_pin_path=str(native_fixture/"scientific.fa"),
            **{"library.0":install["library"], "library.1":str(native/"dependencies")})
        settings.update({key:str(value)for key,value in opt.items()
                         if value is not None and key not in {"columns","verify_acquisition_receipt"}})
        settings.update({"column."+key:value for key,value in opt["columns"].items()})
        input_identity = digest(local(case["input"]["path"]))
        settings.update(acquisition_md5=input_identity["md5"], acquisition_bytes=str(input_identity["bytes"]))
        actual_settings = tsv_pairs(edition.parent/"settings.tsv",("key","value"))
        path_settings = {"input","output","dll","reference","fasta_pin_path","library.0","library.1"}
        verify_call_settings(actual_settings,settings,path_settings,"windows")
        files = {p.name:digest(p) for p in edition.iterdir() if p.is_file()}
        E.match(files, item["edition_files"], "scientific edition complete inventory")
        E.match(E.edition_decision(edition,case,ref,files,assets,local), item["scientific"], "complete full-field scientific decision")
        E.match(E.stable_edition_manifest(edition/"manifest.rds",ref,files,assets), item["stable_manifest"], "native RDS edition manifest")
        metrics = tsv_pairs(edition.parent / "metrics.tsv", ("key","value"))
        E.require(all(metrics.get(key)==value for key,value in tools["runtime"].items())
                  and metrics["dll_md5"] == install["dll_md5"] and metrics["threads"] == "1"
                  and metrics["source_rows"] == metrics["normalized_rows"] == str(count), "Scientific worker selected another package/rows/thread")
        namespaces = E.read_table(edition.parent / "namespaces.tsv", ["package", "version", "path"])
        seen = {}
        base = E.document(evidence / "base-library-before.json")
        for ns in namespaces:
            name, version, physical = ns["package"], ns["version"], npath(ns["path"], "windows")
            E.require(name not in seen, "Duplicate scientific namespace")
            if name == "MungeSumstats":
                wanted_version, wanted_path = install["package_version"], npath(install["package_path"], "windows")
            elif name in install["dependencies"]:
                wanted_version = install["dependencies"][name]
                wanted_path = native / "dependencies" / name
            else:
                E.require(name in base and physical == npath(install["R_home"], "windows") / "library" / name,
                          "Unselected scientific namespace")
                wanted_version, wanted_path = version, physical
            E.require(version.replace("-", ".") == wanted_version.replace("-", ".") and physical == wanted_path,
                      "Scientific namespace version/physical library differs")
            seen[name] = (version, str(physical))
        E.require("MungeSumstats" in seen and "base" in seen, "Missing current/base namespace")
        if base_namespaces is None: base_namespaces = {}
        for key,value in seen.items():
            if key in base_namespaces: E.match(value, base_namespaces[key], "scientific namespace identity")
            base_namespaces[key] = value
        for key in ("command_stdout", "command_stderr", "command_execution"): pin_native(item[key],native,extracted,"windows")
    for name in ("hard-cases", "declared"):
        paths=[evidence/"science"/(name+"-f"+str(n))/"edition/normalized.tsv"for n in(2,3)]
        E.require(paths[0].read_bytes()==paths[1].read_bytes(), "Full cross-format scientific rows differ")
    return {"editions":4,"oracle_fields":19,"cross_format":True,"performance_claim":False}


def lane(root, service, run, jobs, lane, output):
    item = service["artifacts"][lane]
    zip_path = root / (lane + ".zip")
    E.require(item["original_zip"]["artifact_id"] == str(item["service"]["id"]), "Exact ZIP artifact id differs")
    status = read_zip_json(zip_path,"q/artifact-status.json",lane)
    inventory = status["inventory"]["path"]
    native = npath(inventory,lane).parent
    start = read_zip_json(zip_path,"job-start.json",lane)
    matches = [r for r in jobs["jobs"] if str(r["id"]) == start["job_id"]]
    E.require(len(matches)==1, "Original first observation Job not in terminal attempt")
    job = matches[0]
    E.require(start["repository"]==service["repository"] and start["run_id"]==service["run_id"]
              and start["run_attempt"]==service["run_attempt"] and job["name"]==S.JOB_NAMES[lane]
              and job["status"]=="completed" and job["conclusion"]=="success", "Lane terminal service identity/conclusion differs")
    extracted = output / lane
    E.unpack_sealed(zip_path, extracted, item["original_zip"]["bytes"], item["original_zip"]["sha256"], str(native), lane=lane)
    evidence = extracted / "q/evidence"
    checkout, source = verify_checkout(evidence, root, service, run, start)
    runner = E.document(evidence/"runner-provenance.json")
    E.require(runner["checkout_commit"]==checkout["commit"] and runner["job_id"]==start["job_id"]
              and runner["host"]==start["host"] and runner["event"]==start["event"], "Current runner provenance differs")
    logs = command_records(evidence,lane,checkout)
    verify_controls(extracted,evidence,start,job,lane,logs)
    installation, tools = ordinary_install(evidence,source,checkout,lane,native,extracted/"q",logs)
    verify_command_recipes(logs,checkout,installation,tools,native,lane,evidence,extracted/"q")
    candidate = E.document(evidence/"qualification-candidate.json")
    E.require(candidate["schema"]=="mss-current-candidate-v1" and candidate["checks_passed"] is True
              and candidate["checkout_commit"]==checkout["commit"] and candidate["checkout_tree"]==checkout["tree"]
              and candidate["lane"]==lane and candidate["complete_installed_package_terminal_acceptance"] is False
              and candidate["receiver_verified"] is False and candidate["five_legacy_direct_operations"] is True
              and candidate["ordinary_installed_suite"] is (lane == "windows")
              and candidate["typed_native_VCF_and_science"] is (lane == "windows")
              and candidate["same_checkout_Linux_parity"] is (lane == "windows"), "Producer candidate cannot supply terminal acceptance")
    for name in OPERATIONS:
        receipt,_,_ = operation(evidence/"legacy"/name,name,lane,installation)
        legacy_contract(receipt,name)
    checks = {"normal_current_checkout":True,"full_Git_blob_membership":len(checkout["files"]),"native_terminal_closure":True,
              "genuine_five_legacy_operations":True}
    if lane == "windows":
        contract=E.document(CI/"suite-contract.json")
        checks["suite"]=E.verify_installed_suite(evidence/"installed-suite",contract)
        fixture = E.document(CI/"legacy-fixtures.json")["fixtures"]["fixtures/ALSvcf.vcf"]
        checks["retained_ALS_triplets"] = E.verify_als_triplets(evidence/"installed-suite/als-fixtures",
            str(native/"home"),fixture,recorded=True,retained_root=str(native/"evidence/installed-suite/fixtures"))
        E.require(len(contract["source_files"])==93, "Reviewed93 test/helper files required")
        installed_files = E.document(evidence / "installed-before.json")["files"]
        for name, pin in contract["source_files"].items():
            relative = "tests/testthat/" + name
            for recorded in (checkout["files"][relative], installed_files[relative]):
                E.require(recorded["sha256"] == pin["sha256"] and recorded["bytes"] == pin["bytes"],
                          "Installed test/helper bytes do not match checkout and reviewed contract")
        E.require(npath(checks["suite"]["runtime"]["home"],"windows")==native/"home", "Suite did not use private startup home")
        for key in ("package_path","package_version","dll","dll_md5","R_version","R_home"):
            E.require(checks["suite"]["runtime"][key]==installation[key], "Installed suite runtime selected another package")
        checks["native"]=verify_native(evidence,source,checkout,installation,tools,logs,native,extracted/"q")
        checks["science"]=verify_science(evidence,checkout,installation,native,extracted/"q",tools)
        checks["VCF"]=E.verify_vcf_control(evidence/"vcf",npath(installation["library"],"windows"),checkout["commit"],native/"tmp",native_output=native/"evidence/vcf")
    required = {"R-observation","compiler-version","R-config","native-R-runtime", "current-R-CMD-build", "current-R-CMD-INSTALL", "current-installed-identity"}
    required |= {"dependency-"+name for name in installation["dependencies"]}
    required |= {"legacy-"+name for name in OPERATIONS}
    if lane == "windows":
        required |= {"complete-current-installed-suite","scientific-current-references","current-VCF-write-index-query"}
        required |= {"science-"+name+"-f"+str(number)for name in ("hard-cases","declared")for number in (2,3)}
        required |= {"checkout-"+name for name in checks["native"]}
    else: required.add("native-Linux-ownership-controls")
    E.require(set(logs)==required, "Exact complete single-checkout command phase roster differs")
    return installation,evidence,checks


def receive(root, output):
    E.plain_path(root,True)
    E.require(not output.exists() and not output.is_symlink(), "Fresh receiver output required")
    output.mkdir(mode=0o700)
    proof = {"schema":"mss-installed-terminal-decision-v1","qualified":False,"artifact_code_executed":False}
    try:
        service=E.document(root/"service-receipt.json")
        E.require(service["schema"]=="mss-terminal-service-v1" and service["qualified"] is False
                  and service["artifact_code_execution"] is False, "Original independent terminal service receipt required")
        expected_metadata = {"run.json","jobs.json","artifacts.json"}
        for source in service["sources"].values(): expected_metadata |= {source["commit_document"],source["tree_document"]}
        E.require(set(service["terminal_metadata"]) == expected_metadata, "Complete exact independent service metadata census required")
        for name,pin in service["terminal_metadata"].items():
            E.require(Path(name).name==name, "Exact local service metadata leaf required")
            same_pin(pin,root/name)
        run,jobs,artifacts=(E.document(root/name)for name in("run.json","jobs.json","artifacts.json"))
        E.require(run["status"]=="completed" and run["conclusion"]=="success" and str(run["id"])==service["run_id"]
                  and str(run["run_attempt"])==service["run_attempt"] and run["repository"]["full_name"]==service["repository"]
                  and run["path"]==S.WORKFLOW_PATH and len(jobs["jobs"])==jobs["total_count"]==2,
                  "Exact completed two-lane service run required")
        E.require({r["name"]for r in jobs["jobs"]}==set(S.JOB_NAMES.values()) and len({r["id"]for r in jobs["jobs"]})==2
                  and all(str(r["run_id"])==service["run_id"] and r["head_sha"]==run["head_sha"] for r in jobs["jobs"])
                  and artifacts["total_count"]==len(artifacts["artifacts"])<=100 and set(service["artifacts"])==set(S.JOB_NAMES), "Exact producer lanes required")
        for lane_name,retained in service["artifacts"].items():
            matches=[r for r in artifacts["artifacts"]if r["id"]==retained["service"]["id"]]
            E.require(matches==[retained["service"]] and retained["service"]["expired"] is False
                      and str(retained["service"]["workflow_run"]["id"])==service["run_id"]
                      and retained["service"]["workflow_run"]["head_sha"]==run["head_sha"]
                      and retained["service"]["digest"]=="sha256:"+retained["original_zip"]["sha256"]
                      and retained["service"]["name"]=="mss-installed-"+lane_name+"-"+service["run_id"]+"-"+service["run_attempt"], "Exact original selected artifact service identity differs")
        linux,linux_root,linux_checks=lane(root,service,run,jobs,"linux-gold",output)
        windows,windows_root,windows_checks=lane(root,service,run,jobs,"windows",output)
        parity=compare(linux_root/"legacy",windows_root/"legacy",linux,windows,E.document(CI/"closure.json"))
        proof.update(qualified=True,repository=service["repository"],run_id=service["run_id"],run_attempt=service["run_attempt"],
                     checkout_commit=windows["checkout_commit"],checkout_tree=windows["checkout_tree"],
                     coverage={"linux_gold":linux_checks,"windows":windows_checks,"five_legacy_parity":parity},
                     complete_installed_package=True,performance_comparison_claimed=False)
    except Exception as error:
        proof["failure"]={"type":type(error).__name__,"message":str(error)}
        raise
    finally:
        write_json(output/"receiver-result.json",proof)
    return proof


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--service",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    value=receive(args.service.absolute(),args.out.absolute())
    print(json.dumps({key:value[key]for key in("qualified","checkout_commit","run_id","run_attempt")}))


if __name__=="__main__":
    main()
