"""Ordinary checkout installation and portable installed-package producer.

The Windows lane reuses the accepted gated Job controller. Linux gold reuses
the accepted pidfd/subreaper command boundary. Artifacts never grant terminal
acceptance: the separate trusted receiver recomputes every result.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys

CI = Path(__file__).resolve().parent
sys.path.insert(0, str(CI))
import windows_job as w
from checkout import freeze_checkout
from dependencies import acquire, install, reference_payloads
from evidence import document, plain_path, require, verify_installed_suite, verify_vcf_control, verify_als_triplets
from legacy import OPERATIONS, compare, operation, version_maps
from package import build_install, pin
from verify import digest, write_json


def current_time():
    return datetime.now(timezone.utc).isoformat()


def source_framework(source):
    sys.path.insert(0, str(source / "benchmarks/2026-11-19/performance-upgrade/operations"))
    from build import inventory_library
    from run import package_tree
    return inventory_library, package_tree


def observe_runtime(execution, source):
    rscript = Path(shutil.which("Rscript.exe" if os.name == "nt" else "Rscript") or "MISSING").resolve(strict=True)
    log = execution.run("R-observation", [rscript, "--vanilla", "-e",
        'stopifnot(getRversion()=="4.5.3",.Machine$sizeof.pointer==8);writeLines(normalizePath(c(R.home(),R.home("bin")),winslash="/",mustWork=TRUE))'])
    paths = w.small_text(log).splitlines()
    require(len(paths) == 2, "Unexpected current R home/bin observation")
    home, binary = (Path(value).resolve(strict=True) for value in paths)
    r = (binary / "R.exe").resolve(strict=True) if os.name == "nt" else Path(shutil.which("R") or "MISSING").resolve(strict=True)
    compiler = Path("C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe").resolve(strict=True) if os.name == "nt" else Path(shutil.which("g++") or "MISSING").resolve(strict=True)
    compiler_log = execution.run("compiler-version", [compiler, "--version"])
    flags = execution.run("R-config", [r, "CMD", "config", "--all"])
    runtime_path = execution.evidence / "runtime.tsv"
    if os.name == "nt":
        engine = (binary / "Rterm.exe").resolve(strict=True)
        from process_metrics import r_command
        argv = r_command(engine, source / "benchmarks/2026-11-19/performance-upgrade/worker.R", "runtime", str(runtime_path))
    else:
        argv = [rscript, "--vanilla", source / "benchmarks/2026-11-19/performance-upgrade/worker.R", "runtime", runtime_path]
    execution.run("native-R-runtime", argv)
    from verify import tsv_pairs
    runtime = tsv_pairs(runtime_path, ("key", "value"))
    engine = Path(runtime["R_engine"]).resolve(strict=True)
    require(Path(runtime["R_home"]).resolve(strict=True) == home, "R launcher selected a different home")
    native = {name: pin(Path(runtime[name]).resolve(strict=True)) for name in ("R_engine", "R_shared", "BLAS_library", "LAPACK_library")}
    toolchain = dict(schema="mss-current-toolchain-v1",
        Python=dict(command=str(Path(sys.executable).resolve(strict=True)), executable=pin(Path(sys.executable).resolve(strict=True)),
                    version=list(sys.version_info[:3]), implementation=platform.python_implementation()),
        R_launcher=pin(r), Rscript=pin(rscript),
        native_runtime=native, runtime=runtime, compiler=pin(compiler), compiler_version=w.small_text(compiler_log),
        compiler_flags=w.small_text(flags), Rtools="45" if os.name == "nt" else None)
    write_json(execution.evidence / "toolchain.json", toolchain)
    return r, rscript, compiler, engine, home


def stage_legacy_fixtures(execution, source, ci_source):
    folder = execution.evidence / "legacy-inputs"
    folder.mkdir()
    fixtures = folder / "fixtures"
    fixtures.mkdir()
    contract = document(ci_source / "legacy-fixtures.json")
    for name, expected in contract["fixtures"].items():
        original = source / "inst/extdata" / Path(name).name
        require(digest(original) == expected, "Current checkout legacy fixture differs from explicit fixture contract")
        shutil.copyfile(original, folder / name)
        require(digest(folder / name) == expected, "Legacy fixture copy changed")
    shutil.copyfile(ci_source / "legacy-fixtures.json", folder / "INPUTS.json")
    shutil.copyfile(ci_source / "closure.json", folder / "closure.json")
    return folder


def run_legacy(execution, source, installation, dependencies, ci_source, rscript):
    inputs = stage_legacy_fixtures(execution, source, ci_source)
    output = execution.evidence / "legacy"
    output.mkdir()
    installed = execution.evidence / "build/installation-with-payloads.json"
    closure = document(ci_source / "closure.json")
    version_maps(installation, closure)
    results = []
    for name in OPERATIONS:
        execution.run("legacy-" + name, [rscript, "--vanilla", ci_source / "legacy.R", name, installation["lane"],
            installation["library"], dependencies, inputs, output / name, installed],
            timeout=600 if name in OPERATIONS[:3] else 180,
            provenance=dict(checkout_commit=installation["checkout_commit"], driver=pin(ci_source / "legacy.R"),
                            installation=pin(installed)))
        receipt, decoded, schema = operation(output / name, name, installation["lane"], installation)
        results.append(dict(operation=name, checks=receipt["checks"], decoded_files=sorted(decoded), R_type_schema=schema))
    write_json(output / "operations.json", dict(schema="mss-current-legacy-v1",
        checkout_commit=installation["checkout_commit"], lane=installation["lane"], operations=results,
        genuine_Linux_gold=installation["lane"] == "linux-gold", cross_host_parity_still_required=True))
    return output


def qualify(root, evidence, clock):
    request = document(root / "request.json")
    lane = request["lane"]
    if lane == "windows":
        execution = w.Execution(root, evidence, clock)
    else:
        import linux_job
        execution = linux_job.Execution(root, evidence, clock)
    try:
        execution.capacity("current-source-and-metadata", 500_000_000)
        checkout = freeze_checkout(Path(request["checkout"]), evidence / "checkout")
        require(checkout["commit"] == request["checkout_commit"], "Triggering checkout changed before freeze")
        execution.source_pins = {"checkout_commit": checkout["commit"], "checkout_tree": checkout["tree"],
                                 "checkout_receipt": pin(evidence / "checkout/checkout.json")}
        source = Path(checkout["source"])
        ci_source = source / "tools/ci/installed"
        # All ignored-for-R-build CI machinery still belongs to this Git source
        # inventory; installed package code/tests are independently bound by tar.
        closure = document(ci_source / "closure.json")
        if lane == "linux-gold":
            execution.run("native-Linux-ownership-controls", [Path(sys.executable).resolve(strict=True), "-I", "-S", "-B", ci_source / "linux-control.py"],
                          timeout=15, provenance=dict(driver=pin(ci_source / "linux-control.py")))
        write_json(evidence / "runner-provenance.json", request["runner"])
        if lane == "linux-gold":
            system = request["system_packages"]
            actual_system = Path(system["path"])
            require(digest(actual_system) == {key:system[key] for key in ("bytes", "md5", "sha256")}
                    and actual_system.stat().st_size <= 16 * (1 << 20), "Actual Linux system package receipt changed")
            shutil.copyfile(actual_system, evidence / "system-packages.tsv")
        r, rscript, compiler, engine, home = observe_runtime(execution, source)
        inventory_library, package_tree = source_framework(source)
        base_before = inventory_library(home / "library")
        write_json(evidence / "base-library-before.json", base_before)
        rows = acquire(execution, closure, lane, root / "dependency-archives", evidence / "dependency-plan")
        dependencies = root / "dependencies"
        install(execution, rows, root / "dependency-archives", dependencies, r)
        dependency_before = inventory_library(dependencies)
        write_json(evidence / "dependency-library-before.json", dependency_before)
        execution.capacity("ordinary-current-package-build-and-gates", 1_500_000_000)
        installation = build_install(execution, checkout, dependencies, r, rscript, ci_source, evidence)
        installation["reference_payloads"] = reference_payloads(execution, rows, root / "dependency-archives", dependencies)
        installation["reference_payload_inventory_still_required"] = False
        write_json(evidence / "build/installation-with-payloads.json", installation)
        package_before = package_tree(Path(installation["package_path"]))
        write_json(evidence / "installed-before.json", package_before)
        execution.environment["R_LIBS_USER"] = os.pathsep.join([installation["library"], str(dependencies)])
        if lane == "windows":
            suite = evidence / "installed-suite"
            execution.run("complete-current-installed-suite", [rscript, "--vanilla", ci_source / "installed-suite.R",
                installation["library"], dependencies, Path(installation["package_path"]) / "tests/testthat",
                suite, root / "home"], timeout=3600, provenance=dict(driver=pin(ci_source / "installed-suite.R"),
                checkout_commit=checkout["commit"], installed_identity=pin(evidence / "build/installation-with-payloads.json")))
            decision = verify_installed_suite(suite, document(ci_source / "suite-contract.json"))
            fixture = document(ci_source / "legacy-fixtures.json")["fixtures"]["fixtures/ALSvcf.vcf"]
            triplets = verify_als_triplets(suite/"als-fixtures", str(root/"home"), fixture, recorded=False,
                                          retained_root=str(suite/"fixtures"))
            for purpose, files in triplets.items():
                write_json(suite/"als-fixtures"/purpose/"files.json", files)
            decision["retained_ALS_triplets"] = triplets
            write_json(evidence / "installed-suite-decision.json", decision)
            from science import execute as scientific
            scientific(execution, source, installation, dependencies, evidence / "scientific-inputs",
                       evidence / "science", engine, rscript)
            native = evidence / "native"
            native.mkdir()
            from probes import run_current_probes
            run_current_probes(execution, source, installation, dependencies, evidence / "scientific-inputs",
                               native, compiler, rscript)
            vcf = evidence / "vcf"
            execution.run("current-VCF-write-index-query", [rscript, "--vanilla", ci_source / "vcf.R",
                installation["library"], dependencies, checkout["commit"], vcf],
                provenance=dict(checkout_commit=checkout["commit"], driver=pin(ci_source / "vcf.R")))
            verify_vcf_control(vcf, Path(installation["library"]), checkout["commit"], root / "tmp", native_output=vcf)
        legacy = run_legacy(execution, source, installation, dependencies, ci_source, rscript)
        if lane == "windows":
            gold = Path(request["linux_gold"]).resolve(strict=True)
            gold_install = document(gold / "q/evidence/build/installation-with-payloads.json")
            parity = compare(gold / "q/evidence/legacy", legacy, gold_install, installation, closure)
            write_json(evidence / "legacy-parity.json", parity)
        package_after = package_tree(Path(installation["package_path"]))
        dependency_after = inventory_library(dependencies)
        base_after = inventory_library(home / "library")
        write_json(evidence / "installed-after.json", package_after)
        write_json(evidence / "dependency-library-after.json", dependency_after)
        write_json(evidence / "base-library-after.json", base_after)
        require(package_after == package_before and dependency_after == dependency_before and base_after == base_before,
                "Installed package, dependencies or base R changed across gates")
        for row in rows:
            require(digest(root / "dependency-archives" / row["archive"])["sha256"] == row["sha256"],
                    "Frozen public dependency source changed across gates")
        for name, expected in checkout["files"].items():
            actual = digest(source / name)
            require(actual["bytes"] == expected["bytes"] and actual["sha256"] == expected["sha256"],
                    "Frozen checkout changed across build/native/scientific gates")
        clock.remaining("compute")
        write_json(evidence / "qualification-candidate.json", dict(schema="mss-current-candidate-v1",
            checks_passed=True, checkout_commit=checkout["commit"], checkout_tree=checkout["tree"], lane=lane,
            ordinary_installed_suite=lane == "windows", typed_native_VCF_and_science=lane == "windows",
            five_legacy_direct_operations=True, same_checkout_Linux_parity=lane == "windows",
            complete_installed_package_terminal_acceptance=False, receiver_verified=False, finished_utc=current_time()))
    except BaseException as error:
        if not (evidence / "failure.json").exists():
            write_json(evidence / "failure.json", dict(qualified=False, error_type=type(error).__name__,
                message=str(error), finished_utc=current_time(), lane=lane))
        raise


def seal_windows(root, evidence, clock):
    # The qualifier is closed under the accepted owned Job before any reads.
    w.evidence_inventory(evidence, clock, metadata_only=True)
    worker = document(evidence / "qualify-stop.json")
    require(worker["drained"] is True and worker.get("active_after") == 0, "Confirmed complete Job drain required")
    worker_ok = worker["exit_code"] == 0 and worker["reason"] is None and worker["active_before"] == 0
    require(not (evidence / "qualification.json").exists() or document(evidence / "qualification.json")["qualified"] is False,
            "Producer cannot inherit terminal qualification")
    pending = evidence / "qualification.pending"
    write_json(pending, dict(qualified=worker_ok, worker_job_drained=True, worker_qualified=worker_ok,
        scratch_cleanup_passed=worker_ok, deferred_observations={},
        artifact_verification="pending_source_seal_and_receiver_verification"))
    pending.replace(evidence / "qualification.json")
    inventory = w.seal_evidence(evidence, clock)
    write_json(root / "artifact-inventory.json", inventory)
    require(inventory["bytes"] <= w.ARTIFACT_LIMIT, "Evidence exceeds bounded artifact capacity")


def supervise_windows(root, start_file):
    evidence = root / "evidence"
    evidence.mkdir()
    clock = w.job_clock(start_file)
    start = document(start_file)
    request = document(root / "request.json")
    require(w.windows_host() == start["host"], "Current and first-step native Windows host differ")
    clock.attempt = {name: start[name] for name in ("run_id", "run_attempt", "job_id")}
    write_json(root / "deadlines.json", dict(absolute_utc_epoch=clock.wall, upload_reserve_minutes=40))
    control = w.owned_phase(root, clock, "job-control", start_file, Path(__file__))
    if not control["drained"]:
        write_json(root / "UNSEALED.json", dict(qualified=False, reason="Job control drain unconfirmed", control=control))
        return 1
    for suffix in ("stdout", "stderr"):
        (root / ("job-control." + suffix)).rename(evidence / ("job-control." + suffix))
    write_json(evidence / "job-control-stop.json", control)
    if not control.get("control_passed"):
        write_json(evidence / "qualification.json", dict(qualified=False, reason="native Job control failed"))
        write_json(root / "artifact-status.json", dict(qualified=False, native_qualification=False,
            source_sealed=False, receiver_verified=False))
        return 1
    worker = w.owned_phase(root, clock, "qualify", start_file, Path(__file__))
    if not worker["drained"]:
        write_json(root / "UNSEALED.json", dict(qualified=False, reason="qualification Job drain unconfirmed", worker=worker))
        return 1
    for suffix in ("stdout", "stderr"):
        (root / ("qualify." + suffix)).rename(evidence / ("worker." + suffix))
    write_json(evidence / "qualify-stop.json", worker)
    native_ok = worker["exit_code"] == 0 and worker["reason"] is None and worker["active_before"] == 0
    write_json(evidence / "qualification.json", dict(qualified=False, worker_qualified=native_ok,
        worker_job_drained=True, artifact_verification="pending_scratch_cleanup_and_source_seal"))
    sealer = w.owned_phase(root, clock, "seal", start_file, Path(__file__))
    sealed = (sealer["drained"] and sealer.get("exit_code") == 0 and sealer["reason"] is None
              and sealer.get("active_before") == 0 and (root / "artifact-inventory.json").is_file())
    native_ok = native_ok and sealer["drained"] and document(evidence / "qualification.json")["qualified"] is True
    status = dict(qualified=bool(native_ok and sealed), native_qualification=native_ok,
                  source_sealed=bool(sealed), sealer=sealer, receiver_verified=False, lane=request["lane"])
    if sealed:
        with w.physical_directory(root) as fd:
            item, _ = w.transfer_regular(root / "artifact-inventory.json", fd, 16 << 20)
        status["inventory"] = dict(path=str(root / "artifact-inventory.json"), sha256=item["sha256"])
    write_json(root / "artifact-status.json", status)
    return int(not status["qualified"])


def run_linux(root, start_file):
    import io
    import linux_job
    evidence = root / "evidence"
    evidence.mkdir()
    clock = w.job_clock(start_file)
    write_json(root / "deadlines.json", dict(absolute_utc_epoch=clock.wall, native_linux_clock=True))
    error = None
    try:
        qualify(root, evidence, clock)
    except BaseException as caught:
        error = caught
    cleanup = io.StringIO()
    # Last-command receipts already require no cleanup signaling. This final
    # native subreaper closure is independent of the aggregate candidate flag.
    signaled = linux_job.quiesce(None, cleanup)
    closed = not linux_job.owned_processes(linux_job.PROC, os.getpid())
    qualified = error is None and closed and not signaled
    stop = dict(schema="mss-linux-gold-stop-v1", owned_tree_quiescent=closed, cleanup_signaled=signaled,
        final_events=[json.loads(line) for line in cleanup.getvalue().splitlines()], finished_utc=current_time(),
        error=None if error is None else dict(type=type(error).__name__, message=str(error)))
    write_json(root / "linux-stop.json", stop)
    write_json(evidence / "qualification.json", dict(qualified=qualified, lane="linux-gold",
        complete_installed_package=False, terminal_receiver_verified=False))
    inventory = w.seal_evidence(evidence, clock)
    write_json(root / "artifact-inventory.json", inventory)
    write_json(root / "artifact-status.json", dict(qualified=qualified, native_qualification=qualified,
        source_sealed=True, receiver_verified=False, lane="linux-gold",
        inventory=dict(path=str(root / "artifact-inventory.json"), sha256=digest(root / "artifact-inventory.json")["sha256"])))
    return int(not qualified)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--checkout", type=Path)
    parser.add_argument("--job-start", type=Path)
    parser.add_argument("--lane", choices=("windows", "linux-gold"))
    parser.add_argument("--linux-gold", type=Path)
    parser.add_argument("--system-packages", type=Path)
    parser.add_argument("--internal-worker", choices=("qualify", "seal", "job-control"))
    parser.add_argument("--prepare-upload", type=Path)
    args = parser.parse_args()
    require(sys.version_info[:3] == (3, 12, 10), "Pinned native Python3.12.10 required")
    require(platform.machine().lower() in {"amd64", "x86_64"}, "Native x64 required")
    root = args.root.absolute()
    if args.prepare_upload is not None:
        require(not args.internal_worker and args.checkout is None and args.lane is not None,
                "Explicit lane upload-only invocation")
        if not root.exists():
            # Preserve an acquisition/setup failure even if the producer never
            # launched; its real first observation remains separate evidence.
            root.mkdir(mode=0o700)
            (root / "q").mkdir(mode=0o700)
        if not (root / "job-start.json").exists() and args.job_start is not None:
            shutil.copyfile(args.job_start.resolve(strict=True), root / "job-start.json")
        if not (root / "preflight.json").exists():
            write_json(root / "preflight.json", dict(schema="mss-current-preflight-v1",
                scientific_producer_launched=False, lane=args.lane, recorded_utc=current_time()))
        if not (root / "q/request.json").exists() and not (root / "q/UNSEALED.json").exists():
            write_json(root / "q/UNSEALED.json", dict(qualified=False,
                reason="No current-checkout producer request was created; setup/acquisition did not launch"))
        if args.lane == "windows":
            payload, complete = w.prepare_upload_snapshot(root, args.prepare_upload)
        else:
            from snapshot import prepare_linux_snapshot
            payload, complete = prepare_linux_snapshot(root, args.prepare_upload)
        if "GITHUB_OUTPUT" in os.environ:
            with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as stream:
                print("directory=" + str(payload), file=stream)
        else:
            print(json.dumps({"directory": str(payload), "full_source_snapshot": complete}))
        return int(not complete)
    require(args.job_start is not None, "Actual first-step job-start receipt required")
    if args.internal_worker:
        require(os.name == "nt" and args.checkout is None and args.lane is None, "Only native Windows gated workers")
        clock = w.job_clock(args.job_start)
        start = document(args.job_start)
        clock.attempt = {name: start[name] for name in ("run_id", "run_attempt", "job_id")}
        w.wait_for_gate(root, clock, args.internal_worker)
        if args.internal_worker == "qualify":
            qualify(root, root / "evidence", clock)
        elif args.internal_worker == "seal":
            seal_windows(root, root / "evidence", clock)
        else:
            w.job_control_fixture(root, clock)
        return 0
    require(args.checkout is not None and args.lane is not None and not root.exists(),
            "Fresh current-checkout producer root/lane required")
    checkout = args.checkout.resolve(strict=True)
    require(checkout == Path(__file__).resolve().parents[3], "Producer must come from the tested checkout")
    require((args.lane == "windows") == (os.name == "nt"), "Native lane/platform mismatch")
    root.mkdir(mode=0o700)
    q = root / "q"
    q.mkdir(mode=0o700)
    start_file = root / "job-start.json"
    shutil.copyfile(args.job_start.resolve(strict=True), start_file)
    start = document(start_file)
    require(start["lane"] == args.lane, "First-step lane differs")
    from service import bind_running
    runner = bind_running(start, os.environ, checkout)
    commit = subprocess.check_output(["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True, timeout=10).strip()
    request = dict(schema="mss-current-producer-request-v1", checkout=str(checkout), checkout_commit=commit,
        lane=args.lane, runner=runner, linux_gold=str(args.linux_gold.resolve(strict=True)) if args.linux_gold else None)
    if args.lane == "linux-gold":
        require(args.system_packages is not None, "Actual Linux dpkg package receipt required")
        request["system_packages"] = pin(args.system_packages.resolve(strict=True))
    write_json(q / "request.json", request)
    write_json(root / "preflight.json", dict(schema="mss-current-preflight-v1", request=request, started_utc=current_time()))
    require(args.lane != "windows" or args.linux_gold is not None, "Windows requires the same-run Linux gold artifact")
    if args.lane == "linux-gold":
        def interrupted(signum, frame):
            raise KeyboardInterrupt("Native Linux gold interrupted")
        for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            signal.signal(signum, interrupted)
    try:
        if args.lane == "windows":
            return supervise_windows(q, start_file)
        return run_linux(q, start_file)
    except BaseException as error:
        if not (q / "UNSEALED.json").exists():
            write_json(q / "UNSEALED.json", dict(qualified=False, error_type=type(error).__name__, message=str(error)))
        raise


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
