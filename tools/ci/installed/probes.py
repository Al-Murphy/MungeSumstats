"""Current-checkout native/parser/deferral gates; no paired-role benchmark.

The executable probe body is ported from the accepted Windows producer.
All source, compiler, package and driver identities come from the same checkout.
"""
from pathlib import Path
import shutil
from evidence import verify_typed_boundaries, verify_lookup_pressure, verify_parser_gold
from verify import digest, write_json

def pin(path):
    return dict(path=str(Path(path)), **digest(Path(path)))

def run_current_probes(execution, source, installation, library, small, out, compiler, rscript):
    role = "checkout"
    candidate = source
    role_lib = installation["library"]
    scratch = execution.root/"suite-scratch"/"native"
    scratch.mkdir(parents=True)
    phases = {}

    def run(label, argv, driver, driver_role="checkout"):
        if driver_role != "checkout":
            raise ValueError("Current-checkout probe cannot select another role")
        source_pin = dict(role="checkout", commit=installation["checkout_commit"], driver_role="checkout",
                          driver_commit=installation["checkout_commit"], driver=pin(driver))
        log = execution.run(role+"-"+label,argv,provenance=source_pin)
        phases[label] = dict(**source_pin, stdout=pin(log), stderr=pin(log.with_suffix(".stderr")),
                             execution=pin(log.with_suffix(".json")))
        if label.startswith("compile-"):
            phases[label]["compiled_executable"]=pin(Path(argv[-1]))
        return log

    for suite in ("preservation","registered"):
        driver = source/f"tests/standalone/{suite}.R"
        run(suite,[rscript,"--vanilla",driver,source],driver)
    native = [("test_core.cpp","fixtures"),("test_blocks.cpp","blocks"),("projection.cpp","projection.chain"),
              ("test_format3.cpp","format3"),("test_fai.cpp","fai")]
    for filename, arg in native:
        driver = source/"tests/native"/filename
        exe = out/(driver.stem+".exe")
        if not arg.endswith(".chain"):
            (scratch/arg).mkdir()
        run("compile-"+driver.stem,[compiler,"-std=c++17","-Wall","-Wextra","-pedantic","-O2","-UNDEBUG",
                                    driver,"-o",exe],driver)
        run("native-"+driver.stem,[exe,scratch/arg],driver)

    driver = candidate/"tests/standalone/lookup_parity.cpp"
    lookup = out/"lookup.exe"
    options = [compiler,"-std=c++17","-O1","-UNDEBUG","-I",source/"src"]
    run("compile-lookup",[*options,"-DMSS_LOOKUP_COUNTERS",driver,"-o",lookup],driver,"checkout")
    parity = run("lookup-parity",[lookup,"parity",scratch/"lookup-parity"],driver,"checkout")
    pressure = run("lookup-pressure",[lookup,"pressure",scratch/"lookup-pressure"],driver,"checkout")
    for log, name in ((parity,"lookup-parity.stdout"), (pressure,"lookup-pressure.stdout")):
        shutil.copyfile(log,out/name)
    phases["lookup-parity"]["typed_boundary_sha256"] = verify_typed_boundaries(parity)
    phases["lookup-pressure"]["formats"] = verify_lookup_pressure(pressure)

    driver = source/"tests/native/test_normalize_reuse.cpp"
    exe = out/"normalize-reuse.exe"
    run("compile-normalize",[*options,driver,"-o",exe],driver)
    normalize = run("normalize-reuse",[exe,scratch/"normalize-reuse"],driver)
    shutil.copyfile(normalize,out/"normalize-reuse.stdout")

    candidate_ops = source/"benchmarks/2026-11-19/performance-upgrade/operations"
    for mode in ("failure","width","precedence"):
        driver = candidate_ops/("failure_parity.R" if mode == "failure" else "parser_error_order.R")
        dest = out/mode
        extra = [] if mode == "failure" else [mode]
        run("parser-"+mode,[rscript,"--vanilla",driver,role_lib,library,small,dest,*extra],driver,"checkout")
        phases["parser-"+mode]["gold"] = verify_parser_gold(dest/"errors.tsv",mode)

    driver = candidate/"tests/standalone/test_variantannotation_deferral.R"
    for mode in ("startup","bp","parallel-vcf"):
        receipt = out/(mode+".rds")
        run(mode,[rscript,"--vanilla",driver,mode,role_lib,library,receipt],driver,"checkout")
        phases[mode]["receipt"] = pin(receipt)
    write_json(out/"probes.json",phases)
