"""Four current-checkout scientific editions, using existing public R workers.

This is a functional gate: one package installation, no A/B roles, no speedup
summary and no cache-condition claim. Literal golds and every19-field oracle,
source/audit/checksum/RDS manifest decision are retained independently.
"""
import importlib.util
from pathlib import Path
import sys

from evidence import edition_decision, plain_path, require, stable_edition_manifest
from verify import digest, tsv_pairs, write_json

GOLDS = {
    "hard-cases": (32, "ca4f4d2f59d70f6052a5eb149e4e72011fefc175bef074897f8f15578824a565"),
    "declared": (9, "36248164266feb550b628be89212864e80b3d9c0c3b07c71cf61aedbfcf3319e"),
}


def imported(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def pin(path):
    plain_path(path, False)
    return dict(path=str(path), **digest(path))


def prepare(execution, source, installation, dependencies, output, rscript):
    require(not output.exists(), "Fresh scientific fixture root required")
    output.mkdir()
    operations = source / "benchmarks/2026-11-19/performance-upgrade/operations"
    literal = imported("mss_current_literal_science", operations / "scientific.py")
    literal.generate(output, smoke=True)
    for name, (rows, expected) in GOLDS.items():
        require(digest(output / (name + ".gold.tsv"))["sha256"] == expected,
                "Literal full-field oracle changed: " + name)
    execution.run("scientific-current-references", [rscript, "--vanilla", operations / "references.R",
        installation["library"], dependencies, output, output / "unused-matched", "scientific"],
        provenance=dict(checkout_commit=installation["checkout_commit"], driver=pin(operations / "references.R")))
    references = {}
    for number in (2, 3):
        path = output / ("scientific_" + str(number))
        plain_path(path, True)
        require(all(child.is_file() for child in path.iterdir()), "Only complete flat tiny references accepted")
        references[number] = dict(directory=str(path), format=number, fasta=pin(output / "scientific.fa"),
                                  bind_fasta=False, files={child.name: digest(child)["sha256"] for child in sorted(path.iterdir())})
    source_files = {str(path): digest(path) for path in output.rglob("*") if path.is_file()}
    write_json(output / "fixture-provenance.json", dict(schema="mss-current-scientific-v1",
        checkout_commit=installation["checkout_commit"], literal_source=pin(operations / "scientific.py"),
        builder_source=pin(operations / "references.R"), normal_current_install=True,
        references=references, files=source_files, method="Existing literal scientific.py smoke fixtures; all19 gold fields"))
    return references


def execute(execution, source, installation, dependencies, fixtures, output, r_engine, rscript):
    from process_metrics import r_command
    benchmark = source / "benchmarks/2026-11-19/performance-upgrade"
    existing = imported("mss_current_run_settings", benchmark / "run.py")
    literal = imported("mss_current_literal_mappings", benchmark / "operations/scientific.py")
    require(not output.exists(), "Fresh current science output root")
    output.mkdir()
    references = prepare(execution, source, installation, dependencies, fixtures, rscript)
    assets = {str(path): digest(path) for path in fixtures.rglob("*") if path.is_file()}
    def local(value):
        path = Path(value)
        require(path.resolve(strict=True) == path and (fixtures == path or fixtures in path.parents),
                "Scientific asset escaped its sealed root")
        require(digest(path) == assets[str(path)], "Scientific asset changed")
        return path
    before = {str(path): digest(path) for path in fixtures.rglob("*") if path.is_file()}
    results = []
    for name, (rows, expected_gold) in GOLDS.items():
        source_input = fixtures / (name + (".tsv.gz" if name == "hard-cases" else ".tsv"))
        for number in (2, 3):
            reference = references[number]
            options = dict(existing.SETTINGS, frequency_type="EAF", verify_acquisition_receipt=True, columns=literal.DECLARED_MAPPING if name == "declared" else literal.MAPPING)
            case = dict(id=name, rows=rows, input=pin(source_input), oracle=pin(fixtures / (name + ".gold.tsv")), settings=options)
            run = output / (name + "-f" + str(number))
            run.mkdir()
            settings = dict(mode="diagnostic", profile="none", allocation_threshold=10000,
                input=str(source_input), rows=rows, output=str(run), dll=installation["dll"],
                dll_md5=installation["dll_md5"], reference=reference["directory"],
                core_version="mss-preservation-" + str(number), fasta_pin_path=reference["fasta"]["path"],
                **{"library.0": installation["library"], "library.1": str(dependencies)})
            settings.update({key: value for key, value in options.items()
                             if value is not None and key not in {"columns", "verify_acquisition_receipt"}})
            settings.update({"column." + key: value for key, value in options["columns"].items()})
            settings.update(acquisition_md5=case["input"]["md5"], acquisition_bytes=case["input"]["bytes"])
            existing.write_settings(run / "settings.tsv", settings)
            log = execution.run("science-" + name + "-f" + str(number),
                r_command(r_engine, benchmark / "worker.R", "run", str(run / "settings.tsv")),
                provenance=dict(checkout_commit=installation["checkout_commit"], driver=pin(benchmark / "worker.R"),
                                settings=pin(run / "settings.tsv")))
            edition = run / "edition"
            files = {path.name: digest(path) for path in edition.iterdir() if path.is_file()}
            expected_files = {"source-artifact", "source.tsv", "normalized.tsv", "audit.tsv", "metadata.tsv",
                              "columns.tsv", "reference_checksums.tsv", "checksums.tsv", "manifest.rds"}
            require(expected_files <= set(files) <= expected_files | {"reference_contigs.tsv"}
                    and all(child.is_file() for child in edition.iterdir()),
                    "Complete edition file roster differs")
            checked = edition_decision(edition, case, reference, files, assets, local)
            stable = stable_edition_manifest(edition / "manifest.rds", reference, files, assets)
            metrics = tsv_pairs(run / "metrics.tsv", ("key", "value"))
            require(metrics["dll_md5"] == installation["dll_md5"] and metrics["threads"] == "1"
                    and metrics["source_rows"] == metrics["normalized_rows"] == str(rows),
                    "Scientific worker installed package/row/thread identity differs")
            results.append(dict(case=name, format=number, rows=rows, oracle_sha256=expected_gold,
                source_commit=installation["checkout_commit"], command_stdout=pin(log),
                command_stderr=pin(log.with_suffix(".stderr")), command_execution=pin(log.with_suffix(".json")),
                edition_files=files, stable_manifest=stable, scientific=checked, metrics=metrics))
    cross = {}
    for name in GOLDS:
        paths = [output / (name + "-f" + str(number)) / "edition/normalized.tsv" for number in (2, 3)]
        require(paths[0].read_bytes() == paths[1].read_bytes(), "Full cross-format scientific rows differ")
        cross[name] = [pin(path) for path in paths]
    after = {str(path): digest(path) for path in fixtures.rglob("*") if path.is_file()}
    require(before == after, "Immutable literal/reference fixtures changed during execution")
    write_json(output / "scientific-editions.json", dict(schema="mss-current-scientific-v1",
        single_current_install=True, source_commit=installation["checkout_commit"], editions=results,
        complete_independent_oracle_fields=19, cross_format=cross, immutable_fixtures=True,
        performance_comparison_claimed=False))
    return results
