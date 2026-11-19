"""Finite synthetic benchmark-control regressions. Never launch a timed worker."""

from pathlib import Path
import csv
import sys
import unittest
from unittest.mock import patch

BENCH = Path(sys.argv[1]).resolve(strict=True)
OUTPUT = Path(sys.argv[2]).resolve()
OUTPUT.mkdir()
sys.path.insert(0, str(BENCH))
sys.argv = sys.argv[:1]

import run
from verify import digest, verify_normalized, write_json


GOLD_COLUMNS = (
    "source_row", "status", "chr", "pos", "ref", "alt", "id", "effect", "eaf",
    "swapped", "complemented", "candidate_count", "entity_type", "entity_id",
    "normalized_build", "grch38_chr", "grch38_pos", "grch38_ref", "grch38_alt",
)
GOLD = ("1", "swap", "1", "2", "A", "C", "rs42", "0.125", "0.25", "TRUE",
        "FALSE", "1", "small_variant", "rs42", "GRCh37", "NA", "NA", "NA", "NA")


def write_tsv(path, header, values):
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows((header, values))


def pin(path):
    return {"path": str(path), "sha256": digest(path)["sha256"]}


def file_at(path, text="fixture"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def case_at(root, name, rows, *, tags=(), distinct=None, shuffle=True):
    source = pin(file_at(root / f"{name}.input.tsv"))
    receipt = root / f"{name}.keys.json"
    keys = {
        "input_sha256": source["sha256"], "input_rows": rows,
        "distinct_keys": rows if distinct is None else distinct,
        "key_definition": "Synthetic control row key",
        "method": "Deterministic control fixture",
    }
    if shuffle:
        keys["shuffle_derivation"] = "Synthetic permutation, seed 42"
    write_json(receipt, keys)
    case = {"id": name, "rows": rows, "pairs": 1, "tags": list(tags),
            "reference": "fixture2", "input": source, "cardinality_receipt": pin(receipt)}
    if "hard_cases" in tags:
        gold = root / f"{name}.oracle.tsv"
        write_tsv(gold, GOLD_COLUMNS, GOLD)
        case["oracle"] = pin(gold)
    return case


def configuration(root, *, matrix=False):
    r_home = root / "fake-r"
    r_bin = r_home / "bin"
    engine = file_at(r_bin / "exec" / "R")
    launcher = file_at(r_bin / "Rscript")
    shared = file_at(r_home / "lib" / "libR.so")
    toolchain = {"compiler": "synthetic", "compiler_flags": [],
                 "R_version": "synthetic R", "R_platform": "synthetic-platform"}
    config = {
        "schema": 1, "contention_note": "Synthetic control only", "minimum_free_bytes": 0,
        "estimated_output_bytes": 1, "r_executable": pin(launcher), "r_engine": pin(engine),
        "runtime_files": [{"role": "r_shared", **pin(shared)}],
        "dependency_libraries": [], "references": {}, "cases": [],
    }
    for role, source in (("baseline", run.BASELINE), ("candidate", "a" * 40)):
        library = root / role
        package = library / "MungeSumstats"
        file_at(package / "DESCRIPTION", "Package: MungeSumstats\n")
        dll = pin(file_at(package / "libs" / "MungeSumstats.so", role))
        tree = run.package_tree(package)["tree_sha256"]
        build = root / f"{role}-build.json"
        write_json(build, {"source_commit": source, "dll_sha256": dll["sha256"],
                           "package_tree_sha256": tree, "toolchain": toolchain})
        config[role] = {"commit": source, "library": str(library), "dll": dll,
                        "package_tree_sha256": tree, "build_receipt": pin(build)}
    fasta = pin(file_at(root / "reference.fa", ">1\nAC\n"))
    for number in (2, 3):
        refdir = root / f"fixture{number}"
        files = {name: pin(file_at(refdir / name))["sha256"]
                 for name in run.REQUIRED_REFERENCE_FILES}
        config["references"][f"fixture{number}"] = {
            "directory": str(refdir), "format": number, "fasta": fasta,
            "bind_fasta": False, "files": files,
        }
    if matrix:
        for label, rows, pairs in (("1M", 1_000_000, 3), ("10M", 10_000_000, 3), ("40M", 40_000_000, 1)):
            # These are synthetic control labels. The process and scientific
            # boundaries below perform zero source or normalized row work.
            case = case_at(root, label, rows, distinct=1)
            config["cases"].append({**case, "target": True, "pairs": pairs})
        kinds = ("distinct_shuffled", "coordinate_only", "id_only", "negative",
                 "dense", "mixed", "hard_cases")
        case = case_at(root, "supplement-shape", 1, tags=kinds)
        del case["reference"]
        config["cases"].append({**case, "reference_formats": ["fixture2", "fixture3"]})
    else:
        config["cases"].append(case_at(root, "hard", 1, tags=("hard_cases",)))
    return config


class BenchmarkRegressions(unittest.TestCase):
    def setUp(self):
        self.root = OUTPUT / self._testMethodName
        self.root.mkdir()

    def test_hard_case_header_is_rejected_before_any_timed_worker(self):
        config = configuration(self.root)
        gold = self.root / "incomplete-oracle.tsv"
        write_tsv(gold, ("source_row", "normalized_build"), ("1", "GRCh37"))
        config["cases"][0]["oracle"] = pin(gold)
        with patch.object(run, "run_one", side_effect=AssertionError("underspecified gold reached a timed worker")):
            with self.assertRaisesRegex(ValueError, "19|scientific|Hard.case"):
                run.benchmark(config, self.root / "out", suite="diagnostic", profile="none",
                              threshold=65536, configuration_pin={})

    def test_full_oracle_checks_signed_effect_and_eaf(self):
        gold = self.root / "gold.tsv"
        observed = self.root / "normalized.tsv"
        write_tsv(gold, GOLD_COLUMNS, GOLD)
        output_header, output_row = list(GOLD_COLUMNS[1:15]), list(GOLD[1:15])
        output_header += ["source_row", *GOLD_COLUMNS[15:]]
        output_row += ["1", *GOLD[15:]]
        write_tsv(observed, output_header, output_row)
        compared = verify_normalized(observed, 1, gold, encoding="utf-8", max_record_bytes=65536)
        self.assertEqual(compared["oracle_rows_checked"], 1)
        self.assertEqual(compared["oracle_fields"], list(GOLD_COLUMNS))
        for field, value in (("effect", "-0.125"), ("eaf", "0.75")):
            bad = output_row.copy()
            bad[output_header.index(field)] = value
            wrong = self.root / f"{field}.tsv"
            write_tsv(wrong, output_header, bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "oracle mismatch"):
                verify_normalized(wrong, 1, gold, encoding="utf-8", max_record_bytes=65536)

    def test_distinct_case_rejects_repeated_keys(self):
        case = case_at(self.root, "repeated", 2, distinct=1, tags=("distinct_shuffled",))
        with self.assertRaisesRegex(ValueError, "distinct"):
            run.cardinality(case)

    def test_distinct_case_requires_shuffle_provenance(self):
        case = case_at(self.root, "unrecorded", 2, shuffle=False, tags=("distinct_shuffled",))
        with self.assertRaisesRegex(ValueError, "shuffle"):
            run.cardinality(case)

    def test_distinct_case_accepts_unique_keys_with_shuffle_provenance(self):
        case = case_at(self.root, "unique", 2, tags=("distinct_shuffled",))
        self.assertEqual(run.cardinality(case)["shuffle_derivation"], "Synthetic permutation, seed 42")

    def exercise_orchestration(self, path, *, outcome):
        config = configuration(path, matrix=True)
        engine, shared = config["r_engine"]["path"], config["runtime_files"][0]["path"]
        other = file_at(path / "other-r" / "bin" / "exec" / "R")
        other_shared = file_at(path / "other-r" / "lib" / "libR.so")

        def controlled_worker(config, case, directory, role, **kwargs):
            directory.mkdir()
            candidate = role == "candidate"
            file_at(directory / "namespaces.tsv",
                    "package\tversion\tpath\nMungeSumstats\t1.0\t" +
                    str(Path(config[role]["library"]) / "MungeSumstats") + "\n")
            wall = 4 if candidate else 12
            if outcome == "slow_and_unknown" and candidate and case["id"] == "10M":
                wall = 8
            if outcome == "unknown_timing" and candidate and case["id"] == "1M":
                wall = 0
            rss = 200 if outcome == "rss_miss" and candidate else 100
            complete = not (outcome in {"unknown_rss", "slow_and_unknown"} and candidate and case["id"] == "1M")
            observed_engine = other if outcome == "wrong_engine" else Path(engine)
            observed_shared = str(other_shared) if outcome == "wrong_runtime_library" else shared
            stats = {"role": role, "peak_rss_bytes": rss, "wall_seconds": wall,
                     "user_seconds": 1, "system_seconds": 0, "rss_complete": complete}
            metrics = {
                "R_version": "synthetic R", "R_platform": "synthetic-platform", "threads": "1",
                "R_home": str(observed_engine.parent.parent.parent),
                "R_bin": str(observed_engine.parent.parent), "R_engine": str(observed_engine),
                "R_shared": observed_shared, "R_shared_resolution": "mapped",
                "locale": "C", "BLAS": "unavailable", "BLAS_library": "unavailable",
                "LAPACK": "1", "LAPACK_library": "unavailable",
            }
            for phase in ("package_setup", "reference_open", "processing"):
                for unit in ("elapsed", "user.self", "sys.self", "user.child", "sys.child"):
                    metrics[f"{phase}.{unit}"] = str(wall if phase == "processing" and unit == "elapsed" else 0)
            return stats, metrics

        with patch.object(run, "run_one", side_effect=controlled_worker), \
             patch.object(run, "verify_editions", return_value={}), \
             patch.object(run, "verify_manifests"):
            status = run.benchmark(config, path / "out", suite="matrix", profile="none",
                                   threshold=65536, configuration_pin={})
        summary = run.read_json(path / "out" / "summary.json")
        return status, summary

    def test_unknown_rss_leaves_elapsed_pass_indeterminate(self):
        code, summary = self.exercise_orchestration(self.root, outcome="unknown_rss")
        self.assertIsNone(summary["processing_speed_and_rss_target_met"])
        self.assertEqual(code, 3)
        self.assertEqual(summary["target_status"], "indeterminate")

    def test_measured_miss_takes_precedence_over_indeterminate_rss(self):
        code, summary = self.exercise_orchestration(self.root, outcome="slow_and_unknown")
        self.assertIs(summary["processing_speed_and_rss_target_met"], False)
        self.assertEqual(code, 2)

    def test_measured_rss_regression_is_a_miss(self):
        code, summary = self.exercise_orchestration(self.root, outcome="rss_miss")
        self.assertIs(summary["processing_speed_and_rss_target_met"], False)
        self.assertEqual(code, 2)

    def test_unresolved_timing_leaves_target_indeterminate(self):
        code, summary = self.exercise_orchestration(self.root, outcome="unknown_timing")
        self.assertIsNone(summary["processing_speed_and_rss_target_met"])
        self.assertEqual(code, 3)

    def test_speed_and_rss_pass(self):
        code, summary = self.exercise_orchestration(self.root, outcome="pass")
        self.assertIs(summary["processing_speed_and_rss_target_met"], True)
        self.assertEqual(code, 0)

    def test_same_r_version_cannot_substitute_another_engine(self):
        with self.assertRaisesRegex(ValueError, "engine"):
            self.exercise_orchestration(self.root, outcome="wrong_engine")

    def test_runtime_library_must_match_its_declared_role(self):
        with self.assertRaisesRegex(ValueError, "r_shared"):
            self.exercise_orchestration(self.root, outcome="wrong_runtime_library")


if __name__ == "__main__":
    unittest.main(verbosity=2)
