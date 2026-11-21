"""Bind the matched public-API benchmark for the key-page geometry trial."""
import hashlib
import importlib.util
from pathlib import Path
import sys


RUNTIME = Path("/private/tmp/mss-fullfile-20261121/runtime")
PYTHON = Path("/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15")
ORIGINAL = RUNTIME / "namespace-bench.py"
ORIGINAL_PIN = {"path": str(ORIGINAL),
                "sha256": "a61cbab2a7af12e0befd8d4fb6c0e9f81d786160db1ad391c95389717a002e94"}
COUNTER_REVIEW_PIN = {"path": str(RUNTIME / "KEY-GEOMETRY-RESULT-REVIEW.md"),
                      "sha256": "d478ab3a248afb1f94cca78c913037735aadea29d8e6bd9c4fbbd539c4f425d0"}
PRIOR_BASELINE = "14ea7af51a19272c7792fbe582a3d42b291b4a1a"
BASELINE = "f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b"
PRIOR_EAGER_ONLY = frozenset({"AnnotationDbi", "DBI", "GenomicFeatures", "KEGGREST", "RSQLite",
                             "VariantAnnotation", "bit", "bit64", "blob", "cachem", "fastmap", "memoise", "png"})


def load_original():
    if not (Path(sys.executable) == PYTHON and sys.flags.isolated and sys.dont_write_bytecode and sys.flags.no_site):
        raise ValueError("Use the physical pinned Python executable with -I -B -S")
    if ORIGINAL.resolve(strict=True) != ORIGINAL or not ORIGINAL.is_file():
        raise ValueError("Expected the original physical namespace binding source")
    with ORIGINAL.open("rb") as stream:
        data = stream.read(1024**2 + 1)
    if len(data) > 1024**2 or hashlib.sha256(data).hexdigest() != ORIGINAL_PIN["sha256"]:
        raise ValueError("Original namespace binding source digest changed")
    spec = importlib.util.spec_from_file_location("namespace_bench_original", ORIGINAL)
    binding = importlib.util.module_from_spec(spec)
    exec(compile(data, str(ORIGINAL), "exec"), binding.__dict__)
    return binding


def apply_override(binding):
    binding.require(binding.BASELINE == PRIOR_BASELINE and binding.EAGER_ONLY == PRIOR_EAGER_ONLY,
                    "Expected the reviewed original baseline and exact prior 13-name namespace policy")
    original_identity = binding.source_identity

    def source_identity(warm):
        identity = original_identity(warm)
        binding.require(identity["binding"] == ORIGINAL_PIN, "Original namespace binding source digest changed")
        warm.d.check_control(COUNTER_REVIEW_PIN)
        return {**identity, "key_geometry_override": {
            "wrapper": warm.d.pin_control(binding.physical(__file__)),
            "original_baseline": PRIOR_BASELINE,
            "actual_baseline": BASELINE,
            "prior_eager_only": sorted(PRIOR_EAGER_ONLY),
            "revised_eager_only": sorted(binding.EAGER_ONLY),
            "reason": "Both roles defer VariantAnnotation. Key-page geometry requires identical namespace inventories.",
            "counter_review": COUNTER_REVIEW_PIN,
        }}

    binding.BASELINE = BASELINE
    binding.EAGER_ONLY = frozenset()
    binding.source_identity = source_identity


def main():
    binding = load_original()
    apply_override(binding)
    return binding.main()


if __name__ == "__main__":
    sys.exit(main())
