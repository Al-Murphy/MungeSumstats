"""Apply the reviewed OpenTelemetry namespace correction to the frozen binding."""
import hashlib
import importlib.util
from pathlib import Path
import sys


RUNTIME = Path("/private/tmp/mss-fullfile-20261121/runtime")
PYTHON = Path("/opt/homebrew/Cellar/python@3.15/3.15.2/Frameworks/Python.framework/Versions/3.15/bin/python3.15")
ORIGINAL = RUNTIME / "namespace-bench.py"
ORIGINAL_PIN = {"path": str(ORIGINAL),
                "sha256": "a61cbab2a7af12e0befd8d4fb6c0e9f81d786160db1ad391c95389717a002e94"}
AUDIT_PIN = {"path": str(RUNTIME / "OTEL-NAMESPACE-DELTA.md"),
             "sha256": "1baeba297218e4d7fca22cda8706019e9093b20f9530ba9f62997e393d603e79"}
PRIOR_EAGER_ONLY = frozenset({"AnnotationDbi", "DBI", "GenomicFeatures", "KEGGREST", "RSQLite",
                             "VariantAnnotation", "bit", "bit64", "blob", "cachem", "fastmap", "memoise", "png"})
REVISED_EAGER_ONLY = PRIOR_EAGER_ONLY | {"otel"}


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
    binding.require(binding.EAGER_ONLY == PRIOR_EAGER_ONLY, "Expected the reviewed prior 13-name namespace policy")
    original_identity = binding.source_identity

    def source_identity(warm):
        identity = original_identity(warm)
        binding.require(identity["binding"] == ORIGINAL_PIN, "Original namespace binding source digest changed")
        warm.d.check_control(AUDIT_PIN)
        return {**identity, "namespace_policy_override": {
            "wrapper": warm.d.pin_control(binding.physical(__file__)),
            "audit": AUDIT_PIN,
            "prior_eager_only": sorted(PRIOR_EAGER_ONLY),
            "revised_eager_only": sorted(REVISED_EAGER_ONLY),
            "added": sorted(REVISED_EAGER_ONLY - PRIOR_EAGER_ONLY),
            "reason": "DBI loads otel in the eager worker. testthat loaded otel in both roles of the earlier probe.",
        }}

    binding.EAGER_ONLY = REVISED_EAGER_ONLY
    binding.source_identity = source_identity


def main():
    binding = load_original()
    apply_override(binding)
    return binding.main()


if __name__ == "__main__":
    sys.exit(main())
