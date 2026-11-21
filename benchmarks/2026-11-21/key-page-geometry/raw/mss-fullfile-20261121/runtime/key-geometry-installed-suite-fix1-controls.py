"""Exercise post-phase stop handling without spawning a child."""
from contextlib import redirect_stderr
import io
import json
from pathlib import Path
import runpy
import signal
import tempfile
from types import SimpleNamespace
from unittest.mock import patch


RUNTIME = Path("/private/tmp/mss-fullfile-20261121/runtime")
SOURCE = RUNTIME / "key-geometry-installed-suite.py"
PRIOR_SOURCE = RUNTIME / "key-geometry-installed-suite.review0.py"
operator = runpy.run_path(str(SOURCE))
native = operator["load_pinned"](operator["NATIVE"], "_controls_guard")
producer = operator["load_pinned"](operator["PRODUCER"], "_controls_producer")


def exercise(source, scenario):
    namespace = runpy.run_path(str(source))
    run = namespace["run"]
    original_read = namespace["read_file"]
    post_reads, post_publications = [], []
    guard = None

    def simulate_phase(owner, *_args, **_kwargs):
        nonlocal guard
        assert guard is None
        guard = owner
        if scenario == "wall at entry":
            guard.deadline = native.time.monotonic()
            return native.PhaseOutcome("error", "confirmed", -15, native.Failure("wall", "test stop"), (), None)
        if scenario == "shutdown at entry":
            guard._shutdown = signal.SIGTERM
            return native.PhaseOutcome("error", "confirmed", -15, native.Failure("shutdown", "test stop"), (), None)
        return native.PhaseOutcome("error", "confirmed", 1, native.Failure("exit", "test failure"), (), None)

    def read(path, *args, **kwargs):
        data = original_read(path, *args, **kwargs)
        if guard is not None:
            post_reads.append(path)
            if scenario == "shutdown during inventory" and path.is_relative_to(namespace["PACKAGE"]):
                guard._shutdown = signal.SIGTERM
            if scenario == "shutdown during remaining controls" and path == namespace["CONTROL"]:
                guard._shutdown = signal.SIGTERM
        return data

    def write(path, value, owner, **kwargs):
        if guard is not None:
            post_publications.append((path.name, kwargs.get("check_stop", True)))
            if scenario == "shutdown at after publication" and path.name == "after.json":
                guard._shutdown = signal.SIGTERM
        return producer.write_json(path, value, owner, **kwargs)

    def loader(path, _name):
        if path == namespace["NATIVE"]:
            return native
        assert path == namespace["PRODUCER"]
        return SimpleNamespace(write_json=write)

    with tempfile.TemporaryDirectory(prefix="mss-suite-fix1-control-", dir="/private/tmp") as base:
        root = Path(base) / "run"
        with patch.dict(run.__globals__, ROOT=root, __file__=str(SOURCE), read_file=read, load_pinned=loader):
            with patch.object(native.Guard, "run_phase", simulate_phase), redirect_stderr(io.StringIO()):
                assert run() == 1
        result = json.loads((root / "RESULT.json").read_bytes())
        after = root / "after.json"
        assert result["status"] == "failed" and not result["scientific_completion"]
        return result, after.exists(), post_reads, post_publications


old_result, old_after, old_reads, old_writes = exercise(PRIOR_SOURCE, "wall at entry")
assert old_result["outcome"]["primary"]["reason"] == "wall" and old_after and old_reads
assert ("after.json", False) in old_writes
print("Original reproduced the wall-entry defect with post-stop reads and after.json publication.")

for scenario, reason in (("wall at entry", "wall"), ("shutdown at entry", "shutdown")):
    result, after, reads, writes = exercise(SOURCE, scenario)
    assert result["outcome"]["primary"]["reason"] == reason
    assert result["post_verification"] == "incomplete"
    assert not after and not reads and writes == [("RESULT.json", False)]
    print(f"Fixed {scenario}: no post-run reads or after.json, funded failure outcome retained.")

result, after, reads, writes = exercise(SOURCE, "shutdown during inventory")
assert result["post_verification"] == "incomplete" and "Supervisor received signal" in result["error"]
assert len([path for path in reads if path.is_relative_to(operator["PACKAGE"])]) == 1
assert not after and writes == [("RESULT.json", False)]
print("Shutdown during installed inventory stopped after one installed read and before after.json.")

result, after, reads, writes = exercise(SOURCE, "shutdown at after publication")
assert result["post_verification"] == "incomplete" and "Supervisor received signal" in result["error"]
assert not after and writes == [("after.json", True), ("RESULT.json", False)]
print("Shutdown at after.json publication was blocked by the reused stop-aware writer.")

result, after, reads, writes = exercise(SOURCE, "shutdown during remaining controls")
assert result["post_verification"] == "incomplete" and "Supervisor received signal" in result["error"]
assert after and reads[-1] == operator["CONTROL"] and operator["QUALIFY"] not in reads
assert writes == [("after.json", True), ("RESULT.json", False)]
print("Shutdown during the remaining control checks skipped subsequent reads and marked verification incomplete.")

result, after, reads, writes = exercise(SOURCE, "normal closed R failure")
assert result["outcome"]["raw_exit"] == 1 and result["post_verification"] == "complete"
assert after and operator["QUALIFY"] in reads and writes == [("after.json", True), ("RESULT.json", False)]
print("Ordinary closed R failure retained the verified post-run inventory and a failed result.")

guard_for_receipts = object.__new__(native.Guard)
guard_for_receipts._shutdown = None
guard_for_receipts.deadline = native.time.monotonic() + 1
cases, counts, qualified, _ = operator["collect_suite"](
    operator["PRIOR"] / "receipts/deferred", guard_for_receipts._stop_check
)
assert len(cases) == counts["tests"] == 269 and counts["passed"] == 3077 and qualified
print("Retained receipt census still parsed as 269 cases and 3077 passing assertions.")
