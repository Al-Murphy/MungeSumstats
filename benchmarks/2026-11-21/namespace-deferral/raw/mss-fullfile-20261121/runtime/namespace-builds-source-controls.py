"""Exercise only producer metadata parsing and the frozen environment allowlist."""
import ast
from datetime import datetime, timezone
from pathlib import Path
import runpy
import signal
import tempfile

root = Path(__file__).parent
source = root / "namespace-builds.py"
code = compile(ast.parse(source.read_bytes(), filename=str(source)), str(source), "exec")
producer = {"__file__": str(source), "__name__": "namespace_builds_metadata_controls"}
exec(code, producer)


def reject(function, *args):
    try:
        function(*args)
    except ValueError:
        return
    raise AssertionError(f"Accepted invalid metadata in {function.__name__}")


flags = producer["compiler_flags"]
required = producer["REQUIRED_FLAGS"]
config = dict.fromkeys(required, "-O2")
config.update(CPPFLAGS="", CXX17="clang++ -arch arm64", CXX17STD="-std=gnu++17")
capture = "\n".join(f"{name} = {value}" for name, value in config.items())
assert flags(capture) == config
for name in required:
    reject(flags, "\n".join(f"{key} = {value}" for key, value in config.items() if key != name))

parse_runtime = producer["parse_runtime"]
runtime = b'"key"\t"value"\n"R_version"\t"4.5.3"\n"locale"\t"en_US.UTF-8"\n'
expected = {"R_version": "4.5.3", "locale": "en_US.UTF-8"}
assert parse_runtime(runtime, expected) == expected
reject(parse_runtime, runtime, {**expected, "locale": "C"})
reject(parse_runtime, runtime + b'"locale"\t"en_US.UTF-8"\n', expected)
reject(parse_runtime, runtime + b'"wide"\t"row"\t"bad"\n', expected)

guard = Path("/private/tmp/mss-performance-20261119/tools/darwin-full-stage185/guard.py")
guard_source = ast.parse(guard.read_bytes(), filename=str(guard))
allowlist = next(node for node in guard_source.body if isinstance(node, ast.Assign)
                 and isinstance(node.targets[0], ast.Name) and node.targets[0].id == "_ENV_KEYS")
allowed = set(allowlist.value.args[0].func.value.value.split())
env = producer["environment"](Path("/private/tmp/fake-producer"), Path("/private/tmp/fake-common"))
assert {"HOME", "TMPDIR"} <= env.keys() <= allowed
assert env["LC_ALL"] == env["LANG"] == "en_US.UTF-8"
assert env["R_LIBS"] == "/private/tmp/fake-common"
assert env["R_LIBS_USER"] == env["R_LIBS_SITE"] == ""
assert env["MAKEFLAGS"] == "-j1"
assert env.keys().isdisjoint(producer["GIT_SETTINGS"])
assert set(producer["GIT_SETTINGS"]) == {"GIT_CONFIG_NOSYSTEM", "GIT_CONFIG_GLOBAL", "GIT_NO_REPLACE_OBJECTS",
                                        "GIT_NO_LAZY_FETCH"}
print("PASS required and empty compiler flags, strict runtime TSV, and native guard environment compatibility")

native = runpy.run_path(str(guard))
write_json = producer["write_json"]
floor, cap = producer["FLOOR"], producer["CAP"]
with tempfile.TemporaryDirectory(dir=root) as directory:
    output = Path(directory)
    config = native["GuardConfig"](output, datetime.now(timezone.utc), 60, producer["RSS"], cap, floor)
    native_guard = native["Guard"](config)
    sample = {"free_bytes": floor + cap, "block_bytes": 4096, "allocated_bytes": 0, "logical_bytes": 0}
    native_guard._space = lambda *, walk: sample
    result = output / "operator-result.json"
    admission = output / "admission.json"

    def refused(path, value, reason, **kwargs):
        try:
            write_json(path, value, native_guard, **kwargs)
        except native["GuardError"] as error:
            assert error.reason == reason, error
        else:
            raise AssertionError(f"Publication succeeded despite {reason}")
        assert not path.exists()

    sample["allocated_bytes"] = cap + 1
    refused(result, {"status": "failed"}, "disk", check_stop=False)
    sample.update(allocated_bytes=0, free_bytes=floor - 1)
    refused(result, {"status": "failed"}, "filesystem_floor", check_stop=False)
    sample["free_bytes"] = floor
    refused(admission, {"admitted": False}, "evidence_reserve")
    refused(result, {"status": "failed"}, "evidence_reserve", check_stop=False)
    sample["free_bytes"] = floor + 1024**2
    write_json(admission, {"admitted": False}, native_guard)
    assert admission.is_file()

    def expire_after_walk(*, walk):
        assert walk
        native_guard.deadline = 0
        return sample

    native_guard._space = expire_after_walk
    refused(output / "PROVENANCE.json", {}, "wall")
    refused(result, {"status": "completed"}, "wall")
    native_guard.deadline = float("inf")
    native_guard._space = lambda *, walk: sample
    stop_check = native_guard._stop_check

    def stop_during_publication():
        if result.is_file():
            native_guard._signal_handler(signal.SIGTERM, None)
        stop_check()

    native_guard._stop_check = stop_during_publication
    refused(result, {"status": "completed"}, "shutdown")
    write_json(result, {"status": "failed"}, native_guard, check_stop=False)
    assert '"status": "failed"' in result.read_text()
print("PASS guarded failure and rejection funding, deadline after walk, and stop during publication")
