"""Source-only controls for the namespace policy, roster and toolchain metadata."""
import copy
import json
import runpy
from pathlib import Path

source = Path(__file__).with_name("namespace-bench.py")
if not source.exists():
    raise AssertionError("The namespace experiment gate has not been implemented")
module = runpy.run_path(str(source))
NamespacePolicy = module["NamespacePolicy"]
eager_only = (
    "AnnotationDbi", "DBI", "GenomicFeatures", "KEGGREST", "RSQLite", "VariantAnnotation",
    "bit", "bit64", "blob", "cachem", "fastmap", "memoise", "png",
)
common = {"base": ["4.5.3", "/framework/base"], "Rcpp": ["1.0.0", "/common/Rcpp"],
          "MungeSumstats": ["1.0.0", "<selected installation>"]}
installed = {**common, "utils": ["4.5.3", "/framework/utils"],
             **{name: ["2.0.0", f"/common/{name}"] for name in eager_only}}
eager = {**common, **{name: installed[name] for name in eager_only}}
expected = {role: installed for role in ("baseline", "candidate")}
case = "hard-cases-f2"


def changed(state, name, value=None):
    state = copy.deepcopy(state)
    if value is None:
        del state[name]
    else:
        state[name] = value
    return state


def rejects(policy, role, state, reason):
    try:
        policy.observe(case, role, state)
    except ValueError as exc:
        assert reason in str(exc), str(exc)
    else:
        raise AssertionError(f"Accepted {reason}")


policy = NamespacePolicy(expected)
assert policy.observe(case, "baseline", eager) == common
assert policy.observe(case, "candidate", common) == common
assert policy.observe(case, "candidate", common) == common  # pair 2 is candidate-first
assert policy.observe(case, "baseline", eager) == common
assert policy.observe("declared-f3", "baseline", eager) == common
assert policy.observe("declared-f3", "candidate", common) == common
rejects(NamespacePolicy(expected), "baseline", changed(eager, "png"), "baseline missing")
rejects(NamespacePolicy(expected), "candidate", common, "baseline must be observed")
for state, reason in [
    (changed(common, "png", eager["png"]), "candidate loaded eager-only"),
    (changed(common, "utils", installed["utils"]), "candidate-only"),
    (changed(common, "Rcpp"), "additional baseline-only"),
    (changed(common, "base", ["4.5.2", "/framework/base"]), "version or origin"),
    (changed(common, "Rcpp", ["1.0.0", "/unexpected/Rcpp"]), "version or origin"),
    (changed(common, "MungeSumstats", ["1.0.0", "/common/MungeSumstats"]), "version or origin"),
]:
    policy = NamespacePolicy(expected)
    policy.observe(case, "baseline", eager)
    rejects(policy, "candidate", state, reason)
policy = NamespacePolicy(expected)
policy.observe(case, "baseline", eager)
policy.observe(case, "candidate", common)
rejects(policy, "baseline", changed(eager, "utils", installed["utils"]), "changed across repetitions")
rejects(policy, "candidate", changed(common, "Rcpp"), "changed across repetitions")
# Each role can match its own installation, but common metadata must still match across roles.
policy = NamespacePolicy({**expected, "candidate": changed(installed, "MungeSumstats", ["2.0.0", "<selected installation>"])})
policy.observe(case, "baseline", eager)
rejects(policy, "candidate", changed(common, "MungeSumstats", ["2.0.0", "<selected installation>"]), "shared namespace metadata")
print("PASS exact 13-namespace delta and rejection of missing, extra, path, version, and repeat drift")

warm = module["load_frozen"]()
assert len(warm.SOURCES) == 11
assert Path(warm.r.__file__) == warm.BENCH / "run.py"
assert Path(warm.v.__file__) == warm.BENCH / "verify.py"
assert Path(warm.build.__file__) == warm.BUILD
matrix = json.loads(warm.d.check_control(warm.pin(warm.CONFIG, warm.ACCEPTED[warm.CONFIG])))
original = copy.deepcopy(matrix)
diagnostic = module["select_config"](matrix)
assert matrix == original
assert [case["id"] for case in diagnostic["cases"]] == ["hard-cases", "declared", "matched-1M"]
assert [case["pairs"] for case in diagnostic["cases"]] == [3, 3, 3]
assert all(case["target"] is False for case in diagnostic["cases"])
for case in diagnostic["cases"]:
    reference = next(item for item in matrix["cases"] if item["id"] == case["id"])
    assert {k: v for k, v in case.items() if k not in ("pairs", "target")} == {
        k: v for k, v in reference.items() if k not in ("pairs", "target")}
assert diagnostic["references"] == {name: matrix["references"][name]
                                      for name in ("matched_2", "scientific_2", "scientific_3")}
assert 2 * sum(case["pairs"] * len(case.get("reference_formats", [None])) for case in diagnostic["cases"]) == 30
print("PASS authenticated physical frozen modules and preserved five-case, thirty-worker roster")

toolchain_pins = module["toolchain_pins"]
verify_compiler_capture = module["verify_compiler_capture"]
historical_pins = json.loads(warm.d.check_control(matrix["dependency_inventory"]))["build_toolchain_pins"]
sdk = "/Library/Developer/CommandLineTools/SDKs/MacOSX27.2.sdk"
toolchain = {"sdk": sdk, "build_toolchain_pins": historical_pins}
assert list(toolchain_pins(toolchain, historical_pins)) == historical_pins
alternate_sdk = "/Library/Developer/CommandLineTools/SDKs/MacOSX99.0.sdk"
alternate_pin = {"path": f"{alternate_sdk}/SDKSettings.json", "sha256": "a" * 64}
wrong_sdk = {"sdk": alternate_sdk, "build_toolchain_pins": [*historical_pins, alternate_pin]}
wrong_digest = {**toolchain, "build_toolchain_pins": [
    {**pin, "sha256": "a" * 64} if pin["path"] == f"{sdk}/SDKSettings.json" else pin for pin in historical_pins
]}

compiler = ["Apple clang version 21.0.0"]
flags = {"CC": "clang", "CFLAGS": "-O2", "CPPFLAGS": "", "CXX17": "clang++", "CXX17FLAGS": "-O2", "LDFLAGS": ""}
captured = [f"{key} = {value}" for key, value in flags.items()]
declaration = {"compiler": compiler, "compiler_flags": flags}
verify_compiler_capture(declaration, compiler, captured)
missing_capture = [line for line in captured if not line.startswith("LDFLAGS ")]
missing_declaration = {key: value for key, value in flags.items() if key != "LDFLAGS"}

accepted_errors = []
for description, check, args in [
    ("alternate selected SDK with both pins", toolchain_pins, (wrong_sdk, historical_pins)),
    ("changed historical SDK digest", toolchain_pins, (wrong_digest, historical_pins)),
    ("missing required captured flag and null declaration", verify_compiler_capture,
     ({**declaration, "compiler_flags": {**flags, "LDFLAGS": None}}, compiler, missing_capture)),
    ("missing additional captured flag and null declaration", verify_compiler_capture,
     ({**declaration, "compiler_flags": {**flags, "LTO": None}}, compiler, captured)),
    ("null declaration for captured empty flag", verify_compiler_capture,
     ({**declaration, "compiler_flags": {**flags, "LDFLAGS": None}}, compiler, captured)),
    ("missing declared required flag", verify_compiler_capture,
     ({**declaration, "compiler_flags": missing_declaration}, compiler, captured)),
    ("different captured flag value", verify_compiler_capture,
     ({**declaration, "compiler_flags": {**flags, "CFLAGS": "-O3"}}, compiler, captured)),
]:
    try:
        check(*args)
    except ValueError:
        continue
    accepted_errors.append(description)
assert not accepted_errors, f"Accepted invalid toolchain metadata: {', '.join(accepted_errors)}"
print("PASS selected historical SDK and digest, missing/null/mismatched flags, and captured empty flag values")
