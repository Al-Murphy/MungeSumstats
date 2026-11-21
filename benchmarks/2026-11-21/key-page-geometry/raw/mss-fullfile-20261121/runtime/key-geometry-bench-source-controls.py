"""Source controls using the retained deferred-role namespace inventory."""
import json
from pathlib import Path
import runpy

source = Path(__file__).with_name("key-geometry-bench.py")
wrapper = runpy.run_path(str(source))
binding = wrapper["load_original"]()
warm = binding.load_frozen()
prior_identity = binding.source_identity(warm)
case = "hard-cases-f2"
worker = Path("/private/tmp/mss-namespace-benchmark-v2-20261121/diagnostic") / case / "pair-01-candidate"
config = {"candidate": {"library": "/private/tmp/mss-namespace-builds-20261121/candidate/library"}}
warm.d.check_control(warm.pin(worker / "namespaces.tsv", "e6c40018fe19afb387fdb4fe2c0bc2a4f30302208fca91eb56364ac943980c32"))
state = warm.r.dependency_state(worker, config, "candidate")
assert len(state) == 71 and state["MungeSumstats"] == ["1.21.4", "<selected installation>"]
extra = {"testthat": ["3.3.2", "/private/tmp/mss-parallel-support-20261121/common/testthat"]}
installed = {role: {**state, **extra} for role in ("baseline", "candidate")}
without_rcpp = {name: metadata for name, metadata in state.items() if name != "Rcpp"}


def rejects(policy, role, inventory, error):
    try:
        policy.observe(case, role, inventory)
    except ValueError as exc:
        assert error in str(exc), str(exc)
    else:
        raise AssertionError(f"Accepted {error}")


rejects(binding.NamespacePolicy(installed), "baseline", state, "baseline missing")
wrapper["apply_override"](binding)
policy = binding.NamespacePolicy(installed)
assert policy.observe(case, "baseline", state) == state
rejects(policy, "candidate", without_rcpp, "additional baseline-only namespaces ['Rcpp']")
rejects(policy, "candidate", {**state, **extra}, "candidate-only namespaces ['testthat']")
rejects(policy, "candidate", {**state, "Rcpp": [state["Rcpp"][0], "/unexpected/Rcpp"]}, "wrong version or origin")
rejects(policy, "candidate", {**state, "Rcpp": ["0.0", state["Rcpp"][1]]}, "wrong version or origin")
assert policy.observe(case, "candidate", state) == state
assert policy.observe(case, "candidate", state) == state  # candidate-first second pair
assert policy.observe(case, "baseline", state) == state
rejects(policy, "baseline", without_rcpp, "changed across repetitions")
rejects(policy, "candidate", {**state, **extra}, "changed across repetitions")
print("PASS the retained 71-name inventory matches itself and missing, extra, version, path and repeat drift fail")

identity = binding.source_identity(warm)
override = identity["key_geometry_override"]
assert {key: value for key, value in identity.items() if key not in {"key_geometry_override", "baseline_override"}} == {
    key: value for key, value in prior_identity.items() if key != "baseline_override"}
assert Path(binding.__file__) == Path(identity["binding"]["path"])
assert identity["binding"]["sha256"] == "a61cbab2a7af12e0befd8d4fb6c0e9f81d786160db1ad391c95389717a002e94"
assert identity["baseline_override"] == {**prior_identity["baseline_override"],
                                         "actual": "f39ab1b7f0a9fd643d70097e5764ba7e1fcf257b"}
assert override["original_baseline"] == prior_identity["baseline_override"]["actual"]
assert override["actual_baseline"] == identity["baseline_override"]["actual"]
assert override["prior_eager_only"] == sorted({"AnnotationDbi", "DBI", "GenomicFeatures", "KEGGREST", "RSQLite",
                                              "VariantAnnotation", "bit", "bit64", "blob", "cachem", "fastmap", "memoise", "png"})
assert override["revised_eager_only"] == [] and not binding.EAGER_ONLY
assert override["wrapper"] == warm.d.pin_control(source)
assert override["counter_review"] == {"path": str(source.with_name("KEY-GEOMETRY-RESULT-REVIEW.md")),
                                      "sha256": "d478ab3a248afb1f94cca78c913037735aadea29d8e6bd9c4fbbd539c4f425d0"}

matrix = json.loads(warm.d.check_control(identity["matrix_config"]))
diagnostic = binding.select_config(matrix)
expanded = [(f'{case["id"]}-f{diagnostic["references"][ref]["format"]}' if ref else case["id"], case["pairs"])
            for case in diagnostic["cases"] for ref in case.get("reference_formats", [None])]
assert expanded == [(case, 3) for case in ("hard-cases-f2", "hard-cases-f3", "declared-f2", "declared-f3", "matched-1M")]
assert 2 * sum(pairs for _, pairs in expanded) == 30 and all(case["target"] is False for case in diagnostic["cases"])
print("PASS original, wrapper and counter review pins, deferred baseline, and unchanged five-case, thirty-worker roster")
