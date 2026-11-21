"""Check the source-only namespace correction against the two retained worker inventories."""
import json
from pathlib import Path
import runpy

wrapper_path = Path(__file__).with_name("namespace-bench-v2.py")
wrapper = runpy.run_path(str(wrapper_path))
binding = wrapper["load_original"]()
warm = binding.load_frozen()
original_identity = binding.source_identity(warm)
case = "hard-cases-f2"
case_dir = Path("/private/tmp/mss-namespace-benchmark-20261121/diagnostic") / case
builds = Path("/private/tmp/mss-namespace-builds-20261121")
config = {role: {"library": str(builds / role / "library")} for role in ("baseline", "candidate")}
states = {}
for role, digest in (
    ("baseline", "040814dc3654a620b7096ab118e8e72d2d75a9204a77cc849f96b6b501fff255"),
    ("candidate", "e6c40018fe19afb387fdb4fe2c0bc2a4f30302208fca91eb56364ac943980c32"),
):
    root = case_dir / f"pair-01-{role}"
    warm.d.check_control(warm.pin(root / "namespaces.tsv", digest))
    states[role] = warm.r.dependency_state(root, config, role)
baseline, candidate = (states[role] for role in ("baseline", "candidate"))
expected_delta = {"AnnotationDbi", "DBI", "GenomicFeatures", "KEGGREST", "RSQLite", "VariantAnnotation",
                  "bit", "bit64", "blob", "cachem", "fastmap", "memoise", "otel", "png"}
assert (len(baseline), len(candidate)) == (85, 71)
assert baseline.keys() - candidate.keys() == expected_delta
assert candidate.keys() <= baseline.keys()
assert candidate["MungeSumstats"] == baseline["MungeSumstats"] == ["1.21.4", "<selected installation>"]
extra = {"testthat": ["3.3.2", "/private/tmp/mss-parallel-support-20261121/common/testthat"]}
installed = {role: {**baseline, **extra} for role in ("baseline", "candidate")}


def rejects(policy, role, state, error):
    try:
        policy.observe(case, role, state)
    except ValueError as exc:
        assert error in str(exc), str(exc)
    else:
        raise AssertionError(f"Accepted {error}")


policy = binding.NamespacePolicy(installed)
policy.observe(case, "baseline", baseline)
rejects(policy, "candidate", candidate, "additional baseline-only namespaces ['otel']")
print("PASS original policy rejects actual namespace inventories (85 baseline, 71 candidate)")

wrapper["apply_override"](binding)
policy = binding.NamespacePolicy(installed)
assert policy.observe(case, "baseline", baseline) == candidate
assert policy.observe(case, "candidate", candidate) == candidate
assert policy.observe(case, "candidate", candidate) == candidate  # candidate-first second pair
assert policy.observe(case, "baseline", baseline) == candidate
for name in sorted(expected_delta):
    missing = {pkg: identity for pkg, identity in baseline.items() if pkg != name}
    rejects(binding.NamespacePolicy(installed), "baseline", missing, f"baseline missing ['{name}']")
    rejects(policy, "candidate", {**candidate, name: baseline[name]}, "changed across repetitions")
policy = binding.NamespacePolicy(installed)
policy.observe(case, "baseline", baseline)
rejects(policy, "candidate", {**candidate, **extra}, "candidate-only")
rejects(policy, "candidate", {**candidate, "otel": baseline["otel"]}, "candidate loaded eager-only")
rejects(policy, "candidate", {**candidate, "Rcpp": [candidate["Rcpp"][0], "/unexpected/Rcpp"]}, "wrong version or origin")
rejects(policy, "candidate", {**candidate, "Rcpp": ["0.0", candidate["Rcpp"][1]]}, "wrong version or origin")
policy.observe(case, "candidate", candidate)
rejects(policy, "candidate", {pkg: identity for pkg, identity in candidate.items() if pkg != "Rcpp"},
        "changed across repetitions")
policy = binding.NamespacePolicy(installed)
policy.observe(case, "baseline", {**baseline, **extra})
rejects(policy, "candidate", candidate, "additional baseline-only namespaces ['testthat']")
print("PASS exact 14-name policy on retained pair and rejection of omissions, extras, path/version and repeat drift")

identity = binding.source_identity(warm)
override = identity["namespace_policy_override"]
assert {key: value for key, value in identity.items() if key != "namespace_policy_override"} == original_identity
assert Path(binding.__file__) == Path(original_identity["binding"]["path"])
assert original_identity["binding"]["sha256"] == "a61cbab2a7af12e0befd8d4fb6c0e9f81d786160db1ad391c95389717a002e94"
assert override["wrapper"] == warm.d.pin_control(wrapper_path)
assert override["audit"] == {"path": str(wrapper_path.with_name("OTEL-NAMESPACE-DELTA.md")),
                             "sha256": "1baeba297218e4d7fca22cda8706019e9093b20f9530ba9f62997e393d603e79"}
assert override["prior_eager_only"] == sorted(expected_delta - {"otel"})
assert override["revised_eager_only"] == sorted(expected_delta)
assert override["added"] == ["otel"]
matrix = json.loads(warm.d.check_control(original_identity["matrix_config"]))
diagnostic = binding.select_config(matrix)
expanded = [(f'{case["id"]}-f{diagnostic["references"][ref]["format"]}' if ref else case["id"], case["pairs"])
            for case in diagnostic["cases"] for ref in case.get("reference_formats", [None])]
assert expanded == [(case, 3) for case in ("hard-cases-f2", "hard-cases-f3", "declared-f2", "declared-f3", "matched-1M")]
assert 2 * sum(pairs for _, pairs in expanded) == 30
assert all(case["target"] is False for case in diagnostic["cases"])
print("PASS original plus wrapper/audit source identity and unchanged five-case, thirty-worker roster")
