"""Compare genuine same-checkout Linux/Windows legacy outputs as data only.

No artifact-provided code, R, compiler, package loading or shell execution.
The terminal receiver still supplies independently bound install/service/source
identities; this comparison cannot grant complete-package acceptance itself.
"""
import csv
import hashlib
from pathlib import Path
import re

from evidence import document, plain_path, require

OPERATIONS = ("allele-flip", "indels", "multiple-rsids", "tabix", "headers")
USED_REFERENCES = {"BSgenome.Hsapiens.1000genomes.hs37d5", "SNPlocs.Hsapiens.dbSNP144.GRCh37"}


def table(path):
    plain_path(path, False)
    require(path.stat().st_size <= 8 * 1024**2, "Bounded tiny legacy decoded table")
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream, delimiter="\t", strict=True))
    require(rows and len(rows) <= 10000 and rows[0] and len(set(rows[0])) == len(rows[0])
            and all(len(row) == len(rows[0]) for row in rows), "Malformed complete decoded table")
    return rows


def files(root):
    plain_path(root, True)
    found = {}
    for path in root.rglob("*"):
        if path.is_dir():
            plain_path(path, True)
            continue
        plain_path(path, False)
        name = path.relative_to(root).as_posix()
        require(len(found) < 512 and path.stat().st_size <= 8 * 1024**2, "Legacy output roster bound")
        found[name] = path
    return found


def version_maps(installation, closure):
    lane = installation["lane"]
    require(lane in {"linux-gold", "windows"} and installation["schema"] == "mss-current-install-v1"
            and installation["normal_current_checkout"] is True
            and re.fullmatch("[0-9a-f]{40}", installation["checkout_commit"]), "Current actual checkout install required")
    code = {name: row["version"] for name, row in closure["packages"].items() if row["kind"] == "code"}
    refs = closure["lanes"][lane]["reference_packages"]
    wanted = {**code, **{name: closure["packages"][name]["version"] for name in refs}}
    require(len(code) == 148 and set(closure["legacy_used_references"]) == USED_REFERENCES
            and len(wanted) == (150 if lane == "linux-gold" else 154)
            and installation["dependencies"] == wanted and "seqminer" not in wanted, "Explicit150/154 closure differs")
    payloads = installation["reference_payloads"]
    require(set(payloads) == set(refs), "Actual complete reference payload inventory required")
    platform = "windows" if lane == "windows" else "linux"
    for name, payload in payloads.items():
        archive = closure["packages"][name]["artifacts"][platform]
        require(payload["archive"] == {"bytes": int(archive["bytes"]), "sha256": archive["sha256"]}
                and payload["version"] == closure["packages"][name]["version"]
                and payload["installed_payload_matches_source"] is True
                and payload["source_payload"] == payload["installed_payload"]
                and payload["source_payload"], "Reference source/installed payload mismatch: " + name)
        for relative, pin in payload["source_payload"].items():
            require(relative.startswith("extdata/") and ".." not in Path(relative).parts
                    and type(pin["bytes"]) is int and pin["bytes"] > 0
                    and re.fullmatch("[0-9a-f]{64}", pin["sha256"]), "Malformed used reference payload")
    return code


def operation(root, expected_operation, lane, installation):
    plain_path(root, True)
    receipt = document(root / "operation.json")
    require(receipt["status"] == "native_legacy_operation_passed"
            and receipt["operation"] == expected_operation
            and receipt["actual_host"] == ("Windows" if lane == "windows" else "Linux")
            and receipt["declared_source_commit"] == installation["checkout_commit"]
            and receipt["actual_package"] == installation["package_version"]
            and receipt["dll_md5"] == installation["dll_md5"]
            and receipt["pointer_bytes"] == 8 and receipt["R_version"] == "4.5.3"
            and receipt["one_thread"] is True
            and receipt["dependency_code_count"] == 148
            and receipt["installed_dependency_count"] == (154 if lane == "windows" else 150)
            and set(receipt["used_references"]) == USED_REFERENCES
            and receipt["seqminer_installed_in_qualified_library"] is False
            and receipt["error"] is None and receipt["closing_error"] is None,
            "Incomplete current-install native legacy receipt")
    checks = receipt["checks"]
    require(checks and all(value is True for value in checks.values()), "Failed/unperformed legacy scientific check")
    observed = files(root)
    expected = {item["path"]: item for item in receipt["files"]}
    require(len(expected) == len(receipt["files"]) and set(observed) == set(expected) | {"operation.json"},
            "Complete operation output inventory differs")
    for name, item in expected.items():
        require(type(item["bytes"]) in {int, float} and item["bytes"] == observed[name].stat().st_size
                and hashlib.md5(observed[name].read_bytes()).hexdigest() == item["md5"],
                "Legacy retained output changed: " + name)
    decoded = files(root / "decoded")
    require(decoded and "decoded-schema.tsv" in observed, "Missing full decoded results/types")
    schema = table(root / "decoded-schema.tsv")
    require(schema[0] == ["table", "ordinal", "name", "typeof", "class"], "Unexpected R decoded schema fields")
    schema_names = {row[0] + ".tsv" for row in schema[1:]}
    require(schema_names == {name for name in decoded if name.endswith(".tsv")}, "R type schema lacks a decoded table")
    for name, path in decoded.items():
        require(Path(name).name == name and Path(name).suffix in {".tsv", ".txt"}, "Unrecognized decoded output")
        if name.endswith(".tsv"):
            rows = table(path)
            fields = [row for row in schema[1:] if row[0] + ".tsv" == name]
            require([row[1] for row in fields] == list(map(str, range(1, len(rows[0]) + 1)))
                    and [row[2] for row in fields] == rows[0], "Decoded columns/types schema differs")
    return receipt, decoded, schema


def compare(linux_root, windows_root, linux_install, windows_install, closure):
    require(linux_install["lane"] == "linux-gold" and windows_install["lane"] == "windows",
            "Genuine Linux gold then native Windows required")
    code = version_maps(linux_install, closure)
    require(code == version_maps(windows_install, closure), "Code package/version maps differ")
    require(linux_install["checkout_commit"] == windows_install["checkout_commit"]
            and linux_install["checkout_tree"] == windows_install["checkout_tree"]
            and linux_install["source_files"] == windows_install["source_files"],
            "Gold and Windows did not normally install the same checkout")
    for name in USED_REFERENCES:
        require(linux_install["reference_payloads"][name] == windows_install["reference_payloads"][name],
                "Used reference source/payload differs across hosts")
    results = []
    for name in OPERATIONS:
        left, a, schema_a = operation(Path(linux_root) / name, name, "linux-gold", linux_install)
        right, b, schema_b = operation(Path(windows_root) / name, name, "windows", windows_install)
        require(left["checks"] == right["checks"] and set(a) == set(b) and schema_a == schema_b,
                "Operation checks/decoded membership/R types differ: " + name)
        for relative in a:
            if relative.endswith(".tsv"):
                require(table(a[relative]) == table(b[relative]), "Exact complete decoded values differ: " + relative)
            else:
                require(a[relative].read_text(encoding="utf-8").splitlines()
                        == b[relative].read_text(encoding="utf-8").splitlines(),
                        "VCF metadata logical lines differ: " + relative)
        results.append({"operation": name, "complete_decoded_files": sorted(a), "checks": sorted(left["checks"]),
                        "exact_values_and_R_types_match": True})
    return {"schema": "mss-legacy-parity-v1", "same_checkout": linux_install["checkout_commit"],
            "shared_code_packages": 148, "linux_installed_packages": 150, "windows_installed_packages": 154,
            "used_references": sorted(USED_REFERENCES), "operations": results, "five_legacy_operations_verified": True,
            "terminal_service_and_complete_package_verification_still_required": True}
