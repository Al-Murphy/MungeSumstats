#!/usr/bin/env python3
"""Build the frozen runner's one-case 1M configuration and SNP profile."""
import argparse
import copy
import csv
import json
import os
import stat
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

TOP = Path(__file__).resolve().parents[1]
ROOT = TOP / "fullcall-1m-01"
BENCH = Path("/private/tmp/mss-performance-20261119/repo/benchmarks/2026-11-19/performance-upgrade")
FASTA = TOP / "asset-restore-01/GRCh37.ensembl104.chr22.fa"
COMMON = TOP / "restored-01/common"
MANIFEST = "d15ad358b438fef485ddf0980043556a43f3e4d6b72faca35f3429ce4919ab93"
CONTROLS = TOP / "assets/frozen-controls.tsv"
NAMES = ("matched-1M.keys.json", "snp-profile.json", "diagnostic-1m.json")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--admission-note", type=Path, required=True)
    parser.add_argument("--check", action="store_true", help="Validate using private scratch without publishing metadata")
    args = parser.parse_args()
    require(ROOT.resolve(strict=True) == ROOT and ROOT.stat().st_uid == os.getuid()
            and stat.S_IMODE(ROOT.stat().st_mode) == 0o700, "Expected the private diagnostic root")
    require(not any(os.path.lexists(ROOT / name) for name in NAMES), "Existing 1M metadata prevents publication")
    note = args.admission_note.read_text().strip()
    require(note, "Record the actual local diagnostic admission conditions")
    sys.path.insert(0, str(BENCH / "operations"))
    import assets
    import run
    from build import inventory_library
    from verify import digest, read_json, write_json

    assets_table = TOP / "assets/matrix-assets.tsv"
    require(digest(CONTROLS)["sha256"] == "8c32274747fa08ae820c6d3fba748aebacbeed4ed2ad6ca72b0e47f606a55190"
            and digest(assets_table)["sha256"] == "623ecb1a6f4a9fd556a38741c1f72abc7ef5b999787998e700796c7130802ee2",
            "Frozen control or historical asset table changed")
    with CONTROLS.open() as stream:
        controls = {row["control"]: row for row in csv.DictReader(stream, delimiter="\t")}
    require(all(digest(Path(controls[k]["path"]))["sha256"] == controls[k]["sha256"]
                for k in ("assets", "build", "run", "verify", "process_metrics", "worker", "template", "oracle")),
            "Frozen harness or oracle source changed")
    installations_path = ROOT / "installations.json"
    installs = read_json(installations_path)
    require(assets.pin(installations_path)["sha256"] == "8b5b4d9a104f34da2053e21c99f6e235021088ac695741704f815284ba27578b"
            and installs["dependency_libraries"] == [str(COMMON)], "Installations identity changed")
    inventory = installs["dependency_inventory"]
    inventory_path = Path(inventory["path"])
    require(assets.pin(inventory_path) == inventory, "Dependency inventory pin changed")
    trees = read_json(inventory_path)["installed_dependency_tree_sha256"]
    require(all(inventory_library(Path(path)) == packages for path, packages in trees.items()),
            "Common or framework dependency tree changed")

    with assets_table.open() as stream:
        historic = {row["asset"]: row for row in csv.DictReader(stream, delimiter="\t")}
    refdir, input_path = ROOT / "matched-index", Path(historic["1M"]["pinned_path"])
    stable = {"fasta.fai", "keys", "offsets", "records", "rejected.tsv", "source.tsv", "core-version"}
    require({p.name for p in refdir.iterdir()} == stable | {"manifest.rds"}, "Unexpected reference members")
    paths = {**{name: refdir / name for name in stable}, "fasta": FASTA, "1M": input_path}
    entries = {name: digest(path) for name, path in sorted(paths.items())}
    require(all(not paths[name].is_symlink() and observed == {
        "bytes": int(historic[name]["bytes"]), "sha256": historic[name]["sha256"], "md5": historic[name]["md5"]}
        for name, observed in entries.items()), "Historical input or stable reference bytes changed")
    manifest = refdir / "manifest.rds"
    require(not manifest.is_symlink() and manifest.stat().st_size == 712
            and assets.pin(manifest)["sha256"] == MANIFEST, "Rebuilt reference manifest pin changed")
    acceptance = ROOT / "reference-accepted.json"
    accepted = read_json(acceptance)
    require(accepted["new_manifest"] == {**assets.pin(manifest), "bytes": 712}
            and accepted["stable_seven_members_match_historical"] is True, "Reference acceptance disagrees with actual manifest")
    for role in run.ROLES:
        receipt = ROOT / f"open-reference-{role}.jsonl"
        recorded = accepted["opens"][role]
        terminal = json.loads(receipt.read_text().splitlines()[-1])
        require(assets.pin(receipt) == {"path": recorded["receipt"], "sha256": recorded["receipt_sha256"]}
                and assets.pin(receipt.with_suffix(".stdout"))["sha256"] == recorded["stdout_sha256"]
                and terminal["status"] == "completed" and type(terminal["raw_exit"]) is int and terminal["raw_exit"] == 0
                and terminal["closure"] == "confirmed" and terminal["leader_reaped"] is True,
                f"{role}: unsuccessful or changed verified reference open")
    config = read_json(BENCH / "matrix.example.json")
    case = next(item for item in config["cases"] if item["id"] == "matched-1M")
    case.update(target=False, settings={**run.SETTINGS, **case["settings"]})
    ref = dict(directory=str(refdir), format=2, fasta=assets.pin(FASTA), bind_fasta=True,
               files={**{name: entries[name]["sha256"] for name in sorted(stable)}, "manifest.rds": MANIFEST})
    config.update(installs, cases=[case], references={"matched_2": ref}, minimum_free_bytes=186 * 1024**3,
                  estimated_output_bytes=2 * 1024**3, contention_note=note,
                  reference_provenance={"installations": assets.pin(installations_path), "accepted_reference": assets.pin(acceptance),
                                        "stable_payloads": assets.pin(ROOT / "reference-stable-payloads.json")})
    panel, _, _ = assets.load_panel(refdir)
    rows = assets.historical_rows(panel, input_path)
    with TemporaryDirectory(prefix=".prepare-1m-", dir=ROOT / "tmp") as scratch:
        assets.cardinality(Path(scratch), case, input_path, "Compared every raw input record with the independently "
                           "reconstructed alternating panel from historically pinned source.tsv",
                           "Unique catalog (chromosome, position, REF, ALT, current ID); orient both "
                           "input allele orders to the verified panel event", distinct=len(panel), nrows=rows,
                           extras={"reference_source": assets.pin(refdir / "source.tsv")})
        resolved = copy.deepcopy(config)
        hashes, packages, receipts = run.resolve_configuration(resolved, "diagnostic")
        run.read_provenance(resolved, packages, receipts, run.asset_inventory(hashes))
        keys = run.cardinality(resolved["cases"][0])
        profile = {"reference": dict(source=str(refdir / "source.tsv"), fasta=str(FASTA),
                                     fai=str(refdir / "fasta.fai"), build="GRCh37"),
                   "stages": {"1M": dict(rows=rows, input=case["input"])}}
        if not args.check:
            write_json(ROOT / NAMES[0], keys)
            case["cardinality_receipt"] = assets.pin(ROOT / NAMES[0])
            write_json(ROOT / NAMES[1], profile)
            write_json(ROOT / NAMES[2], config)
    print(f"Verified {rows:,} input rows across {len(panel):,} panel events and {case['pairs']} diagnostic pairs")


if __name__ == "__main__":
    os.umask(0o077)
    main()
