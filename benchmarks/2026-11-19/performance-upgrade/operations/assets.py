"""Construct immutable input, reference and matrix receipts on admitted Linux only."""
import argparse
import copy
import csv
import gzip
import hashlib
from pathlib import Path
import subprocess

from build import BENCH, pin
from scientific import DECLARED_MAPPING, MAPPING, generate, write_table
from verify import read_json, write_json

OPS = Path(__file__).resolve().parent
HISTORICAL = BENCH.parents[2] / "inst/verification/fai-safety-stage8/whole-file/profile.json"
PANEL_SIZE = 437458
SEED = "matched-chr22-20261119"
ORDER = ("Sort eligible records by ascending SHA256(UTF8(seed + TAB + current ID)), "
         "breaking digest ties by current ID; seed=" + SEED)
PANEL_HEADER = "SNP CHR BP A1 A2 BETA FRQ SE P".split()


def check_external(path, expected):
    actual = pin(path)
    if actual["sha256"] != expected:
        raise ValueError(f"Historical pin mismatch: {path}")
    return actual


def reference(directory, fasta, number, relocate=False):
    return dict(directory=str(directory), format=number, fasta=pin(fasta), bind_fasta=relocate,
                files={p.name: pin(p)["sha256"] for p in sorted(directory.iterdir()) if p.is_file()})


def cardinality(root, case, path, method, fields, *, distinct=None, nrows=None, extras=None):
    if distinct is None:
        keys = set()
        nrows = 0
        op = gzip.open if path.suffix == ".gz" else open
        with op(path, "rt", encoding="utf-8", newline="") as src:
            source = csv.DictReader(src, delimiter="\t", strict=True)
            for row in source:
                keys.add(tuple(row.get(field) for field in fields))
                nrows += 1
        distinct = len(keys)
    if nrows != case["rows"] or not 1 <= distinct <= nrows:
        raise ValueError(f"Unexpected row or key count for {case['id']}")
    case["input"] = pin(path)
    receipt = {"input_sha256": case["input"]["sha256"], "input_rows": nrows, "distinct_keys": distinct,
               "key_definition": fields if isinstance(fields, str) else
                   "Exact decoded source tuple (" + ", ".join(fields) + "); retain allele order, "
                   "raw missing tokens, contig and ID aliases as supplied; absent fields are null.",
               "method": method, **(extras or {})}
    dest = root / (case["id"] + ".keys.json")
    write_json(dest, receipt)
    case["cardinality_receipt"] = pin(dest)


def prepare_reference(installs, out, inputs, mode):
    cmd = [installs["r_executable"]["path"], "--vanilla", str(OPS / "references.R"),
           installs["baseline"]["library"], *installs["dependency_libraries"], str(out),
           str(inputs / "reference"), mode]
    with (out / (mode + "-prepare.log")).open("xb") as transcript:
        subprocess.run(cmd, stdout=transcript, stderr=subprocess.STDOUT, check=True)


def scientific_cases(installs, cfg, out, inputs, smoke):
    n = generate(out, smoke)
    prepare_reference(installs, out, inputs, "scientific")
    for fmt in [2, 3]:
        name = f"scientific_{fmt}"
        cfg["references"][name] = reference(out / name, out / "scientific.fa", fmt)
    generic = next(c for c in cfg["cases"] if c["id"] == "hard-cases")
    declared = copy.deepcopy(generic)
    declared.update(id="declared", rows=n, source_evidence="Literal synthetic REF/ALT/STRAND "
                    "columns in scientific.py DECLARED; the fixture author supplies the strand.")
    cfg["cases"].append(declared)
    for case in cfg["cases"]:
        name = case["id"]
        if name not in ("hard-cases", "declared", "mixed", "dense"):
            continue
        if smoke and name not in ("hard-cases", "declared"):
            continue
        case["settings"]["columns"] = DECLARED_MAPPING if name == "declared" else MAPPING
        if name in ("hard-cases", "declared"):
            case["oracle"] = pin(out / (name + ".gold.tsv"))
        keys = ["CHR", "BP", "SNP", "A1", "A2", "TYPE"]
        if name == "declared":
            keys += ["REF", "ALT", "STRAND"]
        path = out / (name + (".tsv.gz" if name == "hard-cases" else ".tsv"))
        cardinality(out, case, path, "Exact streaming set census over the generated input", keys)
    write_json(out / "scientific-derivation.json", {
        "authoring_source": pin(OPS / "scientific.py"),
        "reference_source": pin(out / "sites.tsv"), "fasta": pin(out / "scientific.fa"),
        "reference_builder": installs["baseline"]["build_receipt"],
        "method": "Literal source and expected result pairs. Expand explicit null columns, "
                  "source ordinals, source entity identity, the supplied build and absent projections.",
        "cases": {c["id"]: {"input": c["input"], "oracle": c["oracle"],
                            "declarations": c.get("source_evidence", "None")}
                  for c in cfg["cases"] if c["id"] in ("hard-cases", "declared")},
    })


def historical_assets(inputs, history):
    record = history["reference"]
    directory = inputs / "reference/matched-index"
    fasta = inputs / "reference" / Path(record["fasta"]).name
    reference_pins = {k: v["sha256"] for k, v in record["files"].items() if k != "fasta"}
    reference_pins["manifest.rds"] = record["manifest"]["sha256"]
    if {p.name for p in directory.iterdir()} != set(reference_pins):
        raise ValueError("Historical matched directory membership changed")
    for name, wanted in reference_pins.items():
        check_external(directory / name, wanted)
    check_external(fasta, record["files"]["fasta"]["sha256"])
    return directory, fasta


def load_panel(directory):
    with (directory / "source.tsv").open(encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames != "chr pos ref alt id aliases".split():
            raise ValueError("Unexpected reference source")
        panel = []
        all_ids, loci = set(), set()
        for site in reader:
            k = (site["chr"], site["pos"])
            rid, ref, alt = site["id"], site["ref"], site["alt"]
            if (k in loci or rid in all_ids or not rid.startswith("rs") or len(ref) != 1
                    or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt
                    or ref + alt in ("AT", "TA", "CG", "GC")):
                raise ValueError("Expected unique biallelic non-palindromic SNP source panel")
            loci.add(k)
            all_ids.update([rid, *site["aliases"].split(";")])
            panel.append((rid, *k, ref, alt))
    if len(panel) != PANEL_SIZE:
        raise ValueError("Historical source panel has unexpected cardinality")
    return panel, loci, all_ids


def historical_rows(panel, path):
    # Byte-exact generator contract also proves the unique event count for the
    # input itself even when the allele pair is reversed in every other cycle.
    periods = []
    for swapped in (False, True):
        batch = []
        for rid, chrom, pos, ref, alt in panel:
            a, b = (alt, ref) if swapped else (ref, alt)
            fields = (rid, chrom, pos, a, b, "-.125" if swapped else ".125",
                      ".75" if swapped else ".25", ".05", ".01")
            batch.append(("\t".join(fields) + "\n").encode())
        periods.append(b"".join(batch))
    with path.open("rb") as src:
        if src.readline() != ("\t".join(PANEL_HEADER) + "\n").encode():
            raise ValueError("Historical source header changed")
        count = 0
        while True:
            cycle = periods[(count // PANEL_SIZE) % 2]
            data = src.read(len(cycle))
            if not data:
                break
            if not cycle.startswith(data) or not data.endswith(b"\n"):
                raise ValueError("Input deviates from alternating complete reference panel")
            count += data.count(b"\n")
    return count


def matched_cases(installs, cfg, out, inputs, history):
    directory, fasta = historical_assets(inputs, history)
    prepare_reference(installs, out, inputs, "matched")
    cfg["references"]["matched_2"] = reference(directory, fasta, 2, relocate=True)
    cfg["references"]["matched_3"] = reference(out / "matched_3", fasta, 3)
    panel, loci, all_ids = load_panel(directory)
    shuffled = sorted(panel, key=lambda r: (hashlib.sha256((SEED+"\t"+r[0]).encode()).digest(), r[0]))
    profile = {"reference": {"source": str(directory / "source.tsv"), "fasta": str(fasta),
                             "fai": str(directory / "fasta.fai"), "build": "GRCh37"}, "stages": {}}
    for case in cfg["cases"]:
        name = case["id"]
        if name.startswith("matched-"):
            stage = name.split("-")[1]
            path = inputs / "mss" / (name + ".tsv")
            expected = history["stages"][stage]["input"]["sha256"]
            if expected != case["input"]["sha256"]:
                raise ValueError("Template and historical input pin differ")
            check_external(path, expected)
            rows = historical_rows(panel, path)
            cardinality(out, case, path, "Compared every raw input record with the independently "
                        "reconstructed alternating panel from historically pinned source.tsv",
                        "Unique catalog (chromosome, position, REF, ALT, current ID); orient both "
                        "input allele orders to the verified panel event", distinct=len(panel), nrows=rows,
                        extras={"reference_source": pin(directory / "source.tsv")})
            profile["stages"][stage] = dict(rows=rows, input=case["input"])
            continue
        if name not in ("coordinate-only", "id-only", "distinct-shuffled", "negative"):
            continue
        path = out / (name + ".tsv")
        keys = ["SNP", "CHR", "BP"]
        if name == "coordinate-only":
            keys = ["CHR", "BP"]
        elif name == "id-only":
            keys = ["SNP"]
        header = [*keys, *PANEL_HEADER[3:]]
        def rows():
            negative_position = 1
            for i, (rid, chrom, pos, ref, alt) in enumerate(shuffled):
                if name == "negative":
                    while (chrom, str(negative_position)) in loci:
                        negative_position += 1
                    pos = str(negative_position)
                    negative_position += 1
                    rid = f"MSS_ABSENT_{i}"
                    if rid in all_ids:
                        raise ValueError("Synthetic negative ID collides with a source ID or alias")
                values = dict(SNP=rid, CHR=chrom, BP=pos)
                yield [*(values[k] for k in keys), ref, alt, ".125", ".25", ".05", ".01"]
        write_table(path, header, rows())
        cardinality(out, case, path, "Exact set census over generated source columns. " + ORDER,
                    keys, extras={"shuffle_derivation": ORDER, "shuffle_seed": SEED})
        if read_json(Path(case["cardinality_receipt"]["path"]))["distinct_keys"] != PANEL_SIZE:
            raise ValueError(f"Duplicate keys in {name}")
    write_json(out / "snp-oracle.json", profile)


def reduced(config, cases, estimate):
    cfg = copy.deepcopy(config)
    cfg["cases"] = copy.deepcopy(cases)
    used = set()
    for c in cases:
        used.update(c.get("reference_formats", []))
        if "reference" in c:
            used.add(c["reference"])
    cfg["references"] = {k: v for k, v in cfg["references"].items() if k in used}
    cfg["estimated_output_bytes"] = estimate
    return cfg


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["smoke", "full"])
    ap.add_argument("installations", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--inputs", type=Path, default=Path("."))
    ap.add_argument("--admission-note", type=Path, required=True)
    a = ap.parse_args()
    inputs, out = a.inputs.resolve(), a.out.resolve()
    out.mkdir()
    installs = read_json(a.installations)
    cfg = read_json(BENCH / "matrix.example.json")
    cfg.update(installs, contention_note=a.admission_note.read_text().strip(),
               minimum_free_bytes=8_000_000_000, estimated_output_bytes=49_000_000_000)
    smoke = a.mode == "smoke"
    if smoke:
        cfg["cases"] = [c for c in cfg["cases"] if c["id"] == "hard-cases"]
    scientific_cases(installs, cfg, out, inputs, smoke)
    if smoke:
        write_json(out / "smoke.json", reduced(cfg, cfg["cases"], 20_000_000))
        return
    matched_cases(installs, cfg, out, inputs, read_json(HISTORICAL))
    write_json(out / "matrix.json", cfg)
    target = copy.deepcopy(next(c for c in cfg["cases"] if c["id"] == "matched-1M"))
    target.update(pairs=1, target=False)
    write_json(out / "diagnostic.json", reduced(cfg, [target], 1_500_000_000))
    target["legacy_panel"] = pin(inputs / "reference/matched-index/source.tsv")
    write_json(out / "legacy.json", reduced(cfg, [target], 1_500_000_000))


if __name__ == "__main__":
    main()
