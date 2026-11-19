"""Untimed, bounded-memory verification and evidence serialization."""

from collections import Counter
from contextlib import ExitStack
import csv
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path


BLOCK_BYTES = 1024 * 1024
NORMALIZED_COLUMNS = (
    "status", "chr", "pos", "ref", "alt", "id", "effect", "eaf", "swapped",
    "complemented", "candidate_count", "entity_type", "entity_id", "normalized_build",
    "source_row", "grch38_chr", "grch38_pos", "grch38_ref", "grch38_alt",
)
DEFAULT_COLUMNS = {"chr": "CHR", "pos": "BP", "id": "SNP", "other": "A1",
                   "effect_allele": "A2", "effect": "BETA", "eaf": "FRQ"}


def write_json(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def read_json(path):
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def stream_digest(handle):
    sha256 = hashlib.sha256()
    md5 = hashlib.md5(usedforsecurity=False)
    size = 0
    for block in iter(lambda: handle.read(BLOCK_BYTES), b""):
        size += len(block)
        sha256.update(block)
        md5.update(block)
    return {"bytes": size, "sha256": sha256.hexdigest(), "md5": md5.hexdigest()}


def digest(path):
    before = path.stat()
    with path.open("rb") as handle:
        result = stream_digest(handle)
    after = path.stat()
    if (before.st_size != result["bytes"] or before.st_size != after.st_size or
            before.st_mtime_ns != after.st_mtime_ns):
        raise ValueError(f"File changed while hashing: {path}")
    return result


def decoded_digest(path, raw):
    with path.open("rb") as handle:
        gzip_input = handle.read(2) == b"\x1f\x8b"
    if gzip_input:
        with gzip.open(path, "rb") as handle:
            return stream_digest(handle)
    return raw


def tsv_pairs(path, columns, *, encoding="utf-8"):
    with path.open(encoding=encoding, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        if next(reader, None) != list(columns):
            raise ValueError(f"Unexpected TSV header: {path.name}")
        values = {}
        for fields in reader:
            if len(fields) != 2 or fields[0] in values:
                raise ValueError(f"Duplicate key or invalid TSV width: {path.name}")
            values[fields[0]] = fields[1]
        return values


def compare_bytes(paths):
    """Read each byte once, preserving hashes even if editions differ."""
    hashes = [(hashlib.sha256(), hashlib.md5(usedforsecurity=False)) for _ in paths]
    sizes = [0] * len(paths)
    mismatches = {}
    with ExitStack() as stack:
        files = [stack.enter_context(path.open("rb")) for path in paths]
        while True:
            blocks = [handle.read(BLOCK_BYTES) for handle in files]
            if not any(blocks):
                break
            for index, block in enumerate(blocks):
                if index and block != blocks[0]:
                    mismatches.setdefault(str(paths[index]), sizes[0] - len(blocks[0]))
                sizes[index] += len(block)
                for hasher in hashes[index]:
                    hasher.update(block)
    result = [{"bytes": size, "sha256": sha.hexdigest(), "md5": md5.hexdigest()}
              for size, (sha, md5) in zip(sizes, hashes)]
    return result, mismatches


def verify_normalized(path, rows_expected, oracle, *, encoding, max_record_bytes):
    # csv.field_size_limit is process-global. Each case may admit records above
    # the stdlib default, including long specialist identifiers and quoted text.
    csv.field_size_limit(max(131072, max_record_bytes))
    with ExitStack() as stack:
        handle = stack.enter_context(path.open(encoding=encoding, newline=""))
        reader = csv.reader(handle, delimiter="\t", strict=True)
        header = next(reader, None)
        if header != list(NORMALIZED_COLUMNS):
            raise ValueError("Unexpected normalized schema")
        row_column = header.index("source_row")
        status_column = header.index("status")
        fields = []
        if oracle is not None:
            gold = stack.enter_context(oracle.open(encoding="utf-8", newline=""))
            gold_reader = csv.reader(gold, delimiter="\t", strict=True)
            gold_header = next(gold_reader, None)
            if (not gold_header or gold_header[0] != "source_row" or
                    len(gold_header) != len(set(gold_header)) or len(gold_header) < 2 or
                    not set(gold_header).issubset(header)):
                raise ValueError("Oracle must start with source_row and name unique normalized fields")
            fields = [header.index(name) for name in gold_header]
        counts = Counter()
        row_number = 0
        for row_number, row in enumerate(reader, start=1):
            if len(row) != len(header) or row[row_column] != str(row_number):
                raise ValueError(f"Normalized width or source_row mismatch at row {row_number}")
            counts[row[status_column]] += 1
            if fields and [row[index] for index in fields] != next(gold_reader, None):
                raise ValueError(f"Independent scientific oracle mismatch at row {row_number}")
        if row_number != rows_expected:
            raise ValueError(f"Expected {rows_expected} rows, observed {row_number}")
        if fields and next(gold_reader, None) is not None:
            raise ValueError("Independent scientific oracle has extra rows")
    return {"normalized_rows_checked": row_number, "status_counts": dict(counts),
            "oracle_fields": gold_header if fields else None,
            "oracle_rows_checked": row_number if fields else 0}


def verify_editions(case, reference, run_dirs, evidence_dir, assets, encoding):
    editions = [directory / "edition" for directory in run_dirs]
    source_name = "source.tsv" if case["settings"]["delimiter"] == "tab" else "source.csv"
    essential = {"source-artifact", source_name, "normalized.tsv", "audit.tsv",
                 "metadata.tsv", "columns.tsv", "reference_checksums.tsv",
                 "checksums.tsv", "manifest.rds"}
    names = {path.name for path in editions[0].iterdir()}
    if not essential.issubset(names) or names - essential - {"reference_contigs.tsv"}:
        raise ValueError("Unexpected edition contents")
    for edition in editions:
        if {path.name for path in edition.iterdir()} != names:
            raise ValueError("Edition inventories differ")
        if any(not path.is_file() or path.is_symlink() for path in edition.iterdir()):
            raise ValueError("Edition inventory contains an unexpected file type")
    inventory = {str(directory): {} for directory in run_dirs}
    mismatches = {}
    for name in sorted(names):
        paths = [edition / name for edition in editions]
        if name == "manifest.rds":
            entries = [digest(path) for path in paths]
        else:
            entries, offsets = compare_bytes(paths)
            mismatches.update(offsets)
        for directory, entry in zip(run_dirs, entries):
            inventory[str(directory)][name] = entry
    write_json(evidence_dir / "artifact-digests.json", inventory)
    if mismatches:
        write_json(evidence_dir / "byte-deltas.json", mismatches)
        raise ValueError("Stable edition bytes differ. See byte-deltas.json for first differing blocks.")

    input_path = Path(case["input"]["path"])
    source_raw = assets[str(input_path)]
    source_decoded = decoded_digest(input_path, source_raw)
    first = inventory[str(run_dirs[0])]
    if first["source-artifact"] != source_raw or first[source_name] != source_decoded:
        raise ValueError("Edition source does not reproduce input or decompressed input")
    checksums = tsv_pairs(editions[0] / "checksums.tsv", ("filename", "md5"), encoding=encoding)
    expected_checksums = {name: first[name]["md5"] for name in names - {"checksums.tsv", "manifest.rds"}}
    if checksums != expected_checksums:
        raise ValueError("Edition checksums disagree with independently hashed bytes")

    metadata = tsv_pairs(editions[0] / "metadata.tsv", ("key", "value"), encoding=encoding)
    want = {"core_version": f'mss-preservation-{reference["format"]}',
            "reference_path": str(Path(reference["directory"])),
            "source_rows": str(case["rows"]), "normalized_rows": str(case["rows"]),
            "input_checksum": source_raw["md5"], "checksum_algorithm": "MD5",
            "source_file": source_name, "projection": "not_attempted",
            "acquisition_receipt": ("verified" if case["settings"]["verify_acquisition_receipt"]
                                    else "not_supplied")}
    want.update({key: str(case["settings"][key]) for key in (
        "frequency_type", "effect_type", "delimiter", "batch_bytes", "max_record_bytes",
    )})
    if case["settings"].get("build") is not None:
        want["reference_build"] = case["settings"]["build"]
    columns = case["settings"].get("columns", {})
    want["source_strand"] = ("column:" + columns["source_strand"]
                             if "source_strand" in columns else "not_declared")
    if case["settings"].get("source_strand") is not None:
        want["source_strand"] = case["settings"]["source_strand"]
    for name, expected in want.items():
        observed = metadata.get(name)
        if name == "reference_path" and observed is not None:
            observed = str(Path(observed).resolve(strict=True))
        if observed != expected:
            raise ValueError(f"Unexpected metadata receipt field: {name}")

    recorded_columns = tsv_pairs(editions[0] / "columns.tsv", ("field", "source_column"), encoding=encoding)
    if recorded_columns != (columns or DEFAULT_COLUMNS):
        raise ValueError("Column receipt disagrees with the public call")
    ref_checksums = tsv_pairs(editions[0] / "reference_checksums.tsv", ("name", "md5"), encoding=encoding)
    headers = {"fasta", "fasta.fai", "keys", "offsets", "records", "rejected.tsv", "source.tsv", "core-version"}
    if "contig_aliases.tsv" in reference["files"]:
        headers.add("contig_aliases.tsv")
    if not headers.issubset(ref_checksums):
        raise ValueError("Missing reference checksums")
    for name, md5 in ref_checksums.items():
        if name == "fasta":
            path = reference["fasta"]["path"]
        elif Path(name).name == name and name not in {"", ".", ".."} and "\\" not in name:
            path = str(Path(reference["directory"]) / name)
        else:
            raise ValueError("Invalid reference receipt name")
        if path not in assets or assets[path]["md5"] != md5:
            raise ValueError(f"Reference checksum absent from pin or incorrect: {name}")

    result = verify_normalized(
        editions[0] / "normalized.tsv", case["rows"],
        Path(case["oracle"]["path"]) if "oracle" in case else None,
        encoding=encoding, max_record_bytes=case["settings"]["max_record_bytes"],
    )
    audit = tsv_pairs(editions[0] / "audit.tsv", ("status", "count"), encoding=encoding)
    if {name: Decimal(count) for name, count in audit.items()} != result["status_counts"]:
        raise ValueError("Status audit disagrees with complete normalized scan")
    result.update(stable_byte_equality=True, edition_count=len(editions),
                  normalized_columns=list(NORMALIZED_COLUMNS),
                  source_artifact=source_raw, decoded_source=source_decoded)
    return result
