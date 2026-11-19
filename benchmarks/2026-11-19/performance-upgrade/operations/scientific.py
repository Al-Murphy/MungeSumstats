"""Literal artificial association records and their hand-derived truth.

All high-level row assertions are independent of both package installations.
Synthetic contig names and fixture IDs have no human genomic interpretation.
"""
import csv
import gzip
import io


HEADER = "CHR BP SNP A1 A2 BETA FRQ TYPE".split()
ASSERTIONS = "status chr pos ref alt id effect eaf swapped complemented candidate_count".split()
GOLD_HEADER = ["source_row", *ASSERTIONS, "entity_type", "entity_id", "normalized_build",
               "grch38_chr", "grch38_pos", "grch38_ref", "grch38_alt"]
MAPPING = dict(chr="CHR", pos="BP", id="SNP", other="A1", effect_allele="A2",
               effect="BETA", eaf="FRQ", entity_type="TYPE")
DECLARED_MAPPING = dict(MAPPING, source_ref="REF", source_alt="ALT", source_strand="STRAND")

# Columns before | are the source, after | the literal expected scientific row.
# For an unresolved row, '-' expands to nine NA columns (locus, id, numerical
# fields and orientation). The trailing integer is the canonical candidate pool.
GENERIC = """\
S 2 r雪 C T +.1250 .250 small_variant | match S 2 C T r雪 0.125 0.25 FALSE FALSE 1
S 2 r雪 T C -.125 .750 small_variant | swap S 2 C T r雪 0.125 0.25 TRUE FALSE 1
S 2 r雪 G A .125 .25 small_variant | complement S 2 C T r雪 0.125 0.25 FALSE TRUE 1
S 2 r雪 A G -.125 .75 small_variant | swapped_complement S 2 C T r雪 0.125 0.25 TRUE TRUE 1
T 1 tri A G .125 .25 small_variant | match T 1 A G tri 0.125 0.25 FALSE FALSE 2
T 1 tri G A -.125 .75 small_variant | swap T 1 A G tri 0.125 NA TRUE FALSE 2
T 1 tri T C .125 .25 small_variant | complement T 1 A G tri 0.125 0.25 FALSE TRUE 2
T 1 tri C T -.125 .75 small_variant | swapped_complement T 1 A G tri 0.125 NA TRUE TRUE 2
T 1 tri C G .125 .25 small_variant | alt_alt_contrast - 2
Q 1 quad A G .125 .25 small_variant | ambiguous - 3
NA NA shared C T .125 .25 small_variant | ambiguous - 2
chrS 2 ancien_é C T .125 .25 small_variant | match S 2 C T r雪 0.125 0.25 FALSE FALSE 1
NA NA del AA A .125 .25 small_variant | source_coordinates_missing - 1
I 5 NA AA A .125 .25 small_variant | ambiguous - 1
S 6 r雪 C T .125 .25 small_variant | coordinate_id_conflict - 0
S 8 unseen C A .125 .25 small_variant | unmapped - 0
S 9007199254740992 r雪 C T .125 .25 small_variant | invalid_position - 0
S 2 r雪 C C .125 .25 small_variant | unsupported_allele - 1
S 2 hap01 N <DEL> .125 .25 symbolic_sv | specialist_preserved - 0
S 2 fusion C C]X:1] .125 .25 breakend | specialist_preserved - 0
S 2 HLA*01:01 A C .125 .25 hla | specialist_preserved - 0
S 2 STR_雪 A AC .125 .25 repeat | specialist_preserved - 0
S 2 burden"é A C .125 .25 gene_burden | specialist_preserved - 0
S 2 r雪 C T 1e999 2 small_variant | match S 2 C T r雪 NA NA FALSE FALSE 1
S 2 r雪 C T 1e-999 0 small_variant | match S 2 C T r雪 0 0 FALSE FALSE 1
S 2 r雪 C T -2.500 1.0 small_variant | match S 2 C T r雪 -2.5 1 FALSE FALSE 1
S 2 r雪 C T bad 0.5 small_variant | match S 2 C T r雪 NA 0.5 FALSE FALSE 1
S 2 shared C T 2.0 .25 small_variant | match S 2 C T r雪 2 0.25 FALSE FALSE 1
S 6 shared C T .125 .25 small_variant | match S 6 C T second 0.125 0.25 FALSE FALSE 1
NA NA ancien_é C T .125 .25 small_variant | match S 2 C T r雪 0.125 0.25 FALSE FALSE 1
Q 1 quad T C .125 .25 small_variant | ambiguous - 3
T 1 tri A G 0x1p-3 .25 small_variant | match T 1 A G tri 0.125 0.25 FALSE FALSE 2
"""

# Here only the source authors declare REF/ALT/STRAND. The underlying original
# site and allele *pair* are supplied independently of the measured allele.
DECLARED = """\
I 5 del A AA -.125 .75 small_variant AA A + | swap I 1 CA C del 0.125 0.25 TRUE FALSE 1
I 4 del AA A .125 .25 small_variant AA A + | match I 1 CA C del 0.125 0.25 FALSE FALSE 1
Q 1 quad A G .125 .25 small_variant A G + | match Q 1 A G quad 0.125 0.25 FALSE FALSE 1
Q 1 quad G A -.125 .75 small_variant A G + | swap Q 1 A G quad 0.125 NA TRUE FALSE 1
Q 1 quad T C .125 .25 small_variant T C - | complement Q 1 A G quad 0.125 0.25 FALSE TRUE 1
Q 1 quad C T -.125 .75 small_variant T C - | swapped_complement Q 1 A G quad 0.125 NA TRUE TRUE 1
Q 1 quad A G .125 .25 small_variant A G unknown | source_strand_unknown - 0
S 1 unknown A G .125 .25 small_variant A G + | sequence_only S 1 A G NA 0.125 0.25 FALSE FALSE 0
S 1 r雪 A G .125 .25 small_variant A G + | source_id_conflict - 0
"""


def write_table(path, header, rows):
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def fixture_rows(text, declared=False):
    rows, truth = [], []
    for ordinal, line in enumerate(text.splitlines(), 1):
        source, expected = [part.split() for part in line.split("|")]
        if len(source) != len(HEADER) + (3 if declared else 0):
            raise ValueError(line)
        if expected[1] == "-":
            expected = [expected[0], *(["NA"] * 9), expected[2]]
        if len(expected) != len(ASSERTIONS):
            raise ValueError(line)
        # P is preserved verbatim, not mapped. The comment exercises CSV quoting.
        rows.append([*source, "1e-1000", '雪\t"quoted"'])
        truth.append([ordinal, *expected, source[7], source[2], "SYNTHETIC", *(["NA"] * 4)])
    return rows, truth


def generate(root, smoke):
    with (root / "scientific.fa").open("x") as f:
        f.write(">S\nACGTACGT\n>T\nACG\n>Q\nACG\n>I\nCAAAAACGT\n>D\nACG\n")
    sites = [("S", 2, "C", "T", "r雪", "ancien_é;shared"), ("S", 6, "C", "T", "second", "shared"),
             ("T", 1, "A", "C,G", "tri", "."), ("Q", 1, "A", "C,G,T", "quad", "."),
             ("I", 4, "AA", "A", "del", "."), ("I", 5, "AA", "A", "del", ".")]
    # A separate dense locus has 96 actual source records, shared ID/aliases,
    # repeated contexts, and three different canonical events.
    sites.extend(("D", 1, "A", "CGT"[i % 3], "dense", f"alias{i};fan") for i in range(96))
    write_table(root / "sites.tsv", "chr pos ref alt id aliases".split(), sites)
    generic, generic_truth = fixture_rows(GENERIC)
    assert len(generic) == 32
    declared, declared_truth = fixture_rows(DECLARED, declared=True)
    header = [*HEADER, "P", "COMMENT"]
    # Make exactly reproducible gzip bytes with CRLF between records and no final LF.
    raw = io.StringIO(newline="")
    writer = csv.writer(raw, delimiter="\t", lineterminator="\r\n")
    writer.writerow(header)
    writer.writerows(generic)
    with (root / "hard-cases.tsv.gz").open("xb") as dest:
        with gzip.GzipFile(filename="", mode="wb", fileobj=dest, mtime=0) as gz:
            gz.write(raw.getvalue()[:-2].encode("utf-8"))
    write_table(root / "hard-cases.gold.tsv", GOLD_HEADER, generic_truth)
    write_table(root / "declared.tsv", [*HEADER, "REF", "ALT", "STRAND", "P", "COMMENT"], declared)
    write_table(root / "declared.gold.tsv", GOLD_HEADER, declared_truth)
    if smoke:
        return len(declared)
    write_table(root / "mixed.tsv", header, (generic[i % 32] for i in range(100000)))
    dense_ids = ["dense", "fan", *(f"alias{i}" for i in range(96))]
    pairs = [("A", "G"), ("T", "C"), ("C", "G"), ("G", "A")]
    write_table(root / "dense.tsv", header,
                (["D", 1, dense_ids[i % 98], *pairs[i % 4], ".125", ".25", "small_variant",
                  "1e-1000", "dense"] for i in range(10000)))
    return len(declared)
