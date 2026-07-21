import csv
import gzip
import json
from pathlib import Path

import pytest

from trf_pipeline.cli import main
from trf_pipeline.coordinates import reverse_complement_interval
from trf_pipeline.exclusivity import classify_hits, find_hits
from trf_pipeline.fragments import generate_fragments
from trf_pipeline.genome import build_masks, reverse_complement
from trf_pipeline.maturation import extract_and_mature
from trf_pipeline.models import Exclusivity, GenomeHit, TRNALocus
from trf_pipeline.quantification import quantify_fastq
from trf_pipeline.trnascan import parse_trnascan


def locus(**changes):
    values = dict(
        identifier="trna1",
        chromosome="chr1",
        begin=3,
        end=14,
        amino_acid="Ala",
        anticodon="TGC",
        intron_begin=0,
        intron_end=0,
    )
    values.update(changes)
    return TRNALocus(**values)


def test_positive_and_negative_extraction():
    genome = "TTACGTTGCAACGTGG"
    plus = extract_and_mature(locus(), genome)
    minus = extract_and_mature(locus(identifier="minus", begin=14, end=3), genome)
    assert plus.genomic_sequence == genome[2:14]
    assert minus.genomic_sequence == reverse_complement(genome[2:14])


@pytest.mark.parametrize("negative", [False, True])
def test_intron_removal_on_both_strands(negative):
    genome = "TTAAAACCCCGGGGTT"
    item = locus(
        begin=14 if negative else 3, end=3 if negative else 14, intron_begin=7, intron_end=10
    )
    mature = extract_and_mature(item, genome)
    oriented = reverse_complement(genome[2:14]) if negative else genome[2:14]
    assert mature.genomic_sequence == oriented[:4] + oriented[8:]
    assert mature.exon_junction == 4


def test_encoded_and_added_cca():
    encoded = extract_and_mature(locus(begin=1, end=6), "ACGCCA")
    added = extract_and_mature(locus(begin=1, end=6), "ACGTTA")
    assert encoded.sequence == "ACGCCA" and not encoded.added_cca
    assert added.sequence == "ACGTTACCA" and added.added_cca


def test_histidine_minus_one_is_specific():
    his = extract_and_mature(locus(begin=1, end=6, amino_acid="His"), "ACGTTA")
    ala = extract_and_mature(locus(begin=1, end=6, amino_acid="Ala"), "ACGTTA")
    assert his.sequence == "GACGTTACCA" and his.histidine_minus_one
    assert ala.sequence == "ACGTTACCA" and not ala.histidine_minus_one


def test_deterministic_fragments_and_multiple_origins():
    a = extract_and_mature(locus(identifier="b", begin=1, end=8), "ACGTACGT")
    b = extract_and_mature(locus(identifier="a", begin=1, end=8), "ACGTACGT")
    first = generate_fragments([a, b], 4, 4)
    second = generate_fragments([b, a], 4, 4)
    assert [(f.identifier, f.sequence) for f in first] == [
        (f.identifier, f.sequence) for f in second
    ]
    shared = next(f for f in first if f.sequence == "ACGT")
    assert {origin.source_trna_id for origin in shared.origins} == {"a", "b"}


def test_fragment_metadata_for_junction_cca_and_minus_one():
    trna = extract_and_mature(
        locus(begin=1, end=12, amino_acid="His", intron_begin=5, intron_end=6), "ACGTGGACGTTA"
    )
    fragments = generate_fragments([trna], 3, 8)
    assert any(any(o.crosses_spliced_intron for o in f.origins) for f in fragments)
    assert any(any(o.overlaps_added_cca for o in f.origins) for f in fragments)
    assert any(any(o.includes_histidine_minus_one for o in f.origins) for f in fragments)


def test_mask_coordinates_for_both_strands():
    genome = {"chr1": "A" * 20}
    plus = locus(identifier="plus", begin=3, end=6)
    minus = locus(identifier="minus", begin=16, end=13)
    masks = build_masks(genome, [plus, minus])
    assert masks[("chr1", "+")][2:6] == b"\x01" * 4
    assert reverse_complement_interval(12, 16, 20) == (4, 8)
    assert masks[("chr1", "-")][4:8] == b"\x01" * 4


def test_exclusivity_states():
    inside = GenomeHit("chr1", "+", 0, 4, True)
    outside = GenomeHit("chr1", "+", 5, 9, False)
    assert classify_hits([inside]) == Exclusivity.BONA_FIDE
    assert classify_hits([inside, outside]) == Exclusivity.AMBIGUOUS
    assert classify_hits([outside]) == Exclusivity.NON_EXCLUSIVE
    assert classify_hits([]) == Exclusivity.NOT_FOUND
    assert classify_hits([], has_synthetic_bases=True) == Exclusivity.AMBIGUOUS


def test_exact_search_inside_outside_and_absent():
    genome = {"chr1": "AAAACCCCAAAAGGGG"}
    minus_mask = bytearray(16)
    minus_mask[0:4] = b"\x01" * 4
    masks = {("chr1", "+"): bytearray([1] * 8 + [0] * 8), ("chr1", "-"): minus_mask}
    hits = find_hits(["AAAA", "CCCC", "GGGG", "TTTT", "ACAC"], genome, masks)
    assert classify_hits(hits["CCCC"]) == Exclusivity.BONA_FIDE
    assert classify_hits(hits["AAAA"]) == Exclusivity.AMBIGUOUS
    assert classify_hits(hits["GGGG"]) == Exclusivity.NON_EXCLUSIVE
    assert classify_hits(hits["ACAC"]) == Exclusivity.NOT_FOUND


def write_fastq(path: Path, sequences: list[str]):
    text = "".join(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(sequences))
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)


@pytest.mark.parametrize("suffix", [".fq", ".fq.gz"])
def test_quantification_counts_rpm_normalization_and_normalization(tmp_path, suffix):
    path = tmp_path / f"reads{suffix}"
    write_fastq(path, ["acgu", "ACGT", "TTTT", "NNNN"])
    rows = quantify_fastq(path, {"ACGT", "TTTT"})
    assert [(r.sequence, r.raw_count) for r in rows] == [("ACGT", 2), ("TTTT", 1)]
    assert rows[0].rpm_trna_space == pytest.approx(2e6 / 3)
    assert rows[0].rpm_total == 500_000


def test_malformed_inputs(tmp_path):
    bad = tmp_path / "bad.txt"
    bad.write_text("chr1 too few fields\n")
    with pytest.raises(ValueError, match="malformed"):
        parse_trnascan(bad)
    empty = tmp_path / "empty.fq"
    empty.write_text("")
    with pytest.raises(ValueError, match="empty"):
        quantify_fastq(empty, set())
    with pytest.raises(ValueError):
        generate_fragments([], 50, 16)


def test_end_to_end_repeated_output_is_identical(tmp_path):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nACGTACGTACGTACGTACGT\n")
    scan = tmp_path / "scan.txt"
    scan.write_text("chr1 1 1 12 Ala TGC 0 0 50.0\n")
    first, second = tmp_path / "first.tsv", tmp_path / "second.tsv"
    for output in (first, second):
        assert (
            main(
                [
                    "build-lookup",
                    str(genome),
                    str(scan),
                    "--output",
                    str(output),
                    "--min",
                    "4",
                    "--max",
                    "6",
                ]
            )
            == 0
        )
    assert first.read_bytes() == second.read_bytes()
    rows = list(csv.DictReader(first.open(), delimiter="\t"))
    assert rows and all(json.loads(row["origins_json"]) for row in rows)
    fastq = tmp_path / "reads.fq"
    write_fastq(fastq, [rows[0]["sequence"], "TTTT"])
    counts = tmp_path / "counts.tsv"
    assert main(["quantify", str(first), str(fastq), "--output", str(counts)]) == 0
    assert "raw_count" in counts.read_text()
