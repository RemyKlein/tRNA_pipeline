from collections import defaultdict
from typing import Iterable

from .genome import reverse_complement
from .models import Exclusivity, GenomeHit


def find_hits(
    candidates: Iterable[str], genome: dict[str, str], masks: dict[tuple[str, str], bytearray]
) -> dict[str, list[GenomeHit]]:
    """Index candidates by length and scan each genomic window once per length."""
    by_length: dict[int, set[str]] = defaultdict(set)
    for sequence in candidates:
        by_length[len(sequence)].add(sequence)
    hits = {sequence: [] for values in by_length.values() for sequence in values}
    for chrom, forward in genome.items():
        for strand, sequence in (("+", forward), ("-", reverse_complement(forward))):
            mask = masks[(chrom, strand)]
            for length, wanted in by_length.items():
                for start in range(0, len(sequence) - length + 1):
                    value = sequence[start : start + length]
                    if value in wanted:
                        span = mask[start : start + length]
                        hits[value].append(
                            GenomeHit(
                                chrom, strand, start, start + length, bool(span) and 0 not in span
                            )
                        )
    return hits


def classify_hits(hits: list[GenomeHit], *, has_synthetic_bases: bool = False) -> Exclusivity:
    if not hits:
        return Exclusivity.AMBIGUOUS if has_synthetic_bases else Exclusivity.NOT_FOUND
    inside = any(hit.in_trna_space for hit in hits)
    outside = any(not hit.in_trna_space for hit in hits)
    if inside and outside:
        return Exclusivity.AMBIGUOUS
    if inside:
        return Exclusivity.BONA_FIDE
    return Exclusivity.NON_EXCLUSIVE
