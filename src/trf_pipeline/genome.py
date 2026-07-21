from pathlib import Path
from typing import Iterable

from .coordinates import one_based_inclusive_to_half_open, reverse_complement_interval
from .models import TRNALocus

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def normalize_sequence(sequence: str, *, allow_n: bool = True) -> str:
    value = sequence.strip().upper().replace("U", "T")
    allowed = set("ACGTN" if allow_n else "ACGT")
    bad = set(value) - allowed
    if bad:
        raise ValueError(f"unsupported nucleotide(s): {''.join(sorted(bad))}")
    return value


def reverse_complement(sequence: str) -> str:
    return normalize_sequence(sequence).translate(_COMPLEMENT)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    records: dict[str, list[str]] = {}
    current: str | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            current = raw[1:].split()[0]
            if not current or current in records:
                raise ValueError("empty or duplicate FASTA identifier")
            records[current] = []
        elif raw.strip():
            if current is None:
                raise ValueError("FASTA sequence before first header")
            records[current].append(raw.strip())
    if not records:
        raise ValueError("empty FASTA")
    result = {name: normalize_sequence("".join(parts)) for name, parts in records.items()}
    if any(not seq for seq in result.values()):
        raise ValueError("empty FASTA record")
    return result


def build_masks(
    genome: dict[str, str], loci: Iterable[TRNALocus]
) -> dict[tuple[str, str], bytearray]:
    masks = {
        (chrom, strand): bytearray(len(seq))
        for chrom, seq in genome.items()
        for strand in ("+", "-")
    }
    for locus in loci:
        if locus.chromosome not in genome:
            raise ValueError(f"chromosome {locus.chromosome!r} is absent from genome")
        length = len(genome[locus.chromosome])
        low, high = locus.genomic_interval
        start, end = one_based_inclusive_to_half_open(low, high)
        if end > length:
            raise ValueError(f"coordinates outside chromosome bounds for {locus.identifier}")
        if locus.strand == "-":
            start, end = reverse_complement_interval(start, end, length)
        mask = masks[(locus.chromosome, locus.strand)]
        mask[start:end] = b"\x01" * (end - start)
        oriented = (
            genome[locus.chromosome][start:end]
            if locus.strand == "+"
            else reverse_complement(genome[locus.chromosome][low - 1 : high])
        )
        # Synthetic positions used by MINTmap's mask cannot extend beyond a chromosome.
        if start > 0:
            mask[start - 1] = 2
        if not oriented.endswith("CCA") and end + 3 <= length:
            mask[end : end + 3] = b"\x02\x02\x02"
    return masks
