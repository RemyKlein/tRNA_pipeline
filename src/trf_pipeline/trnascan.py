from pathlib import Path
from typing import Iterable

from .models import TRNALocus


def parse_trnascan(
    path: Path, accepted_chromosomes: Iterable[str] | None = None
) -> list[TRNALocus]:
    if not path.is_file():
        raise FileNotFoundError(path)
    accepted = set(accepted_chromosomes) if accepted_chromosomes is not None else None
    loci: list[TRNALocus] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(("Sequence", "Name", "-", "#")):
            continue
        fields = line.split()
        if len(fields) < 9:
            raise ValueError(f"malformed tRNAscan-SE row {line_number}: expected >=9 columns")
        try:
            begin, end, ibegin, iend = map(int, (fields[2], fields[3], fields[6], fields[7]))
        except ValueError as exc:
            raise ValueError(f"malformed coordinates on row {line_number}") from exc
        chrom, aa, anticodon = fields[0], fields[4], fields[5].upper()
        if accepted is not None and chrom not in accepted:
            continue
        if anticodon == "NNN":
            continue
        if (ibegin == 0) != (iend == 0):
            raise ValueError(f"incomplete intron coordinates on row {line_number}")
        low, high = sorted((begin, end))
        if ibegin and not (low <= min(ibegin, iend) <= max(ibegin, iend) <= high):
            raise ValueError(f"intron outside tRNA on row {line_number}")
        loci.append(
            TRNALocus(
                f"tRNA_{aa}_{anticodon}_{chrom}:{begin}-{end}",
                chrom,
                begin,
                end,
                aa,
                anticodon,
                ibegin,
                iend,
            )
        )
    if not loci:
        raise ValueError("no usable tRNAscan-SE records")
    return loci
