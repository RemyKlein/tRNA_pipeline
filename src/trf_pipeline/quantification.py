import gzip
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from .genome import normalize_sequence


@dataclass(frozen=True)
class Count:
    sequence: str
    raw_count: int
    rpm_trna_space: float
    rpm_total: float


def _open_text(path: Path):
    return (
        gzip.open(path, "rt", encoding="utf-8")
        if path.suffix == ".gz"
        else path.open(encoding="utf-8")
    )


def quantify_fastq(path: Path, lookup: set[str]) -> list[Count]:
    if not path.is_file():
        raise FileNotFoundError(path)
    counts: Counter[str] = Counter()
    total = 0
    with _open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
            if (
                not sequence
                or not plus
                or not quality
                or not header.startswith("@")
                or not plus.startswith("+")
            ):
                raise ValueError("malformed FASTQ")
            if len(sequence.strip()) != len(quality.strip()):
                raise ValueError("FASTQ sequence and quality lengths differ")
            value = normalize_sequence(sequence, allow_n=True)
            total += 1
            if value in lookup:
                counts[value] += 1
    if total == 0:
        raise ValueError("empty FASTQ")
    assigned = sum(counts.values())
    return [
        Count(seq, count, count * 1e6 / assigned, count * 1e6 / total)
        for seq, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    ]
