from dataclasses import dataclass, field
from enum import Enum


class Exclusivity(str, Enum):
    BONA_FIDE = "bona_fide"
    AMBIGUOUS = "ambiguous"
    NON_EXCLUSIVE = "non_exclusive"
    NOT_FOUND = "not_found"
    INVALID = "invalid"


@dataclass(frozen=True)
class TRNALocus:
    identifier: str
    chromosome: str
    begin: int
    end: int
    amino_acid: str
    anticodon: str
    intron_begin: int = 0
    intron_end: int = 0

    @property
    def strand(self) -> str:
        return "+" if self.begin <= self.end else "-"

    @property
    def genomic_interval(self) -> tuple[int, int]:
        return min(self.begin, self.end), max(self.begin, self.end)


@dataclass(frozen=True)
class MatureTRNA:
    locus: TRNALocus
    sequence: str
    genomic_sequence: str
    added_cca: bool
    histidine_minus_one: bool
    exon_junction: int | None = None


@dataclass(frozen=True)
class FragmentOrigin:
    source_trna_id: str
    amino_acid: str
    anticodon: str
    chromosome: str
    genomic_locus: str
    strand: str
    mature_start: int
    mature_end: int
    overlaps_added_cca: bool
    includes_histidine_minus_one: bool
    crosses_spliced_intron: bool


@dataclass
class Fragment:
    sequence: str
    origins: list[FragmentOrigin] = field(default_factory=list)

    @property
    def identifier(self) -> str:
        import hashlib

        return f"tRF-{len(self.sequence)}-{hashlib.sha256(self.sequence.encode()).hexdigest()[:12]}"


@dataclass(frozen=True)
class GenomeHit:
    chromosome: str
    strand: str
    start: int
    end: int
    in_trna_space: bool
