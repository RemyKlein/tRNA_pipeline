from collections import defaultdict
from typing import Iterable

from .models import Fragment, FragmentOrigin, MatureTRNA


def generate_fragments(
    trnas: Iterable[MatureTRNA], minimum: int = 16, maximum: int = 50
) -> list[Fragment]:
    if minimum < 1 or minimum > maximum:
        raise ValueError("minimum fragment length must be positive and <= maximum")
    grouped: dict[str, set[FragmentOrigin]] = defaultdict(set)
    for trna in sorted(trnas, key=lambda item: item.locus.identifier):
        has_minus_one = trna.minus_one_base is not None
        genomic_mature_end = (1 if has_minus_one else 0) + len(trna.genomic_sequence)
        cca_start = genomic_mature_end
        for start in range(len(trna.sequence)):
            for length in range(minimum, min(maximum, len(trna.sequence) - start) + 1):
                end = start + length
                grouped[trna.sequence[start:end]].add(
                    FragmentOrigin(
                        trna.locus.identifier,
                        trna.locus.amino_acid,
                        trna.locus.anticodon,
                        trna.locus.chromosome,
                        f"{trna.locus.begin}-{trna.locus.end}",
                        trna.locus.strand,
                        start + (0 if has_minus_one else 1),
                        end - (0 if has_minus_one else 1),
                        trna.added_cca and end > cca_start,
                        trna.minus_one_base if has_minus_one and start == 0 else None,
                        trna.histidine_minus_one and start == 0,
                        trna.exon_junction is not None and start < trna.exon_junction < end,
                    )
                )
    return [
        Fragment(
            sequence,
            sorted(origins, key=lambda x: (x.source_trna_id, x.mature_start, x.mature_end)),
        )
        for sequence, origins in sorted(grouped.items())
    ]
