from .genome import reverse_complement
from .models import MatureTRNA, TRNALocus


def extract_and_mature(locus: TRNALocus, chromosome_sequence: str) -> MatureTRNA:
    low, high = locus.genomic_interval
    if low < 1 or high > len(chromosome_sequence):
        raise ValueError(f"coordinates outside chromosome bounds for {locus.identifier}")
    genomic = chromosome_sequence[low - 1 : high]
    oriented = genomic if locus.strand == "+" else reverse_complement(genomic)
    junction = None
    if locus.intron_begin:
        ilow, ihigh = sorted((locus.intron_begin, locus.intron_end))
        if locus.strand == "+":
            left, right = ilow - low, ihigh - low + 1
        else:
            left, right = high - ihigh, high - ilow + 1
        if not 0 <= left < right <= len(oriented):
            raise ValueError(f"invalid intron for {locus.identifier}")
        oriented = oriented[:left] + oriented[right:]
        junction = left
    added_cca = not oriented.endswith("CCA")
    mature = oriented + ("CCA" if added_cca else "")
    return MatureTRNA(locus, mature, oriented, added_cca, None, junction)


def mintmap_mature_variants(locus: TRNALocus, chromosome_sequence: str) -> list[MatureTRNA]:
    """Return the mature tRNA and all four -1 variants from MINTmap step 4."""
    base = extract_and_mature(locus, chromosome_sequence)
    variants = [base]
    for nucleotide in "ATCG":
        variants.append(
            MatureTRNA(
                locus=base.locus,
                sequence=nucleotide + base.sequence,
                genomic_sequence=base.genomic_sequence,
                added_cca=base.added_cca,
                minus_one_base=nucleotide,
                exon_junction=(base.exon_junction + 1 if base.exon_junction is not None else None),
            )
        )
    return variants
