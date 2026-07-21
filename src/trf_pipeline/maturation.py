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
    is_his = locus.amino_acid.lower() in {"his", "histidine"}
    if is_his:
        mature = "G" + mature
        if junction is not None:
            junction += 1
    return MatureTRNA(locus, mature, oriented, added_cca, is_his, junction)
