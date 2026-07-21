"""Coordinate conversions used throughout the pipeline.

Genomic input is 1-based inclusive. Python and reverse-complement intervals are
0-based half-open.
"""


def one_based_inclusive_to_half_open(start: int, end: int) -> tuple[int, int]:
    if start < 1 or end < start:
        raise ValueError(f"invalid 1-based interval: {start}-{end}")
    return start - 1, end


def reverse_complement_interval(start: int, end: int, sequence_length: int) -> tuple[int, int]:
    """Convert a forward 0-based half-open interval to the reverse complement."""
    if not 0 <= start <= end <= sequence_length:
        raise ValueError("interval is outside the sequence")
    return sequence_length - end, sequence_length - start


def oriented_interval(begin: int, end: int, sequence_length: int) -> tuple[int, int, str]:
    low, high = sorted((begin, end))
    start0, end0 = one_based_inclusive_to_half_open(low, high)
    if begin <= end:
        return start0, end0, "+"
    rc_start, rc_end = reverse_complement_interval(start0, end0, sequence_length)
    return rc_start, rc_end, "-"
