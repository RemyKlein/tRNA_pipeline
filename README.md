# tRF Pipeline

A transparent Python implementation of exact tRNA-derived fragment discovery
and quantification, scientifically inspired by MINTmap 1.0. It is **not an
output-equivalent port**: direct parity testing has not yet been completed.

## Scientific scope

Implemented from the published MINTmap 1.0 method:

- tRNAs are extracted in transcript orientation and introns are spliced;
- a missing 3′ CCA is added, while an encoded terminal CCA is retained once;
- exact candidate sequences of 16–50 nt are enumerated and deduplicated;
- both genomic strands are searched exhaustively for exact matches;
- raw counts and RPM relative to all reads and assigned tRNA-space reads are
  reported.

Following the literal MINTmap 1.0 lookup-building step 4, every mature tRNA is
represented without a −1 extension and with four additional variants prefixed
by A, T, C, and G. The biological discussion specifically associates −1
guanylation with tRNA^His, but the exhaustive four-base enumeration is the
published lookup-generation rule implemented here.

Mitochondrial records use the same configurable maturation rules. No chromosome
name is intrinsically nuclear or mitochondrial, and no chromosome is silently
dropped.

## Coordinate conventions

tRNAscan-SE `Begin` and `End` are interpreted as 1-based inclusive genomic
coordinates. A descending pair denotes the negative strand. Internally all
intervals are 0-based half-open. For a forward interval `[start, end)` on a
chromosome of length `L`, the reverse-complement interval is
`[L-end, L-start)`. Intron coordinates are also 1-based inclusive and are
converted once during maturation.

Every genome sequence, mask, hit, and origin retains explicit chromosome and
strand metadata; FASTA record order has no scientific meaning.

## Exclusivity decision tree

Each occurrence is evaluated separately against the strand-specific tRNA mask.

1. Hits in both tRNA and non-tRNA space → `ambiguous`.
2. Hits only in tRNA space → `bona_fide`.
3. Hits only outside tRNA space → `non_exclusive`.
4. No genomic hit → `not_found`.
5. No contiguous genomic hit for a candidate involving added CCA, histidyl −1
   G, or a splice junction → `ambiguous`, because absence is not evidence of
   exclusivity.

`bona_fide` means genome-exclusive under this reference and mask. It does not
establish biological function.

## Installation

Python 3.9 or newer is supported. Runtime code uses the standard library.

```bash
python -m pip install -e .
python -m pip install -e ".[dev]"  # tests and lint
```

`tRNAscan-SE` is an external dependency only for `run-scan`; absence or a
non-zero subprocess exit is fatal.

## Inputs and species configuration

Genome input is FASTA. Sequence identifiers are preserved verbatim. tRNAscan-SE
rows must contain the conventional columns:

`Sequence, tRNA #, Begin, End, Type, Anticodon, Intron Begin, Intron End, Score`.

Use `--chromosomes chr1,chr2,chrX,chrM` to accept an explicit subset. If omitted,
all chromosomes referenced by valid tRNAscan-SE rows are used; a referenced
chromosome missing from the FASTA is an error. This supports `chr` prefixes,
optional sex chromosomes, arbitrary contigs, and mitochondrial aliases without
mouse-specific assumptions.

FASTQ and `.gz` FASTQ are streamed. Bases are uppercased and `U` is normalized
to `T`; `N` is accepted in reads but will not match the A/C/G/T candidate
lookup. Other symbols, malformed records, and empty files are rejected.

## Commands

```bash
trf-pipeline run-scan genome.fa --output trnascan_out.txt

trf-pipeline build-lookup genome.fa trnascan_out.txt \
  --chromosomes chr1,chr2,chrX,chrY,chrM \
  --min 16 --max 50 --output trf_lookup.tsv

trf-pipeline quantify trf_lookup.tsv sample.fastq.gz \
  --output sample.trf_counts.tsv
```

The compatibility launcher also works after installation:

```bash
python tRF_pipeline.py --help
```

## Output schemas

Lookup TSV:

- `tRF_id`: deterministic SHA-256-derived sequence identifier;
- `sequence`, `length`, `exclusivity`;
- `origins_json`: all source tRNAs with amino acid, anticodon, chromosome,
  locus, strand, mature coordinates, added-CCA overlap, the −1 base when
  present, histidyl −1 inclusion,
  and splice-junction crossing.

Count TSV:

- `sequence`, `raw_count`;
- `RPM_tRNAspace`: raw count / all reads assigned to lookup sequences × 10⁶;
- `RPM_total`: raw count / all FASTQ reads × 10⁶.

## Validation and tests

Synthetic tests cover positive/negative extraction, both intron orientations,
CCA, the four published −1 variants, deterministic multi-origin fragments, masks, all exclusivity
states, special mature-only candidates, gzip FASTQ, both RPM calculations,
malformed inputs, and a repeated end-to-end workflow.

```bash
pytest
ruff check .
ruff format --check .
python tRF_pipeline.py --help
pytest tests/test_pipeline.py::test_end_to_end_repeated_output_is_identical
```

See [VALIDATION.md](VALIDATION.md) for the outstanding machine-readable parity
test plan. Until that comparison is run against the official GRCh37 assets,
this project must not be described as equivalent to MINTmap.

## Performance

Candidates are grouped by length and each genomic window is examined once for
that length, avoiding a complete genome scan for every candidate. This is
deterministic and straightforward to test, but a production mammalian lookup
may still benefit from a compiled multi-pattern index. Masks use one byte per
base per strand, which favors clarity over minimum memory use.

## Limitations and compatibility

- The old 12-stage CLI and unlabeled genome-search-space/mask intermediates were
  removed because they encoded unsafe ordering and coordinate assumptions.
- Structural categories and MINTplates identifiers are not yet generated.
- tRNAscan-SE format variants beyond the documented text table need explicit
  adapters.
- Color-space reads are unsupported.
- Direct MINTmap lookup/output parity remains untested.

The original MINTmap code and lookup are GPL-3.0. This repository retains its
MIT license because this is an independent implementation based on the
published method and does not copy MINTmap source. Copying or distributing
official MINTmap code or lookup assets may impose GPL-3.0 obligations.

## Primary references

- Loher P, Telonis AG, Rigoutsos I. *MINTmap: fast and exhaustive profiling of
  nuclear and mitochondrial tRNA fragments from short RNA-seq data.* Scientific
  Reports 7, 41184 (2017). https://doi.org/10.1038/srep41184
- Official MINTmap 1.0 source and documentation:
  https://github.com/TJU-CMC-Org/MINTmap/
