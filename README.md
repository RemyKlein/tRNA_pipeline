# tRF Pipeline (based on MINTmap)

**Author:** Adapted and implemented by *Rémy Klein*  
**Repository:** [https://github.com/RemyKlein/tRNA_pipeline](https://github.com/RemyKlein/tRNA_pipeline) 

---

**Reference:**  
> Loher, P., Telonis, A. & Rigoutsos, I.  
> MINTmap: fast and exhaustive profiling of nuclear and mitochondrial tRNA fragments from short RNA-seq data. Sci Rep 7, 41184 (2017).
> [https://doi.org/10.1038/srep41184](https://doi.org/10.1038/srep41184)

---

## Overview

This repository provides a **deterministic and reproducible implementation** of the MINTmap tRNA/tRF processing workflow.  
It builds a lookup table of **tRNA-derived fragments (tRFs)** and determines whether each fragment is **exclusive to the tRNA space** within the genome. It designed to perform fast, reproducible, and transparent profiling of **tRNA-derived fragments (tRFs)** from small RNA-seq data.

This pipeline processes a reference genome to:

1. Identify tRNA genes using **tRNAscan-SE**.  
2. Filter tRNAs to retain only canonical chromosomes and valid anticodons.
3. Extract and splice tRNA sequences from the genome.
4. Add 3' **CCA tails** and possible 5′ base extensions.  
5. Generate all possible **k-mers** (candidate tRFs).  
6. Build **genome search space** and generate **exonic masks**.  
7. Split large tRF lookup tables into smaller blocks for memory-safe processing.
8. Check **tRF exclusivity** against tRNA and non-tRNA genomic regions.
9. Count tRF occurrences in trimmed FASTQ reads.
10. Separate **bona fide tRFs** from ambiguous or non-exclusive ones.
11. Add **metadata** to annotated tRF tables (ID, origin, exclusivity).

The final output is a **TSV file** listing all tRFs with exclusivity status (`bona_fide`, `ambiguous`, `non_exclusive`) and associated counts and metadata.

---

## Requirements
- **Python:** ≥ 3.8  
- **tRNAscan-SE:** ≥ 2.0 (must be installed and available in `PATH`)  
- **Dependencies:** lightweight and installable via pip (`biopython`, `pandas`); no external alignment tools required  
- **Recommended system:** Linux or macOS

### Install Python dependencies

```bash
pip install -r requirements.txt
```

Content of `requirements.txt`:
```text
biopython>=1.80
pandas>=1.3
# (Optional: argparse and collections are built-in with Python ≥ 3.6)
```

---

## Genome Input
- Use **Ensembl reference genomes** for consistent chromosome naming (e.g., `1, 2, 3, ..., X, Y, MT`).
- FASTA files **must be decompressed** before running (`.fa` and not `.fa.gz`)

Example for mouse (GRCm39):
```bash
wget https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
```

---

## Installation
Clone the repository:
```bash
git clone https://github.com/remyklein/tRNA_pipeline.git
cd tRNA_pipeline
pip install -r requirements.txt
```

--- 

## Usage Overview
Each stage of the pipeline can be executed independently using subcommands.

| Step | Subcommand          | Description |
|------|-------------------|-------------|
| 1    | `run-scan`           | Identify tRNA genes using tRNAscan-SE |
| 2    | `filter`             | Filter tRNAs for canonical chromosomes and valid anticodons |
| 3    | `extract`            | Extract and splice tRNA sequences |
| 4    | `modification-trna`  | Add 3′ CCA tails and 5′ extensions |
| 5    | `generate-kmers`     | Generate k-mers from mature tRNAs |
| 6    | `genome-search-space`| Prepare genome forward/reverse sequences |
| 7    | `exonic-mask`        | Create exonic mask files for each strand |
| 8    | `split-tsv`          | Split large tRF lookup tables into smaller blocks |
| 9    | `check-exclusivity`  | Classify tRFs (bona fide, ambiguous, non-exclusive) |
| 10   | `trf-count-table`    | Count tRF abundances from FASTQ reads |
| 11   | `split-bona-fide`    | Separate bona fide tRFs |
| 12   | `add-metadata`       | Add tRF metadata (ID, origin, exclusivity) |

---

## Pipeline Steps

### 1. Run tRNAscan-SE
Identify all tRNA genes from the reference genome.
```bash
python tRF_pipeline.py run-scan Mus_musculus.GRCm39.dna.primary_assembly.fa
```
**Output:**
```text
trnascan_out.txt
```
This step detects all putative tRNA genes in the genome using tRNAscan-SE.

### 2. Filter tRNAscan-SE output
Remove entries from non-canonical chromosomes and with unknown anticodons.
```bash
python tRF_pipeline.py filter trnascan_out.txt --output trnascan_out_filtered.txt
```
**Output:**
```text
trnascan_out_filtered.txt
```
To retain only nuclear and mitochondrial tRNAs with valid anticodons (excluding NNN).

### 3. Extract and splice tRNA sequences
Retrieve genomic sequences and remove introns.
```bash
python tRF_pipeline.py extract trnascan_out_filtered.txt Mus_musculus.GRCm39.dna.primary_assembly.fa --output tRNA_spliced.fa
```
**Output:**
```tex
tRNA_spliced.fa
```
Generates clean, intron-free tRNA sequences for further analysis.

### 4. Add CCA tails and 5' extensions
Append "CCA" at the 3' end if missing, and generate 5' N-extended versions (A, T, C, G).
```bash
python tRF_pipeline.py modification-trna tRNA_spliced.fa --output tRNA_all_mature.fa
```
**Output:**
```text
tRNA_all_mature.fa
```
Simulates mature tRNAs and allows detection of fragments originating from the first nucleotide position.

### 5. Generate k-mers
Produce all possible fragments (k-mers) of length 16 - 50 nt.
```bash
python tRF_pipeline.py generate-kmers tRNA_all_mature.fa --prefix trf_lookup --min 16 --max 50
```
**Output:**
```text
trf_lookup_16_50.fa
trf_lookup_16_50.tsv
```
This defines the full search space of possible tRNA-derived fragments (tRF candidates).

### 6. Create genome search space
Concatenate each chromosome's forward and reverse complement strands.
```bash
python tRF_pipeline.py genome-search-space Mus_musculus.GRCm39.dna.primary_assembly.fa --output genome_search_space.txt
```
**Output:**
```text
genome_search_space.txt
```
Used to check whether tRF candidates map exclusively within annotated tRNA regions.

### 7. Generate exonic masks
Mark genomic coordinates corresponding to tRNAs (exons = 1, CCA = 2, other = 0).
```bash
python tRF_pipeline.py exonic-mask genome_search_space.txt trnascan_out_filtered.txt tRNA_spliced.fa
```
**Output directory:**
```text
exonic_masks/
  ├── chr1_plus.mask
  ├── chr1_minus.mask
  ├── ...
```
These binary masks define tRNA-coding regions and are later used for exclusivity checks.

### 8. Split lookup table
Split the tRF TSV table into smaller, memory-safe blocks.
```bash
python tRF_pipeline.py split-tsv trf_lookup_16_50.tsv --output-dir trf_lookup_blocks --lines-per-block 1000000
```
**Output:**
```text
trf_lookup_blocks/
  ├── block_0001.tsv
  ├── block_0002.tsv
  ├── ...
```
Prepares the tRF lookup for large-scale exclusivity testing.

### 9. Check exclusivity
Classify tRFs as **bona fide**, **ambiguous**, or **non-exclusive**.
```bash
python tRF_pipeline.py check-exclusivity Mus_musculus.GRCm39.dna.primary_assembly.fa \
  --block-dir trf_lookup_blocks --mask-dir exonic_masks --output exclusivity_results.tsv
```
**Output:**
```text
exclusivity_results.tsv
```
Ensures each tRF sequence maps uniquely within tRNA regions, following the MINTmap logic.

### 10. Generate tRF count tables
Count tRF occurrences in trimmed FASTQ reads.
```bash
python tRF_pipeline.py trf-count-table trf_lookup_16_50.fa sample1_trimmed.fastq sample2_trimmed.fastq --output tRF_abundance
```
**Output:**
```text
tRF_abundance/
  ├── sample1_tRF_counts.tsv
  ├── sample2_tRF_counts.tsv
```
Quantifies observed tRF fragments directly from sequencing data.

### 11. Split bona fide tRFs
Separate bona fide tRFs from ambiguous/non-exclusive ones.
```bash
python tRF_pipeline.py split-bona-fide exclusivity_results.tsv tRF_abundance --output-dir tRF_bona_fide
```
**Output:**
```text
tRF_bona_fide/
  ├── sample1_bona_fide.tsv
  ├── sample1_ambiguous_or_non_exclusive.tsv
```
Facilitates downstream analyses focused on unique, biologically relevant tRFs.

### 12. Add metadata
Annotate count tables with tRF ID, origins, and exclusivity.
```bash
python tRF_pipeline.py add-metadata exclusivity_results.tsv tRF_abundance trf_lookup_16_50.tsv --output-dir tRF_metadata
```
**Output:**
```text
tRF_metadata/
  ├── sample1_tRF_counts_metadata.tsv
  ├── sample2_tRF_counts_metadata.tsv
```
Generates final, annotated tRF quantification tables ready for statistical analysis.

## Output Structure 
```text
tRNA_pipeline/
├── trnascan_out.txt
├── trnascan_out_filtered.txt
├── tRNA_spliced.fa
├── tRNA_all_mature.fa
├── trf_lookup_16_50.fa
├── trf_lookup_16_50.tsv
├── genome_search_space.txt
├── exonic_masks/
├── trf_lookup_blocks/
├── exclusivity_results.tsv
├── tRF_abundance/
├── tRF_bona_fide/
└── tRF_metadata/
```

---

## Notes & Recommendations
- Use **Ensembl genomes** for consistent chromosome naming.
- Processing all tRFs (16 - 50 nt) may be memory-intensive. The `split-tsv` step avoids memory overflow.
- The exclusivity step (Step 9) is the most time-consuming; parallelization is recommended for large genomes.
- This pipeline was tested for _Mus musculus_ (GRCm39) but can be adapted to other species with minor changes.

---

## Citation
If you use this pipeline, please cite:
>Loher, P., Telonis, A. G., & Rigoutsos, I. (2017). _MINTmap: fast and exhaustive profiling of
>nuclear and mitochondrial tRNA fragments from short RNA-seq data._ Scientific Reports, 7, 41184.
>[https://doi.org/10.1038/srep41184](https://doi.org/10.1038/srep41184).

>Klein, R. (2025). _tRNA_pipeline: a Python implementation of the MINTmap concept for tRF discovery
>and quantification._ GitHub: [https://github.com/RemyKlein/tRNA_pipeline](https://github.com/RemyKlein/tRNA_pipeline) 