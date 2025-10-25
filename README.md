# tRF Pipeline (based on MINTmap)

**Author:** Adapted and implemented by *Rémy Klein*  
**Reference:**  
> Loher, P., Telonis, A. & Rigoutsos, I.  
> MINTmap: fast and exhaustive profiling of nuclear and mitochondrial tRNA fragments from short RNA-seq data. Sci Rep 7, 41184 (2017).
> [https://doi.org/10.1038/srep41184](https://doi.org/10.1038/srep41184)

---

## Overview

This repository provides a **deterministic and reproducible implementation** of the MINTmap tRNA/tRF processing workflow.  
It builds a lookup table of **tRNA-derived fragments (tRFs)** and determines whether each fragment is **exclusive to the tRNA space** within the genome.

This pipeline processes a reference genome to:
1. Identify tRNA genes using **tRNAscan-SE**.  
2. Extract and splice tRNA sequences from the genome.  
3. Add CCA tails and possible 5′ base extensions.  
4. Generate all possible **k-mers** (candidate tRFs).  
5. Build genome search space and exonic masks.  
6. Check **tRF exclusivity** against tRNA and non-tRNA genomic regions.

The final output is a **TSV file** listing all tRFs with exclusivity status (`bona_fide`, `ambiguous`, `non_exclusive`).