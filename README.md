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

This pipeline reproduces **Steps 1–8** from the MINTmap paper and generates a final output file:
