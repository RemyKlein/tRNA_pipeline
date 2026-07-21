# MINTmap parity validation plan

Direct lookup-building parity has not been demonstrated. MINTmap 1.0 ships a
precomputed human GRCh37 lookup but not the scripts used to construct it from
tRNAscan-SE output. This repository therefore does not claim equivalence.

For a future parity run:

1. use the exact 640-entry MINTmap 1.0 reference set and GRCh37 sequence;
2. create candidates and export sequence, origins, and exclusivity as TSV;
3. compare sequences and Y/N lookup values to the official lookup;
4. run both profilers on `ExampleRun/exampleInput.fastq.gz`;
5. compare raw counts, the two documented RPM values, and source annotations;
6. record every difference in a machine-readable JSON report containing tool
   versions, input checksums, totals, and per-sequence differences.

Expected non-parity: this implementation follows the task's biologically
constrained rule of adding only a -1 G to histidyl tRNAs. The literal MINTmap
paper step 4 enumerates A/T/C/G extensions for every reference sequence.

