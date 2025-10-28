import os
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter

def run_trnascan(genome_file):
    genome_file = os.path.abspath(genome_file)
    output_file = "trnascan_out.txt"

    cmd = [
        "tRNAscan-SE",
        "-o", output_file,
        genome_file
    ]
    
    print(f"Running tRNAscan-SE on {genome_file}...")
    try:
        subprocess.run(cmd, check=True)
        print(f"tRNAscan-SE finished successfully!\n")
        print(f"Outputs written to:\n {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error while running tRNAscan-SE: {e}")

def filter_trnascan_output(trnascan_out, trnascan_out_filtered):
    print("=== Step 2: Filtering tRNAscan-SE output ===")
    canonical_chroms = [str(i) for i in range(1, 20)] + ["X", "Y", "MT"]

    with open(trnascan_out, "r") as file, open(trnascan_out_filtered, "w") as out_file:
        for line in file:
            if line.startswith("Sequence") or line.startswith("Name") or line.startswith("-"):
                continue

            parts = line.split()
            chrom = parts[0]
            anti_codon = parts[5]

            if chrom in canonical_chroms and anti_codon != "NNN":
                out_file.write(line)
    
    print(f"File tRNAscan-SE filtered : {trnascan_out_filtered}")
    return trnascan_out_filtered

def extract_trna_sequences(trnascan_out_filtered, trna_spliced, genome_file):
    print("=== Step 3: Extracting and splicing tRNA sequences ===")
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    with open(trnascan_out_filtered, "r") as file, open(trna_spliced, "w") as out_file:

        for line in file:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.split()
            chrom = parts[0]
            start_genomic_pos = int(parts[2])
            end_genomic_pos = int(parts[3])
            start_intron = int(parts[6])
            end_intron = int(parts[7])

            seq_genomic = genome[chrom].seq[min(start_genomic_pos, end_genomic_pos)-1 : max(start_genomic_pos, end_genomic_pos)]

            # Positive Strand 
            if start_genomic_pos < end_genomic_pos:
                tRNA_sequence = seq_genomic
                if start_intron != 0 and end_intron != 0:
                    start_rel = start_intron - start_genomic_pos
                    end_rel = end_intron - start_genomic_pos
                    spliced_seq = tRNA_sequence[:start_rel] + tRNA_sequence[end_rel:]
                else:
                    spliced_seq = tRNA_sequence

            # Negative Strand
            else:
                tRNA_sequence = seq_genomic.reverse_complement()
                if start_intron != 0 and end_intron != 0:
                    gene_len = start_genomic_pos - end_genomic_pos + 1
                    start_rel = gene_len - (start_intron - end_genomic_pos)
                    end_rel   = gene_len - (end_intron - end_genomic_pos)
                    spliced_seq = tRNA_sequence[:start_rel] + tRNA_sequence[end_rel:]
                else:
                    spliced_seq = tRNA_sequence

            tRNA_name = f">tRNA_{parts[4]}_{parts[5]}_{chrom}:{start_genomic_pos}-{end_genomic_pos}"
            out_file.write(f"{tRNA_name}\n{spliced_seq}\n")

    print(f"tRNA sequence extracted and possible introns removed: {trna_spliced}")
    return trna_spliced

def add_cca_and_n1(trna_spliced, trna_all_mature):
    print("=== Step 4: Adding CCA tail and 5'-extensions ===")
    out_cca = "tRNA_CCA.fa"
    
    # Add CCA
    with open(trna_spliced, "r") as file, open(out_cca, "w") as outfile:
        for record in SeqIO.parse(file, "fasta"):
            sequence = str(record.seq)
            if not sequence.endswith("CCA"):
                sequence += "CCA"
            record.seq = Seq(sequence)
            SeqIO.write(record, outfile, "fasta")
    
    # Add N pos n -1
    base = ["A", "T", "C", "G"]
    records_out = []
    with open(out_cca, "r") as file, open(trna_all_mature, "w") as out_file:
        for record in SeqIO.parse(file, "fasta"):
            sequence = str(record.seq)
            records_out.append(record)
            for nt in base:
                new_seq = Seq(nt + sequence)
                new_id = f"{record.id}_5'{nt}"
                new_record = SeqRecord(new_seq, id=new_id, description=f"5' extended with {nt}")
                records_out.append(new_record)
        SeqIO.write(records_out, out_file, "fasta")
    
    print(f"tRNA file modified: {trna_all_mature}")
    return trna_all_mature

def generate_kmers(trna_all_mature, outfile, min_length=16, max_length=50):
    print(f"=== Step 5: Generating k-mers ({min_length}-{max_length}) ===")
    seq_to_origins = defaultdict(set)
    output_fasta = f"{outfile}_{min_length}_{max_length}.fa"
    output_tsv = f"{outfile}_{min_length}_{max_length}.tsv"

    for record in SeqIO.parse(trna_all_mature, "fasta"):
        sequence = str(record.seq)
        tRNA_id = record.id
        L = len(sequence)
        for k in range(min_length, max_length + 1):
            for i in range(0, L - k + 1):
                kmer = sequence[i:i+k]
                seq_to_origins[kmer].add(tRNA_id)

    with open(output_fasta, "w") as f_out, open(output_tsv, "w") as tsv_out:
        tsv_out.write("tRF_id\tsequence\tlength\torigins\n")
        for idx, (kmer, origins) in enumerate(seq_to_origins.items(), start=1):
            trf_id = f"tRF_{idx:07d}"
            origins_str = ";".join(sorted(origins))
            f_out.write(f">{trf_id} {origins_str}\n{kmer}\n")
            tsv_out.write(f"{trf_id}\t{kmer}\t{len(kmer)}\t{origins_str}\n")

    print(f"{len(seq_to_origins)} k-mers generated.")
    return output_fasta, output_tsv

def genome_search_space(genome_file, outfile, num_autosomes):
    print("=== Step 6: Preparing genome search space (forward + reverse complement) ===")
    main_chrom = [str(i) for i in range(1, num_autosomes + 1)] + ["X", "Y", "MT"]

    line_counts = 0
    with open(outfile, "w") as out_file:
        for record in SeqIO.parse(genome_file, "fasta"):
            if record.id in main_chrom:
                sequence = str(record.seq)
                sequence_reverse_complement = str(Seq(sequence).reverse_complement())
                out_file.write(f"{sequence}\n{sequence_reverse_complement}\n")
                line_counts += 2

    print(f"File generated: {outfile}")
    print(f"Expected lines: {len(main_chrom) * 2}, actual lines: {line_counts}")
    return outfile

def generate_exonic_mask(genome_search_space_file, trna_scan_filtered_file, trna_spliced_file):
    print("=== Step 7: Generating exonic masks for chromosomes and strands ===")
    chrom_order = ["1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "2", "3", "4", "5", "6", "7", "8", "9", "MT", "X", "Y"]
    mask = {}

    with open(genome_search_space_file, "r") as infile:
        for index, line in enumerate(infile, start=0):
            strand = "+" if index % 2 == 0 else "-"
            chrom = chrom_order[index // 2]
            if chrom not in mask:
                mask[chrom] = {}
            mask[chrom][strand] = [0] * len(line.strip())

    # Check the presence of CCA before modification
    tRNA_has_CCA = {}
    for record in SeqIO.parse(trna_spliced_file, "fasta"):
        header = record.id
        sequence = str(record.seq)
        chrom_part = header.split("_")[-1]
        chrom, coords = chrom_part.split(":")
        start, end = coords.split("-")
        key = (chrom, int(start), int(end))
        tRNA_has_CCA[key] = sequence.endswith("CCA")


    with open(trna_scan_filtered_file, "r") as infile:
        
        for line in infile:
            line = line.strip()
            parts = line.split()
            chrom, start_position, end_position = parts[0], int(parts[2]), int(parts[3])

            # Positive strand
            if start_position < end_position:
                strand = "+"
            # Negative strand
            else:
                strand = "-"
                start_position, end_position = end_position, start_position

            end_position = min(end_position, len(mask[chrom][strand]))
            start_position = max(1, start_position)

            mask[chrom][strand][start_position - 1:end_position] = [1] * (end_position - start_position + 1)

            key = (chrom, start_position, end_position)
            has_cca = tRNA_has_CCA.get(key, True)

            if not has_cca:
                for pos in range(end_position, end_position + 3):
                    if pos < len(mask[chrom][strand]):
                        mask[chrom][strand][pos] = 2

            minus1_pos = start_position - 2
            if minus1_pos >= 0:
                mask[chrom][strand][minus1_pos] = 2
        
    mask_dir = "exonic_masks"
    os.makedirs(mask_dir, exist_ok=True)

    for chrom in mask:
        for strand in mask[chrom]:
            filename = f"chr{chrom}_{'plus' if strand == '+' else 'minus'}.mask"
            outpath = os.path.join(mask_dir, filename)
            with open(outpath, "w") as f:
                f.write("".join(map(str, mask[chrom][strand])) + "\n")

    print(f"Exonic masks generated in folder: {mask_dir}")
    return mask

def trf_lookup_split(trf_lookup_16_50, lines_per_block, output_dir):
    print("=== Step 8: Splitting large tRF lookup table into smaller TSV blocks ===")
    os.makedirs(output_dir, exist_ok=True)
    block_index = 1
    line_count = 0
    block_file = None

    with open(trf_lookup_16_50, "r") as infile:
        header = next(infile)
        for line in infile:
            if line_count % lines_per_block == 0:
                if block_file:
                    block_file.close()
                block_filename = os.path.join(output_dir, f"block_{block_index:04d}.tsv")
                block_file = open(block_filename, "w")
                block_file.write(header)
                block_index += 1
            block_file.write(line)
            line_count += 1
        if block_file:
            block_file.close()
    print(f"Split complete: {line_count} lines divided into {block_index - 1} blocks in '{output_dir}'")
    return output_dir

def check_exclusivity(genome_file, num_autosomes, block_dir, mask_dir, output_file):
    print("=== Step 9: Checking tRF exclusivity (bona_fide, ambiguous, non_exclusive) ===")
    # Loading genome 
    genome = {}
    main_chrom = [str(i) for i in range(1, num_autosomes + 1)] + ["X", "Y", "MT"]

    print("Loading genome...")
    for record in SeqIO.parse(genome_file, "fasta"):
        if record.id in main_chrom:
            sequence = str(record.seq)
            sequence_reverse_complement = str(Seq(sequence).reverse_complement())
            genome[record.id] = {
                "+":sequence,
                "-":sequence_reverse_complement
            }

    # Load masks
    mask = {}
    print("Loading exonic masks...")
    for fname in os.listdir(mask_dir):
        if fname.endswith(".mask"):
            chrom = fname.split("_")[0].replace("chr","")
            strand = "+" if "plus" in fname else "-"
            with open(os.path.join(mask_dir, fname), "r") as f:
                mask_seq = [int(c) for c in f.read().strip()]
            if chrom not in mask:
                mask[chrom] = {}
            mask[chrom][strand] = mask_seq

    print("Checking tRF exclusivity...")
    with open(output_file, "w") as out_f:
        out_f.write("tRF_id\tsequence\tlength\torigins\texclusivity\n")

        for block_file in sorted(os.listdir(block_dir)):
            if not block_file.endswith(".tsv"):
                continue
            block_path = os.path.join(block_dir, block_file)
            with open(block_path) as f:
                header = next(f)
                for line in f:
                    trf_id, seq, length, origins = line.strip().split("\t")
                    length = int(length)
                    hit_values = set()  

                    for chrom in genome:
                        for strand in ["+", "-"]:
                            seq_genome = genome[chrom][strand]
                            start = 0
                            while True:
                                index = seq_genome.find(seq, start)
                                if index == -1:
                                    break
                                mask_slice = mask[chrom][strand][index:index+length]
                                hit_values.update(mask_slice)
                                start = index + 1

                    if hit_values <= {1,2}:
                        exclusivity = "bona_fide"
                    elif 0 in hit_values and (1 in hit_values or 2 in hit_values):
                        exclusivity = "ambiguous"
                    elif hit_values == {0}:
                        exclusivity = "non_exclusive"
                    else:
                        exclusivity = "ambiguous"

                    out_f.write(f"{trf_id}\t{seq}\t{length}\t{origins}\t{exclusivity}\n")

    print(f"Exclusivity check finished. Results in '{output_file}'")

def run_tRF_abundance_table(lookup, trimmed_fastq, outdir):
    print("=== Step 10: Generating tRF abundance table from FASTQ reads ===")
    os.makedirs(outdir, exist_ok=True)

    lookup_sequences = set()
    for record in SeqIO.parse(lookup, "fasta"):
        seq = str(record.seq)
        lookup_sequences.add(seq)

    print(f"Loaded {len(lookup_sequences):,} sequences from lookup table.")

    for sample in trimmed_fastq:
        sample_name = os.path.basename(sample).replace(".fastq", "")
        output_file = os.path.join(outdir, f"{sample_name}_tRF_counts.tsv")

        read_counts = Counter()
        total_reads = 0
        with open(sample, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                seq = str(record.seq)
                read_counts[seq] += 1
                total_reads += 1

            print(f"Counted {len(read_counts):,} unique sequences from {total_reads:,} reads in FASTQ.")

            # Keeping only reads that appears in lookup_sequence
            filtered_counts = {seq: counts for seq, counts in read_counts.items() if seq in lookup_sequences}

            print(f"{len(filtered_counts):,} reads matched the lookup table sequences (tRF candidates).")

            total_tRF_reads = sum(filtered_counts.values())

            rows = []
            for seq, count in filtered_counts.items():
                rpm_tRNAspace = (count / total_tRF_reads) * 1e6 if total_tRF_reads > 0 else 0
                rpm_total = (count / total_reads) * 1e6 if total_reads > 0 else 0
                rows.append([seq, count, rpm_tRNAspace, rpm_total])

            df = pd.DataFrame(rows, columns=["sequence", "raw_count", "RPM_tRNAspace", "RPM_total"])
            df = df.sort_values(by="raw_count", ascending=False)
            df.to_csv(output_file, sep="\t", index=False)
            print(f"tRF abundance table saved for : {sample_name}-{output_file}")

def split_bona_fide(exclusivity_file, tRF_count_dir, output_dir):
    print("=== Step 11: Splitting bona_fide tRFs from others ===")
    os.makedirs(output_dir, exist_ok=True)

    exclusivity_df = pd.read_csv(exclusivity_file, sep="\t", dtype=str)
    seq_to_status = dict(zip(exclusivity_df["sequence"], exclusivity_df["exclusivity"]))

    for count_file in os.listdir(tRF_count_dir):
        if not count_file.endswith(".tsv"):
            continue
        sample_path = os.path.join(tRF_count_dir, count_file)
        sample_name = os.path.splitext(count_file)[0]

        df = pd.read_csv(sample_path, sep="\t", dtype=str)
        df["exclusivity"] = df["sequence"].map(seq_to_status).fillna("not_in_lookup")

        bona_fide_df = df[df["exclusivity"] == "bona_fide"]
        other_df = df[df["exclusivity"] != "bona_fide"]

        bona_fide_file = os.path.join(output_dir, f"{sample_name}_bona_fide.tsv")
        other_file = os.path.join(output_dir, f"{sample_name}_ambiguous_or_non_exclusive.tsv")

        bona_fide_df.to_csv(bona_fide_file, sep="\t", index=False)
        other_df.to_csv(other_file, sep="\t", index=False)
        print(f"{sample_name}: {len(bona_fide_df)} bona_fide, {len(other_df)} ambiguous/non_exclusive")

def add_metadata(exclusivity_file, count_dir, lookup_tsv, output_dir):
    print("=== Step 12: Adding metadata to tRF count tables ===")
    os.makedirs(output_dir, exist_ok=True)

    exclusivity_df = pd.read_csv(exclusivity_file, sep="\t", dtype=str)
    seq_to_exclusivity = dict(zip(exclusivity_df["sequence"], exclusivity_df["exclusivity"]))

    lookup_df = pd.read_csv(lookup_tsv, sep="\t", dtype=str)
    seq_to_origin = dict(zip(lookup_df["sequence"], lookup_df["origins"]))
    seq_to_id = dict(zip(lookup_df["sequence"], lookup_df["tRF_id"]))

    for count_file in os.listdir(count_dir):
        if not count_file.endswith(".tsv"):
            continue
        sample_path = os.path.join(count_dir, count_file)
        df = pd.read_csv(sample_path, sep="\t", dtype=str)

        df["tRF_id"] = df["sequence"].map(seq_to_id)
        df["origins"] = df["sequence"].map(seq_to_origin)
        df["exclusivity"] = df["sequence"].map(seq_to_exclusivity).fillna("not_in_lookup")

        out_file = os.path.join(output_dir, os.path.basename(count_file).replace(".tsv", "_metadata.tsv"))
        df.to_csv(out_file, sep="\t", index=False)
        print(f"Metadata added for {count_file}: {out_file}")

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "tRNA Processing Pipeline\n"
            "------------------------\n"
            "This pipeline processes tRNA sequences, generates tRF lookup tables, and computes tRF abundances.\n\n"
            "Steps:\n"
            "1. run-scan           : Identify tRNA genes using tRNAscan-SE.\n"
            "2. filter             : Filter tRNAscan-SE output for canonical chromosomes and valid anticodons.\n"
            "3. extract            : Extract tRNA sequences and splice introns.\n"
            "4. modification-trna  : Add CCA tails and 5'-extended sequences.\n"
            "5. generate-kmers     : Create k-mers from mature tRNAs for tRF lookup.\n"
            "6. genome-search-space: Prepare genome sequences for mapping.\n"
            "7. exonic-mask        : Generate exonic masks for each chromosome and strand.\n"
            "8. split-tsv          : Split the large tRF lookup table into smaller blocks for exclusivity check.\n"
            "9. check-exclusivity  : Classify tRFs as bona_fide, ambiguous, or non_exclusive.\n"
            "10. trf-count-table    : Count tRF abundances from FASTQ files.\n"
            "11. split-bona-fide   : Separate bona_fide tRFs from others.\n"
            "12. add-metadata      : Add metadata to tRF count tables.\n\n"
            "Prerequisites:\n"
            "- Genome FASTA must be decompressed and use Ensembl chromosome names.\n"
            "- Biopython installed.\n"
            "- tRNAscan-SE installed and in PATH.\n\n"
            "Example usage:\n"
            "  python tRF_pipeline.py run-scan genome.fa --outdir results\n"
            "  python tRF_pipeline.py generate-kmers tRNA_all_mature.fa --prefix trf_lookup --min 16 --max 50\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    subparser = parser.add_subparsers(dest="cmd", required=True, title="Subcommands",
                                      description="Available pipeline steps:")

    parser_run = subparser.add_parser(
        "run-scan",
        help="Run tRNAscan-SE on a genome FASTA file to identify tRNA genes."
    )
    parser_run.add_argument(
        "genome", type=str,
        help="Path to the genome FASTA file to scan."
    )
    parser_run.set_defaults(func=lambda args: run_trnascan(genome_file=args.genome))

    parser_filter = subparser.add_parser(
        "filter",
        help="Filter tRNAscan-SE output: remove NNN anticodons and non-canonical chromosomes."
    )
    parser_filter.add_argument(
        "input", type=str, help="Input tRNAscan-SE output file."
    )
    parser_filter.add_argument(
        "--output", type=str, default="trnascan_out_filtered.txt",
        help="Path to save the filtered output. Default: trnascan_out_filtered.txt"
    )
    parser_filter.set_defaults(func=lambda args: filter_trnascan_output(
        trnascan_out=args.input, trnascan_out_filtered=args.output
    ))

    parser_extract = subparser.add_parser(
        "extract",
        help="Extract tRNA sequences from genome and splice introns if present."
    )
    parser_extract.add_argument(
        "input", type=str, help="Filtered tRNAscan-SE output file."
    )
    parser_extract.add_argument(
        "genome_file", type=str, help="Reference genome FASTA file."
    )
    parser_extract.add_argument(
        "--output", type=str, default="tRNA_spliced.fa",
        help="Output FASTA file with spliced tRNAs. Default: tRNA_spliced.fa"
    )
    parser_extract.set_defaults(func=lambda args: extract_trna_sequences(
        trnascan_out_filtered=args.input, genome_file=args.genome_file, trna_spliced=args.output
    ))

    parser_mod = subparser.add_parser(
        "modification-trna",
        help="Add CCA tail and generate 5'-extended sequences (A/T/C/G) for tRNAs."
    )
    parser_mod.add_argument(
        "input", type=str, help="FASTA file with spliced tRNAs."
    )
    parser_mod.add_argument(
        "--output", type=str, default="tRNA_all_mature.fa",
        help="Output FASTA file with modifications. Default: tRNA_all_mature.fa"
    )
    parser_mod.set_defaults(func=lambda args: add_cca_and_n1(
        trna_spliced=args.input, trna_all_mature=args.output
    ))

    parser_kmers = subparser.add_parser(
        "generate-kmers",
        help="Generate k-mers from tRNA sequences for tRF lookup (FASTA and TSV outputs)."
    )
    parser_kmers.add_argument("input", type=str, help="FASTA file of mature tRNAs.")
    parser_kmers.add_argument("--prefix", type=str, default="trf_lookup",
                              help="Prefix for output files. Default: trf_lookup")
    parser_kmers.add_argument("--min", type=int, default=16, help="Minimum k-mer length. Default: 16")
    parser_kmers.add_argument("--max", type=int, default=50, help="Maximum k-mer length. Default: 50")
    parser_kmers.set_defaults(func=lambda args: generate_kmers(
        trna_all_mature=args.input, outfile=args.prefix, min_length=args.min, max_length=args.max
    ))

    parser_search_space = subparser.add_parser(
        "genome-search-space",
        help="Prepare genome sequences for tRF mapping (forward and reverse complement)."
    )
    parser_search_space.add_argument("input", type=str, help="Genome FASTA file.")
    parser_search_space.add_argument("--output", type=str, default="genome_search_space.txt",
                                     help="Output search space file. Default: genome_search_space.txt")
    parser_search_space.add_argument("--num-autosomes", type=int, default=19,
                                     help="Number of autosomes (X,Y,MT added automatically). Default: 19")
    parser_search_space.set_defaults(func=lambda args: genome_search_space(
        genome_file=args.input, outfile=args.output, num_autosomes=args.num_autosomes
    ))

    parser_mask = subparser.add_parser(
        "exonic-mask",
        help="Generate exonic masks for chromosomes/strands based on tRNAscan-SE results."
    )
    parser_mask.add_argument("genome-search-space", type=str, help="Genome search space file.")
    parser_mask.add_argument("trnascan-filtered", type=str, help="Filtered tRNAscan-SE output file.")
    parser_mask.add_argument("trna-spliced", type=str, help="FASTA file with spliced tRNA sequences.")
    parser_mask.set_defaults(func=lambda args: generate_exonic_mask(
        genome_search_space_file=getattr(args, "genome-search-space"),
        trna_scan_filtered_file=getattr(args, "trnascan-filtered"),
        trna_spliced_file=getattr(args, "trna-spliced"),
    ))

    parser_split = subparser.add_parser(
        "split-tsv",
        help="Split a large tRF TSV file into smaller blocks for memory-efficient processing."
    )
    parser_split.add_argument(
        "input", type=str,
        help="Path to the large tRF TSV file."
    )
    parser_split.add_argument(
        "--output-dir", type=str, default="trf_lookup_blocks",
        help="Directory to store TSV blocks. Default: trf_lookup_blocks"
    )
    parser_split.add_argument(
        "--lines-per-block", type=int, default=1000000,
        help="Number of lines per block. Default: 1000000"
    )
    parser_split.set_defaults(func=lambda args: trf_lookup_split(
        trf_lookup_16_50=args.input,
        output_dir=args.output_dir,
        lines_per_block=args.lines_per_block
    ))

    parser_exclusivity = subparser.add_parser(
        "check-exclusivity",
        help="Classify tRFs as bona_fide, ambiguous, or non_exclusive using genome and masks."
    )
    parser_exclusivity.add_argument("genome", type=str, help="Genome FASTA file.")
    parser_exclusivity.add_argument("--num-autosomes", type=int, default=19,
                                   help="Number of autosomes (X,Y,MT added automatically). Default: 19")
    parser_exclusivity.add_argument("--block-dir", type=str, default="trf_lookup_blocks",
                                   help="Directory with tRF TSV blocks. Default: trf_lookup_blocks")
    parser_exclusivity.add_argument("--mask-dir", type=str, default="exonic_masks",
                                   help="Directory with exonic mask files. Default: exonic_masks")
    parser_exclusivity.add_argument("--output", type=str, default="exclusivity_results.tsv",
                                   help="Output TSV file. Default: exclusivity_results.tsv")
    parser_exclusivity.set_defaults(func=lambda args: check_exclusivity(
        genome_file=args.genome, num_autosomes=args.num_autosomes,
        block_dir=args.block_dir, mask_dir=args.mask_dir, output_file=args.output
    ))

    parser_tRF_count_table = subparser.add_parser(
        "trf-count-table",
        help="Generate tRF abundance tables from FASTQ reads using a tRF lookup table."
    )
    parser_tRF_count_table.add_argument("lookup", help="FASTA lookup table of tRFs.")
    parser_tRF_count_table.add_argument("file-reads", nargs="+", help="One or more trimmed FASTQ files.")
    parser_tRF_count_table.add_argument("--output", default="tRF_abundance", help="Output directory for TSV files.")
    parser_tRF_count_table.set_defaults(func=lambda args: run_tRF_abundance_table(
        lookup=args.lookup, trimmed_fastq=getattr(args, "file-reads"), outdir=args.output
    ))

    parser_split_bona = subparser.add_parser(
        "split-bona-fide",
        help="Split tRF count tables into bona_fide and others based on exclusivity."
    )
    parser_split_bona.add_argument("exclusivity-file", help="Exclusivity results TSV file.")
    parser_split_bona.add_argument("count-dir", help="Directory with tRF count tables.")
    parser_split_bona.add_argument("--output-dir", default="tRF_bona_fide", help="Output folder for split tables.")
    parser_split_bona.set_defaults(func=lambda args: split_bona_fide(
        exclusivity_file=getattr(args, "exclusivity-file"),
        tRF_count_dir=getattr(args, "count-dir"),
        output_dir=args.output_dir
    ))

    parser_metadata = subparser.add_parser(
        "add-metadata",
        help="Add metadata (tRF ID, origins, exclusivity) to tRF count tables."
    )
    parser_metadata.add_argument("exclusivity-file", help="Exclusivity results TSV file.")
    parser_metadata.add_argument("count-dir", help="Directory with tRF count tables.")
    parser_metadata.add_argument("lookup-tsv", help="tRF lookup TSV table.")
    parser_metadata.add_argument("--output-dir", default="tRF_metadata", help="Directory to save metadata tables.")
    parser_metadata.set_defaults(func=lambda args: add_metadata(
        exclusivity_file=getattr(args, "exclusivity-file"),
        count_dir=getattr(args, "count-dir"),
        lookup_tsv=getattr(args, "lookup-tsv"),
        output_dir=args.output_dir
    ))

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()