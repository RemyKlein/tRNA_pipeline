import os
import argparse
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

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
            trf_id = f"tRF_{idx:07d}_len{len(kmer)}"
            origins_str = ";".join(sorted(origins))
            f_out.write(f">{trf_id} {origins_str}\n{kmer}\n")
            tsv_out.write(f"{trf_id}\t{kmer}\t{len(kmer)}\t{origins_str}\n")

    print(f"{len(seq_to_origins)} k-mers generated.")
    return output_fasta, output_tsv

def genome_search_space(genome_file, outfile, num_autosomes):
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

def main():
    parser = argparse.ArgumentParser(
        description=(
            "tRNA Processing Pipeline\n\n"
            "This pipeline provides multiple steps to process tRNA sequences:\n"
            "1. Run tRNAscan-SE on a genome.\n"
            "2. Filter tRNAscan-SE output.\n"
            "3. Extract tRNA sequences from a genome and remove introns.\n"
            "4. Add CCA tails and 5' nucleotide extensions.\n"
            "5. Generate k-mers from tRNA sequences.\n"
            "6. Generate a genome search space for mapping.\n"
            "7. Generate exonic masks for each chromosome and strand.\n\n"
            "Example usage:\n"
            "  python tRF_pipeline.py run-scan genome.fa --outdir results\n"
            "  python tRF_pipeline.py filter input.txt --output filtered.txt\n"
            "  python tRF_pipeline.py extract filtered.txt genome.fa --output tRNA_spliced.fa\n"
            "  python tRF_pipeline.py modification-trna tRNA_spliced.fa --output tRNA_all_mature.fa\n"
            "  python tRF_pipeline.py generate-kmers tRNA_all_mature.fa --prefix trf_lookup --min 16 --max 50\n"
            "  python tRF_pipeline.py genome-search-space genome.fa --output genome_search_space.txt --num-autosomes 19\n"
            "  python tRF_pipeline.py exonic-mask genome_search_space.txt trnascan_filtered.txt tRNA_spliced.fa\n\n"
            "Note:\n"
            "  - The genome FASTA must be from Ensembl to ensure chromosome naming consistency.\n"
            "  - The genome file must be decompressed before running tRNAscan-SE (e.g., use 'gunzip your_file.fa.gz').\n"
            "  - The exonic mask step uses a fixed output folder: 'exonic_masks/'.\n"
            "  - The pipeline handles + and - strands and inserts CCA/5' extensions as needed."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    subparser = parser.add_subparsers(dest="cmd", required=True, title="Subcommands",
                                      description="Available pipeline steps:")

    parser_run = subparser.add_parser(
        "run-scan",
        help="Run tRNAscan-SE on a genome FASTA file. "
    )
    parser_run.add_argument(
        "genome", type=str, 
        help="Genome FASTA file to scan."
    )
    parser_run.set_defaults(func=lambda args: run_trnascan(
        genome_file=args.genome
    ))

    parser_filter = subparser.add_parser(
        "filter",
        help="Filter tRNAscan-SE output to remove sequences with undetermined anticodons (NNN) and non-canonical chromosomes/contigs."
    )
    parser_filter.add_argument(
        "input", type=str, 
        help="Path to the tRNAscan-SE output file."
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
        help=(
            "Extract tRNA sequences from the genome based on tRNAscan-SE coordinates.\n"
            "Automatically handles + and - strands and removes introns if present."
        )
    )
    parser_extract.add_argument(
        "input", type=str, help="Filtered tRNAscan-SE output file."
    )
    parser_extract.add_argument(
        "genome_file", type=str, help="FASTA file of the reference genome."
    )
    parser_extract.add_argument(
        "--output", type=str, default="tRNA_spliced.fa",
        help="Output FASTA file with spliced tRNA sequences. Default: tRNA_spliced.fa"
    )
    parser_extract.set_defaults(func=lambda args: extract_trna_sequences(
        trnascan_out_filtered=args.input, genome_file=args.genome_file, trna_spliced=args.output
    ))

    parser_mod = subparser.add_parser(
        "modification-trna",
        help=(
            "Add CCA tail to the 3' end of tRNAs that lack it.\n"
            "Additionally, create 4 extra copies for each sequence with an added 5' nucleotide (A,T,C,G)."
        )
    )
    parser_mod.add_argument(
        "input", type=str, help="FASTA file with tRNA sequences to modify."
    )
    parser_mod.add_argument(
        "--output", type=str, default="tRNA_all_mature.fa",
        help="Output FASTA file with CCA added and 5' nucleotide extensions. Default: tRNA_all_mature.fa"
    )
    parser_mod.set_defaults(func=lambda args: add_cca_and_n1(
        trna_spliced=args.input, trna_all_mature=args.output
    ))

    parser_kmers = subparser.add_parser(
        "generate-kmers",
        help=(
            "Generate all k-mers from tRNA sequences within a specified size range.\n"
            "Outputs both a FASTA file of k-mers and a TSV lookup table."
        )
    )
    parser_kmers.add_argument(
        "input", type=str, help="FASTA file of tRNA sequences."
    )
    parser_kmers.add_argument(
        "--prefix", type=str, default="trf_lookup",
        help="Prefix for output FASTA and TSV files. Default: trf_lookup"
    )
    parser_kmers.add_argument(
        "--min", type=int, default=16,
        help="Minimum k-mer length. Default: 16"
    )
    parser_kmers.add_argument(
        "--max", type=int, default=50,
        help="Maximum k-mer length. Default: 50"
    )
    parser_kmers.set_defaults(func=lambda args: generate_kmers(
        trna_all_mature=args.input, outfile=args.prefix, min_length=args.min, max_length=args.max
    ))

    parser_search_space = subparser.add_parser(
        "genome-search-space", 
        help=("Generate a search space for tRNA mapping by writing each chromosome twice:\n"
        "- Forward strand (5'->3')\n"
        "- Reverse complement (5'->3')\n"
        "Only the main chromosomes (1-19, X, Y, MT) are included."
        )
    )
    parser_search_space.add_argument(
        "input", type=str, help="FASTA file of the genome."
    )
    parser_search_space.add_argument(
        "--output", type=str, default="genome_search_space.txt",
        help="Path to save the search space file. Default: genome_search_space.txt"
    )
    parser_search_space.add_argument(
        "--num-autosomes", type=int, default=19,
        help="Number of autosomes for the species. X, Y, and MT chromosomes are added automatically. Default: 19 (mouse)"
    )
    parser_search_space.set_defaults(func=lambda args: genome_search_space(
        genome_file=args.input, outfile=args.output, num_autosomes=args.num_autosomes
    ))

    parser_mask = subparser.add_parser(
        "exonic-mask",
        help=(
            "Generate exonic mask files for each chromosome and strand.\n"
            "Requires genome search space, filtered tRNAscan output, and spliced tRNA sequences.\n"
            "Output files are always written to a folder named 'exonic_masks/'."
        )
    )
    parser_mask.add_argument(
        "genome-search-space", type=str, 
        help="Genome search space file."
    )
    parser_mask.add_argument(
        "trnascan-filtered", type=str, 
        help="Filtered tRNAscan-SE output file."
    )
    parser_mask.add_argument(
        "trna-spliced", type=str, 
        help="Spliced tRNA FASTA file."
    )
    parser_mask.set_defaults(func=lambda args: generate_exonic_mask(
        genome_search_space_file=getattr(args, "genome-search-space"),
        trna_scan_filtered_file=getattr(args, "trnascan-filtered"),
        trna_spliced_file=getattr(args, "trna-spliced"),
    ))

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
