import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

def filter_trnascan_output(infile, outfile):
    with open(infile, "r") as file, open(outfile, "w") as out_file:
        for line in file:
            if line.startswith("Sequence") or line.startswith("Name") or line.startswith("-"):
                continue
            parts = line.split()
            anti_codon = parts[5]

            if anti_codon != "NNN":
                out_file.write(line)
    
    print(f"File tRNAscan-SE filtered : {outfile}")
    return outfile

def extract_trna_sequences(infile, outfile, genome_file):
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    with open(infile, "r") as file, open(outfile, "w") as out_file:

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

    print(f"tRNA sequence extracted: {outfile}")
    return outfile

def add_cca_and_n1(infile, outfile):
    out_cca = "tRNA_CCA.fa"
    
    # Add CCA
    with open(infile, "r") as file, open(out_cca, "w") as outfile:
        for record in SeqIO.parse(file, "fasta"):
            sequence = str(record.seq)
            if not sequence.endswith("CCA"):
                sequence += "CCA"
            record.seq = Seq(sequence)
            SeqIO.write(record, outfile, "fasta")
    
    # Add N pos n -1
    base = ["A", "T", "C", "G"]
    records_out = []
    with open(out_cca, "r") as file, open(outfile, "w") as out_file:
        for record in SeqIO.parse(file, "fasta"):
            sequence = str(record.seq)
            records_out.append(record)
            for nt in base:
                new_seq = Seq(nt + sequence)
                new_id = f"{record.id}_5'{nt}"
                new_record = SeqRecord(new_seq, id=new_id, description=f"5' extended with {nt}")
                records_out.append(new_record)
        SeqIO.write(records_out, out_file, "fasta")
    
    print(f"tRNA file modified: {outfile}")
    return outfile

def generate_kmers(infile, outfile, min_length=16, max_length=50):
    seq_to_origins = defaultdict(set)
    output_fasta = f"{outfile}_{min_length}_{max_length}.fa"
    output_tsv = f"{outfile}_{min_length}_{max_length}.tsv"

    for record in SeqIO.parse(infile, "fasta"):
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

    print(f"OK {len(seq_to_origins)} k-mers generated.")
    return output_fasta, output_tsv

def main():
    parser = argparse.ArgumentParser(description="tRNA processing pipeline: filtering, extraction from the genome target and remove possible introns, CCA addition, and k-mer generation.")
    subparser = parser.add_subparsers(dest="cmd", required=True)

    parser_filter = subparser.add_parser("filter", help="Command to remove all tRNAs not determined by tRNAscan-SE.")
    parser_filter.add_argument("input", type=str, help="Location of files resulting from tRNAscan-SE.")
    parser_filter.add_argument("--output", type=str, default="trnascan_out_filtered.txt", help="Location of the file resulting from this filtering.")
    parser_filter.set_defaults(func=lambda args: filter_trnascan_output(args.input, args.output))

    parser_extract_tRNA = subparser.add_parser("extract", help="Retrieve all tRNA sequences based on the genomic positions given by tRNAscan-SE and remove possible introns, taking into account the + and - strands.")
    parser_extract_tRNA.add_argument("input", type=str, help="Filtered tRNAscan-SE output.")
    parser_extract_tRNA.add_argument("genome_file", type=str, help="FASTA genome file.")
    parser_extract_tRNA.add_argument("--output", type=str, default="tRNA_spliced.fa", help="Output file for spliced tRNAs.")
    parser_extract_tRNA.set_defaults(func=lambda args: extract_trna_sequences(args.input, args.output, args.genome_file))

    parser_tRNA_modification = subparser.add_parser("modification_trna", help="Add the CCA tail sequence to the 3' end for tRNAs that do not already have this sequence. For each sequence, create 4 additional copies to which a nucleotide (A, T, C, G) is added to the 5' end.")
    parser_tRNA_modification.add_argument("input", type=str, help="")
    parser_tRNA_modification.add_argument("--output", type=str, default="tRNA_all_mature.fa", help="")
    parser_tRNA_modification.set_defaults(func=lambda args: add_cca_and_n1(args.input, args.output))

    parser_kmers = subparser.add_parser("generate_kmers", help="")
    parser_kmers.add_argument("input", type=str, help="")
    parser_kmers.add_argument("--prefix", type=str, default="trf_lookup", help="")
    parser_kmers.add_argument("--min", type=int, default=16, help="")
    parser_kmers.add_argument("--max", type=int, default=50, help="")
    parser_kmers.set_defaults(func=lambda args: generate_kmers(
        infile=args.input, 
        outfile=args.prefix, 
        min_length=args.min, 
        max_length=args.max
    ))

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()


