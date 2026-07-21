import argparse
import csv
import json
import logging
import shutil
import subprocess
from pathlib import Path

from .exclusivity import classify_hits, find_hits
from .fragments import generate_fragments
from .genome import build_masks, read_fasta
from .maturation import mintmap_mature_variants
from .quantification import quantify_fastq
from .trnascan import parse_trnascan


def _chromosomes(value: str | None) -> list[str] | None:
    return value.split(",") if value else None


def build_lookup(args: argparse.Namespace) -> None:
    genome = read_fasta(args.genome)
    loci = parse_trnascan(args.trnascan, _chromosomes(args.chromosomes))
    trnas = [
        variant
        for locus in loci
        for variant in mintmap_mature_variants(locus, genome[locus.chromosome])
    ]
    fragments = generate_fragments(trnas, args.minimum, args.maximum)
    masks = build_masks(genome, loci)
    hits = find_hits((item.sequence for item in fragments), genome, masks)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["tRF_id", "sequence", "length", "exclusivity", "origins_json"]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for fragment in fragments:
            synthetic = any(
                o.overlaps_added_cca or o.minus_one_base or o.crosses_spliced_intron
                for o in fragment.origins
            )
            writer.writerow(
                {
                    "tRF_id": fragment.identifier,
                    "sequence": fragment.sequence,
                    "length": len(fragment.sequence),
                    "exclusivity": classify_hits(
                        hits[fragment.sequence], has_synthetic_bases=synthetic
                    ).value,
                    "origins_json": json.dumps(
                        [o.__dict__ for o in fragment.origins],
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                }
            )


def quantify(args: argparse.Namespace) -> None:
    with args.lookup.open(encoding="utf-8") as handle:
        lookup = {row["sequence"] for row in csv.DictReader(handle, delimiter="\t")}
    rows = quantify_fastq(args.fastq, lookup)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sequence", "raw_count", "RPM_tRNAspace", "RPM_total"])
        for row in rows:
            writer.writerow(
                [row.sequence, row.raw_count, f"{row.rpm_trna_space:.6f}", f"{row.rpm_total:.6f}"]
            )


def run_scan(args: argparse.Namespace) -> None:
    executable = shutil.which("tRNAscan-SE")
    if executable is None:
        raise FileNotFoundError("tRNAscan-SE executable not found in PATH")
    subprocess.run([executable, "-o", str(args.output), str(args.genome)], check=True)


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description="Exact tRF discovery inspired by MINTmap 1.0")
    sub = root.add_subparsers(required=True)
    lookup = sub.add_parser("build-lookup", help="build and classify exact candidate tRFs")
    lookup.add_argument("genome", type=Path)
    lookup.add_argument("trnascan", type=Path)
    lookup.add_argument("--output", type=Path, default=Path("trf_lookup.tsv"))
    lookup.add_argument("--chromosomes", help="comma-separated accepted FASTA identifiers")
    lookup.add_argument("--min", dest="minimum", type=int, default=16)
    lookup.add_argument("--max", dest="maximum", type=int, default=50)
    lookup.set_defaults(func=build_lookup)
    count = sub.add_parser("quantify", help="count exact lookup sequences in FASTQ/FASTQ.gz")
    count.add_argument("lookup", type=Path)
    count.add_argument("fastq", type=Path)
    count.add_argument("--output", type=Path, default=Path("trf_counts.tsv"))
    count.set_defaults(func=quantify)
    scan = sub.add_parser("run-scan", help="run tRNAscan-SE and fail on errors")
    scan.add_argument("genome", type=Path)
    scan.add_argument("--output", type=Path, default=Path("trnascan_out.txt"))
    scan.set_defaults(func=run_scan)
    return root


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parser().parse_args(argv)
    args.func(args)
    return 0
