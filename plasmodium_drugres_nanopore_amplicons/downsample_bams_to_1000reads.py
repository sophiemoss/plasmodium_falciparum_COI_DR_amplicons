#!/usr/bin/env python3


##### Run using: #####

# python /mnt/storage11/sophie/gitrepos/smoss_ampseq/nanopore_ampseq_pipeline/plasmodium_capeverde/downsample_to_1000reads.py \
#   --highcov highcov_amplicons.tsv \
#   --bed /mnt/storage11/sophie/plasmodium_cape_verde/CV_PF_DR_amplicons3/pf_drug_resistance/pf3d7_drugresamplicons.bed \
#   --bam-dir bams \
#   --out-dir downsampled_bams \
#   --max-reads 1000 \
#   --seed 42
# 

########################

import argparse
import csv
import os
import random
import subprocess
import sys
from pathlib import Path
import pysam

def run(cmd, check=True):
    """Run a shell command and print it."""
    print("+", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return proc

def read_highcov(tsv_path):
    """Read highcov_amplicons.tsv -> dict[(sample, amplicon)] = coverage"""
    over = {}
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            amp = row["AMPLICON"]
            sample = row["Sample"]
            cov = float(row["Coverage"])
            over[(sample, amp)] = cov
    return over

def read_bed(bed_path):
    """Read BED (chrom, start, end, name) -> dict[name] = (chrom, start, end)"""
    d = {}
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                raise ValueError("BED must have 4 columns: chrom start end name")
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            d[name] = (chrom, start, end)
    return d

def ensure_bai(bam):
    """Ensure BAM has an index."""
    if not Path(bam + ".bai").exists() and not Path(bam[:-4] + ".csi").exists():
        run(["samtools", "index", bam])

def sample_reads_from_region(bam_path, region, n_reads, seed, out_bam):
    """Extract up to n_reads random reads from a region into a new BAM."""
    random.seed(seed)
    bam = pysam.AlignmentFile(bam_path, "rb")
    header = bam.header
    out = pysam.AlignmentFile(out_bam, "wb", header=header)

    reads = list(bam.fetch(region=region))
    total = len(reads)
    if total <= n_reads:
        for r in reads:
            out.write(r)
        out.close()
        bam.close()
        print(f"  {region}: kept all {total} reads (≤ {n_reads})")
        return total
    else:
        chosen = random.sample(reads, n_reads)
        for r in chosen:
            out.write(r)
        out.close()
        bam.close()
        print(f"  {region}: downsampled {total} → {n_reads} reads")
        return n_reads

def main():
    ap = argparse.ArgumentParser(description="Downsample amplicon regions to max 1000 reads using a fixed random seed.")
    ap.add_argument("--highcov", required=True, help="highcov_amplicons.tsv (AMPLICON,Sample,Coverage)")
    ap.add_argument("--bed", required=True, help="BED file with 4 columns: chrom start end AMPLICON")
    ap.add_argument("--bam-dir", required=True, help="Directory with per-sample BAMs named <Sample>.bam")
    ap.add_argument("--out-dir", required=True, help="Output directory for per-sample downsampled BAMs")
    ap.add_argument("--max-reads", type=int, default=1000, help="Maximum reads per amplicon (default 1000)")
    ap.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    args = ap.parse_args()

    over = read_highcov(args.highcov)
    bed = read_bed(args.bed)
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    tmp_root = Path(args.out_dir) / "_tmp_chunks"
    tmp_root.mkdir(parents=True, exist_ok=True)

    samples = sorted({p.stem for p in Path(args.bam_dir).glob("*.bam")})

    for sample in samples:
        bam_in = str(Path(args.bam_dir) / f"{sample}.bam")
        if not Path(bam_in).exists():
            print(f"[!] Missing BAM for {sample}, skipping.")
            continue

        ensure_bai(bam_in)
        sample_tmp = tmp_root / sample
        sample_tmp.mkdir(parents=True, exist_ok=True)
        chunk_bams = []

        for amplicon, (chrom, start, end) in bed.items():
            region = f"{chrom}:{start}-{end}"
            out_chunk = sample_tmp / f"{amplicon}.bam"

            cov = over.get((sample, amplicon), None)
            if cov is not None:
                # Amplicon is in highcov file → cap at max_reads
                sample_reads_from_region(bam_in, region, args.max_reads, args.seed, str(out_chunk))
            else:
                # Amplicon not over-covered → keep all reads
                sample_reads_from_region(bam_in, region, float("inf"), args.seed, str(out_chunk))

            chunk_bams.append(str(out_chunk))

        if not chunk_bams:
            print(f"[!] No highcov amplicons for {sample}, skipping merge.")
            continue

        merged_bam = Path(args.out_dir) / f"{sample}.downsampled.bam"
        if len(chunk_bams) == 1:
            run(["samtools", "sort", "-o", str(merged_bam), chunk_bams[0]])
        else:
            merged_tmp = sample_tmp / "merged.unsorted.bam"
            run(["samtools", "merge", "-f", str(merged_tmp)] + chunk_bams)
            run(["samtools", "sort", "-o", str(merged_bam), str(merged_tmp)])

        run(["samtools", "index", str(merged_bam)])
        print(f"[✔] Wrote {merged_bam}")

if __name__ == "__main__":
    main()
