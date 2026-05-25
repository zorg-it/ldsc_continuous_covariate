#!/usr/bin/env python
import argparse
import gzip
import os
import numpy as np
import pandas as pd
from pybedtools import BedTool


def normalize_chr(x):
    return str(x).replace("chr", "")


def aggregate_values(values, strategy):
    values = [float(v) for v in values]
    if strategy == "mean":
        return float(np.mean(values))
    if strategy == "max":
        return float(np.max(values))
    if strategy == "sum":
        return float(np.sum(values))
    if strategy == "first":
        return float(values[0])
    raise ValueError("Unsupported overlap strategy: {}".format(strategy))


def make_continuous_annot(args):
    df_bim = pd.read_csv(
        args.bimfile,
        delim_whitespace=True,
        usecols=[0, 1, 2, 3],
        names=["CHR", "SNP", "CM", "BP"],
    )
    df_bim["CHR"] = df_bim["CHR"].map(normalize_chr)

    bed_rows = []
    with open(args.bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split()
            if len(fields) < 4:
                raise ValueError("Continuous BED must contain at least 4 columns: chr start end value")
            chrom, start, end, value = fields[:4]
            try:
                value = float(value)
            except ValueError:
                continue
            bed_rows.append([normalize_chr(chrom), int(start), int(end), str(value)])

    if not bed_rows:
        raise ValueError("No valid rows found in BED file: {}".format(args.bed_file))

    bedtool = BedTool(bed_rows).sort()
    bim_rows = [[row.CHR, int(row.BP) - 1, int(row.BP), row.SNP] for row in df_bim.itertuples(index=False)]
    bimbed = BedTool(bim_rows)
    intersect = bimbed.intersect(bedtool, loj=True)

    values_by_snp = {}
    for x in intersect:
        snp = x[3]
        b_value = x[7]
        if b_value == ".":
            continue
        values_by_snp.setdefault(snp, []).append(b_value)

    annot_values = []
    for snp in df_bim["SNP"]:
        vals = values_by_snp.get(snp)
        annot_values.append(aggregate_values(vals, args.overlap_strategy) if vals else 0.0)

    df_out = df_bim[["CHR", "SNP", "CM", "BP"]].copy()
    df_out[args.annot_name] = annot_values

    outdir = os.path.dirname(args.annot_file)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    if args.annot_file.endswith(".gz"):
        with gzip.open(args.annot_file, "wt") as out:
            df_out.to_csv(out, sep="\t", index=False)
    else:
        df_out.to_csv(args.annot_file, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="Create continuous LDSC .annot.gz from BED chr/start/end/value.")
    parser.add_argument("--bed-file", required=True)
    parser.add_argument("--bimfile", required=True)
    parser.add_argument("--annot-file", required=True)
    parser.add_argument("--annot-name", required=True)
    parser.add_argument("--overlap-strategy", default="mean", choices=["mean", "max", "sum", "first"])
    args = parser.parse_args()
    make_continuous_annot(args)


if __name__ == "__main__":
    main()
