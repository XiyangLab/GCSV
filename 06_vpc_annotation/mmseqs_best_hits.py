#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Extract best MMseqs2 hit per query and merge with envHOG annotations."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-e", "--envhog", required=True)
    parser.add_argument("-o", "--output", default="best1.tsv")
    parser.add_argument("--min_identity", type=float, default=0.0)
    parser.add_argument("--max_evalue", type=float, default=1e-5)
    parser.add_argument("--min_qcov", type=float, default=0.8)
    args = parser.parse_args()

    # ----------------- 1. read mmseqs result file -----------------
    colnames = [
        "query","target","evalue","pident","fident",
        "qstart","qend","qlen","tstart","tend","tlen",
        "alnlen","qcov","tcov"
    ]
    df = pd.read_csv(args.input, sep="\t", header=None, names=colnames)

    # ----------------- 2. filtering -----------------
    df = df[
        (df["pident"] >= args.min_identity) &
        (df["evalue"] <= args.max_evalue) &
        (df["qcov"] >= args.min_qcov)
    ]

    if df.empty:
        print("⚠ No hits passed filters.")
        return

    # ----------------- 3. keep best hit per query -----------------
    df_best = (
        df.sort_values(["query", "evalue", "pident"], ascending=[True, True, False])
          .drop_duplicates("query")
    )

    # ----------------- 4. read envhog and merge -----------------
    envhog = pd.read_csv(args.envhog, sep="\t")

    merged = df_best.merge(envhog, how="left",
                           left_on="target", right_on="enVhog")

    # ----------------- 5. save result -----------------
    merged.to_csv(args.output, sep="\t", index=False)
    print(f"✅ Done! Output: {args.output}")
    print(f"   retained queries: {merged['query'].nunique()}")

if __name__ == "__main__":
    main()

