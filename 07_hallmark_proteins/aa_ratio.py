#!/usr/bin/env python3
# coding: utf-8

import sys

if len(sys.argv) < 2:
    print("Usage: python aa_ratio.py input.faa")
    sys.exit(1)

faa = sys.argv[1]

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

seq_id = None
seq = []

print("ID\t" + "\t".join(amino_acids))

with open(faa) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        # fasta header
        if line.startswith(">"):
            if seq_id is not None:
                seq_str = "".join(seq)
                length = len(seq_str)
                ratios = []
                for aa in amino_acids:
                    ratios.append(seq_str.count(aa) / length if length > 0 else 0)
                print(seq_id + "\t" + "\t".join(f"{r:.6f}" for r in ratios))
            
            seq_id = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)

if seq_id is not None:
    seq_str = "".join(seq)
    length = len(seq_str)
    ratios = []
    for aa in amino_acids:
        ratios.append(seq_str.count(aa) / length if length > 0 else 0)
    print(seq_id + "\t" + "\t".join(f"{r:.6f}" for r in ratios))

