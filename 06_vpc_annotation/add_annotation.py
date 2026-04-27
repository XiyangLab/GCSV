#!/usr/bin/env python3

import sys
import pandas as pd

if len(sys.argv) != 4:
    print("Usage: python3 2add_annotation.py db_h result.tsv output.tsv")
    sys.exit(1)

db_h_file = sys.argv[1]
result_file = sys.argv[2]
output_file = sys.argv[3]

annotations = {}
with open(db_h_file, 'r') as f:
    for line in f:
        parts = line.strip().split(maxsplit=1)
        if len(parts) == 2:
            target_id, desc = parts
            annotations[target_id] = desc

df = pd.read_csv(result_file, sep='\t')

df['annotation'] = df['target'].map(annotations)

df.to_csv(output_file, sep='\t', index=False)
