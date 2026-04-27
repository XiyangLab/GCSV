#!/usr/bin/env python3
import sys
import re

if len(sys.argv) != 3:
    print("Usage: python3 1clean_h.py input_file output_tsv")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

def clean_line(line_bytes):

    line = line_bytes.replace(b'\x00', b'')             # remove NULL
    line = re.sub(rb'[^\x20-\x7E\t]+', b'', line)        # Only keep printable characters
    try:
        return line.decode('utf-8').strip()
    except UnicodeDecodeError:
        return None

with open(input_file, 'rb') as infile, open(output_file, 'w', encoding='utf-8') as outfile:
    outfile.write("target\tdescription\n")
    for line_bytes in infile:
        clean = clean_line(line_bytes)
        if clean:
            parts = clean.split(maxsplit=1)
            if len(parts) == 2:
                target, desc = parts
                outfile.write(f"{target}\t{desc}\n")
