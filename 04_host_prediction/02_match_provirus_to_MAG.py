"""
match_provirus_to_MAG.py

Match provirus sequences to their MAG hosts.

Matching logic:
  Virus sequence name -> vOTU id (text before the first space)
                      -> contig id (strip trailing virus-number suffix,
                                    then keep the portion starting at _k/_ctg)
  MAG contig name     -> contig id (keep the portion starting at _k/_ctg)
  Exact contig id match on both sides -> write vOTU id and MAG filename

Usage:
  python match_provirus_to_MAG.py \
      --virus nr_provirus_all_final.fna \
      --mag_dir dereplicated_genomes_rename \
      --output provirus_host_results.tsv
"""

import os
import re
import argparse
from pathlib import Path


# ─────────────────────────────────────────────
# Helper functions
# ─────────────────────────────────────────────

def extract_virus_contig_id(seq_name: str) -> tuple:
    """
    Input : raw FASTA header (without '>')
    Return: (vOTU_id, contig_id)

    Examples:
      'A10_k141_4398338_length_5619_cov_5.8412_1 830-5619/5619'
        -> vOTU_id   = 'A10_k141_4398338_length_5619_cov_5.8412_1'
        -> contig_id = 'k141_4398338_length_5619_cov_5.8412'

      'L_SY621_8-10_ctg1147780_length_12922_coverage_3_circular_no_1 1-11460/12922'
        -> vOTU_id   = 'L_SY621_8-10_ctg1147780_length_12922_coverage_3_circular_no_1'
        -> contig_id = 'ctg1147780_length_12922_coverage_3_circular_no'

      'L_SY631E_8-12_ctg564881_length_3078314_coverage_58_circular_yes_1_3'
        -> vOTU_id   = 'L_SY631E_8-12_ctg564881_length_3078314_coverage_58_circular_yes_1_3'
        -> contig_id = 'ctg564881_length_3078314_coverage_58_circular_yes'
    """
    # Step 1: drop everything after the first space to get the vOTU id
    votu_id = seq_name.split()[0]

    # Step 2: strip one or more trailing virus-number suffixes (e.g. _1, _1_3)
    contig_id_raw = re.sub(r'(_\d+)+$', '', votu_id)

    # Step 3: find the _k<digits> or _ctg boundary and keep the right-hand side
    m = re.search(r'_(k\d|ctg)', contig_id_raw)
    if m:
        contig_id = contig_id_raw[m.start() + 1:]  # +1 to skip the delimiter underscore
    else:
        # No delimiter found; use the full string
        contig_id = contig_id_raw

    return votu_id, contig_id


def extract_mag_contig_id(contig_name: str) -> str:
    """
    Input : contig name as it appears in a MAG FASTA file
    Return: contig_id (the portion starting at _k/_ctg, prefix stripped)

    Example:
      'EGM_E29_sbin_27_k141_1911868_length_21066_cov_6.0000'
        -> 'k141_1911868_length_21066_cov_6.0000'
    """
    m = re.search(r'_(k\d|ctg)', contig_name)
    if m:
        return contig_name[m.start() + 1:]
    return contig_name


# ─────────────────────────────────────────────
# Main pipeline
# ─────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description='Match proviruses to MAG hosts.')
    parser.add_argument('--virus', default='nr_provirus_all_final.fna',
                        help='Provirus FASTA file (default: nr_provirus_all_final.fna)')
    parser.add_argument('--mag_dir', default='dereplicated_genomes_rename',
                        help='Directory containing MAG FASTA files (default: dereplicated_genomes_rename)')
    parser.add_argument('--output', default='provirus_host_results.tsv',
                        help='Output file path (default: provirus_host_results.tsv)')
    return parser.parse_args()


def load_virus_sequences(fna_path: str) -> list:
    """Read virus FASTA; return [(votu_id, contig_id), ...] deduplicated by vOTU id."""
    records = []
    seen_votu = set()
    with open(fna_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                seq_name = line[1:]  # strip leading '>'
                votu_id, contig_id = extract_virus_contig_id(seq_name)
                if votu_id not in seen_votu:
                    seen_votu.add(votu_id)
                    records.append((votu_id, contig_id))
    return records


def build_mag_index(mag_dir: str) -> dict:
    """
    Walk the MAG directory and build a contig_id -> MAG name index.
    Returns: {contig_id: mag_name} or {contig_id: [mag_name1, mag_name2, ...]}
             (list form used when the same contig_id appears in multiple MAGs)
    """
    index = {}
    mag_dir_path = Path(mag_dir)
    fa_files = (list(mag_dir_path.glob('*.fa')) +
                list(mag_dir_path.glob('*.fasta')) +
                list(mag_dir_path.glob('*.fna')))

    if not fa_files:
        raise FileNotFoundError(f"No .fa / .fasta / .fna files found in '{mag_dir}'")

    for fa_file in fa_files:
        mag_name = fa_file.stem  # filename without extension, e.g. EGM_E29_sbin_27
        with open(fa_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    contig_name = line[1:].split()[0]  # first field only
                    contig_id = extract_mag_contig_id(contig_name)
                    if contig_id in index:
                        existing = index[contig_id]
                        if isinstance(existing, list):
                            existing.append(mag_name)
                        else:
                            index[contig_id] = [existing, mag_name]
                    else:
                        index[contig_id] = mag_name
    return index


def match_and_write(virus_records: list, mag_index: dict, output_path: str):
    """Match virus records against the MAG index and write results to file."""
    matched = 0
    unmatched = 0
    with open(output_path, 'w') as out:
        out.write('vOTU_id\thost_MAG\n')
        for votu_id, contig_id in virus_records:
            hit = mag_index.get(contig_id)
            if hit is not None:
                if isinstance(hit, list):
                    for h in hit:
                        out.write(f'{votu_id}\t{h}\n')
                else:
                    out.write(f'{votu_id}\t{hit}\n')
                matched += 1
            else:
                unmatched += 1

    print(f'[Done] Total virus sequences : {len(virus_records)}')
    print(f'       Matched               : {matched}')
    print(f'       Unmatched             : {unmatched}')
    print(f'       Results written to    : {output_path}')


def main():
    args = parse_args()

    print(f'[1/3] Loading virus sequences: {args.virus}')
    virus_records = load_virus_sequences(args.virus)
    print(f'      {len(virus_records)} unique vOTUs loaded')

    print(f'[2/3] Indexing MAG directory: {args.mag_dir}')
    mag_index = build_mag_index(args.mag_dir)
    print(f'      {len(mag_index)} unique contig IDs indexed')

    print(f'[3/3] Matching and writing results...')
    match_and_write(virus_records, mag_index, args.output)


if __name__ == '__main__':
    main()
