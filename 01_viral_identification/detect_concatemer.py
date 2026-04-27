import subprocess
import oxli
from needletail import parse_fastx_file


def get_average_kmer_freq(seq, k=21):
    """Calculate average k-mer frequency."""
    kct = oxli.KmerCountTable(ksize=k)
    seq = seq.upper()
    kct.consume(seq)
    total_kmers = len(seq) - k + 1
    if total_kmers <= 0:
        return 0
    return total_kmers / len(kct)


def repeat_match(record, min_repeat=1000):
    """Run repeat-match to detect internal repeats."""
    p = subprocess.Popen(
        ["repeat-match", "-n", str(min_repeat), "/dev/stdin"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    stdout, stderr = p.communicate(str(record))

    if stdout.startswith("***"):
        return []

    if p.returncode:
        raise RuntimeError(stderr)

    lines = stdout.strip().split("\n")[2:]
    repeats = []

    for line in lines:
        start1, start2, length = line.split()
        strand = -1 if start2.endswith("r") else 1
        start1 = int(start1)
        start2 = int(start2.strip("r"))
        length = int(length)

        repeats.append((start1, start2, length, strand))

    return repeats


def longest_repeat_length(repeats):
    """Return the longest repeat length."""
    if not repeats:
        return 0
    return max(r[2] for r in repeats)


def is_concatemer(genome_len, longest_repeat, kmer_freq):
    """Determine concatemer status."""
    if kmer_freq >= 1.4:
        return True

    if longest_repeat > 0 and longest_repeat < genome_len * 0.9:
        return True

    return False


def main(input_fasta, output_table):

    with open(output_table, "w") as out:

        header = [
            "contig_id",
            "genome_length",
            "longest_repeat",
            "repeat_ratio",
            "kmer_freq",
            "concatemer",
        ]

        out.write("\t".join(header) + "\n")

        for record in parse_fastx_file(input_fasta):

            seq = str(record.seq)
            genome_len = len(seq)

            # kmer frequency
            kmer_freq = get_average_kmer_freq(seq)

            # repeat detection
            repeats = repeat_match(record)
            longest = longest_repeat_length(repeats)

            if genome_len > 0:
                repeat_ratio = longest / genome_len
            else:
                repeat_ratio = 0

            concatemer = is_concatemer(genome_len, longest, kmer_freq)

            row = [
                record.id,
                str(genome_len),
                str(longest),
                f"{repeat_ratio:.4f}",
                f"{kmer_freq:.4f}",
                str(concatemer),
            ]

            out.write("\t".join(row) + "\n")


if __name__ == "__main__":

    import sys

    if len(sys.argv) != 3:
        print("Usage: python detect_concatemer.py input.fasta output.tsv")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_table = sys.argv[2]

    main(input_fasta, output_table)
