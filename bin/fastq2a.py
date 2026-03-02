#!/usr/bin/env python3

import argparse
import sys
import os

VERSION = "0.1"


def mask_fastq_to_fasta(input_file, output_file, threshold):
    """
    Reads a FASTQ file and writes a FASTA file where bases with
    Q-scores <= threshold are replaced with 'N'.
    """

    total_bases = 0
    masked_bases = 0
    total_records = 0

    try:
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:

            while True:
                header = f_in.readline().strip()
                if not header:
                    break

                sequence = f_in.readline().strip()
                plus_line = f_in.readline().strip()
                quality_line = f_in.readline().strip()

                if not quality_line:
                    break

                total_records += 1

                # Convert Phred+33 ASCII to Q scores
                q_scores = [ord(char) - 33 for char in quality_line]

                corrected_seq_list = []

                for i in range(len(sequence)):
                    total_bases += 1

                    if q_scores[i] > threshold:
                        corrected_seq_list.append(sequence[i])
                    else:
                        corrected_seq_list.append('N')
                        masked_bases += 1

                corrected_seq = "".join(corrected_seq_list)

                fasta_header = ">" + header[1:]
                f_out.write(f"{fasta_header}\n{corrected_seq}\n")

        print(f"Masked FASTA written to {output_file}")
        print(f"Masked bases: {masked_bases:,}")

    except FileNotFoundError:
        print(f"Error: Input file not found: {input_file}", file=sys.stderr)
        sys.exit(1)

    except PermissionError:
        print("Error: Permission denied.", file=sys.stderr)
        sys.exit(1)

    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)


def main():

    parser = argparse.ArgumentParser(
        prog="fastq2a.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Convert FASTQ assembly file to FASTA with low-quality bases masked.

Bases with Phred quality scores <= threshold are replaced with 'N'.
Quality scores are assumed to use standard Phred+33 encoding.
""",
        epilog="""
Examples:
  fastq2a.py -i assembly.fastq
  fastq2a.py -i assembly.fastq -o assembly.masked.fasta -q 5

"""
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        metavar="FILE",
        help="Input FASTQ assembly file generated with medaka2 or dorado polish (required)"
    )

    parser.add_argument(
        "-o", "--output",
        metavar="FILE",
        help="Output FASTA file. Default: ./<input_basename>.LQB_masked.fasta"
    )

    parser.add_argument(
        "-q", "--qscore",
        type=int,
        default=10,
        metavar="INT",
        help="Phred quality threshold for masking (default: 10)"
    )

    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"fastq2a.py {VERSION}",
        help="Show program version and exit"
    )

    args = parser.parse_args()

    if not args.output:
        base_name = os.path.splitext(os.path.basename(args.input))[0]
        args.output = f"{base_name}.LQB_masked.fasta"

    mask_fastq_to_fasta(args.input, args.output, args.qscore)


if __name__ == "__main__":
    main()