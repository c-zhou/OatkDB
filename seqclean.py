#!/usr/bin/env python

import argparse
import sys
import re
from collections import defaultdict
from Bio import SeqIO

def calculate_percentile(data, percentile):
    """
    Calculate the percentile of a list of values without using NumPy.

    Args:
    - data (list): The input data.
    - percentile (float): The desired percentile (between 0 and 100).

    Returns:
    - float: The computed percentile value.
    """
    if not 0 <= percentile <= 100:
        raise ValueError("Percentile must be between 0 and 100")

    sorted_data = sorted(data)
    index = (percentile / 100) * (len(sorted_data) - 1)

    if index.is_integer():
        return sorted_data[int(index)]
    else:
        lower_index = int(index)
        upper_index = lower_index + 1
        lower_value = sorted_data[lower_index]
        upper_value = sorted_data[upper_index]
        fraction = index - lower_index
        return lower_value + fraction * (upper_value - lower_value)

# count kmers FASTA records
def count_kmers_in_fasta(fasta_records, k):
    kmer_counts = defaultdict(int)
    for record in fasta_records:
        sequence = str(record.seq)
        for i in range(len(sequence) - k + 1):
            kmer_counts[sequence[i:i + k]] += 1

    return kmer_counts

def sequence_kmer_frequency(fasta_records, k, kmer_counts, len_median):
    kmer_freqs = []
    for record in fasta_records:
        sequence = str(record.seq)
        len_seq = len(sequence)
        n_kmer = len_seq - k + 1
        if n_kmer <= 0:
            kmer_freqs.append(0)
            continue;
        count = 0
        for i in range(n_kmer):
            kmer = sequence[i:i + k]
            if kmer in kmer_counts:
                count += kmer_counts[kmer]
        if len_seq > len_median:
            count *= len_median / len_seq
        kmer_freqs.append(count)

    return kmer_freqs

def kmer_stats(kmer_counts):
    total_num_kmers = sum(kmer_counts.values())
    unique_kmers = len(kmer_counts)
    sys.stderr.write(f'total number of kmers: {total_num_kmers}\n')
    sys.stderr.write(f'unique kmers:          {unique_kmers}\n')
    if unique_kmers:
        max_kmer = max(kmer_counts, key=kmer_counts.get)
        min_kmer = min(kmer_counts, key=kmer_counts.get)
        sys.stderr.write(f'largest kmer count:    [{max_kmer}:{kmer_counts[max_kmer]}]\n')
        sys.stderr.write(f'smallest kmer count:   [{min_kmer}:{kmer_counts[min_kmer]}]\n')

def sequence_filter_by_kmer_frequency(fasta_file, is_protein, kmer_size, rm_iqr_outlier, max_seq, output_file):
    # read sequences from the FASTA file
    if fasta_file and fasta_file != "-":
        fasta_handle = open(fasta_file, 'r')
    else:
        fasta_handle = sys.stdin
    fasta_records = [record for record in SeqIO.parse(fasta_handle, "fasta")]
    n_seq = len(fasta_records)
    lens = [len(record.seq) for record in fasta_records]
    q1 = calculate_percentile(lens, 25)
    q2 = calculate_percentile(lens, 50)
    q3 = calculate_percentile(lens, 75)
    sys.stderr.write(f'==> sequence length Q1: {q1}\n')
    sys.stderr.write(f'==> sequence length Q2: {q2}\n')
    sys.stderr.write(f'==> sequence length Q3: {q3}\n')
    # remove outliers by length IQR
    '''
    if rm_iqr_outlier:
        min_len = q1 - 1.5 * (q3 - q1)
        max_len = q3 + 1.5 * (q3 - q1)
        fasta_records[:] = (record for record in fasta_records if len(record) >= min_len and len(record) <= max_len)
        sys.stderr.write(f'==> number sequences removed by length IQR filtering: {n_seq - len(fasta_records)}\n')
    '''

    # make a kmer dictionary
    kmer_counts = count_kmers_in_fasta(fasta_records, kmer_size)
    sys.stderr.write(f'==> kmer dictionary stats\n')
    kmer_stats(kmer_counts)
    # filter kmers with ambiguous bases
    symbol_regex = re.compile(r'[ACGTUacgtu]+') if not is_protein else re.compile(r'[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]+')
    kmer_counts = {kmer: count for kmer, count in kmer_counts.items() if symbol_regex.match(kmer)}
    sys.stderr.write(f'==> kmer dictionary stats after filtering\n')
    kmer_stats(kmer_counts)
    # calculate average kmer frequency for each sequence
    kmer_freqs = sequence_kmer_frequency(fasta_records, kmer_size, kmer_counts, q2)
    # remove outliers by IQR
    min_freq = min(kmer_freqs)
    sys.stderr.write(f'==> kmer frequency minimum: {min_freq}\n')
    sys.stderr.write(f'==> kmer frequency maximum: {max(kmer_freqs)}\n')
    if rm_iqr_outlier:
        q1 = calculate_percentile(kmer_freqs, 25)
        q3 = calculate_percentile(kmer_freqs, 75)
        min_freq = q1 - 1.5 * (q3 - q1)
        sys.stderr.write(f'==> kmer frequency q1: {q1}\n')
        sys.stderr.write(f'==> kmer frequency q3: {q3}\n')
        sys.stderr.write(f'==> kmer frequency minimum applying IQR: {min_freq}\n')
    # get output sequence indices
    freq_idx = sorted([(i, kmer_freqs[i]) for i in range(len(kmer_freqs))], key=lambda x: x[1], reverse=True)
    max_idx_to = len(freq_idx)
    if max_seq:
        max_idx_to = min(max_seq, max_idx_to)
    outseq_idx = [freq_idx[i][0] for i in range(max_idx_to) if freq_idx[i][1] >= min_freq]
    # output sequences
    if output_file and output_file != "-":
        output_handle = open(output_file, 'w')
    else:
        output_handle = sys.stdout
    for i in range(len(outseq_idx)):
        idx = outseq_idx[i]
        output_handle.write(f">{fasta_records[idx].id} {freq_idx[i][1]} {len(fasta_records[idx].seq)}\n")
        output_handle.write(f"{str(fasta_records[idx].seq)}\n")
    sys.stderr.write(f'==> number sequences after filtering: {len(outseq_idx)}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter sequences based on k-mer similarity.")
    parser.add_argument("fasta_file", type=str, nargs="?", default=None, help="Input FASTA file (default: stdin)")
    parser.add_argument("output_file", type=str, nargs="?", default=None, help="Output file name (default: stdout)")
    parser.add_argument("--protein", action="store_true", help="The input is protein sequences")
    parser.add_argument("--kmer_size", type=int, default=0, help="K-mer size (default: 12 for NT and 6 for AA)")
    parser.add_argument("--rm_iqr_outlier", action="store_true", help="Remove outlier sequences with IQR method")
    parser.add_argument("--max_seq", type=int, default=None, help="Maximum number of sequences to keep (default: 0)")
    args = parser.parse_args()
    
    if not args.kmer_size:
        args.kmer_size = 6 if args.protein else 12

    sequence_filter_by_kmer_frequency(args.fasta_file, args.protein, args.kmer_size, args.rm_iqr_outlier, args.max_seq, args.output_file)

