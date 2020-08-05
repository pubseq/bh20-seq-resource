#!/usr/bin/env python3

# Sort the sequences by quality (percentage of number of N bases not called, descending) and by length (descending).
# The best sequence is the longest one, with no uncalled bases.

import os
import sys
import gzip

def open_gzipsafe(path_file):
    if path_file.endswith('.gz'):
    	return gzip.open(path_file, 'rt')
    else:
        return open(path_file)

path_fasta = sys.argv[1]

header_to_seq_dict = {}
header_percCalledBases_seqLength_list = []

with open_gzipsafe(path_fasta) as f:
    for fasta in f.read().strip('\n>').split('>'):
        header = fasta.strip('\n').split('\n')[0]

        header_to_seq_dict[
            header
        ] = ''.join(fasta.strip('\n').split('\n')[1:])

        seq_len = len(header_to_seq_dict[header])
        header_percCalledBases_seqLength_list.append([
            header, header_to_seq_dict[header].count('N'), (seq_len - header_to_seq_dict[header].count('N'))/seq_len, seq_len
        ])

for header, x, percCalledBases, seqLength_list in sorted(header_percCalledBases_seqLength_list, key=lambda x: (x[-2], x[-1]), reverse = True):
    sys.stdout.write('>{}\n{}\n'.format(header, header_to_seq_dict[header]))
