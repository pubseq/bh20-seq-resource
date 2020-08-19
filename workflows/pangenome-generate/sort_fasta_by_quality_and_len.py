#!/usr/bin/env python3

# Sort the sequences by quality (percentage of number of N bases not called, descending) and by length (descending).
# The best sequence is the longest one, with no uncalled bases.

import os
import sys
import gzip

#import xxhash # Faster library
import hashlib

def open_gzipsafe(path_file):
    if path_file.endswith('.gz'):
    	return gzip.open(path_file, 'rt')
    else:
        return open(path_file)

path_fasta = sys.argv[1]

hash_to_count_and_headers_dict = {}

header_to_seq_dict = {}
header_percCalledBases_seqLength_list = []

with open_gzipsafe(path_fasta) as f:
    for fasta in f.read().strip('\n>').split('>'):
        header = fasta.strip('\n').split('\n')[0]
        sequence = ''.join(fasta.strip('\n').split('\n')[1:])

        ##hash = xxhash.xxh64(sequence).hexdigest() # Faster library
        hash = hashlib.md5(sequence.encode('utf-8')).hexdigest()

        if hash not in hash_to_count_and_headers_dict:
            # New sequence
            hash_to_count_and_headers_dict[hash] = [0, []]

            header_to_seq_dict[header] = sequence

            seq_len = len(sequence)
            header_percCalledBases_seqLength_list.append([header, (seq_len - sequence.count('N'))/seq_len, seq_len])

        hash_to_count_and_headers_dict[hash][0] += 1
        hash_to_count_and_headers_dict[hash][1].append(header)


with open('dups.txt', 'w') as fw:
    for count, header_list in hash_to_count_and_headers_dict.values():
        fw.write('\t'.join([str(count), ', '.join(header_list)]) + '\n')

for header, percCalledBases, seqLength_list in sorted(header_percCalledBases_seqLength_list, key=lambda x: (x[-2], x[-1]), reverse = True):
    sys.stdout.write('>{}\n{}\n'.format(header, header_to_seq_dict[header]))
