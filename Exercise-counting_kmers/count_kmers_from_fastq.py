#!/usr/bin/env python3

import sys, math

from sequence_to_kmer_list import *
from fastq_file_to_sequence_list import *


def count_kmers(kmer_list):

    kmer_count_dict = dict()

    for kmers in kmer_list:
        for kmer in kmers:
            if kmer not in kmer_count_dict:
                kmer_count_dict[kmer] = 1
            else:
                kmer_count_dict[kmer] += 1

    return kmer_count_dict


def main():

    progname = sys.argv[0]

    usage = "\n\n\tusage: {} filename.fastq kmer_length num_top_kmers_show\n\n\n".format(
        progname
    )

    if len(sys.argv) < 4:
        sys.stderr.write(usage)
        sys.exit(1)

    # capture command-line arguments
    fastq_filename = sys.argv[1]
    kmer_length = int(sys.argv[2])
    num_top_kmers_show = int(sys.argv[3])

    seq_list = seq_list_from_fastq_file(fastq_filename)

    all_kmers = list()

    for sequence in seq_list:
        kmers = sequence_to_kmer_list(sequence, kmer_length)
        all_kmers.append(kmers)



    kmer_count_dict = count_kmers(all_kmers)

    unique_kmers = list(kmer_count_dict.keys())
    
    sorted_kmers = sorted(kmer_count_dict, key=kmer_count_dict.get, reverse=True)
    unique_kmers = sorted_kmers
    kmer_entropy = dict()

    for kmer in sorted_kmers:
        probA = kmer.count('A') / len(kmer)
        probG = kmer.count('G') / len(kmer)
        probC = kmer.count('C') / len(kmer)
        probT = kmer.count('T') / len(kmer)
        entropy = 0
        if probA > 0:
            entropy += probA*math.log2(probA)
        if probG > 0:
            entropy += probG*math.log2(probG)
        if probC > 0:
            entropy += probC*math.log2(probC)
        if probT > 0:
            entropy += probT*math.log2(probT)
        if entropy != 0:
            entropy *= -1
        
        kmer_entropy[kmer] = entropy

    ## printing the num top kmers to show
    top_kmers_show = unique_kmers[0:num_top_kmers_show]

    for kmer in top_kmers_show:
        print(f'{kmer}\t{kmer_count_dict[kmer]}\t{kmer_entropy[kmer]:.2f}')

    sys.exit(0)  # always good practice to indicate worked ok!


if __name__ == "__main__":
    main()
