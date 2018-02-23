"""
authors: Kiran Tomlinson and Chae Kim
date: February 16 2018

This file contains methods for simulating noisy reads from a genome.
"""

import sys
from random import randint, choice, random

USAGE = 'Usage: python3 simulate.py fasta coverage read_length error_rate'
NUCLEOTIDES = ['A', 'T', 'G', 'C']


def error(message):
    """
    Print an error message and exit the program.
    :param message: the message to print
    :return: nothing
    """
    print('Error: {}'.format(message))
    print(USAGE)
    sys.exit(1)


def parse_fasta(path):
    """
    Extract a sequence from a .fasta file.
    :param path: the path to the file
    :return: the sequence contained in the file
    """
    try:
        with open(path) as f:
            return ''.join(line.strip().upper() for line in f.readlines()
                           if not line.startswith('>'))
    except FileNotFoundError:
        error('no such file: {}'.format(path))


def simulate_reads(sequence, coverage, read_length, error_rate):
    """
    Simulate a noisy set of reads from a DNA sequence. Reads that go past the
    end of the sequence are cut short.
    :param sequence: the sequence to read from
    :param coverage: the average number of times each base appears in a read
    :param read_length: the (max) length of a read
    :param error_rate: the probability with which to change bases in a read
    :return: a list of the simulated reads
    """

    # When we replace a nucleotide x wih y, we want to make sure x != y
    replacement = {x: [y for y in NUCLEOTIDES if y != x] for x in NUCLEOTIDES}
    reads = []
    g = len(sequence)

    # Perform N random reads
    for i in range(coverage * g // read_length):
        start = randint(0, g - 1)
        read = sequence[start:start+read_length]
        noisy_read = ''

        # Introduce errors in the read
        for j in range(len(read)):
            noisy_read += read[j] if random() > error_rate \
                else choice(replacement[read[j]])

        reads.append(noisy_read)

    return reads


def main(args):
    """
    Take in command line args and create a simulated read file to specification.
    :param args: the command line args passed in
    :return: nothing
    Write to the text file
    """

    if len(args) != 4:
        error('wrong number of arguments')

    sequence = parse_fasta(args[0])

    try:
        coverage = int(args[1])
        read_length = int(args[2])
        error_rate = float(args[3])
    except ValueError as e:
        error('bad input: {}'.format(str(e).split()[-1]))

    with open('reads.txt', 'w') as f:
        f.write('\n'.join(simulate_reads(sequence, coverage, read_length,
                                         error_rate)))

if __name__ == '__main__':
    main(sys.argv[1:])
