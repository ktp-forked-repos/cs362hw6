"""
authors: Kiran Tomlinson and Chae Kim
date: February 20 2018

This file contains methods for assembling contigs from a set of reads.
"""



import sys

USAGE = 'python3 assemble.py reads_file k'


def error(message):
    """
    Print an error message and exit the program.
    :param message: the message to print
    :return: nothing
    """

    print('Error: {}'.format(message))
    print(USAGE)
    exit(1)


def n50(contigs):
    """
    Calculate the N50 score of a list of contigs.
    :param contigs: a list of contigs
    :return: the integer N50 score
    """

    if len(contigs) == 0:
        return 0

    contigs = sorted(contigs, key=len, reverse=True)
    length = sum(len(contig) for contig in contigs)

    fraction = 0
    index = -1
    while fraction < length / 2:
        index += 1
        fraction += len(contigs[index])

    return len(contigs[index])
    

def get_kmers(reads, k):
    """
    Identify all k-mers in a set of reads.
    :param reads: the reads
    :param k: the length of k-mers we want
    :return: the list of all k-mers across the reads
    """
    kmers = []
    for read in reads:
        for i in range(len(read)-k+1):
            kmer = read[i:i+k]
            kmers.append(kmer)

    return kmers


def build_de_bruijn(kmers):
    """
    Build the de Bruijn graph from a list of k-mers.
    :param kmers: the k-mers to build the graph from
    :return: the pair (nodes, edges)
    """

    nodes = set()
    edges = {}
    
    for kmer in kmers:
        left = kmer[:-1]
        nodes.add(left)

        right = kmer[1:]
        nodes.add(right)

        if (left, right) not in edges:
            edges[left, right] = 1
        else:
            edges[left, right] += 1
        
    write_dot(edges)

    return nodes, edges


def remove_tips(nodes, edges, k):
    for x in nodes:
        has_incoming = any((y, x) in edges for y in nodes)
        has_outgoing = any((x, y) in edges for y in nodes)

        if not has_incoming or not has_outgoing:
            if len(x) < 2 * k:
                in_sum = sum(edges[y, x] for y in nodes if (y, x) in edges)
                print('Trim tip with incoming sum: {}'.format(in_sum))

    
def write_dot(edges):
    """
    Write a graph to a dot file.
    :param nodes: a list of nodes in the graph
    :param edges: a dictionary of the edges in the graph
    :return: nothing
    """
    out = 'digraph mygraph {'
    for left, right in edges:
        out += '"{}"->"{}" [label="{}"];\n'.format(left, right,
                                                   edges[left, right])
    out += '}'
    
    with open('test.dot', 'w') as f:
        f.write(out)

        
def assemble(reads, k):
    """
    Assemble a set of reads into a set of contigs using a de Bruijn graph. Write
    the contigs to the file contigs.txt.
    :param reads: The set of reads to assemble
    :param k: the size of k-mer to be used when building the de Bruijn graph
    :return: nothing
    """

    contigs = []
    # TODO: implement
    
    nodes, edges = build_de_bruijn(get_kmers(reads, k))

    remove_tips(nodes, edges, k)
    
    print('N50 score: {}'.format(n50(contigs)))
    
    with open('contigs.txt', 'w') as f:
        f.write('\n'.join(contigs))


def main(args):
    if len(args) != 2:
        error('wrong number of arguments')
        
    try:
        with open(args[0]) as f:
            reads = f.read().splitlines()
    except FileNotFoundError:
        error('no such file: {}'.format(args[0]))
        
    try:
        k = int(args[1])
    except ValueError:
        error('k must be an integer')
        
    assemble(reads, k)
    
if __name__ == '__main__':
    main(sys.argv[1:])