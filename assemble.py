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
        for i in range(len(read)-k):
            kmer = read[i:i+k]
            kmers.append(kmer)

    return kmers


def build_de_bruijn(kmers):
    """
    Build the de Bruijn graph from a list of k-mers.
    :param kmers: the k-mers to build the graph from
    :return: the pair (nodes, edges)
    """

    nodes = []
    edges = {}
    
    for kmer in kmers:
        left = kmer[:-1]
        if left not in nodes:
            nodes.append(left)
        
        right = kmer[1:]
        if right not in nodes:
            nodes.append(right)
    
        left_index = nodes.index(left)
        right_index = nodes.index(right)

        if (left_index, right_index) not in edges:
            edges[left_index, right_index] = 1
        else:
            edges[left_index, right_index] += 1
            
    # for l, r in edges:
    #     print('{} -> {} ({})'.format(nodes[l], nodes[r], edges[l, r]))
        
    write_dot(nodes, edges)

    return nodes, edges
    
    
def write_dot(nodes, edges):
    """
    Write a graph to a dot file.
    :param nodes: a list of nodes in the graph
    :param edges: a dictionary of the edges in the graph
    :return: nothing
    """
    out = 'digraph mygraph {'
    for l, r in edges:
        out += '"{}"->"{}" [label="{}"];\n'.format(nodes[l], nodes[r],
                                                   edges[l, r])
    out += '}'
    
    with open('test.dot', 'w') as f:
        f.write(out)
        
def collapse(nodes, edges):
    colList = []
    newKey = []
    numNodes = len(nodes)
    #newList=[]
    new = ""
    
    for x,y in edges:
        xEdge = sum(1 for k in range numNodes if (x, k) in edges)
        yEdge = sum(1 for k in range numNodes if (k, y) in edges)
        if (xEdge == 1) and (yEdge == 1):
            colList.append((x,y))
            
    for pair in colList:
        new = pair[0]+pair[:-1]
        nodes.remove(pair[0])
        nodes.remove(pair[1])
        nodes.append(new)
        #newList.append(new)
    
    for x,y in edges:
        for k in range numNodes if (x, k) in edges:
            edges[new, k] = edges[x, k]
            del(edges[x,k])
        for k in range numNodes if (k, y) in edges:
            edges[k, new] = edges[k, y]
            del(edges[k,y])
            
    return (nodes, edges)



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
    
    print('N50 score: {}'.format(n50(contigs)))
    
    with open('contigs.txt', 'w') as f:
        f.write('\n'.join(contigs))


def main(args):
    if len(args) != 2:
        error('wrong number of arguments')
        
    try:
        with open(args[0]) as f:
            reads = f.readlines()
    except FileNotFoundError:
        error('no such file: {}'.format(args[0]))
        
    try:
        k = int(args[1])
    except ValueError:
        error('k must be an integer')
        
    assemble(reads, k)
    
if __name__ == '__main__':
    main(sys.argv[1:])