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


def set_map(dictionary, key, values):
    """
    Adds a list of values to a key in dictionary, creating a new dictionary
    entry if necessary
    :param dictionary:
    :param key:
    :param values: a set of values
    :return: nothing
    """
    if key not in dictionary:
        dictionary[key] = values
    else:
        dictionary[key] = dictionary[key].union(values)
    

def get_kmers(reads, k):
    """
    Identify all k-mers in a set of reads.
    :param reads: the reads
    :param k: the length of k-mers we want
    :return: the list of all k-mers across the reads
    """
    kmers = []
    kmer_read_map = {}
    for read in reads:
        for i in range(len(read)-k+1):
            kmer = read[i:i+k]
            kmers.append(kmer)
            set_map(kmer_read_map, kmer, {read})

    return kmers, kmer_read_map


def build_de_bruijn(kmers, kmer_read_map):
    """
    Build the de Bruijn graph from a list of k-mers.
    :param kmer_read_map: a dictionary storing which reads a kmer is in
    :param kmers: the k-mers to build the graph from
    :return: the tuple (nodes, edges, node_read_map), where node_read_map stores which
    reads a node is in
    """

    nodes = set()
    edges = {}
    node_read_map = {}
    
    for kmer in kmers:
        left = kmer[:-1]
        nodes.add(left)
        set_map(node_read_map, left, kmer_read_map[kmer])

        right = kmer[1:]
        nodes.add(right)
        set_map(node_read_map, right, kmer_read_map[kmer])

        if (left, right) not in edges:
            edges[left, right] = 1
        else:
            edges[left, right] += 1

    return nodes, edges, node_read_map


def delete_node(node, nodes, edges):
    """
    Delete a node from the graph.
    :param node: the node to delete
    :param nodes: the set of nodes
    :param edges: the dictionary of edges
    :return: nothing
    """
    nodes.discard(node)
    for x, y in set(edges.keys()):
        if x == node or y == node:
            edges.pop((x, y))


def remove_tips(nodes, edges, k):
    """
    Remove all tips from the de Bruijn graph. A tip is a node that is
    disconnected at one end and whose label is shorter than 2k.
    :param nodes: the set of nodes in the graph
    :param edges: the dictionary of edges
    :param k: the k value used to create the graph
    :return: nothing
    """
    tips = []

    in_neighbors, out_neighbors = get_neighbors(nodes, edges)

    for x in nodes:
        has_incoming = len(in_neighbors[x]) == 0
        has_outgoing = len(out_neighbors[x]) == 0

        if has_incoming != has_outgoing:
            if len(x) < 2 * k:
                tips.append(x)

    for tip in tips:
        delete_node(tip, nodes, edges)

    
def write_dot(nodes, edges, name):
    """
    Write a graph to a dot file.
    :param edges: a dictionary of the edges in the graph
    :param name: the name of the dot file to write
    :return: nothing
    """
    out = 'digraph mygraph {'
    for node in nodes:
        out += '"{}";\n'.format(node)
    for left, right in edges:
        out += '"{}"->"{}" [label="{}"];\n'.format(left, right,
                                                   edges[left, right])
    out += '}'
    
    with open('{}.dot'.format(name), 'w') as f:
        f.write(out)


def collapse(nodes, edges):
    """
    Identify all k-mers in a set of reads.
    :param nodes: a set of nodes obtained by build_de_bruijn
    :param edges: a dictionary of edges in the graph
    :param in_neighbors: dict
    :param out_neighbors: dict
    :param node_read_map: dict of nodes -> reads
    :return: the updated nodes and edges
    """

    in_neighbors, out_neighbors = get_neighbors(nodes, edges)
    chains = []
    
    for x, y in edges:
        xEdge = len(out_neighbors[x])
        yEdge = len(in_neighbors[y])
        # print('--> ', xEdge, yEdge)
        if (xEdge <= 1) and (yEdge <= 1):
            ends_with_x = None
            starts_with_y = None
            for i, chain in enumerate(chains):
                if chain[-1] == x:
                    ends_with_x = i
                elif chain[0] == y:
                    starts_with_y = i
            if ends_with_x is not None and starts_with_y is not None:
                chains[ends_with_x].extend(chains[starts_with_y])
                chains.pop(starts_with_y)
            elif ends_with_x is not None:
                chains[ends_with_x].append(y)
            elif starts_with_y is not None:
                chains[starts_with_y].insert(0, x)
            else:
                chains.append([x, y])

    for chain in chains:
        new = chain[0] + ''.join(node[-1] for node in chain[1:])
        nodes.add(new)
        in_neighbors[new] = set()
        out_neighbors[new] = set()

        # Add in neighbors to the chain node
        for neighbor in in_neighbors[chain[0]]:
            edges[neighbor, new] = edges[neighbor, chain[0]]
            in_neighbors[new].add(neighbor)
            out_neighbors[neighbor].add(new)

        # Add out neighbors to the chain node
        for neighbor in out_neighbors[chain[-1]]:
            edges[new, neighbor] = edges[chain[-1], neighbor]
            out_neighbors[new].add(neighbor)
            in_neighbors[neighbor].add(new)

    for chain in chains:
        for node in chain:
            delete_node(node, nodes, edges)


def trace(read, next_options, node_read_map, remaining_reads):
    """
    Move through the de Bruijn graph according to the current read and available
    neighbors.
    :param read: the current read
    :param next_options: a set of options for the next node to visit
    :param node_read_map: a map from nodes to the reads they are in
    :param remaining_reads: the reads that have not yet been assigned to contigs
    :return: the new state of the trace in a (read, node) tuple
    """

    if len(next_options) == 0:
        return None, None

    # If we can, move to a neighbor on the same read
    # Otherwise, move to a neighbor on an unseen read
    next_node = None
    backup_node = None
    for neighbor in next_options:
        if read in node_read_map[neighbor]:
            next_node = neighbor
        elif len(node_read_map[neighbor].intersection(remaining_reads)) > 0:
            backup_node = neighbor

    # If we aren't staying on the same read, this read is done
    if next_node is None:
        # print('Done with read')

        # If we have no neighbors to visit at all, we're stuck
        if backup_node is None:
            return None, None

        # Go to the back up node, move to new read
        node = backup_node
        read = node_read_map[node].intersection(remaining_reads).pop()
        remaining_reads.remove(read)
    else:
        node = next_node

    return read, node


def get_neighbors(nodes, edges):
    """
    Find all in and out neighbors of the nodes
    :param nodes: the set of nodes
    :param edges: the dict of edges
    :return: the tuple (in, out) of dicts of in and out neighbors
    """
    out_neighbors = {}
    in_neighbors = {}
    for node in nodes:
        out_neighbors[node] = set()
        in_neighbors[node] = set()

    for x, y in edges:
        out_neighbors[x].add(y)
        in_neighbors[y].add(x)

    return in_neighbors, out_neighbors


def assemble(reads, k, dot):
    """
    Assemble a set of reads into a set of contigs using a de Bruijn graph. Write
    the contigs to the file contigs.txt.
    :param reads: The set of reads to assemble
    :param k: the size of k-mer to be used when building the de Bruijn graph
    :return: nothing
    """

    nodes, edges, node_read_map = build_de_bruijn(*get_kmers(reads, k))

    in_neighbors, out_neighbors = get_neighbors(nodes, edges)

    # Map from reads to their nodes
    read_node_map = {}
    for read in reads:
        read_node_map[read] = []

    for node in nodes:
        for read in node_read_map[node]:
            read_node_map[read].append(node)

    # Build contigs from reads
    contigs = []
    remaining_reads = set(reads)
    visited = set()
    while len(remaining_reads) != 0:
        # Choose an arbitrary read to creating a contig around
        read = remaining_reads.pop()

        # Skip reads with no kmers (edge effect)
        if len(read_node_map[read]) == 0:
            continue

        # Skip reads we have covered
        if len(set(read_node_map[read]).difference(visited)) == 0:
            continue

        # Store start of read, put read in contig, mark read nodes as visited
        start = read_node_map[read][0]
        node = start
        contig = start[1:-1]

        # Trace forward until we reach a dead end, building up end of contig
        while node is not None:
            if node not in out_neighbors:
                print('WTF', node in nodes)
            visited.add(node)
            contig += node[-1]
            next_options = out_neighbors[node].difference(visited)
            read, node = trace(read, next_options, node_read_map, remaining_reads)

        # Back up to start and trace backward, building up beginning of contig
        node = start
        while node is not None:
            visited.add(node)
            contig = node[0] + contig
            next_options = in_neighbors[node].difference(visited)
            read, node = trace(read, next_options, node_read_map, remaining_reads)

        contigs.append(contig)

    if dot:
        write_dot(nodes, edges, 'before')

        collapse(nodes, edges)
        remove_tips(nodes, edges, k)
        collapse(nodes, edges)

        write_dot(nodes, edges, 'after')


    print('{}'.format(n50(contigs)))
    with open('contigs.txt', 'w') as f:
        f.write('\n'.join(contigs))


def main(args):
    dot = False
    if '-d' in args:
        args.remove('-d')
        dot = True

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
        
    assemble(reads, k, dot)
    
if __name__ == '__main__':
    main(sys.argv[1:])