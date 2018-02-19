import sys

USAGE = 'python3 assemble.py reads_file k'


def error(message):
    print('Error: {}'.format(message))
    print(USAGE)
    exit(1)


def n50(contigs):
    # TODO: implement
    return -1
    

def get_kmers(reads, k):
    kmers = []
    # TODO: implement
    return kmers


def build_de_bruijn(kmers):
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
            
    for l, r in edges:
        print('{} -> {} ({})'.format(nodes[l], nodes[r], edges[l, r]))
        
    write_dot(nodes, edges)
    
    
def write_dot(nodes, edges):
    out = 'digraph mygraph {'
    for l, r in edges:
        out += '"{}"->"{}"'.format(nodes[l], nodes[r])
    out += '}'
    
    with open('test.dot', 'w') as f:
        f.write(out)

        
def assemble(reads, k):
    contigs = []
    # TODO: implement
    
    build_de_bruijn(['ATG', 'GCG', 'TGG', 'GGC', 'CGT', 'GTG', 'TGC', 'GCA'])
    
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