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
    

def assemble(reads, k):
    contigs = []
    # TODO: implement
    
    
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