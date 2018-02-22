#!/usr/bin/env bash

python3 simulate.py sample_small.fasta 5 10 0.05
python3 assemble.py reads.txt 5
dot -Tpdf trimmed.dot -o trimmed.pdf
open trimmed.pdf
dot -Tpdf untrimmed.dot -o untrimmed.pdf
open untrimmed.pdf