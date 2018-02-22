#!/usr/bin/env bash

python3 simulate.py sample_small.fasta 2 10 0.05
python3 assemble.py reads.txt 5
dot -Tpdf test.dot -o out.pdf
open out.pdf