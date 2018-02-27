#!/usr/bin/env bash

python3 simulate.py sample_small.fasta 10 10 0
python3 assemble.py reads.txt 8
dot -Tpdf before.dot -o before.pdf
open before.pdf
dot -Tpdf after.dot -o after.pdf
open after.pdf