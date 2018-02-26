#!/usr/bin/env bash

python3 simulate.py sample_small.fasta 20 7 0
python3 assemble.py reads.txt 5
dot -Tpdf before.dot -o before.pdf
open before.pdf
dot -Tpdf after.dot -o after.pdf
open after.pdf