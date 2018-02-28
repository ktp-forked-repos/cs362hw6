#!/usr/bin/env bash

python3 simulate.py sample.fasta 12 50 0.01
python3 assemble.py reads.txt 30
dot -Tpdf before.dot -o before.pdf
open before.pdf
dot -Tpdf after.dot -o after.pdf
open after.pdf