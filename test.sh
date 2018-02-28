#!/usr/bin/env bash

python3 simulate.py sample.fasta 12 50 0.01
python3 assemble.py reads.txt 30 -d
./dot.sh