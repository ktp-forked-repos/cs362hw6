#!/usr/bin/env bash

dot -Tpdf before.dot -o before.pdf
open before.pdf
dot -Tpdf after.dot -o after.pdf
open after.pdf