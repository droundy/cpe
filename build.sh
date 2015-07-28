#!/bin/sh

set -ev

(python impedance.py)

(python3 transposed_resistivity.py)

(python poisson.py)

(pdflatex cpe.tex && pdflatex cpe.tex)

