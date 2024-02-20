#!/bin/bash -l
#$ -pe mpi 520
#$ -l mem=4G
#$ -l h_rt=6:00:00
#$ -cwd
#$ -m be
#$ -M ankur.sinha@ucl.ac.uk

source ~/.venv/bin/activate
gerun python3 circuit.py 1234
