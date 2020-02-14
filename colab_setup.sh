#!/bin/bash

# Install Miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -bfp /usr/local
rm Miniconda3-latest-Linux-x86_64.sh

# Install RDKit
conda install -qyc conda-forge rdkit

# Set environment variable
export PYTHONPATH="/usr/local/lib/python3.7/site-packages"

# Install DeepSMILES unofficial GitHub fork
pip install git+git://github.com/mateuszrezler/deepsmiles@master

# Download sample data set and reference script
wget http://www.dna.bio.keio.ac.jp/smiles/TOX21/\
NR-PPAR-gamma_wholetraining.smiles
wget http://www.dna.bio.keio.ac.jp/smiles/feature.py
