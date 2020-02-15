#!/bin/bash

# Install Miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -bfp /usr/local
rm Miniconda3-latest-Linux-x86_64.sh

# Install RDKit
conda install -qyc conda-forge rdkit

# Install DeepSMILES unofficial GitHub fork
pip install git+git://github.com/mateuszrezler/deepsmiles@master
