# SMILearn
A feature engineering package for deep learning models, based on a SMILES 
representation of chemical compounds.  
Enables fast customization and testing of feature matrices.
## Getting started <small>*(no installation required)*</small>
Before installation one can try a demo of this package by opening `demo.ipynb` notebook in Google Colab.<br/><br/>
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mateuszrezler/smilearn/blob/master/demo.ipynb)

This notebook contains complete practical example of usage and basic documentation of the package (section 2.3).
## Installation
The easiest way to recreate environment on a local host is to install Anaconda or Miniconda and buid it from the attached `smilearn.yml` file
```shell
# clone this repository
$ git clone https://github.com/mateuszrezler/smilearn.git

# create environment
$ cd smilearn
$ conda env create -f smilearn.yml

# activate environment
$ conda activate smilearn

# run demo
$ jupyter notebook demo.ipynb
```
## Example
$\gamma$-Cyhalothrin and heatmap of its feature matrix generated using this package.<br/>
![Image gamma-cyhalothrin](img/gamma_cyhalothrin.png)
![Feature matrix](img/matrix.png)
