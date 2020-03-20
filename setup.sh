# clone SMILearn repository
echo 'Cloning SMILearn repository...'
git clone https://github.com/mateuszrezler/smilearn.git

# create environment
echo 'Creating conda environment...'
cd smilearn && conda env create -f smilearn.yml

# activate environment
echo 'Activating conda environment...'
conda activate smilearn

# run demo
echo 'Running demo...'
cd smilearn && jupyter notebook smilearn.ipynb
