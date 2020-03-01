from pandas import read_csv
from src.pipeline import *
from sklearn.pipeline import Pipeline

data = read_csv('datasets/tox21.csv')
af = [lambda mol, index: mol.GetAtomWithIdx(index).GetSymbol()] \
    + [lambda mol, index: int(mol.GetAtomWithIdx(index).GetSymbol() == symbol)
       for symbol in 'HCON'] \
    + [lambda mol, index: int(mol.GetAtomWithIdx(index).GetSymbol()
       not in ('H', 'C', 'O', 'N')),
       lambda mol, index: int(mol.GetAtomWithIdx(index).GetTotalNumHs()),
       lambda mol, index: int(mol.GetAtomWithIdx(index).GetTotalDegree()),
       lambda mol, index: int(mol.GetAtomWithIdx(index).GetFormalCharge()),
       lambda mol, index: int(mol.GetAtomWithIdx(index).GetTotalValence()),
       lambda mol, index: int(mol.GetAtomWithIdx(index).IsInRing()),
       lambda mol, index: int(mol.GetAtomWithIdx(index).GetIsAromatic())] \
    + [lambda mol, index:
       int(str(mol.GetAtomWithIdx(index).GetChiralTag()) == tag)
       for tag in ['CHI_TETRAHEDRAL_CW', 'CHI_TETRAHEDRAL_CCW', 'CHI_OTHER']] \
    + [lambda mol, index:
       int(str(mol.GetAtomWithIdx(index).GetHybridization()) == hyb)
       for hyb in ['S', 'SP', 'SP2', 'SP3', 'SP3D', 'SP3D2', 'OTHER']]
sf = [lambda token: int(token == char)
      for char in '()[].:=#\\/@+-234567<>']

pipeline = Pipeline(
    [('select_smiles_column',
        ColumnSelector(['smiles'])),
     ('rebuild_smiles',
        SmilesRebuilder(columns=['smiles'],
                        allBondsExplicit=True,
                        allHsExplicit=True,
                        canonical=True,
                        isomericSmiles=True)),
     ('translate_to_deepsmiles',
        DeepSmilesEncoder(columns=['smiles'],
                          save_as=['deep_smiles'],
                          branches=True,
                          rings=True)),
     ('tokenize_deepsmiles',
        RegexTokenizer(columns=['deep_smiles'],
                       regex=r'\[.*?\]|%\d{2}|\)+|[\d\(\)\-/\\:=#\$\.]')),
     ('rename_column',
        ColumnRenamer({'deep_smiles': 'tokenized'})),
     ('filter_length',
        LenFilter(columns=['tokenized'],
                  min_len=3,
                  max_len=3,
                  bounded=(True, True))),
     ('drop_nans',
        NanDropper()),
     ('zip_smiles_and_tokenized',
        Zipper(['smiles', 'tokenized'])),
     ('select_zipped_column',
        ColumnSelector(['zipped'])),
     ('vectorize_tokenized',
        SmilesVectorizer(columns=['zipped'],
                         save_as=['smiles_vector'],
                         atom_regex=r'\[.+?\]',
                         atom_functions=af,
                         struct_functions=sf)),
     ('select_smiles_vector_column',
        ColumnSelector(['smiles_vector']))])

df = pipeline.fit_transform(data)
print(df)

