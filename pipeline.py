"""todo"""
from deepsmiles import Converter as DeepSmilesConverter
from numpy import array
from pandas import read_csv, Series
from rdkit import Chem
from re import findall, search
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from src.xrdkit import rsmiles


class PipelineTransformer(BaseEstimator, TransformerMixin):
    """todo"""

    def fit(self, X):
        """todo"""
        return self


class ColumnRenamer(PipelineTransformer):
    """todo"""

    def __init__(self, renamer_dict):
        """todo"""
        self.renamer_dict = renamer_dict

    def transform(self, X):
        """todo"""
        X.columns = [self.renamer_dict[column]
                     if column in self.renamer_dict
                     else column
                     for column in X.columns]
        return X


class ColumnSelector(PipelineTransformer):
    """todo"""

    def __init__(self, columns):
        """todo"""
        self.columns = columns

    def transform(self, X):
        """todo"""
        return X.loc[:, self.columns]


class Zipper(PipelineTransformer):
    """todo"""

    def __init__(self, columns, zip_column='zipped'):
        """todo"""
        self.columns = columns
        self.zip_column = zip_column

    def transform(self, X):
        """todo"""
        values = [X[column].values for column in self.columns]
        X[self.zip_column] = Series(zip(*values))
        return X


class FunctionApplier(PipelineTransformer):
    """todo"""

    def __init__(self, function, columns, save_as=[], **function_kwargs):
        """todo"""
        self.function = function
        self._assign_rest(columns, save_as, **function_kwargs)

    def _assign_rest(self, columns, save_as, **function_kwargs):
        """todo"""
        self.columns = columns
        self.save_as = save_as
        self.function_kwargs = function_kwargs

    def transform(self, X):
        """todo"""
        for idx, column in enumerate(self.columns):
            applied = X[column].apply(
                lambda x: self.function(x, **self.function_kwargs))
            if not self.save_as:
                X[column] = applied
            else:
                X[self.save_as[idx]] = applied
        return X


class DeepSmilesEncoder(FunctionApplier):
    """todo"""

    def __init__(self, columns, save_as=[], branches=True, rings=True):
        """todo"""
        self.function = DeepSmilesConverter(branches, rings).encode
        self._assign_rest(columns, save_as)


class LenFilter(FunctionApplier):
    """todo"""

    def __init__(self, columns, save_as=[], min_len=0, max_len=100,
                 bounded=(True, True)):
        """todo"""
        self.function = self.filter_len
        self._assign_rest(columns, save_as, min_len=min_len, max_len=max_len,
                          bounded=bounded)

    @staticmethod
    def filter_len(obj, min_len, max_len, bounded):
        """todo"""
        signs = ['<'+'='*bounded[i] for i in range(2)]
        return obj if eval(f'min_len {signs[0]} len(obj) {signs[1]} max_len') \
            else None


class NanDropper(PipelineTransformer):
    """todo"""

    def transform(self, X):
        """todo"""
        return X.dropna().reset_index(drop=True)


class SmilesRebuilder(FunctionApplier):
    """todo"""

    def __init__(self, columns, save_as=[], **function_kwargs):
        """todo"""
        self.function = rsmiles
        self._assign_rest(columns, save_as, **function_kwargs)


class SmilesVectorizer(FunctionApplier):
    """todo"""
    def __init__(self,
                 smiles_tuple,
                 save_as=None,
                 atom_regex=r'[A-Z][a-z]|[A-Z]',
                 atom_vect_rules=[],
                 struct_vect_rules=[],
                 pad_len=0):
        if save_as is None:
            save_as = []
        else:
            save_as = [save_as]
        self.function = self.vectorize
        self._assign_rest([smiles_tuple],
                          save_as,
                          atom_regex=atom_regex,
                          atom_vect_rules=atom_vect_rules,
                          struct_vect_rules=struct_vect_rules,
                          pad_len=pad_len)

    @staticmethod
    def vectorize(smiles_tuple, atom_regex, atom_vect_rules,
                  struct_vect_rules, pad_len, max_len=100):

        def append_atom_features(token_vector, mol, atom_index,
                                 atom_vect_rules):
            for atom_featurizer in atom_vect_rules:
                token_vector.append(atom_featurizer(mol, atom_index))

        def append_struct_features(token_vector, mol, token,
                                   struct_vect_rules):
            for struct_featurizer in struct_vect_rules:
                token_vector.append(struct_featurizer(token))

        def calculate_zeros():
            fill_multiplier = len(atom_vect_rules)+len(struct_vect_rules)
            pad_multiplier = pad_len*fill_multiplier
            return fill_multiplier, pad_multiplier

        def fill_with_zeros(mol_vector):
            fill_multiplier, pad_multiplier = calculate_zeros()
            max_mol_vector = max_len*fill_multiplier
            empty_space = max_mol_vector-len(mol_vector)+pad_multiplier
            mol_vector.extend([0]*empty_space)

        def pad_with_zeros(mol_vector):
            pad_multiplier = calculate_zeros()[1]
            mol_vector.extend([0]*pad_multiplier)

        smiles, tokens = smiles_tuple
        mol = Chem.MolFromSmiles(smiles)
        mol_vector = []
        pad_with_zeros(mol_vector)
        atom_index = 0
        for idx, token in enumerate(tokens):
            token_vector = []
            if search(atom_regex, token):
                append_atom_features(token_vector, mol, atom_index,
                                     atom_vect_rules)
                token_vector.extend([0]*len(struct_vect_rules))
                atom_index += 1
            else:
                token_vector.extend([0]*(len(atom_vect_rules)))
                append_struct_features(token_vector, mol, token,
                                       struct_vect_rules)
            mol_vector.extend(token_vector)
        fill_with_zeros(mol_vector)
        pad_with_zeros(mol_vector)
        return array(mol_vector)


class TextTokenizer(FunctionApplier):
    """todo"""

    def __init__(self, columns, save_as=[], regex='.'):
        """todo"""
        self.function = self.tokenize
        self._assign_rest(columns, save_as, regex=regex)

    @staticmethod
    def tokenize(text, regex):
        """todo"""
        return findall(regex, text)


data = read_csv('datasets/tox21.csv')
avr = [lambda mol, i: mol.GetAtomWithIdx(i).GetSymbol()] \
    + [lambda mol, i: int(mol.GetAtomWithIdx(i).GetSymbol() == atom)
       for atom in 'HCON'] \
    + [lambda mol, i: int(mol.GetAtomWithIdx(i).GetSymbol()
       not in ('H', 'C', 'O', 'N')),
       lambda mol, i: int(mol.GetAtomWithIdx(i).GetTotalNumHs()),
       lambda mol, i: int(mol.GetAtomWithIdx(i).GetTotalDegree()),
       lambda mol, i: int(mol.GetAtomWithIdx(i).GetFormalCharge()),
       lambda mol, i: int(mol.GetAtomWithIdx(i).GetTotalValence()),
       lambda mol, i: int(mol.GetAtomWithIdx(i).IsInRing()),
       lambda mol, i: int(mol.GetAtomWithIdx(i).GetIsAromatic())] \
    + [lambda mol, i: int(str(mol.GetAtomWithIdx(i).GetChiralTag()) == tag)
       for tag in ['CHI_TETRAHEDRAL_CW', 'CHI_TETRAHEDRAL_CCW', 'CHI_OTHER']] \
    + [lambda mol, i: int(str(mol.GetAtomWithIdx(i).GetHybridization()) == hyb)
       for hyb in ['S', 'SP2', 'SP3', 'SP3D', 'SP3D2', 'OTHER']]
svr = [lambda x: int(x == char) for char in '()[].:=#\\/@+-234567<>']

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
        TextTokenizer(columns=['deep_smiles'],
                      regex=r'\[.*?\]|%\d{2}|\)+|[\d\(\)\-/\\:=#\$\.]')),
     ('rename_column',
        ColumnRenamer({'deep_smiles': 'tokenized'})),
     ('filter_length',
        LenFilter(columns=['tokenized'],
                  min_len=1,
                  max_len=3,
                  bounded=(False, True))),
     ('drop_nans',
        NanDropper()),
     ('zip_smiles_and_tokenized',
        Zipper(['smiles', 'tokenized'])),
     ('select_zipped_column',
        ColumnSelector(['zipped'])),
     ('vectorize_tokenized',
        SmilesVectorizer(smiles_tuple='zipped',
                         save_as='smiles_vector',
                         atom_regex=r'\[.+?\]',
                         atom_vect_rules=avr,
                         struct_vect_rules=svr,
                         pad_len=0)),
     ('select_smiles_vector_column',
        ColumnSelector(['smiles_vector']))])

print(pipeline.fit_transform(data))
