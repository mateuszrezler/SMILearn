"""todo"""
from deepsmiles import Converter as DeepSmilesConverter
from pandas import Series
from re import findall, search
from rdkit.Chem import MolFromSmiles, MolToSmiles
from sklearn.base import BaseEstimator, TransformerMixin


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


class FunctionApplier(PipelineTransformer):
    """todo"""

    def __init__(self, function, columns, save_as=None, **function_kwargs):
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


class NanDropper(PipelineTransformer):
    """todo"""

    def transform(self, X):
        """todo"""
        return X.dropna().reset_index(drop=True)


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


class DeepSmilesEncoder(FunctionApplier):
    """todo"""

    def __init__(self, columns, save_as=None, branches=True, rings=True):
        """todo"""
        self.function = self.encode
        self._assign_rest(columns, save_as)

    @staticmethod
    def encode(smiles, branches=True, rings=True):
        dsc = DeepSmilesConverter(branches, rings)
        return dsc.encode(smiles)


class LenFilter(FunctionApplier):
    """todo"""

    def __init__(self, columns, save_as=None, min_len=0, max_len=100,
                 bounded=(True, True)):
        """todo"""
        self.function = self.filter_len
        self._assign_rest(columns, save_as, min_len=min_len, max_len=max_len,
                          bounded=bounded)

    @staticmethod
    def filter_len(obj, min_len, max_len, bounded):
        """todo"""
        signs = ['<' + '='*bounded[i] for i in range(2)]
        return obj if eval(f'min_len {signs[0]} len(obj) {signs[1]} max_len') \
            else None


class RegexTokenizer(FunctionApplier):
    """todo"""

    def __init__(self, columns, save_as=None, regex='.'):
        """todo"""
        self.function = self.tokenize
        self._assign_rest(columns, save_as, regex=regex)

    @staticmethod
    def tokenize(text, regex='.'):
        """todo"""
        return findall(regex, text)


class SmilesRebuilder(FunctionApplier):
    """todo"""

    def __init__(self, columns, save_as=None, **mol_to_smiles_kwargs):
        """todo"""
        self.function = self.reb
        self._assign_rest(columns, save_as, **mol_to_smiles_kwargs)

    @staticmethod
    def reb(smiles, **mol_to_smiles_kwargs):
        mol = MolFromSmiles(smiles)
        return MolToSmiles(mol, **mol_to_smiles_kwargs)


class SmilesVectorizer(FunctionApplier):
    """todo"""
    def __init__(self,
                 columns,
                 save_as=None,
                 atom_regex=r'\[.+?\]|Br|Cl|[BCFINOPSbcnops]',
                 atom_functions=[lambda mol, index:
                                 mol.GetAtomWithIdx(index).GetSymbol()],
                 struct_functions=[lambda token: token],
                 max_len=100,
                 pad_len=0):
        self.function = self.vectorize
        self._assign_rest(columns,
                          save_as,
                          atom_regex=atom_regex,
                          atom_functions=atom_functions,
                          struct_functions=struct_functions,
                          max_len=max_len,
                          pad_len=pad_len)

    @staticmethod
    def vectorize(smiles_tuple,
                  atom_regex=r'\[.+?\]|Br|Cl|[BCFINOPSbcnops]',
                  atom_functions=[lambda mol, index:
                                  mol.GetAtomWithIdx(index).GetSymbol()],
                  struct_functions=[lambda token: token],
                  max_len=100,
                  pad_len=0):

        def append_atom_features(token_vector, mol, atom_index,
                                 atom_functions):
            for atom_function in atom_functions:
                token_vector.append(atom_function(mol, atom_index))

        def append_struct_features(token_vector, token, struct_functions):
            for struct_function in struct_functions:
                token_vector.append(struct_function(token))

        def calculate_vec_len_and_pad_zeros():
            vec_len = len(atom_functions) + len(struct_functions)
            return vec_len, pad_len * vec_len

        def fill_with_zeros(mol_vector):
            vec_len, pad_zeros = calculate_vec_len_and_pad_zeros()
            max_zeros = max_len * vec_len
            fill_zeros = max_zeros - len(mol_vector) + pad_zeros
            mol_vector.extend([0] * fill_zeros)

        def pad_with_zeros(mol_vector):
            pad_zeros = calculate_vec_len_and_pad_zeros()[1]
            mol_vector.extend([0] * pad_zeros)

        smiles, tokens = smiles_tuple
        mol = MolFromSmiles(smiles)
        mol_vector = []
        pad_with_zeros(mol_vector)
        atom_index = 0
        for idx, token in enumerate(tokens):
            token_vector = []
            if search(atom_regex, token):
                append_atom_features(token_vector, mol, atom_index,
                                     atom_functions)
                token_vector.extend([0] * len(struct_functions))
                atom_index += 1
            else:
                token_vector.extend([0] * (len(atom_functions)))
                append_struct_features(token_vector, token, struct_functions)
            mol_vector.extend(token_vector)
        fill_with_zeros(mol_vector)
        pad_with_zeros(mol_vector)
        return mol_vector

