"""todo"""
from deepsmiles import Converter as DeepSmilesConverter
from numpy import array
from pandas import Series
from re import findall, search
from rdkit.Chem import MolFromSmiles, MolToSmiles
from sklearn.base import BaseEstimator, TransformerMixin


class PipelineTransformer(BaseEstimator, TransformerMixin):

    def fit(self, X, y=None):

        return self


class ColumnRenamer(PipelineTransformer):

    def __init__(self, renamer_dict):

        self.renamer_dict = renamer_dict

    def transform(self, X, y=None):

        X.columns = [self.renamer_dict[column]
                     if column in self.renamer_dict
                     else column
                     for column in X.columns]
        return X


class ColumnDropper(PipelineTransformer):

    def __init__(self, columns):

        self.columns = columns

    def transform(self, X, y=None):

        return X.drop(self.columns, axis=1)


class ColumnSelector(PipelineTransformer):

    def __init__(self, columns):

        self.columns = columns

    def transform(self, X, y=None):

        return X.loc[:, self.columns]


class FunctionApplier(PipelineTransformer):

    def __init__(self, function, columns, save_as=None, **function_kwargs):

        self.function = function
        self._assign_rest(columns, save_as, **function_kwargs)

    def _assign_rest(self, columns, save_as, **function_kwargs):

        self.columns = columns
        self.save_as = save_as
        self.function_kwargs = function_kwargs

    def transform(self, X, y=None):

        for idx, column in enumerate(self.columns):
            applied = X[column].apply(
                lambda x: self.function(x, **self.function_kwargs))
            if not self.save_as:
                X[column] = applied
            else:
                X[self.save_as[idx]] = applied
        return X


class NanDropper(PipelineTransformer):

    def __init__(self, subset=None):

        self.subset = subset

    def transform(self, X, y=None):

        return X.dropna(subset=self.subset).reset_index(drop=True)


class NanFiller(PipelineTransformer):

    def transform(self, X, y=None):

        return X.fillna(0)


class ToArrayConverter(PipelineTransformer):

    def __init__(self, columns=None, tolist=False):

        self.columns = columns
        self.tolist = tolist

    def transform(self, X, y=None):

        if self.columns:
            values = X[self.columns].values
        else:
            values = X.values
        if self.tolist:
            values = values.tolist()
        return array(values)


class Zipper(PipelineTransformer):

    def __init__(self, columns, zip_column='zipped'):

        self.columns = columns
        self.zip_column = zip_column

    def transform(self, X, y=None):

        values = [X[column].values for column in self.columns]
        X[self.zip_column] = Series(zip(*values))
        return X


class DeepSmilesEncoder(FunctionApplier):

    def __init__(self, columns, save_as=None, branches=True, rings=True):

        self.function = self.encode
        self._assign_rest(columns, save_as)

    @staticmethod
    def encode(smiles, branches=True, rings=True):
        dsc = DeepSmilesConverter(branches, rings)
        return dsc.encode(smiles)


class LenFilter(FunctionApplier):

    def __init__(self, columns, save_as=None, min_len=0, max_len=100,
                 bounded=(True, True)):

        self.function = self.filter_len
        self._assign_rest(columns, save_as, min_len=min_len, max_len=max_len,
                          bounded=bounded)

    @staticmethod
    def filter_len(obj, min_len, max_len, bounded):

        signs = ['<' + '='*bounded[i] for i in range(2)]
        return obj if eval(f'min_len {signs[0]} len(obj) {signs[1]} max_len') \
            else None


class RegexTokenizer(FunctionApplier):

    def __init__(self, columns, save_as=None, regex='.'):

        self.function = self.tokenize
        self._assign_rest(columns, save_as, regex=regex)

    @staticmethod
    def tokenize(text, regex='.'):

        return findall(regex, text)


class SmilesRebuilder(FunctionApplier):

    def __init__(self, columns, save_as=None, **mol_to_smiles_kwargs):

        self.function = self.reb
        self._assign_rest(columns, save_as, **mol_to_smiles_kwargs)

    @staticmethod
    def reb(smiles, **mol_to_smiles_kwargs):
        mol = MolFromSmiles(smiles)
        return MolToSmiles(mol, **mol_to_smiles_kwargs)


class SmilesVectorizer(PipelineTransformer):

    def __init__(self,
                 smiles_column=None,
                 tokens_column=None,
                 ignore_regex=None,
                 atom_regex=r'\[.+?\]|Br|Cl|[BCFINOPSbcnops]',
                 atom_functions=[lambda mol, index:
                                 mol.GetAtomWithIdx(index).GetSymbol()],
                 struct_functions=[lambda token: token],
                 max_len=100,
                 pad_len=0,
                 h_vector=False):

        self.smiles_column = smiles_column
        self.tokens_column = tokens_column
        self.ignore_regex = ignore_regex
        self.atom_regex = atom_regex
        self.atom_functions = atom_functions
        self.struct_functions = struct_functions
        self.max_len = max_len
        self.pad_len = pad_len
        self.h_vector = h_vector
        self._calculate_zeros()

    def _append_atom_features(self, token_vector, mol, atom_index):

        for atom_function in self.atom_functions:
            token_vector.append(atom_function(mol, atom_index))

    def _append_struct_features(self, token_vector, token):

        for struct_function in self.struct_functions:
            token_vector.append(struct_function(token))

    def _calculate_token_vectors(self, mol, tokens):

        token_vectors = []
        atom_index = 0
        for token in tokens:
            token_vector = []
            if self.ignore_regex and search(self.ignore_regex, token):
                continue
            elif self.h_vector and token == 'H':
                token_vector.extend([1] + [0]*(self.n_func-1))
            elif search(self.atom_regex, token):
                self._append_atom_features(token_vector, mol, atom_index)
                token_vector.extend([0]*len(self.struct_functions))
                atom_index += 1
            else:
                token_vector.extend([0]*(len(self.atom_functions)))
                self._append_struct_features(token_vector, token)
            token_vectors.extend(token_vector)
        return token_vectors

    def _calculate_zeros(self):
        self.n_func = len(self.atom_functions) + len(self.struct_functions)
        self.max_zeros = self.max_len*self.n_func
        self.pad_zeros = self.pad_len*self.n_func

    def transform(self, X, y=None):

        smiles_series = X[self.smiles_column].values
        tokens_series = X[self.tokens_column].values
        return self.vectorize_series(smiles_series, tokens_series)

    def vectorize_series(self, smiles_series, tokens_series):

        vectorized_series = []
        for index, tokens in enumerate(tokens_series):
            if len(tokens) > self.max_len:
                print(f'Too many tokens ({len(tokens)})\n{tokens}')
                smiles_vector = [0]*(self.max_zeros + 2*self.pad_zeros)
            else:
                smiles_vector = self.vectorize(smiles_series[index], tokens)
            vectorized_series.append(smiles_vector)
        return array(vectorized_series)#.reshape(len(smiles_series), -1,
                                        #        self.n_func)

    def vectorize(self, smiles, tokens):

        smiles_vector = []
        smiles_vector.extend([0]*self.pad_zeros)
        mol = MolFromSmiles(smiles)
        token_vectors = self._calculate_token_vectors(mol, tokens)
        smiles_vector.extend(token_vectors)
        fill_zeros = self.max_zeros - len(smiles_vector) + self.pad_zeros
        smiles_vector.extend([0]*fill_zeros)
        smiles_vector.extend([0]*self.pad_zeros)
        return array(smiles_vector)


class RingTagInserter(FunctionApplier):

    def __init__(self, columns, save_as=None, **function_kwargs):
        self.function = self.insert
        self._assign_rest(columns, save_as, **function_kwargs)

    @staticmethod
    def insert(smiles, tags='<>', reuse=False):
        numbers = []
        modified_smiles = ''
        smiles_list = findall('\d|%\d+|.+?', smiles)
        for index, element in enumerate(smiles_list):
            if element[-1].isdigit() and smiles[index-1] not in '+-':
                if element not in numbers:
                    modified_smiles += tags[0]
                    numbers.append(element)
                else:
                    modified_smiles += tags[1]
                    if reuse:
                        numbers.remove(element)
            else:
                modified_smiles += element
        return modified_smiles

