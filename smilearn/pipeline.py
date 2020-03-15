from deepsmiles import Converter as DeepSmilesConverter
from numpy import array
from re import findall, search
from rdkit.Chem import MolFromSmiles, MolToSmiles
from sklearn.base import BaseEstimator, TransformerMixin


class PipelineTransformer(BaseEstimator, TransformerMixin):

    def fit(self, X, y=None):
        return self


class FunctionApplier(PipelineTransformer):

    def __init__(self, function, **function_kwargs):
        self.function = function
        self.function_kwargs = function_kwargs

    def transform(self, X, y=None):
        X_new = DataFrame()
        for column in X.columns:
            X_new[column] = X[column].apply(
                lambda x: self.function(x, **self.function_kwargs))
        return X_new


class NanFiller(PipelineTransformer):

    def transform(self, X, y=None):
        return X.fillna(0)


class ToArrayConverter(PipelineTransformer):

    def transform(self, X, y=None):
        return array(X.values.tolist())


class DeepSmilesEncoder(FunctionApplier):

    def __init__(self, **function_kwargs):
        self.function = self.encode
        self.function_kwargs = function_kwargs

    @staticmethod
    def encode(smiles, branches=True, rings=True):
        dsc = DeepSmilesConverter(branches, rings)
        return dsc.encode(smiles)


class RegexTokenizer(FunctionApplier):

    def __init__(self, **function_kwargs):
        self.function = self.tokenize
        self.function_kwargs = function_kwargs

    @staticmethod
    def tokenize(text, regex='.'):
        return findall(regex, text)


class RingTagInserter(FunctionApplier):

    def __init__(self, **function_kwargs):
        self.function = self.insert
        self.function_kwargs = function_kwargs

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


class SmilesFeaturizer(PipelineTransformer):

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
        X_new = DataFrame(X)
        smiles_series = X_new[self.smiles_column].values
        tokens_series = X_new[self.tokens_column].values
        return self.featurize_series(smiles_series, tokens_series)

    def featurize_series(self, smiles_series, tokens_series):
        featurized_series = []
        for index, tokens in enumerate(tokens_series):
            if len(tokens) > self.max_len:
                print(f'Too many tokens ({len(tokens)})\n{tokens}')
                smiles_vector = [0]*(self.max_zeros + 2*self.pad_zeros)
            else:
                smiles_vector = self.featurize(smiles_series[index], tokens)
            featurized_series.append(smiles_vector)
        arr = array(featurized_series)
        return arr.reshape(len(smiles_series), -1, self.n_func)

    def featurize(self, smiles, tokens):
        smiles_vector = []
        smiles_vector.extend([0]*self.pad_zeros)
        mol = MolFromSmiles(smiles)
        token_vectors = self._calculate_token_vectors(mol, tokens)
        smiles_vector.extend(token_vectors)
        fill_zeros = self.max_zeros - len(smiles_vector) + self.pad_zeros
        smiles_vector.extend([0]*fill_zeros)
        smiles_vector.extend([0]*self.pad_zeros)
        return array(smiles_vector)


class SmilesRebuilder(FunctionApplier):

    def __init__(self, **function_kwargs):
        self.function = self.rebuild
        self.function_kwargs = function_kwargs

    @staticmethod
    def rebuild(smiles, **function_kwargs):
        mol = MolFromSmiles(smiles)
        return MolToSmiles(mol, **function_kwargs)

