from deepsmiles import Converter as DeepSmilesConverter
from numpy import array
from pandas import DataFrame, Series
from rdkit.Chem import MolFromSmiles, MolToSmiles
from re import findall, search
from sklearn.base import BaseEstimator, TransformerMixin


class PipelineTransformer(BaseEstimator, TransformerMixin):

    def fit(self, X, y=None):
        return self


class FunctionApplier(PipelineTransformer):

    def __init__(self, function, astype=Series, **function_kwargs):
        self.function = function
        self.astype = astype
        self.function_kwargs = function_kwargs

    def __call__(self, smiles):
        return self.function(smiles)

    def apply(self, smiles):
        return self.function(smiles)

    def transform(self, X, y=None):
        return self.astype(X.apply(self.function))


class SeriesSelector(PipelineTransformer):

    def __init__(self, series):
        self.series = series

    def transform(self, X, y=None):
        return X[self.series]


class ToArrayConverter(PipelineTransformer):

    def transform(self, X, y=None):
        return array(X.values.tolist())


class ZerosFiller(PipelineTransformer):

    def transform(self, X, y=None):
        return X.fillna(0)


class DeepSmilesEncoder(FunctionApplier):

    def __init__(self, branches=True, rings=True):
        super().__init__(self.encode)
        self.encoder = DeepSmilesConverter(branches, rings)

    def encode(self, smiles):
        k = self.function_kwargs
        return self.encoder.encode(smiles)


class RegexTokenizer(FunctionApplier):

    def __init__(self, astype=Series, regex='.'):
        super().__init__(self.tokenize, astype=astype, regex=regex)

    def tokenize(self, text):
        k = self.function_kwargs
        return findall(k['regex'], text)


class RingTagInserter(FunctionApplier):

    def __init__(self, astype=Series, regex=r'%\d+|[\+\-]\d\]|.',
                 reuse=True, tags='<>'):
        super().__init__(self.insert, astype=astype, regex=regex,
                         reuse=reuse, tags=tags)

    def insert(self, smiles):
        k = self.function_kwargs
        modified_smiles = ''
        smiles_list = findall(k['regex'], smiles)
        used_numbers = []
        for element in smiles_list:
            if element[-1].isdigit():
                if element not in used_numbers:
                    modified_smiles += k['tags'][0]
                    used_numbers.append(element)
                else:
                    modified_smiles += k['tags'][1]
                    if k['reuse']:
                        used_numbers.remove(element)
            else:
                modified_smiles += element
        return modified_smiles


class SmilesFeaturizer(PipelineTransformer):

    def __init__(self, atom_functions, struct_functions, smiles_column=0,
                 tokens_column=1, atom_regex=r'\[.+?\]|Br|Cl|.',
                 h_vector=False, ignore_regex=None, max_len=100, pad_len=0):
        self.atom_functions = atom_functions
        self.struct_functions = struct_functions
        self.smiles_column = smiles_column
        self.tokens_column = tokens_column
        self.atom_regex = atom_regex
        self.h_vector = h_vector
        self.ignore_regex = ignore_regex
        self.max_len = max_len
        self.pad_len = pad_len
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
        smiles_array = X_new[self.smiles_column].values
        tokens_array = X_new[self.tokens_column].values
        return self.featurize_series(smiles_array, tokens_array)

    def featurize_series(self, smiles_series, tokens_series):
        featurized_series = []
        for index, tokens in enumerate(tokens_series):
            if len(tokens) > self.max_len:
                print(f'Too many tokens found ({len(tokens)}) at row {index}.')
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

    def __init__(self, astype=Series, allBondsExplicit=False,
                 allHsExplicit=False, canonical=True, doRandom=False,
                 isomericSmiles=True, kekuleSmiles=False, rootedAtAtom=-1):
        super().__init__(self.rebuild, astype=astype,
                         allBondsExplicit=allBondsExplicit,
                         allHsExplicit=allHsExplicit, canonical=canonical,
                         doRandom=doRandom, isomericSmiles=isomericSmiles,
                         kekuleSmiles=kekuleSmiles, rootedAtAtom=rootedAtAtom)

    def rebuild(self, smiles):
        mol = MolFromSmiles(smiles)
        return MolToSmiles(mol, **self.function_kwargs)

