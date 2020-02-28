"""todo"""
from deepsmiles import Converter as DeepSmilesConverter
from pandas import read_csv
from re import findall
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from src.xrdkit import rsmiles


class PipelineTransformer(BaseEstimator, TransformerMixin):
    """todo"""

    def fit(self, X):
        """todo"""
        return self


class ColumnSelector(PipelineTransformer):
    """todo"""

    def __init__(self, column_names):
        """todo"""
        self.column_names = column_names

    def transform(self, X):
        """todo"""
        return X.loc[:, self.column_names]


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

    def _check_types(self):
        """todo"""
        if not (callable(self.function) \
                and isinstance(self.columns, list) \
                and isinstance(self.save_as, list)):
            raise TypeError('todo')

    def transform(self, X):
        """todo"""
        self._check_types()
        for idx, column in enumerate(self.columns):
            applied = X[column].apply(
                lambda x: self.function(x, **self.function_kwargs))
            if self.save_as and self.save_as[idx]:
                X[self.save_as[idx]] = applied
            else:
                X[column] = applied
        return X


class SmilesRebuilder(FunctionApplier):
    """todo"""

    def __init__(self, columns, save_as=[], **function_kwargs):
        """todo"""
        self.function = rsmiles
        self._assign_rest(columns, save_as, **function_kwargs)


class DeepSmilesEncoder(FunctionApplier):
    """todo"""

    def __init__(self, columns, save_as=[], branches=True, rings=True):
        """todo"""
        self.function = DeepSmilesConverter(branches, rings).encode
        self._assign_rest(columns, save_as)


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


# todo: change to `LengthFilter`
class TooLongFilter(FunctionApplier):
    """todo"""

    def __init__(self, columns, save_as=[], max_len=100):
        """todo"""
        self.function = self.filter_too_long
        self._assign_rest(columns, save_as, max_len=max_len)

    @staticmethod
    def filter_too_long(text, max_len):
        """todo"""
        return text if len(text) <= max_len else None


data = read_csv('datasets/tox21.csv')
cs = ColumnSelector(['smiles'])
sr = SmilesRebuilder(columns=['smiles'],
                     allBondsExplicit=True,
                     allHsExplicit=True,
                     canonical=True,
                     isomericSmiles=True)
dse = DeepSmilesEncoder(columns=['smiles'],
                        save_as=['deep_smiles'],
                        branches=True,
                        rings=True)
tt = TextTokenizer(columns=['deep_smiles'],
                   regex=r'\[.*?\]|%\d{2}|\)+|[\d\(\)\-/\\:=#\$\.]')
tlf = TooLongFilter(columns=['deep_smiles'],
                    max_len=100)

pipeline = Pipeline(
    [('select_smiles_column', cs),
     ('rebuild_smiles', sr),
     ('translate_to_deepsmiles', dse),
     ('tokenize_deepsmiles', tt),
     ('filter_too_long_lists', tlf)]
)

print(pipeline.fit_transform(data))

