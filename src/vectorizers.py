from pandas import DataFrame
from numpy import array
from re import findall, search
from rdkit import Chem
from deepsmiles import Converter as DeepSmilesConverter


class VectorizationRules(list):

    def __init__(self, other, **dataframe_args):
        if isinstance(other, list):
            self.extend(other)
            self.as_dataframe = DataFrame(self, **dataframe_args)
        elif isinstance(other, DataFrame):
            self.extend(other.values.tolist())
            self.as_dataframe = other

    @classmethod
    def from_csv(cls, csv_file, **dataframe_args):
        df = read_csv(csv_file, **dataframe_args)
        return cls(df)

    def to_dataframe(self):
        return DataFrame(self)

    def to_csv(self, file_name, index=False, sep='\t', **dataframe_args):
        self.as_dataframe.to_csv(
            file_name, index=index, sep=sep, **dataframe_args
        )


class TextTokenizer(object):
    r"""
    TextTokenizer splits a given string into a list of tokens (n-grams of
    characters) according to a specified regular expression, after
    modifications done by user defined modifier.
    """
    def __init__(self, modifier=lambda x: x, regex='.'):
        """
        Parameters
        ----------
        modifier: function, default lambda x: x
            Function for string modification. Dummy function `lambda x: x`
            is set by default, thus nothing is changed.
        regex: str, default '.'
            Regular expression for string tokenization. Single dot (indicating
            any character) is set by default, thus string is tokenized into
            single characters.
        """
        self.modifier = modifier
        self.regex = regex

    def __call__(self, text):
        """Tokenizes a given string."""
        return self.tokenize(text)

    def tokenize(self, text):
        """Tokenizes a given string."""
        text = self.modifier(text)
        tokens = findall(self.regex, text)
        return tokens


class MolVectorizer(object):
    r"""
    MolVectorizer converts a RDKit Mol object into a customizable mol vector
    in x steps:
    1. The Mol object is turned into a specified SMILES string.
    2. The SMILES string is modified user defined modifier and divided into
       a list of tokens (n-grams of characters, smilist) according to
       a specified regular expression. Smilists larger than a given maximum
       length are not taken into account in the next steps of conversion.
    3. Smilist of appropriate size are vectorized by user defined vectorization
       rules and padded with zeros to the maximum lenght. Finally, there is
       also a possibility to add additional zero-padded margins on both sides
       of the mol vector after length normalization.
    The resulting mol vector serves as the input for machine learning models.
    """
    def __init__(self,
                 all_bonds=False,
                 all_hydrogens=False,
                 kekule=False,
                 smiles_modifier=lambda x: x,
                 smiles_regex=r'[A-Z][a-z]|[A-Zcnops]|as|se|.',
                 atom_regex=r'[A-Z][a-z]|[A-Zcnops]|as|se',
                 atom_vect_rules=[['token', 'str']],
                 struct_vect_rules=[['token', 'str']],
                 max_len=250,
                 pad_len=10):
        """
        Parameters
        ----------
        all_bonds: bool, default False
            Determining whether the new SMILES string will have all bonds
            explicit, including single and aromatic bonds.
        all_hydrogens: bool, default False
            Determining whether the new SMILES string will have all hydrogens
            explicit, hence all atoms in square brackets.
        kekule: bool, default False
            Determining whether the new SMILES string will have Kekule
            representation (all atom symbols with first uppercase letter and
            all aromatic bonds explicit).
        smiles_modifier: function, default lambda x: x
            Function for SMILES string modification. Dummy function
            `lambda x: x` is set by default, thus nothing is changed.
        regex: str, default r'[A-Z][a-z]|[A-Zcnops]|as|se|.'
            Regular expression for SMILES string tokenization. By default,
            atomic symbols and other characters are extracted.
        atom_regex: str, default r'[A-Z][a-z]|[A-Zcnops]|as|se'
            Regular expression for filtering atom tokens from SMILES string.
            By default, bare atomic symbols are extracted.
        """
        self.all_bonds = all_bonds
        self.all_hydrogens = all_hydrogens
        self.kekule = kekule
        self.smiles_modifier = smiles_modifier
        self.smiles_regex = smiles_regex
        self.atom_regex = atom_regex
        self.atom_vect_rules = atom_vect_rules
        self.struct_vect_rules = struct_vect_rules
        self.max_len = max_len
        self.pad_len = pad_len
        self.vectorizers_len = len(atom_vect_rules)+len(struct_vect_rules)
        self.pad_size = pad_len*self.vectorizers_len

    @staticmethod
    def __append_token_scalar(
            token_vector, mol, atom_index, token, token_vectorizers
    ):
        for token_vectorizer in token_vectorizers:
            code, str_astype = token_vectorizer
            code_glob_vars = {
                '__builtins__': {},
                'idx': atom_index,
                'mol': mol,
                'token': token
            }
            str_astype_glob_vars = {
                '__builtins__': {
                    'bool': bool,
                    'float': float,
                    'int': int,
                    'str': str
                }
            }
            astype = eval(str_astype, str_astype_glob_vars)
            token_scalar = astype(eval(code, code_glob_vars))
            token_vector.append(token_scalar)

    def __fill_with_zeros(self, mol_vector):
        max_mol_vector = self.max_len*self.vectorizers_len
        empty_space = max_mol_vector-len(mol_vector)+self.pad_size
        mol_vector.extend([0]*empty_space)

    def __pad_with_zeros(self, mol_vector):
        mol_vector.extend([0]*self.pad_size)

    @staticmethod
    def __is_atom(token, atom_regex):
        return atom_regex and search(atom_regex, token)

    def __vectorize(self, smilist, mol, mol_vector):
        if len(smilist) > self.max_len:
            warn(f'Too long SMILES ({len(smilist)} tokens)'
                  ' is found and ignored.')
            return []
        self.__pad_with_zeros(mol_vector)
        atom_index = 0
        for idx, token in enumerate(smilist):
            token_vector = []
            if self.__is_atom(token, self.atom_regex):
                self.__append_token_scalar(
                    token_vector, mol, atom_index, token, self.atom_vect_rules
                )
                token_vector.extend([0]*len(self.struct_vect_rules))
                atom_index += 1
            else:
                token_vector.extend([0]*(len(self.atom_vect_rules)))
                self.__append_token_scalar(
                    token_vector, mol, None, token, self.struct_vect_rules
                )
            mol_vector.extend(token_vector)
        self.__fill_with_zeros(mol_vector)
        self.__pad_with_zeros(mol_vector)
        return array(mol_vector)

    def vectorize(self, mol):
        mol_vector = []
        smiles = Chem.MolToSmiles(mol, kekuleSmiles=self.kekule,
                                  allBondsExplicit=self.all_bonds,
                                  allHsExplicit=self.all_hydrogens)
        tt = TextTokenizer(self.smiles_modifier, self.smiles_regex)
        smilist = tt.tokenize(smiles)
        mol_vector = self.__vectorize(smilist, mol, mol_vector)
        return mol_vector

avr = VectorizationRules.from_csv(pathjoin('vect_rules', 'avr1.tsv'), sep='\t')
svr = VectorizationRules.from_csv(pathjoin('vect_rules', 'svr1.tsv'), sep='\t')
dsc = DeepSmilesConverter(branches=True, rings=True)
mv = MolVectorizer(atom_vect_rules=avr, struct_vect_rules=svr,
                   smiles_modifier=dsc.encode,
                   smiles_regex=r'\[.*?\]|%\d{2}|\)+|[\d\(\)\-/\\:=#\$\.]',
                   atom_regex=r'\[.+?\]', all_hydrogens=True, all_bonds=True)
molvec = mv.vectorize(Chem.MolFromSmiles(
    'CN1CCC[C@H]1C1=CC=C[N+]([2H])=C1.[Cl-]'))
smimat = molvec.reshape(-1, 43)
df_cols = \
    ('Symbol',) \
    + tuple('HCON') \
    + ('*', 'NHs', 'Deg', 'Chg', 'Val', 'Rng', 'Aro', 'CW', 'CCW', 'Chi*',
       'S_', 'SP', 'SP2', 'SP3', 'SP3D', 'SP3D2', 'Hyb*') \
    + tuple('()[].:=#\\/@+-234567<>')
df = DataFrame(smimat, columns=df_cols)
df[['H', 'C']] = df[['H', 'C']].apply(lambda x: to_numeric(x, errors='coerce',
                                                           downcast='integer'))

print(df[10:15])

