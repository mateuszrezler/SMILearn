from pandas import DataFrame
from numpy import array
from re import findall, search
from rdkit import Chem
from deepsmiles import Converter as DeepSmilesConverter


class MolVectorizer(object):

    def __init__(self, kekule=False, all_bonds=False, all_hydrogens=False,
                 smiles_modifier=None, smiles_regex=None, atom_regex=None,
                 atom_featurizers=None, struct_featurizers=None):
        self.kekule = kekule
        self.all_bonds = all_bonds
        self.all_hydrogens = all_hydrogens
        self.smiles_modifier = smiles_modifier
        self.smiles_regex = smiles_regex
        self.atom_regex = atom_regex
        self.atom_featurizers = atom_featurizers
        self.struct_featurizers = struct_featurizers

    def vectorize(self, mol):

        def is_atom(symbol):
            return self.atom_regex and search(self.atom_regex, symbol)

        smivec = []
        atom_index = 0
        smiles = Chem.MolToSmiles(mol, kekuleSmiles=self.kekule,
                                  allBondsExplicit=self.all_bonds,
                                  allHsExplicit=self.all_hydrogens)
        if self.smiles_modifier:
            smiles = self.smiles_modifier(smiles)
        smilist = smiles
        if self.smiles_regex:
            smilist = findall(self.smiles_regex, smiles)
        for idx, symbol in enumerate(smilist):
            vector = []
            if is_atom(symbol):
                for featurizer in self.atom_featurizers:
                    code = featurizer[0]
                    astype = featurizer[1]
                    feature = astype(eval('mol.' + code))
                    vector.append(feature)
                vector.extend([0]*len(self.struct_featurizers))
                atom_index += 1
            else:
                vector.extend([0]*(len(self.atom_featurizers)))
                for featurizer in self.struct_featurizers:
                    code = featurizer[0]
                    astype = featurizer[1]
                    feature = astype(eval(code))
                    vector.append(feature)
            smivec.extend(vector)
        return smivec


atom_alias = 'GetAtomWithIdx(atom_index).'
af = [[atom_alias + 'GetSymbol()', str],
      [atom_alias + "GetSymbol() == 'H'", int],
      [atom_alias + "GetSymbol() == 'C'", int],
      [atom_alias + "GetSymbol() == 'O'", int],
      [atom_alias + "GetSymbol() == 'N'", int],
      [atom_alias + "GetSymbol() not in ('H', 'C', 'O', 'N')", int],
      [atom_alias + 'GetTotalNumHs()', int],
      [atom_alias + 'GetTotalDegree()', int],
      [atom_alias + 'GetFormalCharge()', int],
      [atom_alias + 'GetTotalValence()', int],
      [atom_alias + 'IsInRing()', int],
      [atom_alias + 'GetIsAromatic()', int],
      [atom_alias + "GetChiralTag() == 'CHI_TETRAHEDRAL_CW'", int],
      [atom_alias + "GetChiralTag() == 'CHI_TETRAHEDRAL_CCW'", int],
      [atom_alias + "GetChiralTag() == 'CHI_OTHER'", int],
      [atom_alias + "GetHybridization() == 'S'", int],
      [atom_alias + "GetHybridization() == 'SP'", int],
      [atom_alias + "GetHybridization() == 'SP2'", int],
      [atom_alias + "GetHybridization() == 'SP3'", int],
      [atom_alias + "GetHybridization() == 'SP3D'", int],
      [atom_alias + "GetHybridization() == 'SP3D2'", int],
      [atom_alias + "GetHybridization() == 'OTHER'", int]]
sf = [[f'symbol == {repr(char[0])}', int] for char in '()[].:=#\\/@+-234567<>']
dsc = DeepSmilesConverter(branches=True, rings=True)
mv = MolVectorizer(atom_featurizers=af, struct_featurizers=sf,
                   smiles_modifier=dsc.encode,
                   smiles_regex=r'\[.*?\]|%\d{2}|\)+|[\d\(\)\-/\\:=#\$\.]',
                   atom_regex=r'\[.+?\]', all_hydrogens=True, all_bonds=True)
smivec = mv.vectorize(Chem.MolFromSmiles(
    'CN1CCC[C@H]1C1=CC=C[N+]([2H])=C1.[Cl-]'))
smimat = array(smivec).reshape(-1, 43)
df_cols = \
    ('Symbol',) \
    + tuple('HCON') \
    + ('*', 'NHs', 'Deg', 'Chg', 'Val', 'Rng', 'Aro', 'CW', 'CCW', 'Chi*',
       'S_', 'SP', 'SP2', 'SP3', 'SP3D', 'SP3D2', 'Hyb*') \
    + tuple('()[].:=#\\/@+-234567<>')
df = DataFrame(smimat, columns=df_cols)
print(df)

