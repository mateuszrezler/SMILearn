import pandas as pd
import re
from rdkit import Chem


def rebuild_smiles(smiles, **options):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), **options)


class FunctionList(list):
    shortcuts = {
        'GetSymbol': 'Sym',
        'GetTotalNumHs': 'NHs',
        'GetTotalDegree': 'Deg',
        'GetFormalCharge': 'Chg',
        'GetTotalValence': 'Val',
        'IsInRing': 'Rng',
        'GetIsAromatic': 'Aro',
        'GetChiralTag': 'Chi',
        'GetHybridization': 'Hyb',
        'CHI_TETRAHEDRAL_CW': 'CW',
        'CHI_TETRAHEDRAL_CCW': 'CCW',
        'CHI_OTHER': 'Chi*',
        'OTHER': 'Hyb*',
    }

    def get_names(self):
        names = []
        for function in self:
            def dict_choice(n):
                return self.shortcuts[function[n]] \
                    if self.shortcuts and function[n] in self.shortcuts \
                    else function[n]
            if len(function) == 1:
                names.append(dict_choice(0))
            else:
                if type(function[1]) is str:
                    names.append(dict_choice(1))
                elif type(function[1]) is list:
                    names.append('['+''.join(function[1])+']')
                elif type(function[1]) is tuple:
                    names.append('!['+''.join(function[1])+']')
                else:
                    names.append(dict_choice(0))
        return names


class SMItrix(list):

    def __init__(self, atom_functions, struct_functions, smiles, **options):
        self.__atom_functions = atom_functions
        self.__struct_functions = struct_functions
        self.__smiles = rebuild_smiles(smiles, **options)
        self.__mol_smiles = rebuild_smiles(smiles, isomericSmiles=True,
                                           kekuleSmiles=False)
        self.__options = options
        self.__build()

    def apply(self, function):
        self.__smiles = function(self.__smiles)
        self.__build()

    def to_smiles(self):
        return self.__smiles

    def __build(self):

        def set_append_type(feature):
            if len(function) == 1:
                vector.append(feature)
            else:
                value = function[1]
                if type(value) is str:
                    vector.append(int(str(feature) == value))
                elif type(value) is list:
                    vector.append(int(feature in value))
                elif type(value) is tuple:
                    vector.append(int(feature not in value))
                elif type(value) in (int, float):
                    vector.append(feature*value)

        self.clear()
        atom_index = 0
        mol = Chem.MolFromSmiles(self.__mol_smiles)
        for idx, char in enumerate(re.findall('.', self.__smiles)):
            vector = []
            if char == 'H':
                all_features_len = len(self.__atom_functions) \
                    + len(self.__struct_functions)
                vector.extend([1] + [0]*(all_features_len-1))
            elif char.isupper():
                for function in self.__atom_functions:
                    e = eval(f'mol.GetAtomWithIdx(atom_index).{function[0]}()')
                    set_append_type(e)
                vector.extend([0]*len(self.__struct_functions))
                atom_index += 1
            elif char.islower():
                continue
            else:
                vector.extend([0]*(len(self.__atom_functions)))
                for function in self.__struct_functions:
                    set_append_type(char)
            self.append(vector)


flist = FunctionList([['GetSymbol', 'H'],
                      ['GetSymbol', 'C'],
                      ['GetSymbol', 'O'],
                      ['GetSymbol', 'N'],
                      ['GetSymbol', ('H', 'C', 'O', 'N')],
                      ['GetTotalNumHs', 1/8],
                      ['GetTotalDegree', 1/4],
                      ['GetFormalCharge', 1/8],
                      ['GetTotalValence', 1/8],
                      ['IsInRing', 1],
                      ['GetIsAromatic', 1],
                      ['GetChiralTag', 'CHI_TETRAHEDRAL_CW'],
                      ['GetChiralTag', 'CHI_TETRAHEDRAL_CCW'],
                      ['GetChiralTag', 'CHI_OTHER'],
                      ['GetHybridization', 'S'],
                      ['GetHybridization', 'SP'],
                      ['GetHybridization', 'SP2'],
                      ['GetHybridization', 'SP3'],
                      ['GetHybridization', 'SP3D'],
                      ['GetHybridization', 'SP3D2'],
                      ['GetHybridization', 'OTHER']])
single_chars = [['char', char] for char in '()[].:=#\\/@+-']
for i in range(2, 8):
    single_chars.append(['digit', str(i)])
single_chars.extend([['opened', '<'],
                     ['closed', '>']])
slist = FunctionList(single_chars)

smitrix = SMItrix(flist, slist, '[Zn+2].c1c[n-]nc1-c1cocc1.[O-]',
                  isomericSmiles=True, kekuleSmiles=True,
                  allBondsExplicit=True)
df = pd.DataFrame(smitrix, columns=flist.get_names()+slist.get_names())
print(df)

