compound = '[2H]C[C@H](I)Br'
mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(compound)))
kek = Chem.MolToSmiles(Chem.MolFromSmiles(compound), kekuleSmiles=True)
function_list = [
                 ['GetSymbol', 'H'],
                 ['GetSymbol', 'C'],
                 ['GetSymbol', 'O'],
                 ['GetSymbol', 'N'],
                 ['GetSymbol', {'H', 'C', 'O', 'N'}]
]
atom_index = 0
matrix = []
for char in kek:
    vector = []
    if char == 'H':
        vector.extend([1] + [0]*len(function_list))
    elif char.isupper():
        features = []
        for function in function_list:
            feature = eval(f'mol.GetAtomWithIdx(atom_index).{function[0]}()')
            if len(function) == 1:
                features.append(feature)
            else:
                if type(function[1]) is str:
                    features.append(int(feature == function[1]))
                elif type(function[1]) is tuple:
                    features.append(int(feature in function[1]))
                elif type(function[1]) is set:
                    features.append(int(feature not in function[1]))
                elif type(function[1]) in (int, float):
                    features.append(feature*function[1])
        vector.extend(features + [0])
        atom_index += 1
    elif char.islower():
        continue
    else:
        vector.extend([0]*len(function_list) + [char])
    matrix.append(vector)
df = pd.DataFrame(matrix)
print(df)

