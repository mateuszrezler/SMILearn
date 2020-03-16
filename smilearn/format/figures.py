from matplotlib.pyplot import show, subplots, tick_params
from rdkit.Chem.Draw import MolToImage


def feature_matrix(X, compound_num, rows, tokens=None, xlabels=None):
    if tokens:
        ylabels = [f'{index}. {token}'
                   for index, token in enumerate(tokens[compound_num])]
    _, axis = subplots(figsize=(20, 100))
    tick_params(bottom=True, labeltop=True, labelbottom=True, top=True)
    heatmap(around(X[compound_num].reshape(rows, len(xlabels)), 3), ax=axis,
            cmap='YlGn', xticklabels=xlabels, yticklabels=ylabels, annot=True)
    show()


def hide_atoms(atoms):
    return {atom: [.75, .75, .75] for atom in atoms}


def molimg(mol, size=(320, 320), **mol_to_image_kwargs):
    display(MolToImage(mol, size, **mol_to_image_kwargs))

