from matplotlib.pyplot import show, subplots, tick_params
from rdkit.Chem.Draw import MolToImage
from seaborn import heatmap


def feature_matrix(X, compound_num, rows, tokens=None,
                   xlabels=None, size=(20, 120)):
    if tokens is not None:
        ylabels = [f'{index}. {token}'
                   for index, token in enumerate(tokens[compound_num])]
    else:
        ylabels = 'auto'
    _, axis = subplots(figsize=size)
    tick_params(bottom=True, labeltop=True, labelbottom=True, top=True)
    heatmap(X[compound_num].reshape(rows, len(xlabels)), annot=True, ax=axis,
            cbar=False, cmap='YlGn', fmt='.1g',
            xticklabels=xlabels, yticklabels=ylabels)
    show()


def hide_atoms(atoms):
    return {atom: [.75, .75, .75] for atom in atoms}


def molimg(mol, size=(320, 320), **mol_to_image_kwargs):
    display(MolToImage(mol, size, **mol_to_image_kwargs))

