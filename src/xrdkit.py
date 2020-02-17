r"""
xrdkit - eXtended RDKit

Module containing extensions and shortcuts for RDKit package.

"""
from rdkit import Chem
from re import findall, sub
from sys import argv


def rsmiles(smiles, *modifiers, **options):
    r"""
    Rebuild SMILES string with given modifiers or options.

    Parameters
    ----------
    smiles : str
        SMILES string to be rebuilt.
    *modifiers
        List of functions that make additional modifications
        after rebuilding process.
    **options
        For keyword-only options,
        see the `rdkit.Chem.MolToSmiles` documentation.

    Returns
    -------
    rsmiles : str
        Rebuilt SMILES string.

    Examples
    --------
    >>> rsmiles('CC(=O)O', allBondsExplicit=True, allHsExplicit=True)
    '[CH3]-[C](=[O])-[OH]'
    >>> rsmiles('C[C@@H](C(=O)O)O', canonical=False, isomericSmiles=False)
    'CC(C(=O)O)O'
    >>> rsmiles('c1ccccc1', kekuleSmiles=True)
    'C1:C:C:C:C:C:1'
    >>> rsmiles('CC(=O)C', rootedAtAtom=2)
    'O=C(C)C'
    >>> from re import sub
    >>> def square_brackets_remover(smiles):
    ...     return sub('[\[\]]', '', smiles)
    ...
    >>> rsmiles('CCO', square_brackets_remover, allHsExplicit=True)
    'CH3CH2OH'

    Call the function without keyword arguments to get a canonical SMILES:

    >>> rsmiles('CC(=O)C')
    'CC(C)=O'

    """
    rsmiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), **options)
    if modifiers:
        for modifier in modifiers:
            rsmiles = modifier(rsmiles)
    return rsmiles


def insert_ring_tags(smiles, tags='<>', opening_before_symbol=True,
                     reuse_closed=True):
    r"""
    Replace numbers in SMILES string with opening and closure tags.

    Parameters
    ----------
    smiles : str
        SMILES string.
    tags : str or iterable, optional
        A pair of opening and closing tags.
    opening_before_symbol : bool, optional
        If this is set to True, an opening tag is inserted
        before the symbol of the first atom in a ring.
    reuse_closed : bool, optional
        If this is set to `True`, a ring closure digit is reused.
        `False` setting is not recommended thus it can be a source of errors.

    Returns
    -------
    modified_smiles : str
        Modified SMILES string.

    Examples
    --------
    >>> insert_ring_tags('c1ccccc1-c1ccccc1')
    '<cccccc>-<cccccc>'
    >>> insert_ring_tags('c1ccccc1-c1ccccc1', opening_before_symbol=False)
    'c<ccccc>-c<ccccc>'
    >>> insert_ring_tags('c1ccccc1-c1ccccc1', reuse_closed=False)
    '<cccccc>-c>ccccc>'

    """
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
                if reuse_closed:
                    numbers.remove(element)
        else:
            modified_smiles += element
    if opening_before_symbol:
        modified_smiles = sub('([A-Z]|[cnops]|as|se|[A-Z][a-z])<', r'<\1',
                              modified_smiles)
    return modified_smiles


if __name__ == '__main__':
    import doctest
    print(__doc__, end='')
    test = doctest.testmod()
    num_errors = test[0]
    if '-v' not in argv and num_errors == 0:
        print('Try `xrdkit.py -v` to run doctest and see docstring examples.')

