def atomf(function, maxval=None, shift=None, eq=None, notin=None, astype=int):

    def _atom_function(mol, index):

        out_function = eval(f'mol.GetAtomWithIdx(index).{function}()')
        if shift:
            out_function += shift
        if maxval:
            out_function /= maxval
        if eq:
            return astype(out_function == eq)
        if notin:
            return astype(out_function not in notin)
        return astype(out_function)

    return _atom_function


def has_chiral_tag(chiral_tag, astype=float):

    def _has_chiral_tag(mol, index):
        atom_chiral_tag = mol.GetAtomWithIdx(index).GetChiralTag()
        return astype(str(atom_chiral_tag) == chiral_tag)

    return _has_chiral_tag


def has_hybridization(hybridization, astype=float):

    def _has_hybridization(mol, index):
        atom_hybridization = mol.GetAtomWithIdx(index).GetHybridization()
        return astype(str(atom_hybridization) == hybridization)

    return _has_hybridization


def has_no_symbols(symbols, astype=float):

    def _has_no_symbols(mol, index):
        atom_symbol = mol.GetAtomWithIdx(index).GetSymbol()
        return astype(atom_symbol not in symbols)

    return _has_no_symbols


def has_symbol(symbol, astype=float):

    def _has_symbol(mol, index):
        atom_symbol = mol.GetAtomWithIdx(index).GetSymbol()
        return astype(atom_symbol == symbol)

    return _has_symbol


def is_aromatic(astype=float):

    def _is_aromatic(mol, index):
        return astype(mol.GetAtomWithIdx(index).GetIsAromatic())

    return _is_aromatic


def is_in_ring(astype=float):

    def _is_in_ring(mol, index):
        return astype(mol.GetAtomWithIdx(index).IsInRing())

    return _is_in_ring


def get_atomic_num(maxval=1, astype=float):

    def _get_atomic_num(mol, index):
        atomic_num = mol.GetAtomWithIdx(index).GetAtomicNum()
        return astype(atomic_num / maxval)

    return _get_atomic_num


def get_charge(maxval=1, shift=0, astype=float):

    def _get_charge(mol, index):
        charge = mol.GetAtomWithIdx(index).GetFormalCharge()
        return astype((charge + shift) / maxval)

    return _get_charge


def get_degree(maxval=1, astype=float):

    def _get_degree(mol, index):
        degree = mol.GetAtomWithIdx(index).GetTotalDegree()
        return astype(degree / maxval)

    return _get_degree


def get_num_hs(maxval=1, astype=float):

    def _get_num_hs(mol, index):
        num_hs = mol.GetAtomWithIdx(index).GetTotalNumHs()
        return astype(num_hs / maxval)

    return _get_num_hs


def get_valence(maxval=1, astype=float):

    def _get_valence(mol, index):
        valence = mol.GetAtomWithIdx(index).GetTotalValence()
        return astype(valence / maxval)

    return _get_valence

