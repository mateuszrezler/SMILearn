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


def has_no_symbol(symbols, astype=float):
    if isinstance(symbols, str):
        symbols = [symbols]

    def _has_no_symbol(mol, index):
        atom_symbol = mol.GetAtomWithIdx(index).GetSymbol()
        return astype(atom_symbol not in symbols)

    return _has_no_symbol


def has_symbol(symbol, astype=float):
    if isinstance(symbols, str):
        symbols = [symbols]

    def _has_symbol(mol, index):
        atom_symbol = mol.GetAtomWithIdx(index).GetSymbol()
        return astype(atom_symbol in symbols)

    return _has_symbol


#def has_symbol(symbol, astype=float):

#    def _has_symbol(mol, index):
#        atom_symbol = mol.GetAtomWithIdx(index).GetSymbol()
#        return astype(atom_symbol == symbol)

#    return _has_symbol


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


def get_group(maxval=1, astype=float,
              group_dict={
                  (1, 3, 11, 19, 37, 55): 1,
                  (4, 12, 20, 38, 56): 2,
                  (21, 39, 57): 3,
                  (22, 40, 72): 4,
                  (23, 41, 73): 5,
                  (24, 42, 74): 6,
                  (25, 43, 75): 7,
                  (26, 44, 76): 8,
                  (27, 45, 77): 9,
                  (28, 46, 78): 10,
                  (29, 47, 79): 11,
                  (30, 48, 80): 12,
                  (5, 13, 31, 49, 81): 13,
                  (6, 14, 32, 50, 82): 14,
                  (7, 15, 33, 51, 83): 15,
                  (8, 16, 34, 52): 16,
                  (9, 17, 35, 53): 17,
                  (2, 10, 18, 36, 54): 18
              }):

    def _get_group(mol, index):
        atomic_num = get_atomic_num()(mol, index)
        for atomic_num_tuple in group_dict:
            if atomic_num in atomic_num_tuple:
                group = group_dict[atomic_num_tuple]
                return astype(group / maxval)
        return astype(0)

    return _get_group


def get_period(maxval=1, astype=float,
               period_dict={
                   (1, 2): 1,
                   range(3, 11): 2,
                   range(11, 19): 3,
                   range(19, 37): 4,
                   range(37, 55): 5,
                   range(55, 84): 6
               }):

    def _get_period(mol, index):
        atomic_num = get_atomic_num()(mol, index)
        for atomic_num_tuple in period_dict:
            if atomic_num in atomic_num_tuple:
                period = period_dict[atomic_num_tuple]
                return astype(period / maxval)
        return astype(0)

    return _get_period


def is_metal(astype=float, atomic_num_tuple=(1, 2, 53, 54, *range(5, 11),
                                             *range(14, 19), *range(32, 37))):

    def _is_metal(mol, index):
        atomic_num = get_atomic_num()(mol, index)
        return astype(atomic_num not in atomic_num_tuple)

    return _is_metal


def is_metalloid(astype=float, atomic_num_tuple=(5, 14, 32, 33, 51, 52)):

    def _is_metalloid(mol, index):
        atomic_num = get_atomic_num()(mol, index)
        return astype(atomic_num in atomic_num_tuple)

    return _is_metalloid

