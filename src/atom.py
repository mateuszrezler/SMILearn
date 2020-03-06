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

