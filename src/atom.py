def atomf(function, maxval=None, shift=None, astype=int):

    def _atom_function(mol, index):
        out_function = eval(f'mol.GetAtomWithIdx(index).{function}()')
        if shift:
            out_function += shift
        if maxval:
            out_function /= maxval
        return astype(out_function)

    return _atom_function

