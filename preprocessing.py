from numpy import array
from re import sub


class FeatureCreator(object):
    r"""
    Creator of feature arrays.

    Methods
    -------
    get_chem_struct_features
    """
    def get_chem_struct_features(self, component, max_branch_len=1,
                                 max_ring_int=9,
                                 struct_features_dict={
                                     '.': (0.5/3, 0, 0, 0, 0),
                                     '-': (1/3, 0, 0, 0, 0),
                                     '/': (1/3, 1, 0, 0, 0),
                                     '\\': (1/3, 0, 1, 0, 0),
                                     ':': (1.5/3, 0, 0, 0, 0),
                                     '=': (2/3, 0, 0, 0, 0),
                                     '#': (3/3, 0, 0, 0, 0),
                                     '(': (0, 0, 0, 1/2, 0),
                                     ')': (0, 0, 0, 2/2, 0)
                                 }):
        r"""
        Translate symbols of structural elements into feature vectors.

        Parameters
        ----------
        component : str
            Input string.
        max_branch_len : int, optional
            Maximum permissible number of closing parentheses describing
            branch length in DeepSMILES notation.
        max_ring_int : int, optional
            Maximum permissible number describing ring connection.
        struct_features_dict : dict, optional
            Dictionary with str keys and tuple values.

        Returns
        -------
        numpy.ndarray
            Structure features (1D array).

        Examples
        --------
        >>> fc = FeatureCreator()
        >>> fc.get_chem_struct_features(':')
        array([0.5, 0. , 0. , 0. , 0. ])
        >>> fc.get_chem_struct_features('))', max_branch_len=10)
        array([0. , 0. , 0. , 0.2, 0. ])
        >>> fc.get_chem_struct_features('%20', max_ring_int=50)
        array([0. , 0. , 0. , 0. , 0.4])
        """
        if component[0] == ')' and len(component) > 1:

            scaled_to_max = len(component)/max_branch_len

            return array((0, 0, 0, scaled_to_max, 0))

        elif component[0].isdigit() or component[0] == '%':

            scaled_to_max = int(sub('%', '', component))/max_ring_int

            return array((0, 0, 0, 0, scaled_to_max))

        return array((struct_features_dict[component]))


if __name__ == '__main__':
    import doctest
    doctest.testmod()

