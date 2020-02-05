from numpy import array


def get_struct_features(component, max_branch=56, max_ring=48,
                        bond_dict={
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

    if component[0] == ')' and len(component) > 1:

        return array((0, 0, 0, len(component)/max_branch, 0))

    elif component[0].isdigit() or component[0] == '%':

        return array((0, 0, 0, 0, int(re.sub('%', '', component))/max_ring))

    return array((bond_dict[component]))

