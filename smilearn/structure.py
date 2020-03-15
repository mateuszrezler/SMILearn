from re import findall, sub


def is_char(char, astype=float):

    def _is_char(token):
        return astype(token == char)

    return _is_char


def get_branch_len(maxval=1, astype=float):

    def _get_branch_len(token):
        branch_len = len(findall('\)', token))
        return astype(branch_len / maxval)

    return _get_branch_len


def get_ring_int(maxval=1, astype=float):

    def _get_ring_int(token):
        if token[-1].isdigit():
            ring_int = int(sub('%?(\d+)', r'\1', token))
            return astype(ring_int / maxval)
        else:
            return astype(0)

    return _get_ring_int

