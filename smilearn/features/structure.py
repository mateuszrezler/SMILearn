from re import findall, sub


def get_branch_len(maxval=1, astype=float, warn=False):

    def _get_branch_len(token):
        branch_len = len(findall('\)', token))
        if branch_len/maxval > 1 and warn:
            print(f'Too long branch found ({branch_len}).')
        return astype(branch_len/maxval)

    return _get_branch_len


def get_ring_int(maxval=1, astype=float, warn=False):

    def _get_ring_int(token):
        if token[-1].isdigit():
            ring_int = int(sub('%?(\d+)', r'\1', token))
            if ring_int/maxval > 1 and warn:
                print(f'Too big ring found ({ring_int}).')
            return astype(ring_int/maxval)
        else:
            return astype(0)

    return _get_ring_int


def is_char(char, astype=float):

    def _is_char(token):
        return astype(token == char)

    return _is_char

