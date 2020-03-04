def is_char(char, astype=int):

    def is_this_char(token):

        return astype(token == char)

    return is_this_char


def get_branch_len(token):

    return len(findall('\)', token))


def get_ring_int(token):

    return int(sub('%?(\d+)', r'\1', token)) if token[-1].isdigit() else 0

