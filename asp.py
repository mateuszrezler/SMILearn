r"""ASP - append site packages
Module appending site packages to `sys.path`.
"""
import sys
from re import search


class PatternError(LookupError):
    r"""An error raised when `/.*/lib/python\d\.\d` pattern
    is not found in `sys.path`.
    """
    def __str__(self):
        """String representation."""
        return 'Pattern: `/.*/lib/python\d\.\d` not found in `sys.path`.'


def main():
    r"""Main function."""
    new_paths = []
    for location in sys.path:
        found = search('/.*/lib/python\d\.\d', location)
        if found:
            new_paths.append(found.group(0) + '/site-packages')
    new_paths_set = set(new_paths)
    sys.path.extend(list(new_paths_set))
    return new_paths_set


if __name__ == '__main__':
    print(__doc__)
    nps = main()
    if nps:
        out_str = '\n'.join(nps)
        print('Generated new site packages path(s):', out_str, sep='\n')
    else:
        raise PatternError
else:
    main()

