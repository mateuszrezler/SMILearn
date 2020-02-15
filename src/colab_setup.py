r"""Module setting up Google Colab Notebook when imported.
Detailed settings can be found and modified in `setup.sh` script.
"""
if __name__ == '__main__':
    print(__doc__)
else:
    import os
    import subprocess
    import sys
    os.chmod('src/colab_setup.sh', 0o700) 
    subprocess.call('src/colab_setup.sh')
    sys.path.append('/usr/local/lib/python3.7/site-packages')

