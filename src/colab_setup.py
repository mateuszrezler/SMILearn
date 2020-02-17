r"""
colab_setup

Module setting up Google Colab Notebook when imported.
Detailed settings can be found and modified in `colab_setup.sh` script.

"""
if __name__ == '__main__':
    print('IMPORT-ONLY MODULE', __doc__, end='', sep='\n')
else:
    import os
    import subprocess
    import sys
    os.chmod('src/colab_setup.sh', 0o700) 
    subprocess.call('src/colab_setup.sh')
    sys.path.append('/usr/local/lib/python3.7/site-packages')

