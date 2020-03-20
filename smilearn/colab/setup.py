r"""
colab_setup

Module setting up Google Colab Notebook when imported.
Detailed settings can be found and modified in `colab_setup.sh` script.

"""
if __name__ != '__main__':
    from os import chmod
    from subprocess import call
    from sys import path
    chmod('smilearn/colab/setup.sh', 0o700) 
    call('smilearn/colab/setup.sh')
    path.append('/usr/local/lib/python3.7/site-packages')
    print('Google Colab mode is on.')
else:
    print('IMPORT-ONLY MODULE', __doc__, end='', sep='\n')

