r"""
clone

Module cloning repositories when imported.
Detailed settings can be found and modified in `clone.sh` script.

"""
if __name__ != '__main__':
    from os import chmod
    from subprocess import call
    chmod('colab/clone.sh', 0o700) 
    call('colab/clone.sh')
else:
    print('IMPORT-ONLY MODULE', __doc__, end='', sep='\n')

