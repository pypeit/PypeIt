"""
Prints the version
"""
from pypeit.scripts import scriptbase

class Version(scriptbase.ScriptBase):

    @staticmethod
    def main(args):
        import pypeit
        print('The version of PypeIt is: {:s}'.format(pypeit.__version__))

