"""
Prints the version
"""
from pypeit.scripts import scriptbase

class Version(scriptbase.ScriptBase):

    @classmethod
    def main(cls, args):
        import pypeit
        print('The version of PypeIt is: {:s}'.format(pypeit.__version__))

