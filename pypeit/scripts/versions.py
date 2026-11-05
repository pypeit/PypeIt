"""
Print PypeIt and dependency versions.
"""
from pypeit.scripts import scriptbase


class Versions(scriptbase.ScriptBase):

    @classmethod
    def main(cls, args):
        from pypeit import io

        for package, _, package_version, _ in io.runtime_versions():
            print(f'{package}: {package_version}')
