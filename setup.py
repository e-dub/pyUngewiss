from distutils.core import setup
import os
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / 'README.md').read_text()


def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths


if __name__ == '__main__':
    extra_files = package_files('pyUngewiss')
    setup(
        name='pyUngewiss',
        version='1.1',
        author='E. J. Wehrle',
        author_email='Erich.Wehrle@unibz.it',
        packages=['pyUngewiss'],
        package_data={'': extra_files},
        description='Python library for UNcertainty analysis in liGhtwEight dsiGn with IntervalS and fuzzy numberS',
        license='GPL3',
        long_description=long_description,
        long_description_content_type='text/markdown',
    )
