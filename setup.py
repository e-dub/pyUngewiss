from distutils.core import setup
import os


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
        version='1.0',
        author='E. J. Wehrle',
        author_email='Erich.Wehrle@unibz.it',
        packages=['pyUngewiss'],
        package_data={'': extra_files},
        description='Python library for UNcertainty analysis in liGhtwEight dsiGn with IntervalS and fuzzy numberS',
        license='GPL3',
    )
