from setuptools import setup
from glob import glob
import os


setup(
    name = 'fshd',
    version = '1.2.4',
    description = 'A Python wrapper for FreeSurfaceHydrodynamics',
    url = 'https://github.com/hamilton8415/FreeSurfaceHydrodynamics',
    author = 'Michael Anderson',
    author_email = 'anderson@mbari.org',
    license = 'Apache 2.0',
    zip_safe = True,
    packages=[''],
    package_dir={'': 'build'},
    package_data={'': [os.path.basename(glob(os.path.join('build', 'fshd.*.so'))[0])]},
)
