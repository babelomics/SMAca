# coding: utf-8

from setuptools import find_packages, setup
from pathlib import Path

# The directory containing this file
HERE = Path(__file__).parent

# The text of the README file
README = (HERE / "README.rst").read_text()

setup(
    name='smaca',
    version='1.0',
    packages=find_packages(),
    url='https://github.com/babelomics/SMAca',
    license='GNU General Public License v.3.0',
    author='Daniel López López, Carlos Loucera',
    author_email=
    'daniel.lopez.lopez@juntadeandalucia.es, carlos.loucera@juntadeandalucia.es',
    long_description=README,
    long_description_content_type="text/x-rst",
    python_requires='>=3.6',
    install_requires=['click', 'cython', 'numpy', 'pysam', 'joblib'],
    include_package_data=True,
    entry_points={'console_scripts': ['smaca = smaca.cli:main']})
