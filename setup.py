from setuptools import setup, find_packages
import pathlib


# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.rst").read_text()

setup(
    name='smaca',
    version='1.0',
    packages=find_packages(),
    url='git@gitlab.cbra.com:dlopez/SMA_test.git',
    license='GNU General Public License v.3.0',
    author='Daniel López López, Carlos Loucera',
    author_email='daniel.lopez.lopez@juntadeandalucia.es, carlos.loucera@juntadeandalucia.es',
    long_description=README,
    long_description_content_type="text/x-rst",
    python_requires='>=3.6',
    install_requires=[
    	'click',
        'cython',
	    'numpy',
	    'pysam',
	    'joblib'
	],
    entry_points={
        'console_scripts': ['smaca = smaca.cli:main']
    }
)
