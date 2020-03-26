from setuptools import find_packages, setup

setup(
    name='smaca',
    version='1.0',
    packages=find_packages(),
    url='git@gitlab.cbra.com:dlopez/SMA_test.git',
    license='GNU General Public License v.3.0',
    author='Daniel López López, Carlos Loucera',
    author_email=
    'daniel.lopez.lopez@juntadeandalucia.es, carlos.loucera@juntadeandalucia.es',
    long_description=
    'A python module for detecting spinal muscular atrophy carriers',
    python_requires='>=3.6',
    install_requires=['click', 'cython', 'numpy', 'pysam', 'joblib'],
    include_package_data=True,
    entry_points={'console_scripts': ['smaca = smaca.cli:main']})
