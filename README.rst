================================
SMAca: SMA Carrier Analysis tool
================================

* `summary`_
* `usage`_
* `output`_
* `interpretation`_
* `installation`_
* `citation`_


summary
-------

Spinal Muscular Atrophy (SMA) is a severe neuromuscular autosomal recessive disorder affecting 1/10,000 live births. Most SMA patients present homozygous deletion of SMN1, while most SMA carriers present only a single SMN1 copy. The sequence similarity between SMN1 and SMN2, and the complexity of the SMN locus, make the estimation of the SMN1 copy-number difficult by next generation sequencing (NGS). 

SMAca is a python tool to detect putative SMA carriers and estimate the absolute SMN1 copy-number in a population. Moreover, SMAca takes advantage of the knowledge of certain variants specific to SMN1 duplication to also identify the so-called “silent carriers” (i.e. individuals with two copies of SMN1 on one chromosome, but none on the other).

This tool is developed with multithreading support to afford high performance and a focus on easy installation. This combination makes it especially attractive to be integrated into production NGS pipelines.





usage
-----

You can run SMAca by typing at the terminal:

::

  $ smaca sample1.bam sample2.bam sample3.bam 



For a large number of samples, the **ncpus** option is recommended:

::

  $ smaca --output results.batch1.csv --ncpus 24 $(cat samplelist.batch1.txt)



For additional options use:

::

  $ smaca --help




output 
------

SMAca outputs a number of statistics for each sample:

:Pi_p: scaled proportion of SMN1 reads for positions *p*.

:cov_x_p: raw coverage of gene *x* at position *p*.

:avg_cov_x: average coverage for the whole gene *x*.

:std_control: standard deviation for the average coverage of the 20 control genes.

:g.27134T>G: consensus sequence at position 27134, as well as counts for "A", "C", "G" and "T".

:g.27706_27707delAT: consensus sequence at positions 27706-27707, as well as counts for "A", "C", "G" and "T".  

:scale_factor: scale factor proportional to the total SMN1 and SMN2 copy-number.




interpretation
--------------

SMA carriers with a single SMN1 copy are expected to have **Pi_b** values under 1/3. However, complex SMN reorganizations may lead to large differences between **Pi_a**, **Pi_b** and **Pi_c**. These cases should be analyzed carefully.

The **scale_factor**, which is proportional to the absolute number of SMN1 and SMN2 copies, and **cov_x_p** can be used to estimate the absolute SMN1:SMN2 copy-number as follows:

+----------+--------------+-----------------------+
| genotype | scale_factor | cov_SMN1_p/cov_SMN2_p |
+==========+==============+=======================+
| 1:3      | 1            | 1/3                   |
+----------+--------------+-----------------------+
| 1:2      | 0.75         | 1/2                   |
+----------+--------------+-----------------------+
| 1:1      | 0.5          | 1                     |
+----------+--------------+-----------------------+

In order to detect the so-called *silent carriers* (i.e. individuals with two copies of SMN1 on one chromosome, but none on the other), the consensus sequence at the two locations should also be taken into account. Depending on the number of SMN2 copies, the expected **scale_factor** should be close to 0.75 (2:1) or 0.5 (2:0) and, in both cases, the scaled proportion of SMN1 reads **Pi_p** should be close to 1/2 in each position.




installation
------------

SMAca is available through PyP:

::

  $ pip install smaca

If you are using the conda packaging manager (e.g. miniconda or anaconda), you can install SMAca from the bioconda channel:

::

  $ conda config --add channels defaults
  $ conda config --add channels conda-forge
  $ conda config --add channels bioconda
  $ conda install smaca

Developers can clone the repository, create a conda/pip environment and install in editable mode:

::

  $ git clone git+https://www.github.com/babelomics/SMAca.git
  $ cd SMAca
  $ python -m venv smaca_venv
  $ source smaca_venv/bin/activate
  $ pip install --editable=.



citation
--------

Daniel Lopez-Lopez, Rosario Carmona, Carlos Loucera, Virginia Aquino, Josefa Salgado, Angel Alonso, Joaquín Dopazo (2020). SMAca: SMN1 copy-number and sequence variant analysis from next generation sequencing data, XXX
