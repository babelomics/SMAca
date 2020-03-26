# coding: utf-8
"""
author: Daniel López
email: daniel.lopez.lopez@juntadeandalucia.es

author: Carlos Loucera
email: carlos.loucera@juntadeandalucia.es

SMA carrier Test utilities module.
"""

import numpy as np
from smaca.constants import SMN


def get_total_depth(sam_file, chrom, start, end):
    """
    Calculate some coverage metrics for a position or a genomic range in a
    given bam file. Genomic ranges are 0-based stop-excluded

    :param sam_file: SAM file already opened for reading
    :param chrom: contig name
    :param start: start position (included) 0-based
    :param end: end position (excluded)
    :return: 4 arrays of counts for each base at each position
      ([[A1,A2,...An],[C1,C2,...,Cn],[G1,G2,...,Gn],[T1,T2,...,Tn]]) and the
      mean coverage along the genomic range
    """

    cov = np.array(sam_file.count_coverage(
        chrom, start, end))  # four arrays of the same length in order A C G T

    return cov, cov.sum() / cov.shape[1]


def get_chr_prefix(sam_file):
    """
    Check if chromosomes names contains "chr".
    :param sam_file: SAM file already opened for reading
    :return: chromosome prefix (either 'chr' or '')
    """
    prefix = ""

    for r in sam_file.fetch():
        if not r.is_unmapped:
            if "chr" in r.reference_name:
                prefix = "chr"
            else:
                prefix = ""

            break

    # safety check
    try:
        c, start, stop = SMN["SMN1"]
        c = f"{prefix}{c}"
        sam_file.count(c, start, stop)

    except ValueError:
        print(
            "The reference genome used is not supported. Supported reference genomes are hg19 and hs37d5"
        )
        raise SystemExit

    return prefix
