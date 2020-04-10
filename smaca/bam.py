# coding: utf-8
"""
author: Daniel López
email: daniel.lopez.lopez@juntadeandalucia.es

author: Carlos Loucera
email: carlos.loucera@juntadeandalucia.es

SMA carrier Test main class.
"""

import numpy as np
import pysam

from smaca.utils import get_chr_prefix, get_total_depth

BASES = np.array(["A", "C", "G", "T"])


class Bam:
    """
    This class contains statistics of coverage calculated on different genomic
    ranges
    """
    def __init__(self, bam_file):
        self.bam_file = bam_file
        self.samfile = pysam.AlignmentFile(bam_file, "rb")
        #TODO:guess reference genome that was used for alignment
        self.chrPrefix = get_chr_prefix(self.samfile)

    def __del__(self):
        self.samfile.close()

    def get_cov_ranges(self, ranges):
        """
        Calculates coverage in a set of ranges

        :param ranges: Array of ranges where coverage should be calculated
        :return: Array of coverages
        """

        mean_cov = np.zeros(len(ranges))

        for i in range(len(ranges)):
            c, start, stop = self.get_genomic_range(list(ranges.values())[i])
            _, mean_cov[i] = get_total_depth(self.samfile, c, start, stop)

        return mean_cov

    def get_consensus_sequence(self, genomic_range):
        """
        calculates consensus sequence from depth of coverage
        :param genomic_range: genomic range in format chr, start, end
        :return: consensus sequence and coverages for each position
        """
        chrom, start, stop = self.get_genomic_range(genomic_range)
        cov, _ = get_total_depth(self.samfile, chrom, start, stop)

        return f'{"".join(BASES[cov.argmax(axis=0).tolist()])} {cov.tolist()}'

    def get_genomic_range(self, genomic_range):
        """
        Add "chr" to the genomic range if necesary
        :param genomic_range:
        """
        c, start, stop = genomic_range

        return f"{self.chrPrefix}{c}", start, stop
