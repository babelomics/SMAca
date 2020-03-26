# coding: utf-8
"""
author: Daniel López
email: daniel.lopez.lopez@juntadeandalucia.es

author: Carlos Loucera
email: carlos.loucera@juntadeandalucia.es

SMA carrier Test main class.
"""

import tempfile

import click
import numpy as np
import pysam
from joblib import Parallel, delayed
from smaca import constants as data
from smaca.utils import get_chr_prefix, get_total_depth

BASES = np.array(["A", "C", "G", "T"])


class SmaCalculator:
    """This class implements some statistics functions to calculate SMN1:SMN2
    proportion in a set of BAMs.
    """
    def __init__(self, bam_list, n_jobs=1):
        """

        :param bam_list: list of bam files (path)
        """
        self.bam_list = np.array(bam_list)
        self.n_bam = len(self.bam_list)
        # number of reads that align to SMN1 at position x
        self.D1_ij = np.zeros((self.n_bam, len(data.SMN1_POS)))
        # number of reads that align to SMN2 at position x
        self.D2_ij = np.zeros((self.n_bam, len(data.SMN2_POS)))
        # total number of reads aligned to the SMN1 region at position j
        # and the analogous SMN2 region
        self.r_ij = np.zeros((self.n_bam, len(data.SMN2_POS)))
        # average coverage for HK gene k
        self.H_ik = np.zeros((self.n_bam, len(data.GENES)))
        # average coverage for the SMNx gene region
        self.c_ix = np.zeros((self.n_bam, len(data.SMN)))
        # scaled coverage for SMN1 and SMN2 for HK gene k
        self.z_ik = np.zeros((self.n_bam, len(data.GENES)))
        # scaled proportion of SMN reads that align to SMN1 in exon 7
        self.pi_ij = np.zeros((self.n_bam, len(data.SMN1_POS)))
        # averaged scaled coverages for SMN1 and SMN2 in the N subjects
        self.zmean_k = np.zeros((len(data.GENES), ))
        # weighted average of the coverage of SMN1 to our K housekeeping genes
        self.theta_i = np.zeros((self.n_bam, ))
        # std of coverage in housekeeping genes in sample i
        self.std_i = np.zeros((self.n_bam, ))
        # std of coverage in housekeeping gene k for each sample
        self.std_k = np.zeros((len(data.GENES), ))
        # get consensus sequence on dup. markers
        self.dup_id = np.empty((self.n_bam, len(data.DUP_MARK)), dtype="S100")

        bam_list = np.array(bam_list)
        n_bam = len(bam_list)

        if n_jobs == 1:
            with click.progressbar(length=self.n_bam,
                                   label='processing BAM files') as bar:
                for i in bar:
                    self.compute(bam_list[i], i, self.D1_ij, self.D2_ij,
                                 self.H_ik, self.c_ix, self.dup_id)
        else:
            from pathlib import Path

            tmp_dir = tempfile.mkdtemp()

            D1_ij_fname_memmap = Path(tmp_dir).joinpath("D1_ij_memmap")
            D1_ij_memmap = np.memmap(D1_ij_fname_memmap.as_posix(),
                                     dtype=np.float,
                                     shape=(n_bam, len(data.SMN1_POS)),
                                     mode='w+')

            D2_ij_fname_memmap = Path(tmp_dir).joinpath("D2_ij_memmap")
            D2_ij_memmap = np.memmap(D2_ij_fname_memmap,
                                     dtype=np.float,
                                     shape=(n_bam, len(data.SMN2_POS)),
                                     mode='w+')

            H_ik_fname_memmap = Path(tmp_dir).joinpath("H_ik_memmap")
            H_ik_memmap = np.memmap(H_ik_fname_memmap,
                                    dtype=np.float,
                                    shape=(n_bam, len(data.GENES)),
                                    mode='w+')

            c_ix_fname_memmap = Path(tmp_dir).joinpath("c_ix_memmap")
            c_ix_memmap = np.memmap(c_ix_fname_memmap,
                                    dtype=np.float,
                                    shape=(n_bam, len(data.SMN)),
                                    mode='w+')

            dup_id_fname_memmap = Path(tmp_dir).joinpath("dup_id_memmap")
            dup_id_memmap = np.memmap(dup_id_fname_memmap,
                                      dtype="S100",
                                      shape=(n_bam, len(data.DUP_MARK)),
                                      mode='w+')

            Parallel(n_jobs=n_jobs)(delayed(
                self.compute)(bam_list[idx], idx, D1_ij_memmap, D2_ij_memmap,
                              H_ik_memmap, c_ix_memmap, dup_id_memmap)
                                    for idx in range(self.n_bam))

            self.D1_ij[:] = D1_ij_memmap[:]
            self.D2_ij[:] = D2_ij_memmap[:]
            self.H_ik[:] = H_ik_memmap[:]
            self.c_ix[:] = c_ix_memmap[:]
            self.dup_id[:] = dup_id_memmap[:]

        self.r_ij = self.D1_ij + self.D2_ij
        #TODO: consider using only "the bests" HK
        self.z_ik = self.c_ix.sum(axis=1).reshape((self.n_bam, 1)) / self.H_ik
        self.std_k = np.std(self.H_ik, axis=0)
        self.std_i = np.std(self.H_ik, axis=1)
        self.zmean_k = self.z_ik.sum(axis=0) / self.n_bam
        self.theta_i = (self.z_ik / self.zmean_k).sum(axis=1) / len(data.GENES)
        self.pi_ij = self.theta_i.reshape(
            (self.n_bam, 1)) * (self.D1_ij / self.r_ij)

    def write_stats(self, output_file):
        """
        Write the calculated stats to and output file in csv format
        :param output_file: output file
        """

        out_array = np.concatenate(
            (self.bam_list.reshape((self.n_bam, 1)), self.pi_ij, self.D1_ij,
             self.D2_ij, self.c_ix, self.theta_i.reshape(
                 (self.n_bam, 1)), self.std_i.reshape(
                     (self.n_bam, 1)), self.dup_id.reshape(
                         (self.n_bam, 2)), self.H_ik),
            axis=1)
        h_quant = np.quantile(self.H_ik, [0.25, 0.5, 0.75])
        footer = f"Control genes min coverage:{np.min(self.H_ik)}\n" \
                 f"Control genes quartile 1,2,3 coverage:{h_quant}\n" \
                 f"Control genes average coverage:{self.H_ik.mean()}\n" \
                 f"Control genes max coverage:{np.max(self.H_ik)}\n" \
                 f"Control genes standard deviation: {self.std_k}"

        np.savetxt(output_file,
                   out_array,
                   fmt='%s',
                   delimiter='|',
                   header=data.HEADER_FILE,
                   footer=footer)

    @staticmethod
    def compute(bam_file, i, D1_ij, D2_ij, H_ik, c_ix, dup_id):
        b = Bam(bam_file)
        D1_ij[i] = b.get_cov_ranges(data.SMN1_POS)
        D2_ij[i] = b.get_cov_ranges(data.SMN2_POS)
        H_ik[i] = b.get_cov_ranges(data.GENES)
        c_ix[i] = b.get_cov_ranges(data.SMN)
        for d in range(len(data.DUP_MARK)):
            dup_id[i][d] = b.get_consensus_sequence(
                list(data.DUP_MARK.values())[d])


class Bam:
    """
    This class contains statistics of coverage calculated on different genomic
    ranges
    """
    def __init__(self, bam_file):
        self.bam_file = bam_file
        self.samfile = pysam.AlignmentFile(bam_file, "rb")
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
