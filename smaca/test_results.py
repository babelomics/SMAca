import filecmp
import unittest
from pathlib import Path

import numpy as np
import pysam

from smaca import constants as C
from smaca import utils
from smaca.sma import Bam, SmaCalculator

DATA_PATH = Path(__file__).parent.joinpath("data")

BAM_hg19 = DATA_PATH.joinpath(
    "HG007.hs37d5.100x.includedRegions.downSampled.bam")
BAM_hg38 = DATA_PATH.joinpath(
    "HG007.GRCh38_full_plus_hs38d1_analysis_"
    "set_minus_alts.100x.includedRegions.downSampled.bam")

BAM_LIST = [
    DATA_PATH.joinpath("151002_7001448_0359_AC7F6GANXX_"
                       "Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008"
                       ".posiSrt.markDup.includedRegions.bam"),
    BAM_hg19
]




class MyTestCase(unittest.TestCase):
    def test_get_total_depth(self):
        sam_file = pysam.AlignmentFile(BAM_LIST[0], "rb")
        c, d = utils.get_total_depth(sam_file, '3', 128598332, 128634910)

        self.assertAlmostEqual(d, 40.88517688227896)

        c, d = utils.get_total_depth(sam_file, '5', 70247900, 70247901)
        np.testing.assert_array_equal(c, [[0], [0], [0], [106]])
        self.assertEqual(d, 106)

    def test_bamclass(self):
        b = Bam(BAM_LIST[0])

        self.assertEqual(
            b.get_genomic_range(C.POSITIONS[C.REF_HG19]["GENES"]['ACAD9']),
            ('3', 128598332, 128634910))
        self.assertEqual(
            b.get_consensus_sequence(
                C.POSITIONS[C.REF_HG19]["SMN1_POS"]['SMN1_a']),
            'G [[0], [0], [178], [0]]')
        self.assertEqual(
            b.get_consensus_sequence(
                C.POSITIONS[C.REF_HG19]["SMN2_POS"]['SMN2_b_e7']),
            'T [[0], [0], [0], [218]]')
        self.assertEqual(
            b.get_consensus_sequence(
                C.POSITIONS[C.REF_HG19]["DUP_MARK"]['g.27134T>G']),
            'T [[0], [0], [0], [106]]')
        self.assertEqual(
            b.get_consensus_sequence(
                C.POSITIONS[C.REF_HG19]["DUP_MARK"]['g.27706_27707delAT']),
            'AT [[13, 0], [0, 0], [0, 0], [0, 13]]')
        self.assertEqual(
            b.get_cov_ranges(C.POSITIONS[C.REF_HG19]["GENES"])[0],
            40.88517688227896)

    def test_sma_stats(self):
        s = SmaCalculator(BAM_LIST, ref=C.REF_HG19)

        np.testing.assert_array_almost_equal(
            s.pi_ij[0], [0.8027337510925755, 0.7296692915157783, 0.6580751624125392])
        self.assertAlmostEqual(s.zmean_k[0], 0.6135924580475719)
        np.testing.assert_array_almost_equal(
            s.std_i, [21.18266345, 5.07109509])
        self.assertAlmostEqual(s.std_k[0], 15.516854393351194)
        self.assertEqual(s.dup_id[0][0], b'T [[0], [0], [0], [106]]')

    def test_get_chr_prefix(self):
        sam_file = pysam.AlignmentFile(BAM_hg19, "rb")
        sam_file_chr = pysam.AlignmentFile(BAM_hg38, "rb")

        self.assertEqual(utils.get_chr_prefix(sam_file), '')
        self.assertEqual(utils.get_chr_prefix(sam_file_chr), 'chr')

    def test_cli(self):
        parallel_out_fpath = DATA_PATH.joinpath("output_par.csv")
        res = SmaCalculator(BAM_LIST, ref=C.REF_HG19, n_jobs=-1)
        res.write_stats(parallel_out_fpath)

        serial_out_fpath = DATA_PATH.joinpath("output_ser.csv")
        res = SmaCalculator(BAM_LIST, ref=C.REF_HG19, n_jobs=1)
        res.write_stats(serial_out_fpath)

        self.assertTrue(parallel_out_fpath.exists())
        self.assertTrue(serial_out_fpath.exists())
        self.assertTrue(filecmp.cmp(parallel_out_fpath, serial_out_fpath))

    def test_hg19_hg38_coverages(self):
        bam_hg19 = Bam(BAM_hg19)
        bam_hg38 = Bam(BAM_hg38)

        for ranges in C.POSITIONS[C.REF_HG19]:
            c_hg19 = bam_hg19.get_cov_ranges(C.POSITIONS[C.REF_HG19][ranges])
            c_hg38 = bam_hg38.get_cov_ranges(C.POSITIONS[C.REF_HG38][ranges])

            np.testing.assert_array_almost_equal(
                c_hg19,
                c_hg38,
                decimal=0,
                err_msg=";".join(C.POSITIONS[C.REF_HG19][ranges]))


if __name__ == '__main__':
    unittest.main()
