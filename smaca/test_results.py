import filecmp
import unittest
from pathlib import Path

import numpy as np
import pysam
from smaca import constants as data
from smaca import utils
from smaca.sma import Bam, SmaCalculator

DATA_PATH = Path(__file__).parent.joinpath("data")
BAM_LIST = [
    DATA_PATH.joinpath("151002_7001448_0359_AC7F6GANXX_"
                       "Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008"
                       ".posiSrt.markDup.includedRegions.bam"),
    DATA_PATH.joinpath("151002_7001448_0359_AC7F6GANXX_"
                       "Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008"
                       ".posiSrt.markDup.includedRegions.bam"),
    DATA_PATH.joinpath("HG004.mate_pair.sorted.includedRegions.bam")
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

        self.assertEqual(b.get_genomic_range(data.GENES['ACAD9']),
                         ('3', 128598332, 128634910))
        self.assertEqual(b.get_consensus_sequence(data.SMN1_POS['SMN1_a']),
                         'G [[0], [0], [178], [0]]')
        self.assertEqual(b.get_consensus_sequence(data.SMN2_POS['SMN2_b_e7']),
                         'T [[0], [0], [0], [218]]')
        self.assertEqual(b.get_consensus_sequence(data.DUP_MARK['g.27134T>G']),
                         'T [[0], [0], [0], [106]]')
        self.assertEqual(
            b.get_consensus_sequence(data.DUP_MARK['g.27706_27707delAT']),
            'AT [[13, 0], [0, 0], [0, 0], [0, 13]]')
        self.assertEqual(b.get_cov_ranges(data.GENES)[0], 40.88517688227896)

    def test_sma_stats(self):
        s = SmaCalculator(BAM_LIST)

        np.testing.assert_array_almost_equal(
            s.pi_ij[0], [0.58264685, 0.5296146, 0.47764956])
        self.assertAlmostEqual(s.zmean_k[0], 1.5592104677167276)
        np.testing.assert_array_almost_equal(
            s.std_i, [21.18266345, 17.4505619, 0.99442209])
        self.assertAlmostEqual(s.std_k[0], 11.03630077)
        self.assertEqual(s.dup_id[0][0], b'T [[0], [0], [0], [106]]')

    def test_get_chr_prefix(self):
        sam_file = pysam.AlignmentFile(BAM_LIST[0], "rb")
        sam_file_chr = pysam.AlignmentFile(BAM_LIST[2], "rb")

        self.assertEqual(utils.get_chr_prefix(sam_file), '')
        self.assertEqual(utils.get_chr_prefix(sam_file_chr), 'chr')

    def test_cli(self):
        parallel_out_fpath = DATA_PATH.joinpath("output_par.csv")
        res = SmaCalculator(BAM_LIST, n_jobs=-1)
        res.write_stats(parallel_out_fpath)

        serial_out_fpath = DATA_PATH.joinpath("output_ser.csv")
        res = SmaCalculator(BAM_LIST, n_jobs=1)
        res.write_stats(serial_out_fpath)

        self.assertTrue(parallel_out_fpath.exists())
        self.assertTrue(serial_out_fpath.exists())
        self.assertTrue(filecmp.cmp(parallel_out_fpath, serial_out_fpath))


if __name__ == '__main__':
    unittest.main()
