# coding: utf-8
"""
author: Daniel López
email: daniel.lopez.lopez@juntadeandalucia.es

author: Carlos Loucera
email: carlos.loucera@juntadeandalucia.es

SMA carrier Test samtools constants file. Genomic ranges are 0-based, stop-excluded
"""

SMN = {"SMN1": [5, 70220767, 70249769], "SMN2": [5, 69345349, 69374349]}

GENES = {
    "ACAD9": [3, 128598332, 128634910],
    "ATR": [3, 142168076, 142297668],
    "CYP11B1": [8, 143953771, 143961262],
    "EDNRB": [13, 78469615, 78493903],
    "FASTKD2": [2, 207630080, 207660913],
    "FOXN1": [17, 26833260, 26865914],
    "HEXB": [5, 73935847, 74018472],
    "IQCB1": [3, 121488609, 121553926],
    "ITGA6": [2, 173292081, 173371181],
    "IVD": [15, 40697685, 40728146],
    "LMNA": [1, 156052363, 156109880],
    "LRPPRC": [2, 44113362, 44223144],
    "NTRK1": [1, 156785431, 156851642],
    "PTEN": [10, 89622869, 89731687],
    "RAB3GAP1": [2, 135809834, 135933964],
    "RAPSN": [11, 47459307, 47470730],
    "SIL1": [5, 138282408, 138629246],
    "SLC22A5": [5, 131705400, 131731306],
    "SLC35D1": [1, 67465014, 67520080],
    "STIM1": [11, 3875756, 4114440]
}

SMN1_POS = {
    "SMN1_a": [5, 70247723, 70247724],
    "SMN1_b_e7": [5, 70247772, 70247773],
    "SMN1_c": [5, 70247920, 70247921]
}

SMN2_POS = {
    "SMN2_a": [5, 69372303, 69372304],
    "SMN2_b_e7": [5, 69372352, 69372353],
    "SMN2_c": [5, 69372500, 69372501]
}

DUP_MARK = {
    "g.27134T>G": [5, 70247900, 70247901],
    "g.27706_27707delAT": [5, 70248472, 70248474]
}

HEADER_FILE = \
 "id," \
 "Pi_a," \
 "Pi_b," \
 "Pi_c," \
 "cov_SMN1_a," \
 "cov_SMN1_b," \
 "cov_SMN1_c," \
 "cov_SMN2_a," \
 "cov_SMN2_b," \
 "cov_SMN2_c," \
 "avg_cov_SMN1," \
 "avg_cov_SMN2," \
 "scale_factor," \
 "std_control," \
 "g.27134T>G," \
 "g.27706_27707delAT," \
 "avg_cov_ACAD9," \
 "avg_cov_ATR," \
 "avg_cov_CYP11B1," \
 "avg_cov_EDNRB," \
 "avg_cov_FASTKD2," \
 "avg_cov_FOXN1," \
 "avg_cov_HEXB," \
 "avg_cov_IQCB1," \
 "avg_cov_ITGA6," \
 "avg_cov_IVD," \
 "avg_cov_LMNA," \
 "avg_cov_LRPPRC," \
 "avg_cov_NTRK1," \
 "avg_cov_PTEN," \
 "avg_cov_RAB3GAP1," \
 "avg_cov_RAPSN," \
 "avg_cov_SIL1," \
 "avg_cov_SLC22A5," \
 "avg_cov_SLC35D1," \
 "avg_cov_STIM1"
