"""Process peptides from the MHC Motif Atlas database."""

import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

project_dir = Path("/path/to/carmen-analysis")

data_dir = project_dir.joinpath("data")
data_pub_dir = data_dir.joinpath("to-be-published")
data_subs_dir = data_dir.joinpath("subsidiary-files")

# MHC Motif Atlas peptides, available from the project website (http://mhcmotifatlas.org)
motif_atlas_peps_src_file = data_dir.joinpath("data_classI_all_peptides.txt")
motif_atlas_peps_col_allele = "Allele"  # column name
motif_atlas_peps_col_peptide = "Peptide"  # column name

motif_atlas_peps_file = data_subs_dir.joinpath("motif-atlas-peptides.csv")

umap_reducer_main_file = data_pub_dir.joinpath("umap-reducer-main.pickle")

main_motif_atlas_emb_file = data_subs_dir.joinpath(
    "main-embeddings_MHCMotifA_pssms.tsv"
)
main_motif_atlas_emb_col_allele = "Allele"  # column name
main_motif_atlas_emb_col_x = "0"  # column name
main_motif_atlas_emb_col_y = "1"  # column name
main_motif_atlas_emb_col_peps = "Peptides"  # column name
main_motif_atlas_emb_list_sep = ","

peptide_len = 9

AMINO_ACIDS = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]


def are_sequence_residues_valid(sequence):
    return all(res in AMINO_ACIDS for res in sequence)


def hobohm(peptides, threshold=0.63):
    peps_n = len(peptides)
    # seq_len = len(peptides[0])
    seq_len = 9

    aa_mat = np.array([list(p) for p in peptides])

    dist_mat = np.zeros((peps_n, peps_n))
    for i in range(peps_n - 1):
        d = 1 - np.sum(aa_mat[i] == aa_mat[i + 1 :], 1) / seq_len
        dist_mat[i, i + 1 :] = d
        dist_mat[i + 1 :, i] = d

    clusters = []
    not_assigned = np.full(peps_n, True)
    for i in range(peps_n):
        if not_assigned[i]:
            c = [i]
            ok_peps = dist_mat[i, i + 1 :] < (1 - threshold)
            c_a = np.nonzero(ok_peps & not_assigned[i + 1 :])[0] + 1 + i
            c += list(c_a)
            clusters.append(c)
            if c_a.shape[0] > 0:
                not_assigned[c_a] = False

    weights = np.zeros(peps_n)
    for c in clusters:
        c_len = len(c)
        for i in c:
            weights[i] = 1 / c_len

    alpha = np.full(seq_len, len(clusters) - 1)

    return alpha, weights


blosum62_matrix = {
    "A": [
        0.2901,
        0.0446,
        0.0427,
        0.041,
        0.065,
        0.0559,
        0.0552,
        0.0783,
        0.042,
        0.0471,
        0.0445,
        0.057,
        0.0522,
        0.0338,
        0.0568,
        0.1099,
        0.073,
        0.0303,
        0.0405,
        0.07,
    ],
    "R": [
        0.031,
        0.345,
        0.0449,
        0.0299,
        0.0163,
        0.0735,
        0.0497,
        0.0229,
        0.0458,
        0.0177,
        0.0243,
        0.1071,
        0.0321,
        0.019,
        0.0258,
        0.0401,
        0.0355,
        0.0227,
        0.028,
        0.0219,
    ],
    "N": [
        0.0256,
        0.0388,
        0.3169,
        0.069,
        0.0163,
        0.0441,
        0.0405,
        0.0391,
        0.0534,
        0.0147,
        0.0142,
        0.0415,
        0.0201,
        0.0169,
        0.0233,
        0.0541,
        0.0434,
        0.0152,
        0.0218,
        0.0165,
    ],
    "D": [
        0.0297,
        0.031,
        0.0831,
        0.3974,
        0.0163,
        0.0471,
        0.0902,
        0.0337,
        0.0382,
        0.0177,
        0.0152,
        0.0415,
        0.0201,
        0.0169,
        0.031,
        0.0489,
        0.0375,
        0.0152,
        0.0187,
        0.0178,
    ],
    "C": [
        0.0216,
        0.0078,
        0.009,
        0.0075,
        0.4837,
        0.0088,
        0.0074,
        0.0108,
        0.0076,
        0.0162,
        0.0162,
        0.0086,
        0.0161,
        0.0106,
        0.0103,
        0.0175,
        0.0178,
        0.0076,
        0.0093,
        0.0192,
    ],
    "Q": [
        0.0256,
        0.0484,
        0.0337,
        0.0299,
        0.0122,
        0.2147,
        0.0645,
        0.0189,
        0.0382,
        0.0133,
        0.0162,
        0.0535,
        0.0281,
        0.0106,
        0.0207,
        0.0332,
        0.0276,
        0.0152,
        0.0218,
        0.0165,
    ],
    "E": [
        0.0405,
        0.0523,
        0.0494,
        0.0914,
        0.0163,
        0.1029,
        0.2965,
        0.0256,
        0.0534,
        0.0177,
        0.0202,
        0.0708,
        0.0281,
        0.019,
        0.0362,
        0.0524,
        0.0394,
        0.0227,
        0.028,
        0.0233,
    ],
    "G": [
        0.0783,
        0.0329,
        0.0652,
        0.0466,
        0.0325,
        0.0412,
        0.035,
        0.5101,
        0.0382,
        0.0206,
        0.0213,
        0.0432,
        0.0281,
        0.0254,
        0.0362,
        0.0663,
        0.0434,
        0.0303,
        0.0249,
        0.0247,
    ],
    "H": [
        0.0148,
        0.0233,
        0.0315,
        0.0187,
        0.0081,
        0.0294,
        0.0258,
        0.0135,
        0.355,
        0.0088,
        0.0101,
        0.0207,
        0.0161,
        0.0169,
        0.0129,
        0.0192,
        0.0138,
        0.0152,
        0.0467,
        0.0082,
    ],
    "I": [
        0.0432,
        0.0233,
        0.0225,
        0.0224,
        0.0447,
        0.0265,
        0.0221,
        0.0189,
        0.0229,
        0.271,
        0.1154,
        0.0276,
        0.1004,
        0.0634,
        0.0258,
        0.0297,
        0.0533,
        0.0303,
        0.0436,
        0.1646,
    ],
    "L": [
        0.0594,
        0.0465,
        0.0315,
        0.028,
        0.065,
        0.0471,
        0.0368,
        0.0283,
        0.0382,
        0.1679,
        0.3755,
        0.0432,
        0.1968,
        0.1142,
        0.0362,
        0.0419,
        0.0651,
        0.053,
        0.0685,
        0.1303,
    ],
    "K": [
        0.0445,
        0.1202,
        0.0539,
        0.0448,
        0.0203,
        0.0912,
        0.0755,
        0.0337,
        0.0458,
        0.0236,
        0.0253,
        0.2781,
        0.0361,
        0.019,
        0.0413,
        0.0541,
        0.0454,
        0.0227,
        0.0312,
        0.0261,
    ],
    "M": [
        0.0175,
        0.0155,
        0.0112,
        0.0093,
        0.0163,
        0.0206,
        0.0129,
        0.0094,
        0.0153,
        0.0368,
        0.0496,
        0.0155,
        0.1606,
        0.0254,
        0.0103,
        0.0157,
        0.0197,
        0.0152,
        0.0187,
        0.0316,
    ],
    "F": [
        0.0216,
        0.0174,
        0.018,
        0.0149,
        0.0203,
        0.0147,
        0.0166,
        0.0162,
        0.0305,
        0.0442,
        0.0547,
        0.0155,
        0.0482,
        0.3869,
        0.0129,
        0.0209,
        0.0237,
        0.0606,
        0.1308,
        0.0357,
    ],
    "P": [
        0.0297,
        0.0194,
        0.0202,
        0.0224,
        0.0163,
        0.0235,
        0.0258,
        0.0189,
        0.0191,
        0.0147,
        0.0142,
        0.0276,
        0.0161,
        0.0106,
        0.4935,
        0.0297,
        0.0276,
        0.0076,
        0.0156,
        0.0165,
    ],
    "S": [
        0.085,
        0.0446,
        0.0697,
        0.0522,
        0.0407,
        0.0559,
        0.0552,
        0.0513,
        0.042,
        0.025,
        0.0243,
        0.0535,
        0.0361,
        0.0254,
        0.0439,
        0.2199,
        0.0927,
        0.0227,
        0.0312,
        0.0329,
    ],
    "T": [
        0.0499,
        0.0349,
        0.0494,
        0.0354,
        0.0366,
        0.0412,
        0.0368,
        0.0297,
        0.0267,
        0.0398,
        0.0334,
        0.0397,
        0.0402,
        0.0254,
        0.0362,
        0.082,
        0.2465,
        0.0227,
        0.028,
        0.0494,
    ],
    "W": [
        0.0054,
        0.0058,
        0.0045,
        0.0037,
        0.0041,
        0.0059,
        0.0055,
        0.0054,
        0.0076,
        0.0059,
        0.0071,
        0.0052,
        0.008,
        0.0169,
        0.0026,
        0.0052,
        0.0059,
        0.4924,
        0.028,
        0.0055,
    ],
    "Y": [
        0.0175,
        0.0174,
        0.0157,
        0.0112,
        0.0122,
        0.0206,
        0.0166,
        0.0108,
        0.0573,
        0.0206,
        0.0223,
        0.0173,
        0.0241,
        0.0888,
        0.0129,
        0.0175,
        0.0178,
        0.0682,
        0.3178,
        0.0206,
    ],
    "V": [
        0.0688,
        0.031,
        0.027,
        0.0243,
        0.0569,
        0.0353,
        0.0313,
        0.0243,
        0.0229,
        0.1767,
        0.0962,
        0.0328,
        0.0924,
        0.055,
        0.031,
        0.0419,
        0.071,
        0.0303,
        0.0467,
        0.2689,
    ],
}
blosum62_matrix = pd.DataFrame(blosum62_matrix).to_numpy()
background = {
    "A": 0.0755236,
    "R": 0.0515842,
    "N": 0.0453131,
    "D": 0.0530344,
    "C": 0.0169811,
    "Q": 0.0402483,
    "E": 0.0632002,
    "G": 0.0684442,
    "H": 0.0224067,
    "I": 0.0573156,
    "L": 0.0934327,
    "K": 0.0594192,
    "M": 0.0235696,
    "F": 0.0407819,
    "P": 0.0492775,
    "S": 0.0722465,
    "T": 0.0574747,
    "W": 0.0125173,
    "Y": 0.0319968,
    "V": 0.0652477,
}
background = pd.DataFrame(pd.Series(background)).T.to_numpy()


def strange_pssm(peptides, halfbits=False, hobohm_cluster=False, kl_logo=False):
    peps_n = len(peptides)
    # seq_len = len(peptides[0])
    seq_len = 9

    # Hobohm algorithm 1 sequence weighting
    if hobohm_cluster:
        alpha, weights = hobohm(peptides)
        alpha = alpha.reshape(-1, 1)
    else:
        alpha = peps_n - 1
        alpha = np.full((seq_len, 1), alpha)
        weights = np.ones(peps_n)
    weights = weights.reshape(-1, 1)

    beta = 200

    # peps = [Seq(p) for p in peptides]
    # m = motifs.create(peps, alphabet=AMINO_ACIDS)
    # # w_obs = pd.DataFrame(m.counts) / peps_n
    # w_obs = pd.DataFrame(m.pwm)
    # aa = w_obs.columns.to_list()
    # w_obs = w_obs.to_numpy()
    aa_mat = np.hsplit(np.array([list(p) for p in peptides]), seq_len)
    w_obs = []
    for i in range(seq_len):
        aa_w = np.sum((AMINO_ACIDS == aa_mat[i]) * weights, axis=0)
        w_obs.append(aa_w)
    w_obs /= np.sum(w_obs, 1).reshape(-1, 1)

    p_c = np.dot(w_obs, blosum62_matrix)
    probs = (alpha * w_obs + beta * p_c) / (alpha + beta)  # self.p
    p_not_0 = probs > 0
    g_ref = np.ones(probs.shape) * background
    if np.any(g_ref == 0):
        for i in np.nonzero(g_ref == 0)[0]:
            g_ref[i, :] = 1 / probs.shape[1]

    pssm = np.zeros(probs.shape)  # self.w
    pssm[probs == 0] = -99.999
    if halfbits:
        pssm[p_not_0] = np.log2(probs[p_not_0] / g_ref[p_not_0]) * 2
        if np.any(pssm < -99.999):
            pssm[pssm < -99.999] = -99.999
    else:
        pssm[p_not_0] = np.log2(probs[p_not_0] / g_ref[p_not_0])
        if np.any(pssm < -50):
            pssm[pssm < -50] = -50

    # Kullback-Leibler logo
    if kl_logo:
        pssm = np.sum(pssm * probs, axis=1).reshape(-1, 1) * probs * np.sign(pssm)

    return pd.DataFrame(pssm, columns=AMINO_ACIDS)


# Filter out the source data and save only clean peptides of one length

motif_atlas_peps = pd.read_csv(motif_atlas_peps_src_file, sep="\t")
motif_atlas_peps[motif_atlas_peps_col_peptide] = motif_atlas_peps[
    motif_atlas_peps_col_peptide
].str.upper()

# Filter out mouse alleles
motif_atlas_peps = motif_atlas_peps.loc[
    ~motif_atlas_peps[motif_atlas_peps_col_allele].str.startswith("H2-")
]

# Get only clean peptides of one length
motif_atlas_peps = motif_atlas_peps.loc[
    motif_atlas_peps[motif_atlas_peps_col_peptide].str.len() == peptide_len
]
motif_atlas_peps = motif_atlas_peps.loc[
    motif_atlas_peps[motif_atlas_peps_col_peptide].apply(are_sequence_residues_valid)
]

motif_atlas_peps.to_csv(motif_atlas_peps_file, index=False)

# Code each allele's peptides into PSSMs

motif_atlas_alleles = motif_atlas_peps.groupby(motif_atlas_peps_col_allele).agg(list)

motif_atlas_alleles = motif_atlas_alleles.rename(
    columns={
        motif_atlas_peps_col_allele: main_motif_atlas_emb_col_allele,
        motif_atlas_peps_col_peptide: main_motif_atlas_emb_col_peps,
    }
)
motif_atlas_alleles[main_motif_atlas_emb_col_peps] = motif_atlas_alleles[
    main_motif_atlas_emb_col_peps
].apply(lambda x: sorted(set(x)))

motif_atlas_alleles = motif_atlas_alleles.reset_index()

# Select PSSM type
halfbits = True
hobohm_cluster = True

alleles_peptides = motif_atlas_alleles[main_motif_atlas_emb_col_peps]

alleles_pssms = []

for peptides in tqdm(alleles_peptides):
    pssm = strange_pssm(
        peptides, halfbits=halfbits, hobohm_cluster=hobohm_cluster
    ).to_numpy()

    alleles_pssms.append(pssm.flatten().tolist())

alleles_pssms = pd.DataFrame(alleles_pssms)

# Load the UMAP model and create embeddings

with open(umap_reducer_main_file, "rb") as handle:
    reducer = pickle.load(handle)

embedding = reducer.transform(alleles_pssms.to_numpy())

embedding = pd.DataFrame(embedding)
embedding.columns = [main_motif_atlas_emb_col_x, main_motif_atlas_emb_col_y]

# Gather everything and save the table

motif_atlas_alleles[main_motif_atlas_emb_col_peps] = motif_atlas_alleles[
    main_motif_atlas_emb_col_peps
].apply(lambda x: main_motif_atlas_emb_list_sep.join(x))

motif_atlas_alleles = pd.concat([motif_atlas_alleles, embedding], axis=1)

motif_atlas_alleles.to_csv(main_motif_atlas_emb_file, sep="\t", index=False)
