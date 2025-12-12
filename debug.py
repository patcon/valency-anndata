import valenci_anndata as val

adata = val.datasets.polis("https://pol.is/report/r2dfw8eambusb8buvecjt")
# adata = val.datasets.polis("https://pol.is/3ntrtcehas")

val.preprocessing.rebuild_vote_matrix(adata)
print(adata)