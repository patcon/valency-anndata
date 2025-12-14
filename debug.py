import valency_anndata as val

adata = val.datasets.polis("https://pol.is/report/r2dfw8eambusb8buvecjt")
# adata = val.datasets.polis("https://pol.is/3ntrtcehas")

print(adata)
print(adata.X.shape)
print(adata.X)
print(adata.uns["schema"])