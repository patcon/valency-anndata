import valency_anndata as val

# adata = val.datasets.polis("https://pol.is/report/r2dfw8eambusb8buvecjt") # small
adata = val.datasets.polis.load("https://pol.is/report/r29kkytnipymd3exbynkd", translate_to="en") # Chile
# adata = val.datasets.polis("https://pol.is/3ntrtcehas")

print(adata)
if adata.X is not None:
    print(adata.X.shape)
print(adata.X)
print(adata.uns["schema"])
print(adata.var)

val.datasets.polis.translate_statements(adata, translate_to="ja")
print(adata.var)