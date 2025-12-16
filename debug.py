import valency_anndata as val
import numpy as np

# adata = val.datasets.polis.load("https://pol.is/report/r2dfw8eambusb8buvecjt") # small
adata = val.datasets.polis.load("https://pol.is/report/r29kkytnipymd3exbynkd", translate_to=None) # Chile
# adata = val.datasets.polis.load("https://pol.is/3ntrtcehas")

print(adata)
if adata.X is not None:
    print(adata.X.shape)
print(adata.X)
print(adata.uns["schema"])
print(adata.var)

# val.datasets.polis.translate_statements(adata, translate_to="ja")
if not isinstance(adata.X, np.ndarray):
    raise TypeError("valency-anndata Polis pipeline assumes dense adata.X")

val.tools.recipe_polis(adata, key_added="X_polis")
val.scanpy.pl.embedding(adata, basis="polis")