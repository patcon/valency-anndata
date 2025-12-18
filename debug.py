import valency_anndata as val

# adata = val.datasets.polis.load("https://pol.is/report/r2dfw8eambusb8buvecjt") # small
adata = val.datasets.polis.load("https://pol.is/report/r29kkytnipymd3exbynkd", translate_to=None) # Chile
# adata = val.datasets.polis.load("https://pol.is/3ntrtcehas")
# adata = val.datasets.polis.load("https://pol.is/89euwydkce") # australia solar
adata = val.datasets.polis.load("https://pol.is/2demo") # minimum wage

print(adata)
if adata.X is not None:
    print(adata.X.shape)
print(adata.X)
print(adata.uns["schema"])
print(adata.var)

val.tools.recipe_polis(adata, key_added_pca="X_pca", key_added_kmeans="kmeans_polis")
val.scanpy.pl.pca(adata, color="kmeans_polis")