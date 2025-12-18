import valency_anndata as val
from anndata import AnnData
import numpy as np
import pandas as pd

if False:
    adata = val.datasets.polis.load("https://pol.is/report/r2dfw8eambusb8buvecjt") # small
    # adata = val.datasets.polis.load("https://pol.is/report/r29kkytnipymd3exbynkd", translate_to=None) # Chile
    # adata = val.datasets.polis.load("https://pol.is/3ntrtcehas")
    # adata = val.datasets.polis.load("https://pol.is/89euwydkce") # australia solar
    # adata = val.datasets.polis.load("https://pol.is/2demo") # minimum wage

def make_fake_polis_adata(
    n_participants: int = 400,
    n_statements: int = 60,
) -> AnnData:
    """
    Create a fake Polis-like AnnData object with:
      - X: vote matrix (participants × statements)
      - obs: participant metadata
      - var: statement metadata
    """

    rng = np.random.default_rng(42)

    # -------------------
    # Fake vote matrix
    # Values roughly match Polis: -1 disagree, 0 pass, 1 agree
    # -------------------
    X = rng.choice(
        [-1.0, 0.0, 1.0],
        size=(n_participants, n_statements),
        p=[0.35, 0.30, 0.35],
    )

    # -------------------
    # obs: participants
    # -------------------
    obs = pd.DataFrame(
        {
            "participant_id": [f"p{i:04d}" for i in range(n_participants)],
            "n_votes": (X != 0).sum(axis=1),
        }
    ).set_index("participant_id")

    # -------------------
    # var: statements
    # -------------------
    var = pd.DataFrame(
        {
            "statement_id": [f"s{i:03d}" for i in range(n_statements)],
            "text_length": rng.integers(20, 200, size=n_statements),
        }
    ).set_index("statement_id")

    return AnnData(X=X, obs=obs, var=var)

def fake_recipe_polis(
    adata: AnnData,
    *,
    key_added_pca: str = "X_pca",
    key_added_kmeans: str = "kmeans_polis",
    n_pca: int = 2,
    n_clusters: int = 5,
):
    rng = np.random.default_rng(0)

    n_obs = adata.n_obs

    # -------------------
    # Fake PCA (obsm)
    # -------------------
    adata.obsm[key_added_pca] = rng.normal(size=(n_obs, n_pca))

    # -------------------
    # Fake KMeans labels (obs)
    # -------------------
    labels = rng.integers(0, n_clusters, size=n_obs)
    adata.obs[key_added_kmeans] = pd.Categorical(labels)

    # -------------------
    # Fake layers
    # -------------------
    # Layer: X_masked (simulate some missing values)
    X_masked = adata.X.copy()
    mask = rng.random(size=X_masked.shape) < 0.2  # 20% missing
    X_masked[mask] = 0
    adata.layers["X_masked"] = X_masked

    # Layer: X_masked_imputed_mean (replace missing with column mean)
    X_imputed = X_masked.copy()
    col_means = np.where(X_masked != 0, X_masked, np.nan).mean(axis=0)
    inds = np.where(X_imputed == 0)
    X_imputed[inds] = np.take(col_means, inds[1])
    adata.layers["X_imputed_mean"] = X_imputed

    # -------------------
    # Provenance marker (uns)
    # -------------------
    adata.uns.setdefault("polis", {})
    adata.uns["polis"]["fake"] = True

adata = make_fake_polis_adata()

adata_snap = adata.copy()

with val.viz.schematic_diagram(diff_from=adata):
    fake_recipe_polis(
        adata,
        key_added_pca="X_pca",
        key_added_kmeans="kmeans_polis",
    )
    # val.tools.recipe_polis(adata, key_added_pca="X_pca", key_added_kmeans="kmeans_polis")

val.viz.schematic_diagram(adata)

print(adata_snap)
print("=========")
print(adata)

adata_snap_layers = adata.copy()
del adata.obs["n_votes"]
del adata.layers["X_masked"]
val.viz.schematic_diagram(adata, diff_from=adata_snap_layers)

# val.scanpy.pl.pca(adata, color="kmeans_polis")