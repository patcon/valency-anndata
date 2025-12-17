import numpy as np
from reddwarf.sklearn.transformers import calculate_scaling_factors
from anndata import AnnData
import valency_anndata as val

def recipe_polis(
    adata: AnnData,
    *,
    key_added_pca: str = "X_polis",
    key_added_kmeans: str = "kmeans_polis",
    inplace: bool = True,
):
    if not inplace:
        adata = adata.copy()

    # Preconditions
    assert isinstance(adata.X, np.ndarray)

    # 1. Pre-reducer mask
    adata.var["pre_reducer_mask"] = adata.var.eval("~is_meta and moderation_state > -1")

    mask = adata.var["pre_reducer_mask"].to_numpy()

    # 2. Mask votes
    X_masked = adata.X.copy()
    X_masked[:, ~mask] = 0
    adata.layers["X_masked"] = X_masked

    # 3. Impute
    val.preprocessing.impute(
        adata,
        strategy="mean",
        source_layer="X_masked",
        target_layer="X_masked_imputed_mean",
    )

    # 4. PCA (unscaled)
    val.tools.pca(
        adata,
        layer="X_masked_imputed_mean",
        key_added="X_pca_masked_unscaled",
    )

    # 5. Scale PCA using participation structure
    scaling_factors = calculate_scaling_factors(adata.X)
    X_pca_unscaled = adata.obsm["X_pca_masked_unscaled"]
    adata.obsm["X_pca_masked_scaled"] = X_pca_unscaled / scaling_factors[:, None]

    # Set a recognizable key
    adata.obsm[key_added_pca] = adata.obsm["X_pca_masked_scaled"]

    val.tools.kmeans(
        adata,
        use_rep=key_added_pca,
        k_bounds=(2, 5),
        init="polis",
        key_added=key_added_kmeans,
        inplace=inplace,
    )

    if not inplace:
        return adata