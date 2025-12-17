import numpy as np
from anndata import AnnData
from reddwarf.sklearn.cluster import BestPolisKMeans
from typing import Optional, Tuple


def kmeans(
    adata: AnnData,
    use_rep: Optional[str] = None,
    k_bounds: Optional[Tuple[int, int]] = None,
    init: str = "kmeans++",
    init_centers: Optional[np.ndarray] = None,
    random_state: Optional[int] = None,
    key_added: str = "kmeans",
    inplace: bool = True,
) -> AnnData | None:
    """
    Apply BestPolisKMeans clustering to an AnnData object.

    Parameters
    ----------
    adata : AnnData
        Input data. Must have `.X` as a numpy array.
    k_bounds : tuple[int, int] or list[int], optional
        Minimum and maximum number of clusters to try. Defaults to [2, 5].
    init : {'k-means++', 'random', 'polis'}
        Initialization method for KMeans. Defaults to 'polis'.
    init_centers : array-like, optional
        Initial cluster centers to use.
    random_state : int, optional
        Random seed for reproducibility.
    key_added : str
        Name of the column to store cluster labels in `adata.obs`.
    inplace : bool
        If True, modify `adata` in place and return None.
        If False, return a copy with the clustering added.

    Returns
    -------
    AnnData or None
        Returns a copy if `inplace=False`, otherwise modifies in place.
    """
    if use_rep is None:
        X = adata.X
    else:
        if use_rep not in adata.obsm:
            raise KeyError(f"use_rep='{use_rep}' not found in adata.obsm")
        X = adata.obsm[use_rep]

    if k_bounds is None:
        k_bounds_list = [2, 5]
    else:
        k_bounds_list = list(k_bounds)

    if not isinstance(X, np.ndarray):
        raise ValueError("adata.X must be a numpy array.")

    best_kmeans = BestPolisKMeans(
        k_bounds=k_bounds_list,
        init=init,
        init_centers=init_centers,
        random_state=random_state,
    )
    best_kmeans.fit(X)

    if not best_kmeans.best_estimator_:
        raise RuntimeError("BestPolisKMeans did not find a valid estimator.")

    labels = best_kmeans.best_estimator_.labels_

    if inplace:
        adata.obs[key_added] = labels
        adata.uns[key_added] = {
            "best_k": best_kmeans.best_k_,
            "best_score": best_kmeans.best_score_,
            "init": init,
        }
        return None
    else:
        adata_copy = adata.copy()
        adata_copy.obs[key_added] = labels
        adata_copy.uns[key_added] = {
            "best_k": best_kmeans.best_k_,
            "best_score": best_kmeans.best_score_,
            "init": init,
        }
        return adata_copy