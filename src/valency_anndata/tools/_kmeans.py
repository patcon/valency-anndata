import numpy as np
import pandas as pd
from anndata import AnnData
from reddwarf.sklearn.cluster import BestPolisKMeans
from typing import Optional, Tuple
from scanpy.get import _check_mask
from numpy.typing import NDArray


def kmeans(
    adata: AnnData,
    use_rep: Optional[str] = None,
    k_bounds: Optional[Tuple[int, int]] = None,
    init: str = "kmeans++",
    init_centers: Optional[np.ndarray] = None,
    random_state: Optional[int] = None,
    mask_obs: NDArray[np.bool_] | str | None = None,
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
    mask_obs
        Restrict clustering to a certain set of observations. The mask is
        specified as a boolean array or a string referring to an array in
        :attr:`~anndata.AnnData.obs`.
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

    if X is None:
        raise ValueError("No data matrix found for clustering.")

    if k_bounds is None:
        k_bounds_list = [2, 5]
    else:
        k_bounds_list = list(k_bounds)

    mask = _check_mask(adata, mask_obs, "obs")
    if mask is None:
        X_cluster = X
    else:
        X_cluster = X[mask]
        if X_cluster.shape[0] == 0:
            raise ValueError("mask_obs excludes all observations.")

    if not isinstance(X, np.ndarray):
        raise ValueError("adata.X must be a numpy array.")

    best_kmeans = BestPolisKMeans(
        k_bounds=k_bounds_list,
        init=init,
        init_centers=init_centers,
        random_state=random_state,
    )
    best_kmeans.fit(X_cluster)

    if not best_kmeans.best_estimator_:
        raise RuntimeError("BestPolisKMeans did not find a valid estimator.")

    raw_labels = best_kmeans.best_estimator_.labels_

    if mask is None:
        full_labels = raw_labels
    else:
        # dtype=object keeps labels from casting to float.
        full_labels = np.full(adata.n_obs, np.nan, dtype=object)
        full_labels[mask] = raw_labels

    labels = pd.Categorical(full_labels)

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