from typing import Optional
from anndata import AnnData
from scanpy import logging as logg

def localmap(
    adata: AnnData,
    *,
    layer: str = "X_imputed",
    n_neighbors: Optional[int] = None,
    n_components: int = 2,
    mask_var: str | None = None,
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | None:
    """
    Compute LocalMAP dimensionality reduction.

    Parameters
    ----------
    adata
        AnnData object.
    layer
        Layer to use for computation. Default is "X_imputed".
    n_neighbors
        Number of neighbors for LocalMAP.
    n_components
        Number of dimensions for the embedding. Default is 2.
    mask_var
        Column name in `adata.var` to use for masking variables.
        If provided, only variables where `mask_var` is True will be used.
    key_added
        Key under which to store the embedding in `adata.obsm`.
        Default is "X_localmap".
    copy
        Return a copy instead of modifying adata in place.

    Returns
    -------
    AnnData | None
        Returns AnnData if `copy=True`, otherwise returns None.
    """
    adata = adata.copy() if copy else adata

    key_obsm, key_uns = ("X_localmap", "localmap") if key_added is None else [key_added] * 2

    start = logg.info("computing LocalMAP")

    from pacmap import LocalMAP

    import inspect as _inspect
    _localmap_params = _inspect.signature(LocalMAP.__init__).parameters
    _knn_kwargs = {"knn_backend": "voyager"} if "knn_backend" in _localmap_params else {}

    estimator = LocalMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        **_knn_kwargs,
    )

    # Get data from layer, optionally filtering by mask_var
    X = adata.layers[layer]
    if mask_var is not None:
        mask = adata.var[mask_var].values
        X = X[:, mask]

    X_reduced = estimator.fit_transform(X)

    adata.obsm[key_obsm] = X_reduced

    return adata if copy else None

def pacmap(
    adata: AnnData,
    *,
    layer: str = "X_imputed",
    n_neighbors: Optional[int] = None,
    n_components: int = 2,
    mask_var: str | None = None,
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | None:
    """
    Compute PaCMAP dimensionality reduction.

    Parameters
    ----------
    adata
        AnnData object.
    layer
        Layer to use for computation. Default is "X_imputed".
    n_neighbors
        Number of neighbors for PaCMAP.
    n_components
        Number of dimensions for the embedding. Default is 2.
    mask_var
        Column name in `adata.var` to use for masking variables.
        If provided, only variables where `mask_var` is True will be used.
    key_added
        Key under which to store the embedding in `adata.obsm`.
        Default is "X_pacmap".
    copy
        Return a copy instead of modifying adata in place.

    Returns
    -------
    AnnData | None
        Returns AnnData if `copy=True`, otherwise returns None.
    """
    adata = adata.copy() if copy else adata

    key_obsm, key_uns = ("X_pacmap", "pacmap") if key_added is None else [key_added] * 2

    start = logg.info("computing PaCMAP")

    from pacmap import PaCMAP

    import inspect as _inspect
    _pacmap_params = _inspect.signature(PaCMAP.__init__).parameters
    _knn_kwargs = {"knn_backend": "voyager"} if "knn_backend" in _pacmap_params else {}

    estimator = PaCMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
        **_knn_kwargs,
    )

    # Get data from layer, optionally filtering by mask_var
    X = adata.layers[layer]
    if mask_var is not None:
        mask = adata.var[mask_var].values
        X = X[:, mask]

    X_reduced = estimator.fit_transform(X)

    adata.obsm[key_obsm] = X_reduced

    return adata if copy else None