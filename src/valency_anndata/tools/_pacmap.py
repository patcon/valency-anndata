from typing import Optional
from anndata import AnnData
from scanpy import logging as logg

def localmap(
    adata: AnnData,
    *,
    layer: str = "X_imputed",
    n_neighbors: Optional[int] = None,
    n_components: int = 2,
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | None:
    """
    """
    adata = adata.copy() if copy else adata

    key_obsm, key_uns = ("X_localmap", "localmap") if key_added is None else [key_added] * 2

    start = logg.info("computing LocalMAP")

    from pacmap import LocalMAP

    estimator = LocalMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
    )

    X_reduced = estimator.fit_transform(adata.layers[layer])

    adata.obsm[key_obsm] = X_reduced

    return adata if copy else None

def pacmap(
    adata: AnnData,
    *,
    layer: str = "X_imputed",
    n_neighbors: Optional[int] = None,
    n_components: int = 2,
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | None:
    """
    """
    adata = adata.copy() if copy else adata

    key_obsm, key_uns = ("X_pacmap", "pacmap") if key_added is None else [key_added] * 2

    start = logg.info("computing PaCMAP")

    from pacmap import PaCMAP

    estimator = PaCMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
    )

    X_reduced = estimator.fit_transform(adata.layers[layer])

    adata.obsm[key_obsm] = X_reduced

    return adata if copy else None