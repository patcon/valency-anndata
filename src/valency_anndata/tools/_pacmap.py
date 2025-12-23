from typing import Optional
from anndata import AnnData
from scanpy import logging as logg
from scanpy._utils import NeighborsView

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

    pacmap = PaCMAP(
        n_components=n_components,
        n_neighbors=n_neighbors,
    )

    X_pacmap = pacmap.fit_transform(adata.layers[layer])

    adata.obsm[key_obsm] = X_pacmap

    return adata if copy else None