from scanpy.preprocessing._pca import pca
from scanpy.tools._tsne import tsne
from scanpy.tools._umap import umap
from scanpy.tools._leiden import leiden
from ._kmeans import kmeans
from ._polis import recipe_polis
from ..utils import _reexport_with_doc

pca  = _reexport_with_doc(pca)
tsne = _reexport_with_doc(tsne)
umap = _reexport_with_doc(umap)
leiden = _reexport_with_doc(leiden)

__all__ = [
    "kmeans",
    "recipe_polis",

    # Simple re-export of scanpy.
    "pca",
    "tsne",
    "umap",
    "leiden",
]
