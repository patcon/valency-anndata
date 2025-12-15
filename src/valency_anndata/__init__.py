from . import datasets, preprocessing, tools
from . import scanpy

pp = preprocessing
tl = tools

__all__ = [
    "datasets",
    "preprocessing",
    "tools",
    # Backward-compat with scanpy.
    "pp",
    "tl",
    # Make all of scanpy accessible within valency_anndata
    "scanpy",
]