from scanpy.plotting._tools.scatterplots import (
    pca,
    umap,
)
from ..utils import _reexport_with_doc

from ._langevitour import langevitour
from .schematic_diagram import schematic_diagram


pca = _reexport_with_doc(pca)
umap = _reexport_with_doc(umap)

__all__ = [
    "langevitour",
    "schematic_diagram",

    # Simple re-export of scanpy.
    "pca",
    "umap",
]