from ._schematic import adata_structure_svg
from ._context import _SchematicDiagramContext
from ._utils import _show_svg
from typing import Literal, Optional, overload
from anndata import AnnData

@overload
def schematic_diagram(adata: AnnData) -> None: ...
@overload
def schematic_diagram(adata: AnnData, *, diff_from: Optional[AnnData] | Literal[False]) -> None: ...
@overload
def schematic_diagram(*, diff_from: AnnData) -> _SchematicDiagramContext: ...

def schematic_diagram(
    adata: Optional[AnnData] = None,
    *,
    diff_from: Optional[AnnData] | Literal[False] = False,
):
    """
    Render a schematic diagram of an AnnData object, optionally highlighting
    structural differences relative to a snapshot.

    This function supports two usage modes: **render mode** and **context-manager mode**.

    -------------------------------------------------------------------------
    1. Render mode
    -------------------------------------------------------------------------
    Render a diagram of `adata` immediately.

    Examples
    --------
    >>> schematic_diagram(adata)
    >>> schematic_diagram(adata, diff_from=adata_snapshot)

    Behavior
    --------
    - Visualizes `adata` structure (`X`, `obs`, `var`, `layers`, `obsm`).
    - If `diff_from` is provided:
        - Highlights additions and removals relative to `diff_from`.
    - If `diff_from` is `None`:
        - Highlights all entries as additions (diff from empty AnnData).
    - If `diff_from` is `False`:
        - No diff highlighting is applied.
    - The diagram is displayed inline (notebooks) or in a browser (script).

    -------------------------------------------------------------------------
    2. Context-manager mode
    -------------------------------------------------------------------------
    Capture a snapshot on entering a `with` block, rendering a diff on exit.

    Example
    -------
    >>> with schematic_diagram(diff_from=adata):
    ...     some_transformation(adata)

    Behavior
    --------
    - `diff_from` must be provided; `adata` must be omitted.
    - On entry, a snapshot of `diff_from` is recorded.
    - On exit, a diff diagram between the snapshot and current `adata` is rendered.
    - Exceptions inside the `with` block prevent rendering.

    Parameters
    ----------
    adata : AnnData, optional
        The AnnData object to visualize (required in render mode, must be omitted in
        context-manager mode).
    diff_from : AnnData, None, or False, default False
        Determines the snapshot to diff against: (must be AnnData in context-manager mode)
        - `AnnData` instance: highlights differences from the snapshot.
        - `None`: highlights all entries as additions (diff from empty).
        - `False`: disables diff highlighting.

    Returns
    -------
    None
        In render mode, the diagram is displayed; nothing is returned.
    _SchematicDiagramContext
        In context-manager mode, a context manager for automatic diff rendering.

    Notes
    -----
    - Explicit diff rendering always takes precedence over context-manager snapshots.
    - Snapshots are stored internally to allow nested diff scopes.
    - This function does not mutate `adata`.
    """
    if adata is None:
        # ------------------
        # Context mode
        # ------------------
        if isinstance(diff_from, AnnData):
            return _SchematicDiagramContext(diff_from)
    else:
        # ------------------
        # Render mode (explicit or implicit diff)
        # ------------------
        if diff_from is False:
            base = None
        elif diff_from is None:
            base = AnnData()
        else:
            base = diff_from

        dwg = adata_structure_svg(adata, diff_from=base)
        _show_svg(dwg)
        return None

    raise TypeError("Invalid schematic_diagram() call")
