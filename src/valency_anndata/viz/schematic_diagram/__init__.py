from ._schematic import adata_structure_svg
from ._context import _SchematicDiagramContext
from ._utils import _show_svg
from typing import Optional, overload
from anndata import AnnData

@overload
def schematic_diagram(adata: AnnData) -> None: ...
@overload
def schematic_diagram(adata: AnnData, *, diff_from: AnnData) -> None: ...
@overload
def schematic_diagram(*, diff_from: AnnData) -> _SchematicDiagramContext: ...

def schematic_diagram(
    adata: Optional[AnnData] = None,
    *,
    diff_from: Optional[AnnData] = None,
):
    """
    Render a schematic SVG diagram of an AnnData object, optionally highlighting
    structural differences relative to a snapshot.

    This function supports two usage modes:

    -------------------------------------------------------------------------
    1. Render mode
    -------------------------------------------------------------------------
    Render a schematic diagram of ``adata``.

    Examples
    --------
    >>> schematic_diagram(adata)
    >>> schematic_diagram(adata, diff_from=adata_snapshot)

    In render mode:
    - The diagram visualizes the structure of ``adata`` (X, obs, var).
    - If ``diff_from`` is provided, labels in ``obs``, ``var``, ``layers``, and ``obsm``
      are highlighted to indicate additions and removals relative to ``diff_from``.
    - If ``diff_from`` is not provided, **no diff highlighting is applied**; the
      diagram is rendered normally.
    - The diagram is displayed inline (in notebooks) or opened in a browser
      (in script mode).

    Parameters
    ----------
    adata : AnnData, optional
        The AnnData object to visualize. Must be provided in render mode.
    diff_from : AnnData, optional
        An optional snapshot to diff against.

    Returns
    -------
    None
        The diagram is displayed; nothing is returned.

    -------------------------------------------------------------------------
    2. Context-manager mode
    -------------------------------------------------------------------------
    Capture a snapshot of an AnnData object at the beginning of a ``with`` block,
    and automatically render a diff schematic on exit.

    Examples
    --------
    >>> with schematic_diagram(diff_from=adata):
    ...     val.tools.some_transformation(adata, inplace=True)
    >>> with schematic_diagram():
    ...     adata = val.datasets.polis.load(...)

    In context-manager mode:
    - If ``diff_from`` is provided, a snapshot of it is taken on entry and used
      as the baseline for highlighting differences when the block exits.
    - If ``diff_from`` is not provided, a fresh empty ``AnnData()`` is used as the
      baseline, so all entries in ``adata`` at the end of the block will appear
      as “added”.
    - ``adata`` should be omitted when using the context manager.
    - On exit, a schematic diff between the snapshot and the current state of
      ``adata`` is rendered automatically.
    - Exceptions inside the ``with`` block will prevent rendering.

    Returns
    -------
    _SchematicDiagramContext
        A context manager that records a diff snapshot and renders automatically.

    Notes
    -----
    - Explicit diffs (``schematic_diagram(adata, diff_from=...)``) always take
      precedence over implicit diffs captured via a context manager.
    - Snapshots are stored internally to allow nested diff scopes.
    - This function does not mutate ``adata``.
    """
    # ------------------
    # Context mode
    # ------------------
    if adata is None:
        if diff_from is None:
            # Use empty AnnData for automatic diffs
            diff_from = AnnData()
        return _SchematicDiagramContext(diff_from)

    # ------------------
    # Render mode (explicit or implicit diff)
    # ------------------
    if adata is not None:
        dwg = adata_structure_svg(adata, diff_from=diff_from)
        _show_svg(dwg)
        return None
    raise TypeError("Invalid schematic_diagram() call")
