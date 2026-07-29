from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import colormaps as _mpl_colormaps, pyplot as _plt


def _resolve_mask(
    adata: AnnData, mask: "np.ndarray | str | None", dim: str
) -> "np.ndarray | None":
    """Resolve ``mask_obs``/``mask_var`` to a boolean array, per scanpy's mask convention."""
    if mask is None:
        return None
    annot = getattr(adata, dim)
    n = adata.n_obs if dim == "obs" else adata.n_vars
    if isinstance(mask, str):
        if mask not in annot.columns:
            msg = f"Did not find `adata.{dim}[{mask!r}]`."
            raise ValueError(msg)
        mask_array = annot[mask].to_numpy()
    else:
        mask_array = np.asarray(mask)
        if len(mask_array) != n:
            msg = f"The shape of the mask does not match adata.n_{dim}."
            raise ValueError(msg)
    if not pd.api.types.is_bool_dtype(mask_array.dtype):
        msg = f"mask_{dim} array must be boolean."
        raise ValueError(msg)
    return mask_array


# Register a brighter discrete-friendly variant of RdYlGn.
# Takes the red and green endpoints directly from RdYlGn, but overrides
# the muted midpoint yellow with a fully saturated yellow (#ffff00).
_rdylgn = _plt.get_cmap("RdYlGn")
_RdYlGnBright = ListedColormap(
    [_rdylgn(0.0), "#ffff00", _rdylgn(1.0)],
    name="RdYlGnBright",
)
try:
    _mpl_colormaps.register(_RdYlGnBright)
except ValueError:
    pass  # already registered

# Primary-colour variant: fully saturated red, yellow, and green.
_RdYlGnPrimary = ListedColormap(["#d73027", "#ffff00", "#1a9850"], name="RdYlGnPrimary")
try:
    _mpl_colormaps.register(_RdYlGnPrimary)
except ValueError:
    pass  # already registered


def heatmap(
    adata: AnnData,
    groupby: str | None = None,
    mask_obs: np.ndarray | str | None = None,
    mask_var: np.ndarray | str | None = None,
    cmap: str = "RdYlGn",
    discrete: bool = False,
    show_labels: bool = True,
    max_tick_labels: int | None = 50,
    show: bool = True,
    **kwargs,
):
    """
    Plot a vote-matrix heatmap with Polis-friendly defaults.

    A thin wrapper around :func:`scanpy.pl.heatmap` that adds optional discrete
    colorbar labelling and removes the requirement for a ``groupby`` column.

    Parameters
    ----------
    adata
        AnnData object with participants as observations and statements as variables.
    groupby
        Column in ``adata.obs`` to group participants by. When ``None``,
        participants are shown in their current index order with no grouping.
    mask_obs
        Boolean mask (or name of a boolean ``adata.obs`` column) selecting which
        participants to include, following scanpy's ``mask_obs`` convention
        (see e.g. :func:`scanpy.pl.pca`). Useful for excluding participants
        with no ``groupby`` assignment (e.g. unclustered rows).
    mask_var
        Boolean mask (or name of a boolean ``adata.var`` column) selecting
        which statements to include, following scanpy's ``mask_var``
        convention (see e.g. :func:`scanpy.pp.pca`). Useful for excluding
        moderated-out statements, e.g. ``mask_var=(adata.var["moderation_state"] >= 1).to_numpy()``.
    cmap
        Colormap name. Defaults to ``"RdYlGn"``. Also accepts ``"RdYlGnBright"``,
        a custom :class:`~matplotlib.colors.ListedColormap` with fully saturated
        red, yellow, and green (``["#d73027", "#ffff00", "#1a9850"]``).
    discrete
        When ``True``, renders a segmented colorbar with labelled ticks
        (``"disagree (-1)"``, ``"pass (0)"``, ``"agree (+1)"``), using a
        :class:`~matplotlib.colors.BoundaryNorm` with boundaries at
        ``[-1.5, -0.5, 0.5, 1.5]``.
    show_labels
        Whether to show participant (row) and statement (column) tick labels.
        Defaults to ``True``.
    max_tick_labels
        Maximum number of tick labels to show on each axis when
        ``show_labels=True``. Labels are thinned by a uniform stride so at
        most this many appear. Set to ``None`` to show all labels. Defaults
        to ``50``.
    show
        Whether to call ``plt.show()`` at the end. Set to ``False`` to get back
        the axes dict for further customisation.
    **kwargs
        Additional keyword arguments forwarded to :func:`scanpy.pl.heatmap`.

    Returns
    -------
    None or dict
        When ``show=True`` (default) returns ``None``. When ``show=False``,
        returns the axes dictionary from :func:`scanpy.pl.heatmap`.

    Examples
    --------

    ```py
    adata = val.datasets.polis.load("https://pol.is/report/r29kkytnipymd3exbynkd")

    val.viz.heatmap(adata, discrete=True)
    ```
    """
    import scanpy as sc
    import matplotlib.pyplot as plt

    obs_mask_array = _resolve_mask(adata, mask_obs, "obs")
    var_mask_array = _resolve_mask(adata, mask_var, "var")
    if obs_mask_array is not None or var_mask_array is not None:
        adata = adata[
            obs_mask_array if obs_mask_array is not None else slice(None),
            var_mask_array if var_mask_array is not None else slice(None),
        ].copy()

    _dummy_col = None
    if groupby is None:
        _dummy_col = "__heatmap_dummy__"
        adata.obs[_dummy_col] = "all"
        adata.obs[_dummy_col] = adata.obs[_dummy_col].astype("category")
        groupby = _dummy_col

    if discrete:
        ticks = np.array([-1.0, 0.0, 1.0])
        boundaries = np.array([-1.5, -0.5, 0.5, 1.5])
        norm = BoundaryNorm(boundaries, ncolors=plt.get_cmap(cmap).N)
        # scanpy forbids passing norm alongside vmin/vmax/vcenter
        kwargs.pop("vmin", None)
        kwargs.pop("vmax", None)
        kwargs.pop("vcenter", None)
        kwargs["norm"] = norm

    kwargs.pop("show_gene_labels", None)
    if show_labels and max_tick_labels is None:
        # Let scanpy render x-axis labels natively — it resizes the figure to
        # fit them properly. y-axis (obs) labels are never set by scanpy.
        kwargs["show_gene_labels"] = True

    # Suppress scanpy's "Gene labels are not shown when more than 50 genes"
    # warning — only fires when show_gene_labels is not set, i.e. when we
    # will manage labels ourselves via post-hoc FixedLocator.
    _prev_verbosity = sc.settings.verbosity
    sc.settings.verbosity = 0  # errors only
    try:
        axes_dict = sc.pl.heatmap(
            adata,
            var_names=adata.var_names,
            groupby=groupby,
            cmap=cmap,
            show=False,
            **kwargs,
        )
    finally:
        sc.settings.verbosity = _prev_verbosity

    if _dummy_col is not None:
        del adata.obs[_dummy_col]
        # Hide the uninformative groupby axis strip
        if axes_dict and "groupby_ax" in axes_dict:
            axes_dict["groupby_ax"].set_visible(False)

    if show_labels and max_tick_labels is not None:
        # Post-hoc thinning: scanpy's native rendering is bypassed, so we set
        # both axes manually via FixedLocator/FixedFormatter with a stride.
        heatmap_ax = axes_dict.get("heatmap_ax") if axes_dict else None
        if heatmap_ax is not None:
            import matplotlib.ticker as ticker

            def _strided(names, max_n):
                stride = max(1, len(names) // max_n)
                indices = list(range(0, len(names), stride))
                return indices, [names[i] for i in indices]

            # x-axis
            var_names = list(adata.var_names)
            x_indices, x_labels = _strided(var_names, max_tick_labels)
            heatmap_ax.xaxis.set_major_locator(ticker.FixedLocator(x_indices))
            heatmap_ax.xaxis.set_major_formatter(ticker.FixedFormatter(x_labels))
            heatmap_ax.tick_params(axis="x", labelbottom=True, labelsize=8, rotation=90)

            # y-axis: scanpy always sets labelleft=False and leaves a
            # FuncFormatter returning empty strings, so replace both.
            obs_names = list(adata.obs_names)
            y_indices, y_labels = _strided(obs_names, max_tick_labels)
            heatmap_ax.yaxis.set_major_locator(ticker.FixedLocator(y_indices))
            heatmap_ax.yaxis.set_major_formatter(ticker.FixedFormatter(y_labels))
            heatmap_ax.tick_params(axis="y", labelleft=True, labelsize=8)

    if discrete:
        ticklabels = ["disagree (-1)", "pass (0)", "agree (+1)"]
        # scanpy renders the colorbar via plt.colorbar(image, cax=...), so the
        # Colorbar object is stored on the image artist, not on the QuadMesh.
        cbar = None
        for ax in plt.gcf().axes:
            for img in ax.images:
                cb = getattr(img, "colorbar", None)
                if cb is not None:
                    cbar = cb
                    break
            if cbar is not None:
                break
        if cbar is not None:
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(ticklabels)

    if show:
        plt.show()
        return None

    return axes_dict
