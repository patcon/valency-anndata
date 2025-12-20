import pandas as pd
from anndata import AnnData
from typing import Sequence, Tuple, Optional, Any
import re


# Use Any for now to handle all possible AnnData matrix types
Matrix = Any

def langevitour(
    adata: AnnData,
    *,
    use_reps: Optional[Sequence[str]] = None,
    color: Optional[str] = None,
    scale: Optional[str] = None,
    initial_axes: Optional[list[str]] = None,
    point_size: int = 2,
    **kwargs,
):
    """
    Interactive Langevitour visualization over one or more representations.

    Parameters
    ----------
    adata
        AnnData object.
    use_reps
        Representations to include, e.g. ["X_pca[:10]", "X_umap[:2]"].
    color
        obs column for grouping / coloring.
    scale
        obs column for point scaling.
    initial_axes
        Set up to 3 initial axes xyz, e.g. ["U_umap[:2]", "X_pca[0]"]
    point_size
        Base point size.
    **kwargs
        Passed through to `Langevitour`.
        See R docs: https://logarithmic.net/langevitour/reference/langevitour.html#arguments
    """
    import warnings

    with warnings.catch_warnings():
        # Prevent setuptools from showing a warning about
        # Langevitour using `import pkg_resources`.
        warnings.filterwarnings(
            "ignore",
            message="pkg_resources is deprecated as an API",
            category=UserWarning,
        )
        from langevitour import Langevitour

    X_df = resolve_use_reps(adata, use_reps)

    group = adata.obs[color].tolist() if color is not None else None

    state = {}
    if initial_axes:
        # default positions in Y, X, Z plane
        default_positions = [
            [0.85, 0],  # pseudo-X axis
            [0, 0.95],  # pseudo-Y axis
            [0.6, 0.6], # pseudo-Z axis
        ]

        labelPos = {}
        for i, rep_str in enumerate(initial_axes):
            if i >= 3:
                break  # only support 3 initial axes
            key, dim, _ = parse_rep(rep_str)
            if dim is None:
                dim = 0  # default to first dimension if not specified
            col_name = format_rep_column(key, dim + 1)
            labelPos[col_name] = default_positions[i]

        state["labelPos"] = labelPos

    if scale is None:
        s = X_df.std() * 4
        scale_factors = [s] if isinstance(s, (float, int)) else s.tolist()
    else:
        scale_factors = scale

    return Langevitour(
        X_df,
        group=group,
        scale=scale_factors,
        point_size=point_size,
        state=state,
        **kwargs,
    )

def resolve_use_reps(
    adata: AnnData,
    use_reps: Optional[Sequence[str]],
    *,
    default_rep: str = "X_pca",
    default_n_dims: int = 10,
) -> pd.DataFrame:
    """
    Resolve and concatenate representations from `adata.obsm`
    into a single DataFrame.

    See: https://logarithmic.net/langevitour/

    Additional example: https://colab.research.google.com/github/pfh/langevitour/blob/main/py/examples/langevitour.ipynb

    Parameters
    ----------
    adata
        AnnData object.
    use_reps
        Sequence like ["X_pca[:10]", "X_umap[:2]", "X_pca[0]"].
        If None, defaults to [f"{default_rep}[:{default_n_dims}]"].
    default_rep
        Representation to use if `use_reps` is None.
    default_n_dims
        Number of dimensions for default representation.

    Returns
    -------
    pd.DataFrame
        Concatenated feature space with named columns.
    """
    if use_reps is None:
        use_reps = [f"{default_rep}[:{default_n_dims}]"]

    dfs = []

    for rep in use_reps:
        key, n_dims, single = parse_rep(rep)

        if key not in adata.obsm:
            raise KeyError(
                f"Representation '{key}' not found in adata.obsm. "
                f"Available: {list(adata.obsm.keys())}"
            )

        X: Matrix = adata.obsm[key]
        if n_dims is not None:
            if single:
                # pick single column
                X = X[:, n_dims:n_dims + 1]
            else:
                # [:N] slice
                X = X[:, :n_dims]

        # n_dims is None → take all columns

        n_actual = X.shape[1]

        cols = [format_rep_column(key, i + 1) for i in range(n_actual)]
        df = pd.DataFrame(X, index=adata.obs_names, columns=pd.Index(cols))
        dfs.append(df)

    return pd.concat(dfs, axis=1)

def format_rep_column(rep_key: str, i: int) -> str:
    """
    Format a representation column name for display.

    Examples
    --------
    X_pca, 1   -> PC1   (special case)
    X_umap, 2  -> UMAP2 (uppercased)
    X_foo, 3   -> FOO3  (uppercased)
    """
    if rep_key.startswith("X_"):
        rep_key = rep_key[2:]

    rep_label = rep_key.upper()

    if rep_label == "PCA":
        rep_label = "PC"

    return f"{rep_label}{i}"

_REP_RE = re.compile(
    r"""
    ^
    (?P<key>[A-Za-z0-9_]+)          # obsm key, e.g. X_pca
    (?:\[\s*(?P<start>\d*)          # optional [start
        (?::(?P<stop>\d+))?         # optional :stop
    \s*\])?
    $
    """,
    re.VERBOSE,
)


def parse_rep(rep: str) -> Tuple[str, Optional[int], bool]:
    """
    Parse representation spec like:
      - 'X_umap'   -> all columns
      - 'X_pca[:2]' -> first 2 columns
      - 'X_pca[0]'  -> single column

    Returns
    -------
    key: str
        obsm key
    n: Optional[int]
        number of columns to take if a slice (None = take all)
    single: bool
        True if a single column, False if slice or full
    """
    m = _REP_RE.match(rep)
    if not m:
        raise ValueError(
            f"Invalid rep specification '{rep}'. "
            "Expected e.g. 'X_pca[:10]', 'X_umap', or 'X_pca[0]'."
        )

    key = m.group("key")
    start = m.group("start")
    stop = m.group("stop")

    if stop is not None:
        # [:N] slice case
        return key, int(stop), False
    elif start:
        # single column
        return key, int(start), True
    else:
        # no slice, take all
        return key, None, False
