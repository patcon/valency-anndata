import re
from typing import Optional, Sequence

import scanpy as sc


_COLOR_SPEC_RE = re.compile(
    r"""
    ^
    (?P<key>[A-Za-z_][A-Za-z0-9_]*)
    \[
        \s*
        (?:
            (?P<index>\d+)
            |
            (?P<start>\d*)\s*:\s*(?P<stop>\d*)
        )
        \s*
    \]
    $
    """,
    re.VERBOSE,
)


def _parse_color_spec(color: str) -> Optional[tuple[str, int, Optional[int]]]:
    if "[" not in color and "]" not in color:
        return None

    m = _COLOR_SPEC_RE.match(color)
    if not m:
        raise ValueError(
            f"Invalid embedding color spec '{color}'. "
            "Expected 'foo[i]', 'foo[a:b]', or 'foo[:b]'."
        )

    key = m.group("key")
    index = m.group("index")

    if index is not None:
        return key, int(index), None

    start = int(m.group("start") or 0)
    stop_text = m.group("stop")

    if stop_text == "":
        raise ValueError(
            f"Invalid embedding color slice '{color}'. "
            "Expected 'foo[a:b]' or 'foo[:b]' with stop > start."
        )

    stop = int(stop_text)
    if stop <= start:
        raise ValueError(
            f"Invalid embedding color slice '{color}'. "
            "Expected 'foo[a:b]' or 'foo[:b]' with stop > start."
        )

    return key, start, stop


def _expand_color_spec(color: str):
    parsed = _parse_color_spec(color)
    if parsed is None:
        return color

    key, start, stop = parsed
    if stop is None:
        return color

    return [f"{key}[{i}]" for i in range(start, stop)]


def _expand_color_specs(color: Optional[str | Sequence[str]]):
    if color is None:
        return None

    if isinstance(color, str):
        return _expand_color_spec(color)

    expanded = []
    for item in color:
        item_expanded = _expand_color_spec(item)
        if isinstance(item_expanded, list):
            expanded.extend(item_expanded)
        else:
            expanded.append(item_expanded)

    return expanded


def _rewrite_color(adata, color):
    expanded = _expand_color_specs(color)
    if expanded is None:
        return adata, None

    colors = [expanded] if isinstance(expanded, str) else list(expanded)
    forwarded = []
    adata_plot = adata

    for color_name in colors:
        parsed = _parse_color_spec(color_name)
        if parsed is None:
            forwarded.append(color_name)
            continue

        key, index, stop = parsed
        if stop is not None:
            raise ValueError(
                f"Expected expanded single-column color spec, got '{color_name}'."
            )

        if key not in adata.obsm:
            raise KeyError(
                f"Color spec '{color_name}' references missing adata.obsm['{key}']."
            )

        X = adata.obsm[key]
        if index >= X.shape[1]:
            raise IndexError(
                f"Color spec '{color_name}' is out of bounds for "
                f"adata.obsm['{key}'] with {X.shape[1]} columns."
            )

        if adata_plot is adata:
            adata_plot = adata.copy()

        adata_plot.obs[color_name] = adata_plot.obsm[key][:, index]
        forwarded.append(color_name)

    if isinstance(expanded, str):
        return adata_plot, forwarded[0]

    return adata_plot, forwarded


def embedding(adata, *args, color=None, **kwargs):
    adata_plot, color = _rewrite_color(adata, color)
    return sc.pl.embedding(adata_plot, *args, color=color, **kwargs)
