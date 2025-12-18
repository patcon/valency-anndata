import math
import platform
import subprocess
import svgwrite
import tempfile
from anndata import AnnData
from typing import TYPE_CHECKING, Any, Optional, overload


def get_default_browser_name() -> str | None:
    """
    Best-effort detection of the system default web browser.

    Returns
    -------
    str | None
        A human-readable browser identifier (e.g. 'Chrome', 'Firefox', 'Safari'),
        or None if it cannot be determined.

    Notes
    -----
    - This is inherently OS-specific and fragile.
    - Do NOT rely on this for logic or correctness.
    - Intended for diagnostics, logging, or user-facing messages only.
    """

    system = platform.system()

    # ------------------
    # macOS
    # ------------------
    if system == "Darwin":
        try:
            # Ask LaunchServices for the default HTTP handler
            cmd = [
                "plutil",
                "-p",
                str(
                    subprocess.check_output(
                        [
                            "defaults",
                            "read",
                            "com.apple.LaunchServices/com.apple.launchservices.secure",
                            "LSHandlers",
                        ]
                    )
                ),
            ]
        except Exception:
            pass

        try:
            output = subprocess.check_output(
                [
                    "defaults",
                    "read",
                    "com.apple.LaunchServices/com.apple.launchservices.secure",
                    "LSHandlers",
                ],
                stderr=subprocess.DEVNULL,
            ).decode()

            for block in output.split("},"):
                is_url_handler = (
                    'LSHandlerURLScheme = http' in block
                    or 'LSHandlerURLScheme = https' in block
                )
                if is_url_handler:
                    for line in block.splitlines():
                        if "LSHandlerRoleAll" in line or "LSHandlerRoleViewer" in line:
                            value = line.split("=", 1)[-1].strip().strip('";')

                            # Skip version numbers like "6533.100"
                            if "." not in value or value.replace(".", "").isdigit():
                                continue

                            bundle = value.lower()

                            if "chrome" in bundle:
                                return "Chrome"
                            if "firefox" in bundle:
                                return "Firefox"
                            if "safari" in bundle:
                                return "Safari"
                            if "edge" in bundle:
                                return "Edge"

                            return bundle

        except Exception:
            pass

    # ------------------
    # Windows
    # ------------------
    if system == "Windows":
        try:
            import winreg

            # Suppress type warnings on non-windows platforms
            if TYPE_CHECKING:
                winreg: Any

            with winreg.OpenKey(
                winreg.HKEY_CURRENT_USER,
                r"Software\Microsoft\Windows\Shell\Associations\UrlAssociations\http\UserChoice",
            ) as key:
                progid, _ = winreg.QueryValueEx(key, "ProgId")

            progid = progid.lower()
            if "chrome" in progid:
                return "Chrome"
            if "firefox" in progid:
                return "Firefox"
            if "edge" in progid:
                return "Edge"
            if "safari" in progid:
                return "Safari"

            return progid

        except Exception:
            pass

    # ------------------
    # Linux
    # ------------------
    if system == "Linux":
        try:
            output = subprocess.check_output(
                ["xdg-settings", "get", "default-web-browser"],
                stderr=subprocess.DEVNULL,
            ).decode().strip().lower()

            if "chrome" in output:
                return "Chrome"
            if "firefox" in output:
                return "Firefox"
            if "edge" in output:
                return "Edge"
            if "safari" in output:
                return "Safari"

            return output or None

        except Exception:
            pass

    return None

def diff_text_style(status: str | None) -> dict[str, str]:
    if status == "added":
        return {
            "fill": "green",
            "font_weight": "bold",
        }
    if status == "removed":
        return {
            "fill": "red",
            "font_weight": "bold",
        }
    return {
        "fill": "black",
        "font_weight": "normal",
    }

# ------------------------------------------------------------
# SVG primitives
# ------------------------------------------------------------
def draw_grid_block(
    dwg,
    *,
    x,
    y,
    width,
    height,
    rows,
    cols,
    label,
    stroke="#333",
    grid_stroke="#ccc",
):
    group = dwg.g()

    # Outer rectangle
    group.add(
        dwg.rect(
            insert=(x, y),
            size=(width, height),
            fill="none",
            stroke=stroke,
            stroke_width=2,
        )
    )

    # Horizontal grid
    if rows > 1:
        row_h = height / rows
        for i in range(1, rows):
            group.add(
                dwg.line(
                    start=(x, y + i * row_h),
                    end=(x + width, y + i * row_h),
                    stroke=grid_stroke,
                )
            )

    # Vertical grid
    if cols > 1:
        col_w = width / cols
        for j in range(1, cols):
            group.add(
                dwg.line(
                    start=(x + j * col_w, y),
                    end=(x + j * col_w, y + height),
                    stroke=grid_stroke,
                )
            )

    # Label (centered)
    lines = label.split("\n")
    line_height = 14
    start_y = y + height / 2 - (len(lines) - 1) * line_height / 2

    for i, line in enumerate(lines):
        group.add(
            dwg.text(
                line,
                insert=(x + width / 2, start_y + i * line_height),
                text_anchor="middle",
                font_size=12,
                font_family="sans-serif",
            )
        )

    dwg.add(group)


# ------------------------------------------------------------
# AnnData → SVG with diff
# ------------------------------------------------------------
def adata_structure_svg(adata: AnnData, diff_from: AnnData | None = None):
    cell = 18
    max_cells = 10
    pad = 40
    line_height = 14
    obs_key_spacing = 15  # horizontal spacing between rotated keys

    if diff_from is not None:
        # Start with the original order from diff_from
        obs_keys = list(diff_from.obs.keys())
        # Add any new keys that appear in adata but not in diff_from
        obs_keys += [k for k in adata.obs.keys() if k not in obs_keys]

        # Similarly for var keys
        var_keys = list(diff_from.var.keys())
        var_keys += [k for k in adata.var.keys() if k not in var_keys]
    else:
        obs_keys = list(adata.obs.keys())
        var_keys = list(adata.var.keys())

    # -------------------
    # Determine diff sets
    # -------------------
    obs_status: dict[str, str] = {}
    var_status: dict[str, str] = {}

    if diff_from is not None:
        obs_prev = set(diff_from.obs.keys())
        obs_now = set(adata.obs.keys())
        for key in obs_now - obs_prev:
            obs_status[key] = "added"
        for key in obs_prev - obs_now:
            obs_status[key] = "removed"

        var_prev = set(diff_from.var.keys())
        var_now = set(adata.var.keys())
        for key in var_now - var_prev:
            var_status[key] = "added"
        for key in var_prev - var_now:
            var_status[key] = "removed"

    # -------------------
    # Determine matrix size
    # -------------------
    obs_cells = min(max_cells, math.ceil(math.sqrt(adata.n_obs)))
    var_cells = min(max_cells, math.ceil(math.sqrt(adata.n_vars)))

    X_width = var_cells * cell
    X_height = obs_cells * cell

    # -------------------
    # Var block height
    # -------------------
    var_block_height = max(60, len(var_keys) * line_height)

    # -------------------
    # Obs block width
    # -------------------
    min_obs_width = 60
    needed_obs_width = len(obs_keys) * obs_key_spacing
    obs_width = max(min_obs_width, needed_obs_width)

    font_size = 12
    tilt_factor = 0.707  # sin/cos 45°
    last_key_extra = len(obs_keys[-1]) * (font_size * 0.5) * tilt_factor if obs_keys else 0
    extra_canvas_padding = last_key_extra + 10

    x0 = pad + 120
    y0 = pad + var_block_height + 30
    canvas_width = x0 + X_width + 30 + obs_width + extra_canvas_padding
    canvas_height = X_height + var_block_height + 150

    dwg = svgwrite.Drawing(
        size=(canvas_width, canvas_height),
        profile="full",
    )

    # -------------------
    # X block
    # -------------------
    draw_grid_block(
        dwg,
        x=x0,
        y=y0,
        width=X_width,
        height=X_height,
        rows=obs_cells,
        cols=var_cells,
        label=f"X\n{adata.n_obs} x {adata.n_vars}",
        stroke="#2ecc71",
    )

    # -------------------
    # Obs block (right)
    # -------------------
    draw_grid_block(
        dwg,
        x=x0 + X_width + 30,
        y=y0,
        width=obs_width,
        height=X_height,
        rows=obs_cells,
        cols=1,
        label=f"obs\n{adata.n_obs} x {adata.obs.shape[1]}",
        stroke="#3498db",
    )

    # -------------------
    # Obs keys (rotated 45° above obs block)
    # -------------------
    baseline_y = y0 - 7
    for i, key in enumerate(obs_keys):
        x = x0 + X_width + 30 + 10 + i * obs_key_spacing
        style = diff_text_style(obs_status.get(key))

        dwg.add(
            dwg.text(
                key,
                insert=(x, baseline_y),
                font_size=font_size,
                font_family="sans-serif",
                text_anchor="start",
                transform=f"rotate(-45,{x},{baseline_y})",
                **style,
            )
        )

    # -------------------
    # Var block (top)
    # -------------------
    draw_grid_block(
        dwg,
        x=x0,
        y=y0 - var_block_height - 30,
        width=X_width,
        height=var_block_height,
        rows=1,
        cols=var_cells,
        label=f"var\n{adata.n_vars} x {adata.var.shape[1]}",
        stroke="#9b59b6",
    )

    # -------------------
    # Var keys (left of var block)
    # -------------------
    for i, key in enumerate(var_keys):
        y = y0 - 27 - (len(var_keys) - i) * line_height + line_height / 2
        style = diff_text_style(var_status.get(key))

        dwg.add(
            dwg.text(
                key,
                insert=(x0 - 10, y),
                font_size=12,
                font_family="sans-serif",
                text_anchor="end",
                **style,
            )
        )

    return dwg

def _display_svg_in_notebook(svg_text: str) -> bool:
    """Return True if displayed successfully, False otherwise."""
    try:
        from IPython.display import SVG, display  # type: ignore[reportMissingImports]
        display(SVG(svg_text))
        return True
    except ImportError:
        return False

def _show_svg(dwg):
    svg_text = dwg.tostring()
    if _display_svg_in_notebook(svg_text):
        return

    # Otherwise, assume script mode → open in default browser
    import webbrowser
    with tempfile.NamedTemporaryFile(suffix=".svg", delete=False, mode="w") as f:
        f.write(svg_text)
        temp_svg_path = f.name
    browser = get_default_browser_name()
    if browser:
        webbrowser.get(browser).open(f"file://{temp_svg_path}")
        print(f"Opened in your default browser ({browser})")
    else:
        webbrowser.open(f"file://{temp_svg_path}")
        print("Opened in your browser")

class _SchematicDiagramContext:
    """
    Context manager that captures a snapshot of an AnnData object on entry,
    and automatically renders a diff schematic on exit.
    The snapshot only exists for the duration of the block.
    """

    def __init__(self, adata: AnnData):
        self.adata = adata
        self.snapshot = None

    def __enter__(self):
        # Take a snapshot of the current state
        self.snapshot = self.adata.copy()
        return self.adata

    def __exit__(self, exc_type, exc, tb):
        # Render diff if no exception
        if exc_type is None:
            dwg = adata_structure_svg(self.adata, diff_from=self.snapshot)
            _show_svg(dwg)

        # Forget the snapshot: no global registration
        self.snapshot = None

        # Do not suppress exceptions
        return False

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
    - If ``diff_from`` is provided, labels in ``obs`` and ``var`` are highlighted
      to indicate additions and removals relative to ``diff_from``.
    - If ``diff_from`` is not provided, no diff highlighting will be applied.
    - The diagram is displayed inline (in notebooks) or opened in a browser
      (in script mode).

    Parameters
    ----------
    adata : AnnData, optional
        The AnnData object to visualize.
    diff_from : AnnData, optional
        An optional snapshot to diff against. Explicit diffs always take
        precedence over any snapshot captured via a context manager.

    Returns
    -------
    None
        The diagram is displayed; nothing is returned.

    -------------------------------------------------------------------------
    2. Context-manager mode
    -------------------------------------------------------------------------
    Capture a snapshot of an AnnData object at the beginning of a ``with`` block,
    and automatically render a diff schematic on exit.

    Example
    -------
    >>> with schematic_diagram(diff_from=adata):
    ...     some_transformation(adata)

    In context-manager mode:
    - ``diff_from`` must be provided, ``adata`` must be omitted.
    - A snapshot of ``diff_from`` is taken on entry to the ``with`` block.
    - On exit, a schematic diff between the original snapshot and the
      current state of ``adata`` is automatically rendered.
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
    if adata is None and diff_from is not None:
        return _SchematicDiagramContext(diff_from)

    # ------------------
    # Render mode (explicit or implicit diff)
    # ------------------
    if adata is not None:
        dwg = adata_structure_svg(adata, diff_from=diff_from)
        _show_svg(dwg)
        return None

    raise TypeError("Invalid schematic_diagram() call")