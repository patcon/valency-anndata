import math
import platform
import subprocess
import svgwrite
import tempfile
from anndata import AnnData
from typing import TYPE_CHECKING, Any


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

# ------------------------------------------------------------
# Environment detection
# ------------------------------------------------------------
def in_notebook() -> bool:
    import sys
    return "ipykernel" in sys.modules


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

    # -------------------
    # Determine diff sets
    # -------------------
    obs_color_map = {}
    var_color_map = {}
    if diff_from is not None:
        obs_prev = set(diff_from.obs.keys())
        obs_now = set(adata.obs.keys())
        for key in obs_now:
            if key not in obs_prev:
                obs_color_map[key] = "green"  # added
        for key in obs_prev:
            if key not in obs_now:
                obs_color_map[key] = "red"    # removed

        var_prev = set(diff_from.var.keys())
        var_now = set(adata.var.keys())
        for key in var_now:
            if key not in var_prev:
                var_color_map[key] = "green"
        for key in var_prev:
            if key not in var_now:
                var_color_map[key] = "red"

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
    var_keys = [k for k in adata.var]
    var_block_height = max(60, len(var_keys) * line_height)

    # -------------------
    # Obs block width
    # -------------------
    obs_keys = [k for k in adata.obs]
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
        color = obs_color_map.get(key, "black")
        dwg.add(
            dwg.text(
                key,
                insert=(x, baseline_y),
                font_size=font_size,
                font_family="sans-serif",
                text_anchor="start",
                fill=color,
                transform=f"rotate(-45,{x},{baseline_y})",
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
        color = var_color_map.get(key, "black")
        dwg.add(
            dwg.text(
                key,
                insert=(x0 - 10, y),
                font_size=12,
                font_family="sans-serif",
                text_anchor="end",
                fill=color,
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

def schematic_diagram(adata: AnnData, diff_from: AnnData | None = None):
    dwg = adata_structure_svg(adata, diff_from=diff_from)
    _show_svg(dwg)