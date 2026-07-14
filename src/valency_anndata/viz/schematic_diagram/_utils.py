import tempfile
from typing import Optional
import webbrowser

from ._browser import get_default_browser_name


def _display_svg_in_jupyter(svg_text: str, meta_filename: Optional[str] = None) -> bool:
    """
    Attempt to display an SVG string in an IPython/Jupyter notebook.

    Parameters
    ----------
    svg_text : str
        The raw SVG content to display.
    meta_filename : str | None, optional
        Optional filename to attach as metadata in the notebook output cell.

    Returns
    -------
    bool
        True if the SVG was successfully displayed in a notebook, False otherwise.

    Notes
    -----
    - Requires IPython/Jupyter environment. Returns False if not available.
    - Metadata is passed as {"filename": meta_filename} only if meta_filename is not None.
    """
    try:
        from IPython import get_ipython # type: ignore[reportMissingImports]
        from IPython.display import SVG, display # type: ignore[reportMissingImports]
    except ImportError:
        return False

    ip = get_ipython()
    if ip is None:
        return False

    metadata = {"filename": meta_filename} if meta_filename else {}
    display(SVG(svg_text), metadata=metadata)
    return True

def _display_svg_in_marimo(svg_text: str) -> bool:
    """
    Attempt to display an SVG string inline in a marimo notebook.

    Parameters
    ----------
    svg_text : str
        The raw SVG content to display.

    Returns
    -------
    bool
        True if the SVG was successfully displayed in a marimo notebook,
        False otherwise.

    Notes
    -----
    - Requires a running marimo kernel. Returns False if marimo is not
      installed or not currently running a notebook.
    """
    try:
        import marimo
    except ImportError:
        return False

    if not marimo.running_in_notebook():
        return False

    marimo.output.replace(marimo.Html(svg_text))
    return True

def _show_svg(dwg, filename: Optional[str] = None) -> Optional[str]:
    """
    Display an SVG drawing in the best available environment.

    The function first attempts to render the SVG inline in a marimo or
    IPython/Jupyter notebook. If that fails, it falls back to writing the
    SVG to a temporary file and opening it in a web browser.

    Parameters
    ----------
    dwg
        An object providing a ``tostring() -> str`` method returning SVG text
        (e.g. an svgwrite.Drawing).
    filename : str | None, optional
        Optional filename metadata for notebook display, or a hint for users
        when opened externally.

    Returns
    -------
    str | None
        The path to the temporary SVG file if opened in a browser, otherwise
        None if displayed inline.
    """
    svg_text = dwg.tostring()

    # Try inline notebook display first (marimo, then IPython/Jupyter)
    if _display_svg_in_marimo(svg_text):
        return None
    if _display_svg_in_jupyter(svg_text, meta_filename=filename):
        return None

    # Derive a friendly prefix from the provided filename
    if filename:
        import os
        base = os.path.basename(filename)
        stem, _ = os.path.splitext(base)
        prefix = f"{stem}-"
    else:
        prefix = ""

    # Otherwise, assume script mode → open in default browser
    with tempfile.NamedTemporaryFile(prefix=prefix, suffix=".svg", delete=False, mode="w") as f:
        f.write(svg_text)
        temp_svg_path = f.name

    browser = get_default_browser_name()
    if browser:
        webbrowser.get(browser).open(f"file://{temp_svg_path}")
        print(f"Opened in your default browser ({browser})")
    else:
        webbrowser.open(f"file://{temp_svg_path}")
        print("Opened in your browser")

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
