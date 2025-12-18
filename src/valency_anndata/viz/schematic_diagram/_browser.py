import platform
import subprocess
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
