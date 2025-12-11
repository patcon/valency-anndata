import re
from types import SimpleNamespace
from urllib.parse import urlparse

DEFAULT_BASE = "https://pol.is"

REPORT_RE = re.compile(r"^r[a-z0-9]{15,}$")     # e.g. r4zdxrdscmukmkakmbz3k
CONVO_RE  = re.compile(r"^[0-9][a-z0-9]{8,}$")  # e.g. 4asymkcrjf (starts with digit)

def parse_polis_source(source: str):
    """
    Returns a SimpleNamespace with:
        base_url
        report_id
        conversation_id
    """
    source = source.strip()

    # ───────────────────────────────────────────
    # 1. URL case
    # ───────────────────────────────────────────
    if source.startswith("http://") or source.startswith("https://"):
        url = urlparse(source)
        base_url = f"{url.scheme}://{url.netloc}"

        # normalize path segments
        parts = [p for p in url.path.split("/") if p]

        report_id = None
        conversation_id = None

        if len(parts) == 2 and parts[0] == "report":
            # /report/<report_id>
            report_id = parts[1]
        elif len(parts) == 1:
            # /<conversation_id>
            conversation_id = parts[0]

        return SimpleNamespace(
            base_url=base_url,
            report_id=report_id,
            conversation_id=conversation_id,
        )

    # ───────────────────────────────────────────
    # 2. Bare IDs (conversation or report)
    # ───────────────────────────────────────────
    # Starts with digit → conversation_id
    if CONVO_RE.match(source):
        return SimpleNamespace(
            base_url=DEFAULT_BASE,
            report_id=None,
            conversation_id=source,
        )

    # Starts with "r" → report_id
    if REPORT_RE.match(source):
        return SimpleNamespace(
            base_url=DEFAULT_BASE,
            report_id=source,
            conversation_id=None,
        )

    raise ValueError(f"Unrecognized Polis source format: {source}")
