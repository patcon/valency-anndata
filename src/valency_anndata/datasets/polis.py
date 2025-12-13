from io import StringIO
import anndata as ad
import pandas as pd
from polis_client import PolisClient

import re
from types import SimpleNamespace
from urllib.parse import urlparse

DEFAULT_BASE = "https://pol.is"

REPORT_RE = re.compile(r"^r[a-z0-9]{15,}$")     # e.g. r4zdxrdscmukmkakmbz3k
CONVO_RE  = re.compile(r"^[0-9][a-z0-9]{8,}$")  # e.g. 4asymkcrjf (starts with digit)

def _parse_polis_source(source: str):
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

def polis(source: str):
    """
    Accepts any of the following:
    - https://pol.is/report/<report_id>
    - https://pol.is/<conversation_id>
    - https://<host>/report/<report_id>
    - https://<host>/<conversation_id>
    - <conversation_id> (starts with digit; base_url=pol.is)
    - <report_id>       (starts with 'r'; base_url=pol.is)

    Returns:
        AnnData object
    """
    adata = ad.AnnData(X=None)

    convo_meta = _parse_polis_source(source)
    client = PolisClient(base_url=convo_meta.base_url)
    # client = PolisClient(base_url=convo_meta.base_url, xid="foobar")

    if convo_meta.report_id:
        votes_csv_text = client.get_export_file(
            filename="votes.csv",
            report_id=convo_meta.report_id,
        )
        votes_csv_df = pd.read_csv(StringIO(votes_csv_text))
        votes_csv_df.sort_values(by="timestamp", ascending=True, inplace=True)

        adata.uns["votes_csv"] = votes_csv_df

        report = client.get_report(report_id=convo_meta.report_id)
        assert report is not None
        convo_meta.conversation_id = report.conversation_id

    elif convo_meta.conversation_id:
        votes_list = client.get_all_votes_slow(conversation_id=convo_meta.conversation_id)
        votes_api_df = pd.DataFrame(votes_list)
        votes_api_df.sort_values(by="modified", ascending=True, inplace=True)

        adata.uns["votes_api"] = votes_api_df

    else:
        raise

    statements = client.get_comments(conversation_id=convo_meta.conversation_id)
    assert statements is not None
    adata.uns["statements"] = pd.DataFrame([s.to_dict() for s in statements])

    # if convo_meta.conversation_id:
    #     xids = client.get_xids(conversation_id=convo_meta.conversation_id)
    #     adata.uns["xids"] = pd.DataFrame(xids)

    return adata
