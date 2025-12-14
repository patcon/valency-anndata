import re
import anndata as ad
import pandas as pd
from dataclasses import dataclass
from io import StringIO
from polis_client import PolisClient
from urllib.parse import urlparse
from ..preprocessing import rebuild_vote_matrix


DEFAULT_BASE = "https://pol.is"

REPORT_RE = re.compile(r"^r[a-z0-9]{15,}$")     # e.g. r4zdxrdscmukmkakmbz3k
CONVO_RE  = re.compile(r"^[0-9][a-z0-9]{8,}$")  # e.g. 4asymkcrjf (starts with digit)

@dataclass
class PolisSource:
    base_url: str
    conversation_id: str | None = None
    report_id: str | None = None

def _parse_polis_source(source: str):
    """
    Returns a PolisSource with:
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

        return PolisSource(
            base_url=base_url,
            report_id=report_id,
            conversation_id=conversation_id,
        )

    # ───────────────────────────────────────────
    # 2. Bare IDs (conversation or report)
    # ───────────────────────────────────────────
    # Starts with digit → conversation_id
    if CONVO_RE.match(source):
        return PolisSource(
            base_url=DEFAULT_BASE,
            report_id=None,
            conversation_id=source,
        )

    # Starts with "r" → report_id
    if REPORT_RE.match(source):
        return PolisSource(
            base_url=DEFAULT_BASE,
            report_id=source,
            conversation_id=None,
        )

    raise ValueError(f"Unrecognized Polis source format: {source}")


def polis(source: str, *, build_X: bool = True) -> ad.AnnData:
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
    adata = _load_raw_polis_data(source)

    if build_X:
        rebuild_vote_matrix(adata, trim_rule=1.0, inplace=True)

    # if convo_meta.conversation_id:
    #     xids = client.get_xids(conversation_id=convo_meta.conversation_id)
    #     adata.uns["xids"] = pd.DataFrame(xids)

    return adata


def _load_raw_polis_data(source):
    adata = ad.AnnData()

    convo_src = _parse_polis_source(source)
    client = PolisClient(base_url=convo_src.base_url)
    # client = PolisClient(base_url=convo_meta.base_url, xid="foobar")

    vote_frames = []
    vote_sources = {}

    # ───────────────────────────────────────────
    # Load votes
    # ───────────────────────────────────────────
    if convo_src.report_id:
        votes_csv_text = client.get_export_file(
            filename="votes.csv",
            report_id=convo_src.report_id,
        )
        df = pd.read_csv(StringIO(votes_csv_text))
        df["source"] = "csv"
        df["source_id"] = convo_src.report_id
        df.sort_values(by="timestamp", ascending=True, inplace=True)

        vote_frames.append(df)

        vote_sources["csv"] = {
            "via": "live_csv",
            "report_id": convo_src.report_id,
            "base_url": convo_src.base_url,
            "retrieved_at": pd.Timestamp.utcnow().isoformat(),
        }

        report = client.get_report(report_id=convo_src.report_id)
        assert report is not None
        convo_src.conversation_id = report.conversation_id or None

    elif convo_src.conversation_id:
        votes_list = client.get_all_votes_slow(conversation_id=convo_src.conversation_id)
        df = pd.DataFrame(votes_list)
        df["source"] = "api"
        df["source_id"] = convo_src.conversation_id
        df.rename(columns={"modified": "timestamp"}, inplace=True)

        vote_frames.append(df)

        vote_sources["api"] = {
            "via": "live_api",
            "conversation_id": convo_src.conversation_id,
            "base_url": convo_src.base_url,
            "retrieved_at": pd.Timestamp.utcnow().isoformat(),
        }
        df.sort_values(by="timestamp", ascending=True, inplace=True)

    if not vote_frames:
        raise ValueError("No votes could be loaded")

    votes = (
        pd.concat(vote_frames, ignore_index=True)
          .sort_values("timestamp", kind="stable")
          .reset_index(drop=True)
    )

    adata.uns["votes"] = votes
    adata.uns["votes_meta"] = {
        "sources": vote_sources,
        "sorted_by": "timestamp",
    }

    statements = client.get_comments(conversation_id=convo_src.conversation_id)
    assert statements is not None

    adata.uns["statements"] = (
        pd.DataFrame([s.to_dict() for s in statements])
            .set_index("tid", drop=False)
    )

    adata.uns["statements_meta"] = {
        "source": {
            "via": "live_api",
            "conversation_id": convo_src.conversation_id,
            "base_url": convo_src.base_url,
            "retrieved_at": pd.Timestamp.utcnow().isoformat(),
        },
    }

    pd.Timestamp.utcnow().isoformat()

    adata.uns["source"] = {
        "base_url": convo_src.base_url,
        "conversation_id": convo_src.conversation_id,
        "report_id": convo_src.report_id,
    }

    adata.uns["schema"] = {
        "X": "participant × statement vote matrix (derived)",
        "votes": "raw vote events",
    }

    return adata