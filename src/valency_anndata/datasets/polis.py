import asyncio
import re
import anndata as ad
import pandas as pd
from dataclasses import dataclass
from googletrans import Translator
from io import StringIO
from polis_client import PolisClient
from typing import List, Optional
from urllib.parse import urlparse

from ..preprocessing import rebuild_vote_matrix
from ..utils import run_async


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


def load(source: str, *, translate_to: Optional[str] = None, build_X: bool = True) -> ad.AnnData:
    """
    Load a Polis conversation or report into an AnnData object.

    This function accepts either a URL or an ID for a Polis conversation or report,
    fetches raw vote events and statements via the Polis API or CSV export, and
    optionally constructs a participant × statement vote matrix in `adata.X`.

    Parameters
    ----------
    source : str
        The Polis source to load. Supported formats include:
        - Full report URL: `https://pol.is/report/<report_id>`
        - Conversation URL: `https://pol.is/<conversation_id>`
        - Custom host URLs: `https://<host>/report/<report_id>` or `https://<host>/<conversation_id>`
        - Bare IDs:
            - Conversation ID (starts with a digit), e.g., `4asymkcrjf`
            - Report ID (starts with 'r'), e.g., `r4zdxrdscmukmkakmbz3k`
        The function will automatically parse the source to determine whether
        it refers to a conversation or report and fetch the appropriate data.


    translate_to : str or None, optional
        Target language code (e.g., "en", "fr", "es") for translating statement text.
        If provided, the original statement text in `adata.uns["statements"]["txt"]`
        is translated and stored in `adata.var["content"]`. The `adata.var["language_current"]`
        field is updated to the target language, and `adata.var["is_translated"]` is set to True.
        Defaults to None (no translation).

    build_X : bool, default True
        If True, constructs a participant × statement vote matrix from the raw votes
        using `rebuild_vote_matrix()`. This populates `adata.X`, `adata.obs`, and `adata.var`.
        After the first build, a snapshot of this initial matrix is stored in `adata.raw`.

    Returns
    -------
    adata : anndata.AnnData
        An AnnData object containing the loaded Polis data.

        Attributes populated:
        - `adata.uns["votes"]` : pd.DataFrame
            Raw vote events fetched from the API or CSV export.
        - `adata.uns["votes_meta"]` : dict
            Metadata about the sources of votes, e.g., API vs CSV.
        - `adata.uns["statements"]` : pd.DataFrame
            Raw statements/comments for the conversation.
        - `adata.uns["statements_meta"]` : dict
            Metadata about the statements source.
        - `adata.uns["source"]` : dict
            Basic information about the Polis source (base URL, conversation ID, report ID).
        - `adata.uns["schema"]` : dict
            High-level description of `X` and `votes`.
        - `adata.X` : np.ndarray (if `build_X=True`)
            Participant × statement vote matrix (rows = participants, columns = statements).
        - `adata.obs` : pd.DataFrame (if `build_X=True`)
            Participant metadata (index = voter IDs).
        - `adata.var` : pd.DataFrame (if `build_X=True`)
            Statement metadata (index = statement IDs).
        - `adata.raw` : anndata.AnnData (if `build_X=True`)
            Snapshot of the first vote matrix and associated metadata. This allows
            downstream filtering or processing without losing the original vote matrix.

    Notes
    -----
    - If `build_X=False`, only `adata.uns` will be populated, containing the raw
      votes and statements, and `.X`, `.obs`, `.var`, and `.raw` will remain empty.
    - `adata.raw` is assigned only after the first vote matrix build and is intended
      to be immutable.
    - If `translate_to` is provided, `adata.var["content"]` is updated with translated
    text and `adata.var["language_current"]` is set to the target language.
    - The vote matrix is derived from the most recent votes per participant per statement,
      sorted by timestamp.
    """
    adata = _load_raw_polis_data(source)

    if build_X:
        rebuild_vote_matrix(adata, trim_rule=1.0, inplace=True)
        adata.raw = adata.copy()

    _populate_var_statements(adata, translate_to=translate_to)

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

def _populate_var_statements(adata, translate_to: Optional[str] = None):
    statements_aligned = adata.uns["statements"].copy()
    statements_aligned.index = statements_aligned.index.astype(str)
    statements_aligned = statements_aligned.reindex(adata.var_names)

    adata.var["content"] = statements_aligned["txt"]
    adata.var["participant_id_authored"] = statements_aligned["pid"]
    adata.var["created_date"] = statements_aligned["created"]
    adata.var["is_seed"] = statements_aligned["is_seed"]
    adata.var["is_meta"] = statements_aligned["is_meta"]
    adata.var["moderation_state"] = statements_aligned["mod"]
    adata.var["language_original"] = statements_aligned["lang"]

    adata.var["language_current"] = adata.var["language_original"]
    adata.var["is_translated"] = False

    if translate_to is not None:
        translate_statements(adata, translate_to=translate_to, inplace=True)

# private async function
async def _translate_texts_async(texts: List[str], dest_lang: str) -> List[str]:
    async with Translator() as translator:
        results = await asyncio.gather(
            *(translator.translate(t, dest=dest_lang) for t in texts)
        )
    return [r.text for r in results]

# public synchronous wrapper
def translate_texts(texts: List[str], dest_lang: str) -> List[str]:
    return run_async(_translate_texts_async(texts, dest_lang))

def translate_statements(
    adata: ad.AnnData,
    translate_to: str,
    inplace: bool = True
) -> Optional[list[str]]:
    """
    Translate statements in `adata.uns['statements']['txt']` into another language,
    and optionally update `adata.var['content']` and `adata.var['language_current']`.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing `uns['statements']` and `var_names`.
    translate_to : str
        Target language code (e.g., "en", "fr", "es").
    inplace : bool, default True
        If True, updates `adata.var['content']` and `adata.var['language_current']`.
        If False, returns a list of translated strings without modifying `adata`.

    Returns
    -------
    translated_texts : list[str] | None
        List of translated texts if `inplace=False`, else None.
    """
    statements_aligned = adata.uns["statements"].copy()
    statements_aligned.index = statements_aligned.index.astype(str)
    statements_aligned = statements_aligned.reindex(adata.var_names)

    original_texts = statements_aligned["txt"].tolist()
    translated_texts = run_async(_translate_texts_async(original_texts, translate_to))

    if inplace:
        adata.var["content"] = translated_texts
        adata.var["language_current"] = translate_to
        adata.var["is_translated"] = True
        return None
    else:
        return translated_texts

__all__ = [
    "load",
    "translate_statements",
]