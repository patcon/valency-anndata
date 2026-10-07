"""anywidget-based voter vignette browser (``variant="v2"``).

The frontend lives in ``js/src/voter_vignette.ts`` and is built with Vite into
``static/`` (run ``make js``).
"""

import random
from pathlib import Path

import anywidget
import numpy as np
import pandas as pd
import traitlets
from anndata import AnnData

from ._voter_vignette import _prepare_votes_df, _user_lists

_STATIC = Path(__file__).parent / "static"


def _to_epoch_ms(values: pd.Series) -> pd.Series:
    """Coerce epoch timestamps to milliseconds, detecting seconds vs ms by magnitude."""
    numeric = pd.to_numeric(values, errors="coerce")
    median = numeric.median()
    if pd.notna(median) and median < 1e11:  # looks like seconds
        numeric = numeric * 1000
    return numeric


def _none_if_na(value):
    return None if pd.isna(value) else value


def _user_payload(
    votes_df: pd.DataFrame,
    statements: pd.DataFrame,
    user_id: str,
) -> dict:
    """
    Build the JSON-serializable data the frontend renders for one participant.

    Parameters
    ----------
    votes_df:
        Votes for (at least) this user, with columns ``voter-id`` (str),
        ``comment-id``, ``vote`` and ``t_ms`` (epoch milliseconds).
    statements:
        Statement table indexed by statement id (str), with columns ``content``,
        ``participant_id_authored``, ``moderation_state`` and ``t_ms``.
    user_id:
        Participant to build the payload for.
    """
    user_id = str(user_id)
    user_votes = votes_df[votes_df["voter-id"] == user_id].sort_values("t_ms")
    voted_ids = user_votes["comment-id"].astype(str)
    voted_statements = statements.reindex(voted_ids)

    votes = [
        {
            "t": int(t),
            "vote": int(vote),
            "statement_id": sid,
            "content": _none_if_na(content),
            "moderation_state": None if pd.isna(mod) else int(mod),
        }
        for t, vote, sid, content, mod in zip(
            user_votes["t_ms"],
            user_votes["vote"],
            voted_ids,
            voted_statements["content"],
            voted_statements["moderation_state"],
        )
    ]

    authored = statements[statements["participant_id_authored"] == user_id]
    authored = authored[authored["t_ms"].notna()].sort_values("t_ms")
    authored_list = [
        {
            "t": int(t),
            "statement_id": str(sid),
            "content": _none_if_na(content),
            "moderation_state": None if pd.isna(mod) else int(mod),
        }
        for sid, t, content, mod in zip(
            authored.index,
            authored["t_ms"],
            authored["content"],
            authored["moderation_state"],
        )
    ]

    return {
        "user_id": user_id,
        "votes": votes,
        "statements": authored_list,
    }


def _prepare_statements(adata: AnnData) -> pd.DataFrame:
    var = adata.var
    statements = pd.DataFrame(index=var.index.astype(str))
    statements["content"] = var["content"].to_numpy() if "content" in var else None
    statements["participant_id_authored"] = (
        var["participant_id_authored"].astype("string").to_numpy()
    )
    if "moderation_state" in var:
        statements["moderation_state"] = pd.to_numeric(
            var["moderation_state"], errors="coerce"
        ).to_numpy()
    else:
        statements["moderation_state"] = np.nan
    statements["t_ms"] = _to_epoch_ms(var["created_date"]).to_numpy()
    return statements


def _id_sort_key(user_id: str):
    """Sort numeric ids numerically ("2" before "10"), then any others as text."""
    return (0, int(user_id), "") if user_id.isdigit() else (1, 0, user_id)


class VoterVignetteWidget(anywidget.AnyWidget):
    """Zoomable voter timeline. Created via ``voter_vignette_browser(adata, variant="v2")``."""

    _esm = _STATIC / "voter_vignette.js"
    _css = _STATIC / "voter_vignette.css"

    # Parallel arrays (one entry per participant, ids sorted numerically) feeding
    # the user picker's counts, filters and sorting.
    all_users = traitlets.List(traitlets.Unicode()).tag(sync=True)
    user_vote_counts = traitlets.List(traitlets.Int()).tag(sync=True)
    user_statement_counts = traitlets.List(traitlets.Int()).tag(sync=True)
    user_id = traitlets.Unicode("").tag(sync=True)
    user_data = traitlets.Dict().tag(sync=True)

    def __init__(self, adata: AnnData, user_id: str | None = None, **kwargs):
        votes_df = _prepare_votes_df(adata)
        all_voters, all_commenters, all_users = _user_lists(adata, votes_df)
        users = sorted(all_users.tolist(), key=_id_sort_key)
        vote_counts = votes_df["voter-id"].value_counts()
        statement_counts = (
            adata.var["participant_id_authored"].dropna().astype(str).value_counts()
        )

        votes_df["t_ms"] = _to_epoch_ms(votes_df["timestamp"])
        self._votes_df = votes_df[["voter-id", "comment-id", "vote", "t_ms"]]
        # Row positions per voter, so switching users doesn't rescan every vote.
        self._vote_rows = self._votes_df.groupby("voter-id").indices
        self._statements = _prepare_statements(adata)

        if user_id is None:
            pool = all_commenters if len(all_commenters) else all_voters
            user_id = random.choice(pool.tolist()) if len(pool) else ""

        super().__init__(
            all_users=users,
            user_vote_counts=[int(vote_counts.get(u, 0)) for u in users],
            user_statement_counts=[int(statement_counts.get(u, 0)) for u in users],
            **kwargs,
        )
        # Set after init so the observer below fills user_data.
        self.user_id = str(user_id)

    @property
    def voters(self) -> list[str]:
        """Participants with at least one vote."""
        return [u for u, n in zip(self.all_users, self.user_vote_counts) if n]

    @property
    def commenters(self) -> list[str]:
        """Participants who authored at least one statement."""
        return [u for u, n in zip(self.all_users, self.user_statement_counts) if n]

    @traitlets.observe("user_id")
    def _on_user_id(self, change):
        user_id = change["new"]
        rows = self._vote_rows.get(user_id, [])
        self.user_data = _user_payload(
            self._votes_df.iloc[rows], self._statements, user_id
        )
