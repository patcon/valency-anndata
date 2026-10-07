"""Tests for val.viz.voter_vignette_browser (data prep for the v2 anywidget)."""

from pathlib import Path

import pandas as pd
import pytest

import valency_anndata as val
from valency_anndata.viz._voter_vignette_v2 import (
    VoterVignetteWidget,
    _prepare_statements,
    _to_epoch_ms,
    _user_payload,
)

_REAL_FIXTURE = Path(__file__).parent / "fixtures" / "polis_real"


@pytest.fixture(scope="module")
def real_adata():
    return val.datasets.load(str(_REAL_FIXTURE))


@pytest.fixture(scope="module")
def widget(real_adata):
    return VoterVignetteWidget(real_adata, user_id="0")


def test_to_epoch_ms_detects_seconds_and_ms():
    seconds = pd.Series([1_632_834_136, 1_632_834_200])
    assert _to_epoch_ms(seconds).tolist() == [1_632_834_136_000, 1_632_834_200_000]
    ms = pd.Series([1_632_834_136_000, 1_632_834_200_000])
    assert _to_epoch_ms(ms).tolist() == ms.tolist()


def test_bundle_is_built():
    static = Path(val.viz.__file__).parent / "static"
    for name in ("voter_vignette.js", "voter_vignette.css"):
        assert (static / name).stat().st_size > 0, f"run `make js` to build {name}"


def test_payload_joins_statement_content(real_adata, widget):
    payload = widget.user_data
    assert payload["user_id"] == "0"
    assert payload["votes"], "participant 0 should have votes in the fixture"

    var = real_adata.var
    for vote in payload["votes"]:
        assert vote["vote"] in (-1, 0, 1)
        assert vote["content"] == var.loc[vote["statement_id"], "content"]

    times = [v["t"] for v in payload["votes"]]
    assert times == sorted(times)
    # Epoch milliseconds, not seconds.
    assert all(t > 1e12 for t in times)


def test_payload_authored_statements(real_adata, widget):
    authored = real_adata.var[
        real_adata.var["participant_id_authored"].astype(str) == "0"
    ]
    payload = widget.user_data
    assert {s["statement_id"] for s in payload["statements"]} == set(authored.index)
    # Fixture created_date is in seconds; payload must be in ms.
    assert all(s["t"] > 1e12 for s in payload["statements"])


def test_payload_unknown_user_is_empty(real_adata):
    votes_df = real_adata.uns["votes"].assign(
        **{"voter-id": lambda d: d["voter-id"].astype(str), "t_ms": 0}
    )
    payload = _user_payload(votes_df, _prepare_statements(real_adata), "no-such-user")
    assert payload == {"user_id": "no-such-user", "votes": [], "statements": []}


def test_changing_user_id_updates_payload(widget):
    other = next(u for u in widget.voters if u != widget.user_id)
    widget.user_id = other
    assert widget.user_data["user_id"] == other
    widget.user_id = "0"


def test_user_ids_sorted_numerically(widget):
    ids = [int(u) for u in widget.all_users]
    assert ids == sorted(ids)


def test_user_counts_align_with_users(real_adata, widget):
    n = len(widget.all_users)
    assert len(widget.user_vote_counts) == len(widget.user_statement_counts) == n
    assert sum(widget.user_vote_counts) == len(real_adata.uns["votes"])
    assert sum(widget.user_statement_counts) == (
        real_adata.var["participant_id_authored"].notna().sum()
    )


def test_user_lists_are_synced(widget):
    assert widget.commenters and set(widget.commenters) <= set(widget.all_users)
    assert set(widget.voters) <= set(widget.all_users)


def test_variant_v2_returns_widget(real_adata):
    w = val.viz.voter_vignette_browser(real_adata, variant="v2")
    assert isinstance(w, VoterVignetteWidget)
    assert w.user_id in w.commenters


def test_default_variant_is_v2(real_adata):
    assert isinstance(val.viz.voter_vignette_browser(real_adata), VoterVignetteWidget)


def test_unknown_variant_raises(real_adata):
    with pytest.raises(ValueError, match="variant"):
        val.viz.voter_vignette_browser(real_adata, variant="v3")
