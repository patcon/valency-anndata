from typing import Literal

from anndata import AnnData


def _prepare_votes_df(adata: AnnData):
    """Copy of ``adata.uns["votes"]`` with string voter ids and a unit sanity check."""
    import warnings

    votes_df = adata.uns["votes"].copy()
    votes_df["voter-id"] = votes_df["voter-id"].astype(str)

    # Heuristic check: votes should be seconds, median ~1e9–1e10
    votes_median = votes_df["timestamp"].median()
    if votes_median > 1e11:  # looks too large for seconds
        warnings.warn(
            f"Median timestamp in votes is {votes_median}, which seems too large. "
            "Expected seconds. If these are milliseconds, divide by 1000."
        )
    return votes_df


def _user_lists(adata: AnnData, votes_df):
    """Return (all_voters, all_commenters, all_users) as arrays of string ids."""
    import numpy as np

    all_voters = votes_df["voter-id"].unique()
    all_commenters = adata.var["participant_id_authored"].dropna().astype(str).unique()
    all_users = np.unique(np.concatenate([all_voters, all_commenters]))
    return all_voters, all_commenters, all_users


def voter_vignette_browser(
    adata: AnnData,
    variant: Literal["v1", "v2"] = "v1",
):
    """
    Interactive browser for quickly surveying many voting timelines of random
    participants alongside statements they authored.

    Parameters
    ----------
    adata:
        An AnnData object loaded from a Polis conversation.<br/>
        (See Assumptions below)
    variant:
        Which implementation to render.

        - `"v1"` (default): ipywidgets controls with a static matplotlib plot.
        - `"v2"`: an [anywidget](https://anywidget.dev)-based timeline that you can
          zoom and pan, with hover tooltips showing the statement behind each vote.

    Returns
    -------
    None for `"v1"`. For `"v2"`, the widget object, which Jupyter displays when it
    is the last expression in a cell. Setting `widget.user_id = "<id>"` switches
    the selected participant.

    Assumptions
    -----------

    - Votes are stored in `adata.uns["votes"]` with columns:
        - `voter-id`
        - `vote` (-1, 0, 1)
        - `timestamp` (seconds since epoch)

    - Statements are stored in `adata.var` with columns:
        - `participant_id_authored`
        - `created_date` (milliseconds since epoch for `"v1"`; `"v2"` accepts
          seconds or milliseconds)
        - `content`
        - `moderation_state` (optional, -1/0/1)

    Behavior
    --------

    - Renders a dropdown to select a user, with buttons for random voter or commenter.
    - Plots votes over time with colors (red/neutral/green).
    - Draws vertical bars for authored statements with moderation-state coloring.
    - Displays statements below the plot in submission order.
    - Warns if vote or statement timestamps appear out of expected ranges.

    With `variant="v2"`, additionally:

    - Scroll to zoom and drag to pan the timeline, down to second-level resolution.
      Double-click or press "Reset zoom" to return to the full range.
    - Hover a vote marker to see the statement voted on, the vote, and its time.
    - Click an authored statement in the list to zoom to the votes around it.

    Examples
    --------

    ```py
    adata = val.datasets.polis.load("https://pol.is/report/r29kkytnipymd3exbynkd", translate_to="en")

    val.viz.voter_vignette_browser(adata)

    # Zoomable anywidget version
    val.viz.voter_vignette_browser(adata, variant="v2")
    ```
    <img src="../../assets/documentation-examples/viz--voter-vignette-browser.png">
    """
    if variant == "v2":
        from ._voter_vignette_v2 import VoterVignetteWidget

        return VoterVignetteWidget(adata)
    if variant != "v1":
        raise ValueError(f"Unknown variant {variant!r}; expected 'v1' or 'v2'.")

    import random
    import pandas as pd
    import matplotlib.pyplot as plt
    from ipywidgets import widgets
    from IPython.display import display, Markdown
    import warnings

    # -----------------------------
    # Prepare votes dataframe
    # -----------------------------
    votes_df = _prepare_votes_df(adata)
    votes_df["timestamp"] = pd.to_datetime(votes_df["timestamp"], unit="s")

    # -----------------------------
    # Core plotting function
    # -----------------------------
    def plot_user_activity(user_id: str):
        user_id = str(user_id)

        # --- Votes ---
        user_votes = votes_df[votes_df["voter-id"] == user_id]
        n_votes = len(user_votes)

        if user_votes.empty:
            first_vote = last_vote = None
            delta = pd.Timedelta(0)
        else:
            first_vote = user_votes["timestamp"].min()
            last_vote = user_votes["timestamp"].max()
            delta = last_vote - first_vote

        # Adaptive duration string
        if delta < pd.Timedelta(hours=1):
            duration_str = f"{delta.total_seconds()/60:.1f} minutes"
        elif delta < pd.Timedelta(days=1):
            duration_str = f"{delta.total_seconds()/3600:.1f} hours"
        else:
            duration_str = f"{delta.days} days"

        plt.figure(figsize=(12, 4))

        if not user_votes.empty:
            vote_colors = {-1: "red", 0: "gold", 1: "green"}
            colors = user_votes["vote"].map(vote_colors)
            plt.scatter(user_votes["timestamp"], user_votes["vote"], c=colors, s=50)

        plt.yticks([-1, 0, 1], ["Disagree", "Pass", "Agree"])
        plt.xlabel("Time")
        plt.ylabel("Vote")

        # --- Statements ---
        statements = adata.var
        user_statements = statements[
            statements["participant_id_authored"].astype(str) == user_id
        ]
        n_statements = len(user_statements)

        if not user_statements.empty:
            created_ms = pd.to_numeric(
                user_statements["created_date"], errors="coerce"
            )

            # Heuristic check: statements in milliseconds
            statements_median = created_ms.median()
            if statements_median < 1e11 or statements_median > 1e14:
                warnings.warn(
                    f"Median created_date in statements is {statements_median}. "
                    "Expected milliseconds."
                )

            statement_times = pd.to_datetime(created_ms, unit="ms")

            # Map moderation_state for plotting vertical bars
            mod_colors = {1: "green", 0: "gray", -1: "red"}
            moderation_states = (
                user_statements.get("moderation_state", 0)
                .fillna(0)
                .astype(int)
            )

            for t, mod in zip(statement_times, moderation_states):
                plt.axvline(x=t, color=mod_colors.get(mod, "gray"), lw=2, alpha=0.7)

        # --- Legend proxies ---
        vote_proxy = plt.Line2D(
            [0], [0],
            marker="o",
            color="black",
            markersize=8,
            linestyle="None",
            label=f"Votes ({n_votes})"
        )
        statement_proxy = plt.Line2D(
            [0], [0],
            color="black",
            lw=2,
            label=f"Statements ({n_statements})"
        )
        plt.legend(handles=[vote_proxy, statement_proxy], loc="center left", bbox_to_anchor=(1, 0.5))

        # --- Title ---
        if not user_votes.empty:
            plt.title(
                f"User {user_id} activity | {first_vote} → {last_vote} ({duration_str})"
            )
        else:
            plt.title(f"User {user_id} activity | No votes")

        plt.ylim(-1.5, 1.5)
        plt.tight_layout()
        plt.show()

        # --- Statements text ---
        if not user_statements.empty:
            user_statements_sorted = (
                user_statements.assign(
                    created_dt=pd.to_datetime(
                        pd.to_numeric(user_statements["created_date"], errors="coerce"),
                        unit="ms"
                    )
                )
                .sort_values("created_dt")
            )

            md_text = f"**Statements by {user_id} in submission order:**\n\n"
            for t, s in zip(user_statements_sorted["created_dt"], user_statements_sorted["content"]):
                md_text += f"- {t}: {s}\n"
            display(Markdown(md_text))

    # -----------------------------
    # User selection widgets
    # -----------------------------
    all_voters, all_commenters, all_users = _user_lists(adata, votes_df)
    initial_user = random.choice(all_commenters.tolist())

    user_dropdown = widgets.Dropdown(
        options=sorted(all_users),
        value=initial_user,
        description="User ID:"
    )

    random_voter_btn = widgets.Button(description="Random voter")
    random_commenter_btn = widgets.Button(description="Random commenter")

    def pick_random_voter(_):
        user_dropdown.value = random.choice(all_voters)

    def pick_random_commenter(_):
        user_dropdown.value = random.choice(all_commenters)

    random_voter_btn.on_click(pick_random_voter)
    random_commenter_btn.on_click(pick_random_commenter)

    display(
        widgets.VBox([
            widgets.HBox([user_dropdown, random_voter_btn, random_commenter_btn]),
            widgets.interactive_output(plot_user_activity, {"user_id": user_dropdown})
        ])
    )
