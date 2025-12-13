from datetime import datetime
from anndata import AnnData
from .utils import trim_by_time
import pandas as pd

def rebuild_vote_matrix(
    data: AnnData,
    trim_rule: int | float | str | datetime = 1.0,
    time_col: str = "timestamp",
    inplace: bool = True,
) -> AnnData | None:
    """
    Rebuild a vote matrix from votes CSV, trimming by time, deduplicating votes,
    and returning a new AnnData with vote matrix. Preserves metadata from original.
    """

    # Load votes CSV
    votes_df = data.uns.get("votes_csv")
    if votes_df is None:
        raise KeyError("`uns['votes_csv']` not found in AnnData")
    votes_df = votes_df.copy()

    # Trim by time
    votes_df = votes_df.pipe(trim_by_time, rule=trim_rule, col=time_col)

    # Sort & deduplicate
    votes_df = votes_df.sort_values(time_col)
    votes_df = votes_df.drop_duplicates(
        subset=["voter-id", "comment-id"], keep="last"
    )

    # Pivot into voter × comment
    vote_matrix_df = votes_df.pivot(
        index="voter-id", columns="comment-id", values="vote"
    )

    # Build a new AnnData
    new_adata = AnnData(
        X=vote_matrix_df.to_numpy(dtype=float),
        obs=pd.DataFrame(index=vote_matrix_df.index.astype(str)),
        var=pd.DataFrame(index=vote_matrix_df.columns.astype(str))
    )

    # Copy over other metadata
    new_adata.uns.update(data.uns)
    new_adata.obsm.update(data.obsm)
    new_adata.layers.update(data.layers)

    if inplace:
        # Replace all internal state of the original AnnData
        data._init_as_actual(new_adata)
        return None
    else:
        return new_adata
