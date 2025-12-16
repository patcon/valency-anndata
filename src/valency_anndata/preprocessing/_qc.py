import numpy as np
import pandas as pd
from typing import Optional, Tuple, Any
import anndata as ad
from scanpy import _compat

# Use Any for now to handle all possible AnnData matrix types
Matrix = Any

def describe_obs(X: Matrix, obs_names: Optional[pd.Index] = None) -> pd.DataFrame:
    """Participant-level QC metrics for Polis-style vote matrix."""
    if isinstance(X, _compat.CSBase) :
        X = X.toarray()
    n_obs, n_var = X.shape
    mask_seen = ~np.isnan(X)

    n_votes = mask_seen.sum(axis=1)
    n_agree = ((X == 1) & mask_seen).sum(axis=1)
    n_disagree = ((X == -1) & mask_seen).sum(axis=1)
    n_pass = ((X == 0) & mask_seen).sum(axis=1)
    n_engaged = n_agree + n_disagree

    pct_agree = np.divide(n_agree, n_votes, out=np.zeros_like(n_agree, dtype=float), where=n_votes>0)
    pct_disagree = np.divide(n_disagree, n_votes, out=np.zeros_like(n_disagree, dtype=float), where=n_votes>0)
    pct_pass = np.divide(n_pass, n_votes, out=np.zeros_like(n_pass, dtype=float), where=n_votes>0)
    pct_engaged = np.divide(n_engaged, n_votes, out=np.zeros_like(n_engaged, dtype=float), where=n_votes>0)
    pct_agree_engaged = np.divide(n_agree, n_engaged, out=np.zeros_like(n_agree, dtype=float), where=n_engaged>0)
    pct_seen = n_votes / n_var
    mean_vote = np.divide(np.nansum(X, axis=1), n_votes, out=np.zeros_like(n_votes, dtype=float), where=n_votes>0)

    return pd.DataFrame(
        {
            "n_votes": n_votes,
            "n_agree": n_agree,
            "n_disagree": n_disagree,
            "n_pass": n_pass,
            "n_engaged": n_engaged,
            "pct_agree": pct_agree,
            "pct_disagree": pct_disagree,
            "pct_pass": pct_pass,
            "pct_engaged": pct_engaged,
            "pct_agree_engaged": pct_agree_engaged,
            "pct_seen": pct_seen,
            "mean_vote": mean_vote,
        },
        index=obs_names,
    )


def describe_var(X: Matrix, var_names: Optional[pd.Index] = None) -> pd.DataFrame:
    """Statement-level QC metrics for Polis-style vote matrix."""
    if isinstance(X, _compat.CSBase):
        X = X.toarray()
    X = X.T  # transpose for statements
    n_var, n_obs = X.shape
    mask_seen = ~np.isnan(X)

    n_votes = mask_seen.sum(axis=1)
    n_agree = ((X == 1) & mask_seen).sum(axis=1)
    n_disagree = ((X == -1) & mask_seen).sum(axis=1)
    n_pass = ((X == 0) & mask_seen).sum(axis=1)
    n_engaged = n_agree + n_disagree

    pct_agree = np.divide(n_agree, n_votes, out=np.zeros_like(n_agree, dtype=float), where=n_votes>0)
    pct_disagree = np.divide(n_disagree, n_votes, out=np.zeros_like(n_disagree, dtype=float), where=n_votes>0)
    pct_pass = np.divide(n_pass, n_votes, out=np.zeros_like(n_pass, dtype=float), where=n_votes>0)
    pct_engaged = np.divide(n_engaged, n_votes, out=np.zeros_like(n_engaged, dtype=float), where=n_votes>0)
    pct_agree_engaged = np.divide(n_agree, n_engaged, out=np.zeros_like(n_agree, dtype=float), where=n_engaged>0)
    pct_seen = n_votes / n_obs
    mean_vote = np.divide(np.nansum(X, axis=1), n_votes, out=np.zeros_like(n_votes, dtype=float), where=n_votes>0)

    return pd.DataFrame(
        {
            "n_votes": n_votes,
            "n_agree": n_agree,
            "n_disagree": n_disagree,
            "n_pass": n_pass,
            "n_engaged": n_engaged,
            "pct_agree": pct_agree,
            "pct_disagree": pct_disagree,
            "pct_pass": pct_pass,
            "pct_engaged": pct_engaged,
            "pct_agree_engaged": pct_agree_engaged,
            "pct_seen": pct_seen,
            "mean_vote": mean_vote,
        },
        index=var_names,
    )


def calculate_qc_metrics(
    adata: ad.AnnData,
    *,
    inplace: bool = False,
) -> Optional[Tuple[pd.DataFrame, pd.DataFrame]]:
    """Compute participant- and statement-level metrics using describe_obs and describe_var."""
    X = adata.X
    if X is None:
        raise ValueError("adata.X is None")
    obs_metrics = describe_obs(X, obs_names=adata.obs_names)
    var_metrics = describe_var(X, var_names=adata.var_names)

    if inplace:
        adata.obs[obs_metrics.columns] = obs_metrics
        adata.var[var_metrics.columns] = var_metrics
        return None

    return obs_metrics, var_metrics
