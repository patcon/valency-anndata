#!/usr/bin/env python3
"""
Convert ESS human values block (PVQ-21 items) to a polis-like AnnData object.

Downloads ESS round data via the ESS API, extracts the 21 Portrait Values
Questionnaire items, and maps 1-6 Likert responses to {+1, 0, -1} votes.

Usage:
    uv run scripts/ess_to_adata.py --user-id YOUR_USER_ID
    uv run scripts/ess_to_adata.py --user-id YOUR_USER_ID --rubric strict

Find your ESS user ID at https://ess.sikt.no/en/api after logging in.
"""

import argparse
import os
import sys
import urllib.request
import urllib.error
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from pandas.io.stata import StataReader
from platformdirs import user_cache_dir


ESS_API_BASE = "https://api.ess.sikt.no"
ESS_DOI_PREFIX = "10.21338"

_ROUND_DOIS = {
    11: "ess11e01_0",
}

# Portrait Values Questionnaire items (ESS human values module)
# All use a 1-6 scale: 1=Very much like me, 6=Not like me at all
HUMAN_VALUES_VARS = [
    "ipcrtiva", "impricha", "ipeqopta", "ipshabta", "impsafea", "impdiffa",
    "ipfrulea", "ipudrsta", "ipmodsta", "ipgdtima", "impfreea", "iphlppla",
    "ipsucesa", "ipstrgva", "ipadvnta", "ipbhprpa", "iprspota", "iplylfra",
    "impenva",  "imptrada", "impfuna",
]

# Vote mapping rubrics for the 1-6 scale
# 1=Very much like me (+1), 6=Not like me at all (-1)
# Missing codes 66/77/88/99 are recoded to NaN by the API before download.
RUBRICS = {
    "no_pass":            {1: 1,    2: 1,     3: 1,     4: -1,    5: -1,    6: -1},
    "narrow_pass":        {1: 1,    2: 1,     3: 0,     4: 0,     5: -1,    6: -1},
    "wide_pass":          {1: 1,    2: 0,     3: 0,     4: 0,     5: 0,     6: -1},
    "fractional_linear":  {1: 1.0,  2: 0.67,  3: 0.33,  4: -0.33, 5: -0.67, 6: -1.0},
    # tanh(2x) / tanh(2) applied to fractional_linear values
    "fractional_sigmoid": {1: 1.0,  2: 0.91,  3: 0.60,  4: -0.60, 5: -0.91, 6: -1.0},
}


def _get_cache_path(round_num: int) -> Path:
    cache_dir = Path(user_cache_dir("valency-anndata")) / "ess"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"round{round_num}_data.dta"


def _download_ess_file(round_num: int, user_id: str) -> bytes:
    doi_suffix = _ROUND_DOIS.get(round_num)
    if doi_suffix is None:
        supported = sorted(_ROUND_DOIS.keys())
        raise ValueError(f"Round {round_num} not supported. Supported rounds: {supported}")

    url = (
        f"{ESS_API_BASE}/v1/data/dataFile/{ESS_DOI_PREFIX}/{doi_suffix}"
        f"?userId={user_id}&fileFormat=dta&recodeMissingValues"
    )

    print(f"Downloading ESS round {round_num} data from API...")
    try:
        with urllib.request.urlopen(url) as response:
            return response.read()
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise ValueError(
                f"ESS round {round_num} data file not found (HTTP 404). "
                f"Check that DOI suffix '{doi_suffix}' is correct."
            ) from e
        raise


def _get_or_download_dta(round_num: int, user_id: str, skip_cache: bool) -> Path:
    cache_path = _get_cache_path(round_num)

    if not skip_cache and cache_path.exists():
        print(f"Using cached file: {cache_path}")
        return cache_path

    data = _download_ess_file(round_num, user_id)
    cache_path.write_bytes(data)
    print(f"Cached to {cache_path}")
    return cache_path


def _load_variable_labels(path: Path) -> dict:
    with StataReader(str(path)) as reader:
        return reader.variable_labels()  # {varname: human-readable label}


def _apply_rubric(series: pd.Series, mapping: dict) -> pd.Series:
    def _map(x):
        if pd.isna(x):
            return np.nan
        try:
            return mapping.get(int(x), np.nan)
        except (ValueError, TypeError):
            return np.nan
    return series.map(_map)


def build_adata(dta_path: Path, rubric: str, round_num: int) -> AnnData:
    mapping = RUBRICS[rubric]

    df = pd.read_stata(str(dta_path), convert_categoricals=False)
    var_labels = _load_variable_labels(dta_path)

    missing = [v for v in HUMAN_VALUES_VARS if v not in df.columns]
    if missing:
        raise ValueError(f"Variables not found in data: {missing}")

    votes = df[HUMAN_VALUES_VARS].apply(lambda col: _apply_rubric(col, mapping))

    obs_index = (
        df["cntry"].astype(str)
        + "_"
        + df["idno"].fillna(0).astype(int).astype(str)
    )
    votes.index = obs_index

    obs = pd.DataFrame(
        {"country": df["cntry"].values, "idno": df["idno"].values},
        index=obs_index,
    )
    var = pd.DataFrame(
        {
            "content":       [var_labels.get(v, v) for v in HUMAN_VALUES_VARS],
            "spectrum_low":  ["(1) Very much like me"] * len(HUMAN_VALUES_VARS),
            "spectrum_high": ["(6) Not like me at all"] * len(HUMAN_VALUES_VARS),
        },
        index=HUMAN_VALUES_VARS,
    )

    adata = AnnData(X=votes.values.astype(float), obs=obs, var=var)
    adata.uns["source"] = {
        "kind": "ess",
        "round": round_num,
        "doi": f"{ESS_DOI_PREFIX}/{_ROUND_DOIS[round_num]}",
        "rubric": rubric,
    }

    # Long-format votes table expected by downstream apps (e.g. cartography viewer)
    votes_long = (
        votes.reset_index()
        .rename(columns={"index": "voter-id"})
        .melt(id_vars="voter-id", var_name="comment-id", value_name="vote")
        .dropna(subset=["vote"])
    )
    adata.uns["votes"] = votes_long

    return adata


def main():
    parser = argparse.ArgumentParser(
        description="Convert ESS human values (PVQ-21) to a polis-like AnnData .h5ad file."
    )
    parser.add_argument(
        "--user-id",
        default=os.environ.get("ESS_USER_ID"),
        help="ESS user ID for API tracking. Falls back to $ESS_USER_ID env var.",
    )
    parser.add_argument(
        "--round",
        type=int,
        default=11,
        dest="round_num",
        help="ESS round number (default: 11).",
    )
    parser.add_argument(
        "--rubric",
        choices=list(RUBRICS),
        default="narrow_pass",
        help="Vote mapping rubric (default: narrow_pass). "
             "no_pass: 1-3=+1, 4-6=-1. narrow_pass: 1/2=+1, 3/4=0, 5/6=-1. "
             "wide_pass: 1=+1, 2-5=0, 6=-1. fractional_linear: 1.0/0.67/0.33/-0.33/-0.67/-1.0. "
             "fractional_sigmoid: 1.0/0.91/0.60/-0.60/-0.91/-1.0.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output .h5ad path (default: exports/ess_round{round}_human_values.h5ad).",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Re-download even if a cached file exists.",
    )
    args = parser.parse_args()

    if not args.user_id:
        print(
            "Error: ESS user ID required. Pass --user-id or set $ESS_USER_ID.\n"
            "Find your user ID at https://ess.sikt.no/en/api after logging in.",
            file=sys.stderr,
        )
        sys.exit(1)

    output_path = args.output or Path(f"exports/ess_round{args.round_num}_human_values.h5ad")

    dta_path = _get_or_download_dta(args.round_num, args.user_id, args.no_cache)
    adata = build_adata(dta_path, args.rubric, args.round_num)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    print(f"Wrote {adata.n_obs} participants × {adata.n_vars} variables to {output_path}")


if __name__ == "__main__":
    main()
