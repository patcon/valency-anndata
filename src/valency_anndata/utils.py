import pandas as pd
import numpy as np
from datetime import datetime
import asyncio
from typing import Any, Coroutine


def trim_by_time(
    df: pd.DataFrame,
    rule: int | float | str | datetime,
    col: str = "timestamp",
) -> pd.DataFrame:
    """
    Trim a timestamp-sorted DataFrame using a flexible one-argument rule.

    All rules can be negative to indicate slicing from the end instead of the start.

    Trimming strategies based on the type of `rule`:

    1. **Float (0 < |rule| < 1)**: proportional fraction of rows.
       - Positive → keep first fraction.
       - Example: 0.8 keeps first 80%.
       - Example: -0.2 keeps last 20%.

    2. **Percent integer (1–100)**: fraction of rows as a percentage.
       - Example: 50 keeps first 50%.
       - Example: -10 keeps last 10%.

    3. **Row count integer (>100)**: keep N rows.
       - Example: 500 keeps first 500 rows.
       - Example: -700 keeps last 700 rows.

    4. **Absolute timestamp cutoff (int ≥ 1_200_000_000)**: keeps rows where df[col] <= rule.

    5. **Datetime cutoff (datetime.datetime)**: converted to Unix timestamp and used as cutoff.

    6. **Statistical rule (str)**: "mean-2std" or "mean+3std" computes mean ± k*std and trims.

    7. **Callable (function)**: accepts a timestamp value and returns True/False.
       - Example: lambda t: t < df[col].quantile(0.9)

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing the timestamp column.
    rule : float | int | datetime | str | callable
        Single argument determining the trimming strategy.
    col : str, default "timestamp"
        Name of the timestamp column.

    Returns
    -------
    pd.DataFrame
        A filtered DataFrame with the same columns as the input.

    Raises
    ------
    KeyError
        If the requested column does not exist.
    ValueError
        If rule is zero or an invalid string/float format.
    TypeError
        If `rule` has an unsupported type.

    Notes
    -----
    - Designed to be `.pipe()`-friendly:
          df.pipe(trim_by_time, 0.8)
          df.pipe(trim_by_time, "mean+2std", col="created_at")
    - Requires timestamps to be sorted ascending for proportional trimming.

    Examples
    --------
    >>> df.pipe(trim_by_time, 0.75)          # first 75%
    >>> df.pipe(trim_by_time, -0.25)         # last 25%
    >>> df.pipe(trim_by_time, 50)            # first 50%
    >>> df.pipe(trim_by_time, -10)           # last 10%
    >>> df.pipe(trim_by_time, 500)           # first 500 rows
    >>> df.pipe(trim_by_time, -700)          # last 700 rows
    >>> df.pipe(trim_by_time, 1_500_000_000) # absolute timestamp cutoff
    >>> df.pipe(trim_by_time, "mean-2std")   # statistical cutoff
    >>> df.pipe(trim_by_time, lambda t: t < df.timestamp.quantile(0.9))  # custom
    """

    if col not in df.columns:
        raise KeyError(f"Column '{col}' not found in DataFrame.")

    ts = df[col].to_numpy()  # fast NumPy vector
    n = len(ts)

    # ---------------------------------------------------------------
    # 1. Float proportions (positive = start, negative = end)
    # ---------------------------------------------------------------
    if isinstance(rule, float):
        if rule == 0:
            raise ValueError("Proportional trim cannot be 0")

        # special case: include all rows
        if rule == 1.0:
            return df.copy()
        if rule == -1.0:
            return df.copy()

        if 0 < rule < 1:
            cutoff_idx = int(n * rule)
            cutoff_ts = ts[cutoff_idx]
            return df.loc[df[col] <= cutoff_ts]
        if -1 < rule < 0:
            cutoff_idx = int(n * abs(rule))
            cutoff_ts = ts[-cutoff_idx]
            return df.loc[df[col] >= cutoff_ts]

    # ---------------------------------------------------------------
    # 2. Percent integers (1–100 or -1 to -100)
    # ---------------------------------------------------------------
    if isinstance(rule, (int, np.integer)):
        if rule == 0:
            raise ValueError("Percent trim cannot be 0")

        # special case: 100% or -100%
        if rule == 100 or rule == -100:
            return df.copy()

        if 1 <= rule < 100:
            fraction = rule / 100.0
            cutoff_idx = int(n * fraction)
            cutoff_ts = ts[cutoff_idx]
            return df.loc[df[col] <= cutoff_ts]
        if -100 < rule <= -1:
            fraction = abs(rule) / 100.0
            cutoff_idx = int(n * fraction)
            cutoff_ts = ts[-cutoff_idx]
            return df.loc[df[col] >= cutoff_ts]

        # ---------------------------------------------------------------
        # 3. Row count integers (>100 or <-100)
        # ---------------------------------------------------------------
        if rule > 100:
            return df.iloc[:rule]
        if rule < -100:
            return df.iloc[rule:]

        # ---------------------------------------------------------------
        # 4. Absolute timestamp cutoff (large int)
        # ---------------------------------------------------------------
        if rule >= 1_200_000_000:
            return df.loc[df[col] <= rule]

    # ---------------------------------------------------------------
    # 5. Datetime cutoff
    # ---------------------------------------------------------------
    if isinstance(rule, datetime):
        cutoff = int(rule.timestamp())
        return df.loc[df[col] <= cutoff]

    # ---------------------------------------------------------------
    # 6. Statistical rule: e.g. "mean-2std", "mean+3std"
    # ---------------------------------------------------------------
    if isinstance(rule, str):
        if "mean" in rule and "std" in rule:
            mu = ts.mean()
            sd = ts.std()

            if "-" in rule:
                k = float(rule.split("-")[1].replace("std", ""))
                cutoff = mu - k * sd
            elif "+" in rule:
                k = float(rule.split("+")[1].replace("std", ""))
                cutoff = mu + k * sd
            else:
                raise ValueError(f"Could not parse std rule: {rule}")

            return df.loc[df[col] <= cutoff]

        raise ValueError(f"Unknown trim rule string: {rule}")

    # ---------------------------------------------------------------
    # 7. Callable rule
    # ---------------------------------------------------------------
    if callable(rule):
        mask = df[col].map(rule)
        return df.loc[mask]

    # ---------------------------------------------------------------
    # Unsupported type
    # ---------------------------------------------------------------
    raise TypeError(f"Unsupported rule type: {type(rule)}")

def run_async(coro: Coroutine[Any, Any, Any]) -> Any:
    """
    Run an async coroutine safely in scripts or Jupyter notebooks.

    Usage:
        result = run_async(my_async_function(...))
    """
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = None

    if loop and loop.is_running():
        # Jupyter notebook or other running loop
        import nest_asyncio
        nest_asyncio.apply()
        return loop.run_until_complete(coro)
    else:
        # normal script
        return asyncio.run(coro)