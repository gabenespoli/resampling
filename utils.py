from typing import Union

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pandas import Series
from scipy.stats import norm


def get_ci(
    x: Union[list, Series],
    y: float = 0.95,
    method: str = "normal_range",
) -> tuple:
    """Get confidence interval for a vector of values.

    Args:
        x: Vector of values to get CI.
        y: Numeric input to function being called.
        method: Which function to call. See those functions for detail about their
            input. Possible values are "normal_range" and "percentile".

    Returns:
        tuple: A tuple with the lower and upper bounds of the confidence interval.

    """
    if method in [1, "normal_range"]:
        return get_ci_normal_range(x, confidence=y)
    elif method in [2, "percentile"]:
        return get_ci_percentile(x, percentile=y)
    raise ValueError(f"Unknown method for getting CI: {method}")


def get_ci_normal_range(x: Union[list, Series], confidence: float = 0.95) -> tuple:
    """Get a confidence interval (CI) for a vector.

    Uses the following formula:
      mean +- ((z * sd) / sqrt(n)), where z is the critical value.

    Args:
        x: Vector from which to calculate a CI.
        confidence: How much confidence the interval should have as a number between 0
            and 1.

    Returns:
        tuple: A tuple with the lower and upper bounds of the confidence interval.

    """
    z = norm.ppf((1 + confidence) / 2)
    n = len(x)
    mean = sum(x) / n
    sd = (sum((xi - mean) ** 2 for xi in x) / n) ** 0.5
    delta = z * sd / n**0.5
    return mean - delta, mean + delta


def get_ci_percentile(x, percentile: float = 0.95) -> tuple:
    n = len(x)
    lower = int(n * (1 - percentile) / 2)
    upper = n - lower
    x_sorted = sorted(x)
    return x_sorted[lower], x_sorted[upper]


def plot(
    df_resampled: pd.DataFrame,
    title: str = "",
    ci_line_style: dict = dict(mode="lines", line=dict(dash="dash", color="black")),
):
    ci_lines = df_resampled.groupby("n")[["ci_upper", "ci_lower"]].mean().reset_index()

    fig = px.scatter(
        df_resampled,
        x="n",
        y="mean",
        # color="iter",
        title=title,
    )

    fig.add_trace(
        go.Scatter(
            x=ci_lines["n"],
            y=ci_lines["ci_upper"],
            **ci_line_style,
        ),
    )

    fig.add_trace(
        go.Scatter(
            x=ci_lines["n"],
            y=ci_lines["ci_lower"],
            **ci_line_style,
        ),
    )

    fig.update_layout(
        xaxis_title="n Sampled",
        yaxis_title="Centered Rating",
        showlegend=False,
    )

    return fig
