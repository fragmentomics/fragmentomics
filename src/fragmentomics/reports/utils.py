"""
Utility functions for report generation.
"""

from __future__ import annotations

import base64
from io import BytesIO
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from matplotlib.figure import Figure


def fig_to_base64(fig: Figure, format: str = "png", dpi: int = 150) -> str:
    """
    Convert matplotlib figure to base64-encoded image.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure
    format : str, default "png"
        Image format
    dpi : int, default 150
        Resolution

    Returns
    -------
    str
        Base64-encoded image data URL
    """
    buf = BytesIO()
    fig.savefig(buf, format=format, dpi=dpi, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    return f"data:image/{format};base64,{b64}"
