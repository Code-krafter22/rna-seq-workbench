"""Per-method gene expression-profile figure for the gene-of-interest view.

Builds one boxplot+jitter panel per DE method (each in its own normalization
space) for a single gene, with a user-selectable diverging color scheme and
large fonts, and exports the figure to PNG / SVG / JPEG bytes for download.
"""
from __future__ import annotations

import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# two-color schemes (first color = first condition, second = second condition)
PALETTES: dict[str, list[str]] = {
    "Red–Blue":    ["#d62728", "#1f77b4"],
    "Yellow–Blue": ["#e1b12c", "#0097e6"],
    "Blue–Orange": ["#0984e3", "#e17055"],
    "Pink–Green":  ["#e84393", "#00b894"],
}


def colors_for(conditions: list[str], palette_name: str) -> dict[str, str]:
    pal = PALETTES.get(palette_name, PALETTES["Red–Blue"])
    return {c: pal[i % len(pal)] for i, c in enumerate(conditions)}


def build_gene_figure(gene: str, per_method_long: dict, conditions: list[str],
                      palette_name: str, font_size: int = 16):
    """per_method_long: {method: DataFrame[sample, expression, condition]}."""
    methods = list(per_method_long.keys())
    colors = colors_for(conditions, palette_name)
    n = max(1, len(methods))
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5.4), squeeze=False)
    axes = axes[0]

    for ax, m in zip(axes, methods):
        df = per_method_long[m]
        data = [df.loc[df.condition == c, "expression"].to_numpy() for c in conditions]
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False,
                        medianprops=dict(color="black", linewidth=2))
        for patch, c in zip(bp["boxes"], conditions):
            patch.set_facecolor(colors[c])
            patch.set_alpha(0.55)
        rng = np.random.default_rng(0)
        for i, c in enumerate(conditions, start=1):
            y = df.loc[df.condition == c, "expression"].to_numpy()
            x = rng.normal(i, 0.06, size=len(y))
            ax.scatter(x, y, color=colors[c], edgecolor="black", linewidth=0.3,
                       s=30, alpha=0.85, zorder=3)
        ax.set_xticks(range(1, len(conditions) + 1))
        ax.set_xticklabels(conditions, fontsize=font_size)
        ax.set_title(m, fontsize=font_size + 2, fontweight="bold")
        ax.tick_params(axis="y", labelsize=font_size)
        ax.set_ylabel("normalized expression", fontsize=font_size)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    fig.suptitle(f"{gene}: expression by group (per method)",
                 fontsize=font_size + 5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


def figure_bytes(fig, fmt: str) -> bytes:
    """Render the figure to bytes in 'png', 'svg', or 'jpeg'."""
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=200, bbox_inches="tight",
                facecolor="white")
    buf.seek(0)
    return buf.getvalue()
