"""Draw time_slicing_empty_slices.png for docs/user/slicing.rst, section "Empty and partially filled slices".

A 60 s run is sliced every 10 s from the run start, and only the events after 25 s are loaded
(``FilterByTimeStart = 25``). The windows 0-10 s and 10-20 s hold no events and are dropped. The window
20-30 s keeps the events from 25 s onward. The remaining slices are labelled with their output suffix.

Usage: python time_slicing_empty_slices.py <output directory>
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

INK, DROPPED, BOX_FACE, BOX_EDGE = "#1f1f1f", "#c8c8c8", "#f4f4f4", "#555555"
SLICE_COLORS = ["#8fb8f6", "#f4736c", "#fde49b", "#8ab56a"]  # colors of the other time-slicing figures
RUN_END, INTERVAL, FILTER_START = 60.0, 10.0, 25.0  # seconds
Y_LINE = 0.0  # height of the time axis


def draw(output_directory: str) -> str:
    """Draw the figure and return the path of the saved PNG."""
    fig, ax = plt.subplots(figsize=(10.8, 3.5))

    # time axis, run start and run end
    ax.annotate("", xy=(RUN_END + 8.0, Y_LINE), xytext=(-4.0, Y_LINE), arrowprops=dict(arrowstyle="-|>", color=INK))
    ax.text(RUN_END + 8.0, Y_LINE - 0.25, "Time (s)", fontsize=14, fontweight="bold", ha="center", va="top")
    for t, label in [(0.0, "Run\nstart"), (RUN_END, "Run\nend")]:
        ax.plot([t, t], [Y_LINE - 0.35, Y_LINE + 1.2], color=INK, lw=1.2)
        ax.text(t, Y_LINE + 1.75, label, fontsize=14, fontweight="bold", ha="center", va="center")

    # load-time cutoff, as in the figure of Example 4
    ax.add_patch(Rectangle((FILTER_START, 0.7), RUN_END - FILTER_START - 0.6, 0.4, fc=BOX_FACE, ec=BOX_EDGE, lw=1.0))
    ax.text((FILTER_START + RUN_END) / 2, 0.9, "FilterByTimeStart = 25 s", fontsize=13, ha="center", va="center")

    # time windows: the first two are empty and dropped, the others become output suffixes _0 to _3
    n_windows = int(RUN_END / INTERVAL)
    for index in range(n_windows):
        start, end = index * INTERVAL, (index + 1) * INTERVAL
        dropped = end <= FILTER_START
        color = DROPPED if dropped else SLICE_COLORS[(index - 2) % len(SLICE_COLORS)]
        ax.plot([start, end], [Y_LINE, Y_LINE], color=color, lw=8, ls=(0, (2, 1)) if dropped else "-",
                solid_capstyle="butt")  # fmt: skip
        label = "dropped" if dropped else f"_{index - 2}"
        ax.text((start + end) / 2, Y_LINE + 0.35, label, fontsize=13, ha="center", va="center",
                color="#777777" if dropped else INK)  # fmt: skip
    for t in range(0, int(RUN_END) + 1, int(INTERVAL)):
        ax.plot([t, t], [Y_LINE - 0.12, Y_LINE + 0.12], color=INK, lw=1.2)
        ax.text(t, Y_LINE - 0.45, f"{t:.0f}", fontsize=12, ha="center", va="top")

    # annotations below the axis
    for x, text in [(10.0, "no events"), (27.5, "events from 25 s kept")]:
        ax.annotate(text, xy=(x, Y_LINE - 0.85), xytext=(x, Y_LINE - 1.6), fontsize=12, ha="center", va="center",
                    arrowprops=dict(arrowstyle="-|>", color=INK))  # fmt: skip

    ax.set_xlim(-5.0, RUN_END + 11.0)
    ax.set_ylim(-1.9, 2.2)
    ax.axis("off")
    path = os.path.join(output_directory, "time_slicing_empty_slices.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_directory", help="directory where the PNG is written")
    draw(parser.parse_args().output_directory)
