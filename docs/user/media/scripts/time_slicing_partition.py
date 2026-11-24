"""Draw the time-slicing partition figures for docs/user/slicing.rst.

- ``time_slicing_offset.png`` (Example 2): a 10 s run sliced every 2 s from an offset of 1 s. The final slice
  absorbs the 1 s leftover of the run and lasts 3 s, from 7 s to 10 s.
- ``time_slicing_interval_list.png`` (Example 3): a 13 s run sliced with the list ``[2, 4]``. The last value, 4 s, is
  repeated, and the final slice absorbs the 3 s leftover and lasts 7 s, from 6 s to 13 s.

Usage: python time_slicing_partition.py <output directory>
"""

import argparse
import os
from typing import List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK = "#1f1f1f"
SLICE_COLORS = ["#8fb8f6", "#f4736c", "#fde49b", "#8ab56a", "#cfa0c8"]  # colors of the other time-slicing figures
Y_LINE = 0.0  # height of the time axis


def draw(path: str, run_end: float, offset: float, windows: List[Tuple[float, float]], notes: List[Tuple]) -> str:
    """Draw one time axis with its slices, and return the path of the saved PNG.

    Parameters
    ----------
    path
        Output PNG file
    run_end
        Run end, in seconds from the run start
    offset
        Start of the first slice, in seconds
    windows
        ``(start, end)`` of each slice, in seconds
    notes
        ``(x, text)`` annotations drawn below the axis, pointing at time ``x``
    """
    fig, ax = plt.subplots(figsize=(10.8, 3.5))

    # time axis, run start and run end
    ax.annotate("", xy=(run_end * 1.12, Y_LINE), xytext=(-0.06 * run_end, Y_LINE),
                arrowprops=dict(arrowstyle="-|>", color=INK))  # fmt: skip
    ax.text(run_end * 1.12, Y_LINE - 0.25, "Time (s)", fontsize=14, fontweight="bold", ha="center", va="top")
    for t, label in [(0.0, "Run\nstart"), (run_end, "Run\nend")]:
        ax.plot([t, t], [Y_LINE - 0.35, Y_LINE + 1.2], color=INK, lw=1.2)
        ax.text(t, Y_LINE + 1.75, label, fontsize=14, fontweight="bold", ha="center", va="center")

    # offset, as a double arrow between the run start and the first slice
    if offset > 0:
        ax.annotate("", xy=(0.0, Y_LINE), xytext=(offset, Y_LINE),
                    arrowprops=dict(arrowstyle="<|-|>", color=INK, lw=1.5))  # fmt: skip
        ax.text(offset / 2, Y_LINE + 0.35, "offset", fontsize=13, ha="center", va="center")

    # slices, labelled with their index and length
    boundaries = sorted({0.0, run_end, *[t for window in windows for t in window]})
    for index, (start, end) in enumerate(windows):
        ax.plot([start, end], [Y_LINE, Y_LINE], color=SLICE_COLORS[index % len(SLICE_COLORS)], lw=8,
                solid_capstyle="butt", zorder=3)  # fmt: skip
        ax.text((start + end) / 2, Y_LINE + 0.35, f"slice {index}", fontsize=13, ha="center", va="center")
        ax.text((start + end) / 2, Y_LINE + 0.75, f"{end - start:g} s", fontsize=11, ha="center", va="center",
                color="#555555")  # fmt: skip
    for t in boundaries:
        ax.plot([t, t], [Y_LINE - 0.12, Y_LINE + 0.12], color=INK, lw=1.2)
        ax.text(t, Y_LINE - 0.45, f"{t:.1f}", fontsize=12, ha="center", va="top")

    # annotations below the axis
    for x, text in notes:
        ax.annotate(text, xy=(x, Y_LINE - 0.85), xytext=(x, Y_LINE - 1.6), fontsize=12, ha="center", va="center",
                    arrowprops=dict(arrowstyle="-|>", color=INK))  # fmt: skip

    ax.set_xlim(-0.1 * run_end, run_end * 1.2)
    ax.set_ylim(-1.9, 2.2)
    ax.axis("off")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def draw_all(output_directory: str) -> List[str]:
    """Draw both figures in `output_directory` and return the paths of the saved PNGs."""
    offset_figure = draw(
        os.path.join(output_directory, "time_slicing_offset.png"),
        run_end=10.0,
        offset=1.0,
        windows=[(1.0, 3.0), (3.0, 5.0), (5.0, 7.0), (7.0, 10.0)],
        notes=[(8.5, "final slice absorbs the 1 s leftover")],
    )
    list_figure = draw(
        os.path.join(output_directory, "time_slicing_interval_list.png"),
        run_end=13.0,
        offset=0.0,
        windows=[(0.0, 2.0), (2.0, 6.0), (6.0, 13.0)],
        notes=[(9.5, "last value (4 s) repeated, plus the 3 s leftover")],
    )
    return [offset_figure, list_figure]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_directory", help="directory where the PNG files are written")
    draw_all(parser.parse_args().output_directory)
