"""Draw auto_wedge_criteria.png for docs/user/binning.rst, section "How Automatic Wedge Finding Works".

The figure shows the contrast H/A0 and the significance H/sigma_H of fictitious Q-rings against Q, with the
default thresholds of autoWedgeContrastMin and autoWedgeSignalToNoiseMin and the status of each ring. Most rejected
rings at high Q are significant but have a low contrast, the case where significance alone would accept
systematic artifacts.

Usage: python auto_wedge_criteria.py <output.png>
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CONTRAST_MIN, SIGNIFICANCE_MIN = 0.2, 10.0  # defaults of autoWedgeContrastMin and autoWedgeSignalToNoiseMin
# fictitious rings: Q (1/A), contrast H/A0 and significance H/sigma_H of the weaker peak, and status
RINGS = [
    (0.02, 0.45, 40.0, "seed"),
    (0.04, 0.70, 65.0, "seed"),
    (0.06, 0.85, 90.0, "seed"),
    (0.08, 0.55, 70.0, "seed"),
    (0.10, 0.30, 45.0, "seed"),
    (0.12, 0.12, 22.0, "rescued"),
    (0.14, 0.050, 25.0, "rejected"),
    (0.16, 0.040, 30.0, "rejected"),
    (0.18, 0.060, 18.0, "rejected"),
    (0.20, 0.020, 6.0, "rejected"),
    (0.22, 0.045, 14.0, "rejected"),
    (0.24, 0.010, 1.5, "rejected"),
]
STYLES = {  # marker shape and color, so that the status does not rely on color alone
    "seed": dict(marker="o", color="#2a6fdb", label="seed"),
    "rescued": dict(marker="^", color="#d9822b", label="rescued: matches the seed orientation"),
    "rejected": dict(marker="x", color="#6b6b6b", label="rejected"),
}
INK = "#1f1f1f"


def draw(output: str) -> None:
    """Draw the figure and save it to ``output``."""
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8.0, 6.4), dpi=150, sharex=True)
    for axis, column, threshold, name in [
        (top, 1, CONTRAST_MIN, "contrast $H/A_0$"),
        (bottom, 2, SIGNIFICANCE_MIN, r"significance $H/\sigma_H$"),
    ]:
        for status, style in STYLES.items():
            rings = [ring for ring in RINGS if ring[3] == status]
            axis.plot([ring[0] for ring in rings], [ring[column] for ring in rings], ls="none", ms=8, mew=2,
                      **style)  # fmt: skip
        axis.axhline(threshold, color=INK, lw=1, ls="--")
        axis.text(0.003, threshold * 1.12, f"seed threshold {threshold:g}", ha="left", va="bottom", color=INK,
                  fontsize=10)  # fmt: skip
        axis.set_yscale("log")
        axis.set_ylabel(name, fontsize=12)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", which="major", color="#e5e5e5", lw=0.8)
        axis.set_axisbelow(True)
    bottom.text(0.18, 45.0, "significant, but low contrast", ha="center", color=INK, fontsize=10)
    bottom.set_xlabel("Q (1/Å)", fontsize=12)
    bottom.set_xlim(0.0, 0.26)
    top.set_title("Weaker peak of each Q-ring", fontsize=12)
    top.legend(loc="lower left", frameon=False, fontsize=10)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output", help="path of the PNG file to write")
    draw(parser.parse_args().output)
