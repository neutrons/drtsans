"""Draw wedge_binning_7.png for docs/user/binning.rst, section "Automatic Finding of the Wedge".

The left half is the BIOSANS main detector with the peak wedges w1, w2 and the background wedges b1, b2 found for a
real sample. The right half is a fictitious azimuthal profile I(phi | q_i) with the same wedge orientation, showing
how the wedge limits follow from the fitted centers c1, c2 and widths FWHM1, FWHM2.

Usage: python auto_wedge_biosans.py <detector.png> <output.png>
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SIGMA_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))
BACKGROUND = 0.25  # fictitious flat background A0
PEAKS = [(345.0, 30.0, 0.50), (165.0, 20.0, 0.75)]  # fictitious center (deg), FWHM (deg) and height of peaks 1, 2
PEAK_WIDTH, BACKGROUND_WIDTH = 0.25, 1.5  # autoWedgePeakWidth and autoWedgeBackgroundWidth of the JSON example
INK, MUTED, CURVE, PEAK_SHADE, BACKGROUND_SHADE = "#1f1f1f", "#6b6b6b", "#2a6fdb", "#f6c9a8", "#dcdcf0"


def profile(phi: np.ndarray) -> np.ndarray:
    """Flat background plus the two Gaussian peaks, periodic in phi (deg)."""
    total = np.full_like(phi, BACKGROUND)
    for center, fwhm, height in PEAKS:
        distance = (phi - center + 180.0) % 360.0 - 180.0
        total += height * np.exp(-0.5 * (distance * SIGMA_TO_FWHM / fwhm) ** 2)
    return total


def wedges(ax: plt.Axes) -> None:
    """Shade the peak wedges w1, w2 and the background wedges b1, b2 and label their limits."""
    (c1, fwhm1, _), (c2, fwhm2, _) = PEAKS
    b1 = ((c1 + BACKGROUND_WIDTH * fwhm1) % 360.0, c2 - BACKGROUND_WIDTH * fwhm2)
    b2 = (c2 + BACKGROUND_WIDTH * fwhm2, c1 - BACKGROUND_WIDTH * fwhm1)
    for name, (start, stop) in [("$b_1$", b1), ("$b_2$", b2)]:
        ax.axvspan(start, stop, color=BACKGROUND_SHADE, lw=0, zorder=0)
        ax.text((start + stop) / 2, 0.06, name, fontsize=18, ha="center", color=INK)
    for index, (center, fwhm, height) in enumerate(PEAKS):
        ax.axvspan(center - PEAK_WIDTH * fwhm, center + PEAK_WIDTH * fwhm, color=PEAK_SHADE, lw=0, zorder=0)
        ax.text(center, 0.06, f"$w_{index + 1}$", fontsize=18, ha="center", color=INK)
        half = BACKGROUND + height / 2
        ax.annotate("", xy=(center - fwhm / 2, half), xytext=(center + fwhm / 2, half),
                    arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2))  # fmt: skip
        side = 1.0 if center < 270.0 else -1.0  # keep the label inside the plot
        ax.text(center + side * (fwhm / 2 + 3), half, rf"$\mathrm{{FWHM}}_{index + 1}$", fontsize=12, va="center",
                ha="left" if side > 0 else "right", color=INK)  # fmt: skip
        ax.plot([center, center], [0.0, BACKGROUND + height], color=MUTED, lw=0.8, ls=":")
        ax.text(center, BACKGROUND + height + 0.03, f"$c_{index + 1}$", fontsize=14, ha="center", color=INK)
    ax.text(255.0, 1.25, r"$w_i$ spans $c_i \pm$ autoWedgePeakWidth $\times\,\mathrm{FWHM}_i$", fontsize=12,
            ha="center", va="bottom", color=INK)  # fmt: skip
    for edge, label in [
        (b1[0], r"$c_1 +$ autoWedgeBackgroundWidth $\times\,\mathrm{FWHM}_1$"),
        (b1[1], r"$c_2 -$ autoWedgeBackgroundWidth $\times\,\mathrm{FWHM}_2$"),
    ]:
        ax.axvline(edge, color=MUTED, lw=1, ls="--", zorder=1)
        ax.text(edge + 2.5, 0.2, label, rotation=90, fontsize=11, va="bottom", color=INK)


def draw(detector: str, output: str) -> None:
    """Draw the figure from the detector image ``detector`` and save it to ``output``."""
    fig, (left, right) = plt.subplots(1, 2, figsize=(14.0, 6.6), dpi=170, gridspec_kw=dict(width_ratios=[1, 1.25]))
    left.imshow(plt.imread(detector))
    left.axis("off")
    phi = np.linspace(0.0, 360.0, 1441)
    wedges(right)
    right.plot(phi, profile(phi), color=CURVE, lw=2.2, zorder=2)
    right.set_xlim(0.0, 360.0)
    right.set_ylim(0.0, 1.42)
    right.set_xticks(np.arange(0, 361, 90))
    right.set_yticks([])
    right.set_xlabel("azimuthal angle φ (deg)", fontsize=13)
    right.set_ylabel(r"$I(\phi\,|\,q_i)$", fontsize=14)
    right.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("detector", help="path of the cropped BIOSANS detector image")
    parser.add_argument("output", help="path of the PNG file to write")
    arguments = parser.parse_args()
    draw(arguments.detector, arguments.output)
