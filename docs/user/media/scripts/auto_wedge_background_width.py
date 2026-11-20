"""Draw the auto_wedge_background_*.png figures for docs/user/binning.rst, section "From peaks to wedges".

Each panel shows a fictitious azimuthal profile I(phi) with two Gaussian peaks, unrolled from 0 to 360 deg, the two
peak wedges and the two background wedges, and the angular width of each wedge. The background width is applied as
in drtsans.auto_wedge.getWedgeSelection: kept, reduced, or rejected.

Usage: python auto_wedge_background_width.py <output directory>
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SIGMA_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))
PEAK_WIDTH = 0.25  # default of autoWedgePeakWidth
INK, CURVE, PEAK_COLOR, BACKGROUND_COLOR, FAIL_COLOR = "#1f1f1f", "#3d3d3d", "#2a6fdb", "#e8955a", "#c62828"
# per figure, two panels: title, peak centers (deg), FWHM (deg), requested autoWedgeBackgroundWidth
FIGURES = {
    "kept_reduced": [
        ("(a) FWHM 30°: the requested 1.5 is kept", (90.0, 270.0), 30.0, 1.5),
        ("(b) FWHM 95°: the requested 1.5 is reduced to 0.70", (90.0, 270.0), 95.0, 1.5),
    ],
    "no_room_overlap": [
        ("(a) peaks 140° apart, FWHM 150°: no room for background wedges, fails", (90.0, 230.0), 150.0, 1.5),
        ("(b) FWHM 80°, requested 0.1 ≤ autoWedgePeakWidth: kept, wedges overlap", (90.0, 270.0), 80.0, 0.1),
    ],
    "not_antipodal": [
        ("(a) peaks 160° apart, FWHM 20°: background wedges of 100° and 140°", (90.0, 250.0), 20.0, 1.5),
        ("(b) peaks 50° apart, FWHM 20°: a background wedge of 270°, fails", (90.0, 140.0), 20.0, 1.5),
    ],
}


def applied_background_width(centers: tuple[float, float], fwhm: float, requested: float) -> tuple[float, bool]:
    """Background width applied by getWedgeSelection, and whether the peaks leave room for background wedges."""
    if requested <= PEAK_WIDTH:  # kept as requested, the wedges overlap
        return requested, True
    separation = min(abs(centers[1] - centers[0]), 360.0 - abs(centers[1] - centers[0]))
    if (requested + PEAK_WIDTH) * 2 * fwhm <= separation:
        return requested, True
    reduced = separation / (2 * fwhm) - PEAK_WIDTH
    return reduced, reduced > PEAK_WIDTH


def spans(start: float, width: float) -> list[tuple[float, float]]:
    """Split the angular range [start, start + width] (deg) into pieces within [0, 360]."""
    start %= 360.0
    if start + width <= 360.0:
        return [(start, start + width)]
    return [(start, 360.0), (0.0, start + width - 360.0)]


def shade(ax: plt.Axes, start: float, width: float, color: str, label: str, hatch: str | None = None) -> None:
    """Shade a wedge of ``width`` deg starting at ``start`` deg and write its width above it."""
    for low, high in spans(start, width):
        ax.axvspan(low, high, color=color, alpha=0.35 if hatch is None else 0.15, hatch=hatch, lw=0)
    widest = max(spans(start, width), key=lambda piece: piece[1] - piece[0])
    ax.text((widest[0] + widest[1]) / 2, 0.385, label, color=INK, ha="center", va="bottom", fontsize=10)


def panel(ax: plt.Axes, title: str, centers: tuple[float, float], fwhm: float, requested: float) -> None:
    """Draw one panel: the profile, the peak wedges and the background wedges."""
    width, room = applied_background_width(centers, fwhm, requested)
    phi = np.linspace(0.0, 360.0, 1441)
    profile = np.full_like(phi, 0.2)
    for center in centers:
        distance = (phi - center + 180.0) % 360.0 - 180.0
        profile += 0.15 * np.exp(-0.5 * (distance * SIGMA_TO_FWHM / fwhm) ** 2)
    ax.plot(phi, profile, color=CURVE, lw=1.6)
    for center in centers:
        shade(ax, center - PEAK_WIDTH * fwhm, 2 * PEAK_WIDTH * fwhm, PEAK_COLOR, f"{2 * PEAK_WIDTH * fwhm:.1f}°")
    distance = (centers[1] - centers[0]) % 360.0
    for start, extent in [(centers[0], distance), (centers[1], 360.0 - distance)]:
        background = extent - 2 * width * fwhm  # angular width of this background wedge
        failed = not room or not 0.0 < background < 180.0
        shade(ax, start + width * fwhm, background, FAIL_COLOR if failed else BACKGROUND_COLOR,
              f"{background:.1f}°", hatch="//" if failed else None)  # fmt: skip
    if width < PEAK_WIDTH:  # mark where the background wedges enter the peak wedges, on both sides of each peak
        for center in centers:
            for sign in (1.0, -1.0):
                start = center + sign * (width if sign > 0 else PEAK_WIDTH) * fwhm
                for low, high in spans(start, (PEAK_WIDTH - width) * fwhm):
                    ax.axvspan(low, high, fc="none", ec=INK, hatch="xx", lw=0)
    ax.set_title(title, loc="left", fontsize=11, color=FAIL_COLOR if not room or "fails" in title else INK)
    ax.set_xlim(0.0, 360.0)
    ax.set_ylim(0.15, 0.43)
    ax.set_xticks(np.arange(0, 361, 45))
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)


def draw(name: str, output: str) -> None:
    """Draw the figure ``name`` and save it to ``output``."""
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 5.6), dpi=150, sharex=True)
    for ax, (title, centers, fwhm, requested) in zip(axes, FIGURES[name]):
        panel(ax, title, centers, fwhm, requested)
    axes[-1].set_xlabel("azimuthal angle φ (deg)", fontsize=11)
    handles = [plt.Rectangle((0, 0), 1, 1, color=PEAK_COLOR, alpha=0.35),
               plt.Rectangle((0, 0), 1, 1, color=BACKGROUND_COLOR, alpha=0.35),
               plt.Rectangle((0, 0), 1, 1, fc=FAIL_COLOR, alpha=0.15, hatch="//")]  # fmt: skip
    labels = ["peak wedge", "background wedge", "invalid background wedge"]
    if name == "no_room_overlap":
        handles.append(plt.Rectangle((0, 0), 1, 1, fc="none", ec=INK, hatch="xx", lw=0))
        labels.append("overlap of background and peak wedges")
    fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False, fontsize=10)
    fig.tight_layout(rect=(0.0, 0.06, 1.0, 1.0))
    fig.savefig(output)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_dir", help="directory where to write the PNG files")
    directory = parser.parse_args().output_dir
    for figure_name in FIGURES:
        draw(figure_name, os.path.join(directory, f"auto_wedge_background_{figure_name}.png"))
