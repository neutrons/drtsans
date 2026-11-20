"""Draw auto_wedge_fit_model.png for docs/user/binning.rst, section "How Automatic Wedge Finding Works".

The figure shows a fictitious azimuthal profile I(phi) of one Q-ring, the auto-wedge fit model (a flat background
plus two Gaussian peaks) and the fit parameters: background A0, height H with its uncertainty sigma_H, centers c1
and c2, and FWHM.

Usage: python auto_wedge_fit_model.py <output.png>
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

SIGMA_TO_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))
BACKGROUND, HEIGHT, FWHM, CENTERS = 0.20, 0.15, 20.0, (180.0, 0.0)  # 1/cm, 1/cm, deg, deg
INK, MUTED, FIT_COLOR, DATA_COLOR = "#1f1f1f", "#6b6b6b", "#2a6fdb", "#3d3d3d"


def model(phi: np.ndarray) -> np.ndarray:
    """Flat background plus two Gaussian peaks, periodic in phi (deg)."""
    sigma = FWHM / SIGMA_TO_FWHM
    total = np.full_like(phi, BACKGROUND)
    for center in CENTERS:
        distance = (phi - center + 180.0) % 360.0 - 180.0
        total += HEIGHT * np.exp(-0.5 * (distance / sigma) ** 2)
    return total


def draw(output: str) -> None:
    """Draw the figure and save it to ``output``."""
    rng = np.random.default_rng(7)
    phi = np.arange(-89.5, 270.0, 2.0)  # 2 deg bins over the 360 deg fit window
    error = np.full_like(phi, 0.008)
    data = model(phi) + rng.normal(0.0, error)
    fine = np.linspace(-90.0, 270.0, 2000)

    fig, ax = plt.subplots(figsize=(9.0, 4.8), dpi=150)
    ax.errorbar(phi, data, yerr=error, fmt="o", ms=3.5, color=DATA_COLOR, ecolor="#9a9a9a", elinewidth=1,
                label="I(φ) of one Q-ring")  # fmt: skip
    ax.plot(fine, model(fine), color=FIT_COLOR, lw=2, label="fit: A₀ + two Gaussians")
    ax.axhline(BACKGROUND, color=MUTED, lw=1.2, ls="--")
    ax.annotate(r"flat background $A_0$", xy=(95, BACKGROUND), xytext=(95, 0.158), color=INK, ha="center",
                fontsize=11, arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))  # fmt: skip

    center, top, half = CENTERS[0], BACKGROUND + HEIGHT, BACKGROUND + HEIGHT / 2
    # height H, from the background to the top of the peak
    ax.annotate("", xy=(center + 26, top), xytext=(center + 26, BACKGROUND),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2))  # fmt: skip
    ax.plot([center, center + 30], [top, top], color=MUTED, lw=0.8, ls=":")
    ax.text(center + 30, half, r"height $H \pm \sigma_H$", color=INK, va="center", fontsize=11)
    # center and FWHM
    ax.axvline(center, ymax=0.81, color=MUTED, lw=0.8, ls=":")
    ax.text(center, top + 0.022, r"center $c_1$", color=INK, ha="center", fontsize=11)
    ax.annotate("", xy=(center - FWHM / 2, half), xytext=(center + FWHM / 2, half),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2))  # fmt: skip
    ax.text(center - FWHM / 2 - 3, half, "FWHM", color=INK, ha="right", va="center", fontsize=11)
    # second peak, 180 deg away
    ax.text(CENTERS[1], top + 0.022, r"center $c_2$" + "\n" + r"($c_1 + 180°$ in symmetric mode)", color=INK,
            ha="center", fontsize=10)  # fmt: skip

    ax.set_xlim(-90, 270)
    ax.set_ylim(0.14, 0.42)
    ax.set_xticks(np.arange(-90, 271, 45))
    ax.set_xlabel("azimuthal angle φ (deg)", fontsize=12)
    ax.set_ylabel("intensity (1/cm)", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e5e5e5", lw=0.8)
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", frameon=False, fontsize=10)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output", help="path of the PNG file to write")
    draw(parser.parse_args().output)
