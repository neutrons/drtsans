"""Draw wedge_binning_1.png, _2, _3, _5 and _6 for docs/user/binning.rst, section "Wedge Binning".

Each figure shows the detector plane on the left, with the wedges in azimuthal angle phi and the Q-bin boundaries
between q_m and q_M, and a fictitious intensity profile on the right. The symmetric wedge is the wedge rotated by
180 degrees, consistent with I(q, phi + pi) = I(q, phi).

Usage: python wedge_schematics.py <output directory>
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Arc, Circle, Polygon, Rectangle  # noqa: E402

INK, FRAME, SHADE = "#1f1f1f", "#1f2a44", "#fdf0c4"
PHI_MAX_COLOR, PHI_MIN_COLOR, QMIN_COLOR, QMAX_COLOR = "#e8955a", "#e6b800", "#1f6fbf", "#e0201a"
PHI_MIN, PHI_MAX = -22.0, 22.0  # deg, the wedge along X
BOUNDARIES = np.linspace(0.16, 0.94, 10)  # Q-bin boundaries in units of the detector half-width, q_m to q_M
Q_LABELS = [f"{q:.2f}" for q in np.arange(0.03, 0.125, 0.01)]  # 1/A, the boundaries of figures 5 and 6
ARROW = dict(arrowstyle="->", color=INK, lw=1.2)


def _polar(radius: float, angle: float) -> tuple[float, float]:
    """Cartesian coordinates of a point given its radius and its azimuthal angle (deg)."""
    return radius * np.cos(np.radians(angle)), radius * np.sin(np.radians(angle))


def detector(ax: plt.Axes) -> Rectangle:
    """Draw the detector frame and the X and Y axes, and return the frame to clip the drawings inside it."""
    frame = Rectangle((-1.0, -0.5), 2.0, 1.0, fill=False, ec=FRAME, lw=2.0, zorder=5)
    ax.add_patch(frame)
    ax.annotate("", xy=(1.1, 0.0), xytext=(-1.0, 0.0), arrowprops=ARROW)
    ax.annotate("", xy=(0.0, 0.6), xytext=(0.0, -0.5), arrowprops=ARROW)
    ax.text(1.1, -0.04, "X", fontsize=20, va="top", ha="center", color=INK)
    ax.text(0.03, 0.6, "Y", fontsize=20, va="center", color=INK)
    ax.set_xlim(-1.02, 1.16)
    ax.set_ylim(-0.52, 0.64)
    ax.set_aspect("equal")
    ax.axis("off")
    return frame


def wedge_edges(ax: plt.Axes, rotation: float = 0.0) -> None:
    """Draw the dashed edges of the wedge rotated by ``rotation`` (deg) and the Q-bin boundaries inside it."""
    for angle, color in [(PHI_MIN + rotation, PHI_MIN_COLOR), (PHI_MAX + rotation, PHI_MAX_COLOR)]:
        ax.plot(*zip((0.0, 0.0), _polar(1.08, angle)), color=color, ls="--", lw=1.5)
    for index, radius in enumerate(BOUNDARIES):
        color = QMIN_COLOR if index == 0 else QMAX_COLOR if index == len(BOUNDARIES) - 1 else INK
        ax.add_patch(Arc((0, 0), 2 * radius, 2 * radius, theta1=PHI_MIN + rotation, theta2=PHI_MAX + rotation,
                         color=color, lw=1.4))  # fmt: skip
    ax.text(*_polar(BOUNDARIES[0] + 0.02, -42 + rotation), "$q_m$", color=QMIN_COLOR, fontsize=16, ha="center",
            va="center")  # fmt: skip
    ax.text(*_polar(BOUNDARIES[-1] + 0.02, -25 + rotation), "$q_M$", color=QMAX_COLOR, fontsize=16, ha="center",
            va="center")  # fmt: skip
    sign = np.cos(np.radians(rotation))
    for index in range(3):  # q_i at the bin centers
        center = sign * (BOUNDARIES[index] + BOUNDARIES[index + 1]) / 2
        ax.plot([center, center], [-0.015, 0.015], color=INK, lw=1.2)
        ax.text(center, -0.03, f"$q_{index + 1}$", fontsize=14, ha="center", va="top", color=INK)


def phi_labels(ax: plt.Axes, rotation: float = 0.0, suffix: str = "", radius: float = 0.72) -> None:
    """Label the edges of the wedge rotated by ``rotation`` (deg)."""
    side = 1.0 if np.cos(np.radians(rotation)) > 0 else -1.0  # the labels sit outside the wedge
    for angle, name, color, shift in [(PHI_MAX, "M", PHI_MAX_COLOR, 6), (PHI_MIN, "m", PHI_MIN_COLOR, -6)]:
        ax.text(*_polar(radius, angle + rotation + side * shift), rf"$\phi_{name}${suffix}", color=color, fontsize=18,
                ha="center", va="center", bbox=dict(fc="white", ec="none", pad=1))  # fmt: skip


def shaded(ax: plt.Axes, frame: Rectangle, rotations: list[float], q_angle: float) -> None:
    """Draw the dotted Q-ring boundaries of figures 5 and 6 and the shaded wedges, clipped to the frame."""
    for rotation in rotations:
        corners = [(0.0, 0.0), _polar(2.0, PHI_MIN + rotation), _polar(2.0, PHI_MAX + rotation)]
        ax.add_patch(Polygon(corners, color=SHADE, lw=0, zorder=0)).set_clip_path(frame)
    for index, radius in enumerate(BOUNDARIES):
        color = QMIN_COLOR if index == 0 else QMAX_COLOR if index == len(BOUNDARIES) - 1 else INK
        ax.add_patch(Circle((0, 0), radius, fill=False, ec=color, ls=":", lw=1.2)).set_clip_path(frame)
        if index < 4:
            ax.text(*_polar(radius, q_angle), Q_LABELS[index], color=color, fontsize=13, va="center",
                    ha="left" if np.cos(np.radians(q_angle)) > 0 else "right",
                    bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.8))  # fmt: skip
    corner = (0.92, -0.44) if np.cos(np.radians(q_angle)) > 0 else (-0.92, -0.44)
    ax.text(*corner, Q_LABELS[-1], color=QMAX_COLOR, fontsize=13, ha="center", va="center")


def profile_panel(ax: plt.Axes, kind: str) -> None:
    """Draw a fictitious I(phi | q_i) profile for each q_i (``kind="phi"``) or a fictitious I(q) (``kind="q"``)."""
    ax.annotate("", xy=(1.12, 0.0), xytext=(-0.05, 0.0), arrowprops=ARROW)
    ax.set_xlim(-0.08, 1.15)
    ax.set_ylim(-0.25, 1.15)
    ax.axis("off")
    if kind == "phi":
        ax.annotate("", xy=(0.5, 1.08), xytext=(0.5, -0.04), arrowprops=ARROW)
        phi = np.linspace(0.0, 1.0, 200)
        for height in (0.65, 0.75, 0.85):
            ax.plot(phi, 0.02 + height * np.exp(-0.5 * ((phi - 0.5) / 0.14) ** 2), color=INK, lw=1.4)
        ax.text(0.55, 1.06, r"$I(\phi\,|\,q_i)$", fontsize=18, va="center", color=INK)
        for position, name in [(0.0, "m"), (1.0, "M")]:
            ax.plot([position, position], [-0.025, 0.025], color=INK, lw=1.4)
            ax.text(position, -0.07, rf"$\phi_{name}$", fontsize=18, ha="center", va="top", color=INK)
        return
    ax.annotate("", xy=(0.0, 1.08), xytext=(0.0, -0.04), arrowprops=ARROW)
    ax.text(0.08, 1.06, "$I(q)$", fontsize=18, va="center", color=INK)
    edges = np.linspace(0.1, 1.0, 10)
    q = np.linspace(edges[0] + 0.05, edges[-1] - 0.03, 200)
    intensity = 0.08 + 0.82 / (1.0 + np.exp((q - 0.72) / 0.06))  # fictitious I(q) in the wedge
    ax.plot(q, intensity, color="#b8860b", lw=2.2)
    for index in range(3):
        center = edges[0] + 0.05 + 0.13 * index
        ax.plot([center, center], [0.0, np.interp(center, q, intensity)], color=INK, ls=":", lw=1)
        ax.text(center, -0.04, f"$q_{index + 1}$", fontsize=15, ha="center", va="top", color=INK)
    for position, name, color in [(edges[0], "m", QMIN_COLOR), (edges[-1], "M", QMAX_COLOR)]:
        ax.plot([position, position], [-0.04, 0.04], color=color, lw=3)
        ax.text(position, -0.11, f"$q_{name}$", fontsize=15, ha="center", va="top", color=color)


def draw(number: int, output: str) -> None:
    """Draw figure ``number`` and save it to ``output``."""
    fig, (left, right) = plt.subplots(1, 2, figsize=(13.0, 4.4), dpi=200, gridspec_kw=dict(width_ratios=[2.2, 1]))
    frame = detector(left)
    if number in (1, 2, 3):
        wedge_edges(left)
        phi_labels(left)
    if number == 1:  # the azimuthal angle phi of a point in the wedge
        left.plot(*zip((0.0, 0.0), _polar(0.82, 12.0)), color="#4c9a2a", lw=1.4)
        left.add_patch(Arc((0, 0), 1.7, 1.7, theta1=0, theta2=12, color="#4c9a2a", lw=1.4))
        left.annotate("", xy=_polar(0.85, 12.5), xytext=_polar(0.85, 10.0),
                      arrowprops=dict(arrowstyle="->", color="#4c9a2a", lw=1.4))  # fmt: skip
        left.text(*_polar(0.91, 6.0), r"$\phi$", color="#4c9a2a", fontsize=18, ha="center", va="center")
    if number == 3:  # the symmetric wedge, rotated by 180 degrees
        wedge_edges(left, rotation=180.0)
        phi_labels(left, rotation=180.0, suffix=" + 180°", radius=0.68)
    if number == 5:
        shaded(left, frame, [0.0, 180.0], q_angle=-55.0)
        for angle, name in [(PHI_MAX + 4, "M"), (PHI_MIN - 4, "m")]:
            left.text(*_polar(0.6, angle), rf"$\phi_{name}$?", fontsize=18, ha="center", va="center", color=INK)
    if number == 6:
        shaded(left, frame, [180.0], q_angle=-120.0)
        for angle, name in [(180.0 + PHI_MIN - 6, "m"), (180.0 + PHI_MAX + 6, "M")]:
            left.text(*_polar(0.6, angle), rf"$\phi_{name}$?", fontsize=18, ha="center", va="center", color=INK)
        left.add_patch(Rectangle((0.0, -0.5), 1.0, 1.0, fc="#d9d9d9", ec="#9a9a9a", hatch="///", alpha=0.9,
                                 lw=0, zorder=3))  # fmt: skip
        left.text(0.5, 0.0, "masked region", fontsize=22, style="italic", ha="center", va="center", color=INK,
                  zorder=4, bbox=dict(fc="#d9d9d9", ec="none"))  # fmt: skip
    profile_panel(right, "phi" if number == 1 else "q")
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_dir", help="directory where to write the PNG files")
    directory = parser.parse_args().output_dir
    for figure_number in (1, 2, 3, 5, 6):
        draw(figure_number, os.path.join(directory, f"wedge_binning_{figure_number}.png"))
