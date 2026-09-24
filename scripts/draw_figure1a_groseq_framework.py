#!/usr/bin/env python3
"""Draw the project-specific Figure 1a workflow at publication resolution.

The design is original but preserves the broad panel proportions of the
reference article: a full-height input panel, three upper workflow panels,
and a lower model-architecture strip.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle
import numpy as np


GREEN = "#188A68"
GREEN_LIGHT = "#74CFAE"
PURPLE = "#6F5BD3"
PURPLE_LIGHT = "#A899F2"
INK = "#222B35"
MUTED = "#66717D"
ORANGE = "#E39B3A"
BLUE = "#4B84C4"
PANEL = "#F8FAFC"
LINE = "#44505C"


def rounded_box(ax, x, y, w, h, facecolor="white", edgecolor=LINE, lw=1.05, radius=0.012, z=1):
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.005,rounding_size={radius}",
        linewidth=lw, edgecolor=edgecolor, facecolor=facecolor, zorder=z,
    )
    ax.add_patch(box)
    return box


def arrow(ax, x1, y1, x2, y2, color=INK, lw=1.5, scale=10):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=scale,
        linewidth=lw, color=color, shrinkA=0, shrinkB=0, zorder=8,
    ))


def section_title(ax, x, y, w, text):
    ax.text(x + w / 2, y, text, ha="center", va="top", fontsize=6.6,
            fontweight="bold", color=INK, linespacing=.92, zorder=10)


def draw_database(ax, cx, cy, w=0.045, h=0.075):
    x = cx - w / 2
    ax.add_patch(Rectangle((x, cy - h / 2), w, h, facecolor="#DCE3EA", edgecolor="#7E8B98", lw=0.7))
    for yy in (cy - h / 2, cy, cy + h / 2):
        ax.add_patch(Ellipse((cx, yy), w, h * 0.20, facecolor="#EAF0F4", edgecolor="#7E8B98", lw=0.7))
    ax.add_patch(Ellipse((cx, cy - h / 2), w, h * 0.20, facecolor="#D3DCE4", edgecolor="#7E8B98", lw=0.7))


def draw_plant(ax, cx, cy, scale=1.0, kind="arabidopsis"):
    if kind == "arabidopsis":
        ax.plot([cx, cx], [cy - 0.022 * scale, cy + 0.026 * scale], color="#477A3D", lw=1.0, zorder=5)
        for dx, dy, angle in [(-.014, -.012, 25), (.014, -.008, -25), (-.012, .006, 35), (.012, .012, -35)]:
            ax.add_patch(Ellipse((cx + dx * scale, cy + dy * scale), .020 * scale, .009 * scale,
                                 angle=angle, facecolor="#6BAE55", edgecolor="#477A3D", lw=.4, zorder=5))
        for dx in (-.008, 0, .008):
            ax.add_patch(Circle((cx + dx * scale, cy + .030 * scale), .0035 * scale,
                                facecolor="#F2C94C", edgecolor="#AD812A", lw=.3, zorder=6))
    elif kind == "wheat":
        ax.plot([cx, cx], [cy - .030 * scale, cy + .032 * scale], color="#8A7836", lw=1.0, zorder=5)
        for i in range(5):
            yy = cy + (-.004 + i * .008) * scale
            side = -1 if i % 2 == 0 else 1
            ax.add_patch(Ellipse((cx + side * .006 * scale, yy), .015 * scale, .0055 * scale,
                                 angle=side * 25, facecolor="#D3AE42", edgecolor="#8A7836", lw=.35, zorder=5))
        ax.plot([cx, cx - .016 * scale], [cy - .013 * scale, cy + .002 * scale], color="#7B9846", lw=.8)
        ax.plot([cx, cx + .016 * scale], [cy - .004 * scale, cy + .010 * scale], color="#7B9846", lw=.8)
    else:
        ax.plot([cx, cx], [cy - .034 * scale, cy + .030 * scale], color="#3C7C43", lw=1.4, zorder=5)
        for side, yy, angle in [(-1, -.018, 25), (1, -.008, -28), (-1, .003, 32), (1, .012, -32)]:
            ax.add_patch(Ellipse((cx + side * .012 * scale, cy + yy * scale), .032 * scale, .008 * scale,
                                 angle=side * angle, facecolor="#59A85C", edgecolor="#3C7C43", lw=.4, zorder=5))
        for dx in (-.008, -.004, 0, .004, .008):
            ax.plot([cx, cx + dx * scale], [cy + .028 * scale, cy + .040 * scale], color="#7C6E37", lw=.45)


def signal_profile(n=90, seed=1):
    rng = np.random.default_rng(seed)
    x = np.linspace(0, 1, n)
    y = np.zeros_like(x)
    for center, height, width in zip(rng.uniform(.05, .95, 7), rng.uniform(.25, 1, 7), rng.uniform(.012, .06, 7)):
        y += height * np.exp(-0.5 * ((x - center) / width) ** 2)
    y += rng.uniform(0, .035, n)
    return x, y / max(y.max(), 1e-6)


def draw_track(ax, x, y, w, h, color, seed, baseline=True, signed=False):
    xx, yy = signal_profile(100, seed)
    xx = x + xx * w
    yy = yy * h
    if signed:
        yy = -yy
    base = y
    ax.fill_between(xx, base, base + yy, color=color, alpha=.88, lw=0, zorder=4)
    if baseline:
        ax.plot([x, x + w], [base, base], color="#87919C", lw=.35, zorder=5)


def draw_one_hot(ax, x, y, w, h):
    colors = ["#EB6B62", "#4FA9D8", "#66B96B", "#E4B647"]
    seq = [0, 2, 2, 1, 3, 0, 1, 2, 3, 3, 0, 2]
    cw, ch = w / len(seq), h / 4
    for col, base in enumerate(seq):
        for row in range(4):
            fc = colors[row] if row == base else "#E8EDF1"
            ax.add_patch(Rectangle((x + col * cw, y + (3 - row) * ch), cw * .93, ch * .88,
                                   facecolor=fc, edgecolor="white", lw=.2, zorder=4))
    for row, label in enumerate("ATCG"):
        ax.text(x - .006, y + (3.5 - row) * ch, label, ha="right", va="center",
                fontsize=4.2, color=MUTED)


def small_step(ax, x, y, w, text, color):
    rounded_box(ax, x, y, w, .050, facecolor=color, edgecolor="white", lw=.5, radius=.006, z=3)
    ax.text(x + w / 2, y + .025, text, ha="center", va="center", fontsize=5.0,
            color=INK, fontweight="bold", zorder=5)


def draw_arch_box(ax, x, y, w, h, title, lines, color):
    rounded_box(ax, x, y, w, h, facecolor=color, edgecolor="#74808C", lw=.75, radius=.007, z=2)
    ax.text(x + w / 2, y + h - .015, title, ha="center", va="top", fontsize=4.9,
            fontweight="bold", color=INK, linespacing=.92, zorder=5)
    ax.text(x + w / 2, y + h * .39, lines, ha="center", va="center", fontsize=4.6,
            color="#3B4650", linespacing=1.12, zorder=5)


def build_figure():
    fig = plt.figure(figsize=(7.2, 3.9), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Panel layout follows the approximate aspect and occupancy of the source Fig. 1a.
    left = (0.030, 0.055, 0.205, 0.885)
    prep = (0.258, 0.390, 0.242, 0.550)
    train = (0.523, 0.390, 0.205, 0.550)
    predict = (0.751, 0.390, 0.219, 0.550)
    arch = (0.258, 0.055, 0.712, 0.295)

    ax.text(.012, .968, "a", fontsize=10, fontweight="bold", color=INK, ha="left", va="top")
    for panel in (left, prep, train, predict, arch):
        rounded_box(ax, *panel, facecolor="white", edgecolor=LINE, lw=1.05, radius=.014)

    # 1. Input data.
    x, y, w, h = left
    section_title(ax, x, y + h - .018, w, "Public GRO-seq\ndata")
    draw_database(ax, x + .043, y + h - .105, .050, .075)
    ax.text(x + .078, y + h - .095, "Wild-type / untreated", ha="left", va="center",
            fontsize=5.4, color=MUTED)
    ax.text(x + .078, y + h - .121, "nascent RNA profiles", ha="left", va="center",
            fontsize=5.4, color=MUTED)
    species = [
        ("A. thaliana", "arabidopsis", GREEN, 10),
        ("T. aestivum", "wheat", BLUE, 20),
        ("Z. mays", "maize", ORANGE, 30),
    ]
    row_y = [y + .625, y + .385, y + .145]
    for (label, kind, accent, seed), cy in zip(species, row_y):
        rounded_box(ax, x + .014, cy - .090, w - .028, .185, facecolor="#FBFCFD",
                    edgecolor="#D7DEE5", lw=.55, radius=.008)
        draw_plant(ax, x + .046, cy + .035, .88, kind)
        ax.text(x + .076, cy + .044, label, fontsize=5.6, fontstyle="italic", fontweight="bold",
                color=INK, va="center")
        draw_track(ax, x + .025, cy - .008, w - .050, .040, GREEN, seed)
        draw_track(ax, x + .025, cy - .035, w - .050, .034, PURPLE, seed + 1, signed=True)
        ax.text(x + .023, cy - .078, "+ / − strand GRO-seq", fontsize=4.4, color=MUTED, va="center")

    # 2. Strand-resolved preprocessing.
    x, y, w, h = prep
    section_title(ax, x, y + h - .018, w, "Strand-resolved\npreprocessing")
    sx, sy, sw = x + .014, y + h - .120, (w - .046) / 3
    small_step(ax, sx, sy, sw, "QC + trim", "#E8F1F7")
    arrow(ax, sx + sw, sy + .025, sx + sw + .012, sy + .025, color="#7B8792", lw=.8, scale=6)
    small_step(ax, sx + sw + .014, sy, sw, "rRNA filter", "#FCEBD8")
    arrow(ax, sx + 2 * sw + .014, sy + .025, sx + 2 * sw + .026, sy + .025, color="#7B8792", lw=.8, scale=6)
    small_step(ax, sx + 2 * sw + .028, sy, sw, "Bowtie2", "#E7F4EC")
    ax.text(x + w / 2, sy - .020, "SRR runs combined within SRX", ha="center", va="top",
            fontsize=4.7, color=MUTED)
    # Replicate tracks and consensus interval.
    base_y = y + .298
    for i, seed in enumerate((41, 44, 47)):
        draw_track(ax, x + .026, base_y + i * .037, w - .052, .025, GREEN_LIGHT, seed)
    ax.text(x + .017, base_y + .045, "replicates", rotation=90, fontsize=4.1, color=MUTED,
            va="center", ha="center")
    ax.add_patch(Rectangle((x + .092, base_y - .022), .075, .013, facecolor=GREEN,
                           edgecolor="none", alpha=.95, zorder=5))
    ax.text(x + .1295, base_y - .035, "strand-specific consensus", ha="center", va="top",
            fontsize=4.4, color=INK)
    # Window and one-hot encoding.
    win_y = y + .085
    ax.plot([x + .035, x + w - .035], [win_y + .095, win_y + .095], color="#8A949E", lw=.7)
    ax.add_patch(Rectangle((x + .075, win_y + .080), .093, .030, facecolor="#FFF2C7",
                           edgecolor=ORANGE, lw=.7, zorder=4))
    ax.text(x + .1215, win_y + .119, "1,024-bp window", ha="center", va="bottom",
            fontsize=4.6, color=INK, fontweight="bold")
    arrow(ax, x + .1215, win_y + .076, x + .1215, win_y + .060, color="#7B8792", lw=.8, scale=6)
    draw_one_hot(ax, x + .045, win_y - .005, w - .090, .052)
    ax.text(x + w / 2, win_y - .019, "one-hot encoded DNA", ha="center", va="top",
            fontsize=4.7, color=MUTED)

    # 3. Species-specific training.
    x, y, w, h = train
    section_title(ax, x, y + h - .018, w, "Species-specific\nlearning")
    ax.text(x + w / 2, y + h - .090, "chromosome-level train / validation / test split",
            ha="center", va="center", fontsize=4.5, color=MUTED)
    centers = [y + .365, y + .235, y + .105]
    labels = [("Arabidopsis", "arabidopsis", "Model A"), ("Wheat", "wheat", "Model W"), ("Maize", "maize", "Model M")]
    for cy, (label, kind, model_label) in zip(centers, labels):
        draw_plant(ax, x + .040, cy, .65, kind)
        arrow(ax, x + .060, cy, x + .092, cy, color="#7B8792", lw=.9, scale=7)
        rounded_box(ax, x + .096, cy - .040, .082, .080, facecolor="#EDF3F8",
                    edgecolor=BLUE, lw=.7, radius=.007)
        ax.text(x + .137, cy + .010, model_label, ha="center", va="center", fontsize=5.4,
                fontweight="bold", color=INK)
        ax.text(x + .137, cy - .014, label, ha="center", va="center", fontsize=4.2,
                color=MUTED)
    ax.text(x + w / 2, y + .026, "independent two-output regression", ha="center", va="center",
            fontsize=4.7, color=INK, fontweight="bold")

    # 4. Genome-wide predictions.
    x, y, w, h = predict
    section_title(ax, x, y + h - .018, w, "Genome-wide strand\nprediction")
    ax.text(x + w / 2, y + h - .090, "DNA sequence only", ha="center", va="center",
            fontsize=4.8, color=MUTED)
    ax.text(x + .020, y + h - .115, "…ACGTTGCACTGATC…", ha="left", va="center",
            fontsize=4.7, family="monospace", color=INK)
    arrow(ax, x + .155, y + h - .115, x + .185, y + h - .115, color="#7B8792", lw=.8, scale=6)
    track_starts = [y + .340, y + .215, y + .090]
    for i, (cy, label) in enumerate(zip(track_starts, ("Arabidopsis", "Wheat", "Maize"))):
        ax.text(x + .014, cy + .036, label, ha="left", va="center", fontsize=4.5,
                color=INK, fontweight="bold")
        draw_track(ax, x + .018, cy, w - .036, .055, GREEN, 70 + i * 4)
        draw_track(ax, x + .018, cy, w - .036, .050, PURPLE, 72 + i * 4, signed=True)
    ax.text(x + w - .014, y + .024, "continuous y+ / y−", ha="right", va="center",
            fontsize=4.7, color=MUTED)

    # Main process arrows.
    arrow(ax, left[0] + left[2] + .004, .675, prep[0] - .006, .675, lw=1.8, scale=11)
    arrow(ax, prep[0] + prep[2] + .004, .675, train[0] - .006, .675, lw=1.8, scale=11)
    arrow(ax, train[0] + train[2] + .004, .675, predict[0] - .006, .675, lw=1.8, scale=11)

    # Lower architecture strip.
    x, y, w, h = arch
    ax.text(x + w / 2, y + h - .018, "Prediction process and model architecture",
            ha="center", va="top", fontsize=7.2, fontweight="bold", color=INK)
    inner_y, inner_h = y + .030, h - .078
    gap = .012
    widths = [.088, .090, .144, .126, .088, .098, .094]
    usable = w - .030
    scale = (usable - gap * (len(widths) - 1)) / sum(widths)
    widths = [value * scale for value in widths]
    xx = x + .015
    specs = [
        ("DNA sequence", "1,024 bp\nA / T / C / G", "#F4F6F8"),
        ("One-hot", "4 × 1,024\nchannels", "#E8F2FA"),
        ("Residual Conv", "3 blocks\nk=9\n480→640→960", "#EAF1FB"),
        ("Dilated Conv", "5 residual blocks\nd=2,4,8,16,25", "#FFF1D8"),
        ("B-spline", "position basis\n16 features", "#E8F5E8"),
        ("Regression\nhead", "Dense 256\nSigmoid", "#F9E6EC"),
        ("Outputs", "y+   y−\ncontinuous", "#EEE9FC"),
    ]
    for idx, ((title, lines, color), ww) in enumerate(zip(specs, widths)):
        draw_arch_box(ax, xx, inner_y, ww, inner_h, title, lines, color)
        if title == "Outputs":
            ax.add_patch(Rectangle((xx + .012, inner_y + .021), ww * .32, .020,
                                   facecolor=GREEN, edgecolor="none", zorder=5))
            ax.add_patch(Rectangle((xx + ww * .58, inner_y + .021), ww * .32, .020,
                                   facecolor=PURPLE, edgecolor="none", zorder=5))
        if idx < len(specs) - 1:
            arrow(ax, xx + ww + .002, inner_y + inner_h / 2, xx + ww + gap - .002,
                  inner_y + inner_h / 2, color="#65717D", lw=.9, scale=7)
        xx += ww + gap

    return fig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-prefix", default="outputs/figure1a_groseq_framework")
    args = parser.parse_args()
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    fig = build_figure()
    fig.savefig(prefix.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(prefix.with_suffix(".pdf"), facecolor="white")
    fig.savefig(prefix.with_suffix(".svg"), facecolor="white")
    tiff_path = prefix.with_suffix(".tiff")
    fig.savefig(tiff_path, dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    from PIL import Image
    # Matplotlib writes RGBA TIFFs. Journals generally expect an opaque RGB
    # figure, so flatten the already-white canvas while preserving 600-dpi
    # metadata and lossless LZW compression.
    with Image.open(tiff_path) as image:
        rgb = image.convert("RGB")
    rgb.save(tiff_path, dpi=(600, 600), compression="tiff_lzw")
    with Image.open(tiff_path) as image:
        print({
            "tiff": str(tiff_path),
            "pixels": image.size,
            "dpi": image.info.get("dpi"),
            "mode": image.mode,
            "compression": image.info.get("compression"),
        })


if __name__ == "__main__":
    main()
