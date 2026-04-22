#!/usr/bin/env python3
"""
Recreate the LOESS animated GIFs for Week 8/9 lectures.

Generates two GIFs:
  1. week_08_S121_loess_animation.gif
     - Single-panel animation showing a focal point sweeping across x,
       with tri-cube weighted points and a local regression line.

  2. week_08_S122_loess_multi_span_animation.gif
     - 2x2 faceted animation comparing four span values (0.1, 0.15, 0.25, 0.66)
       as the focal point sweeps simultaneously across all panels.

Data: Simulated non-linear trend (sinusoidal + noise), matching the
      Week 9 lecture code: margin = 0.05*sin(day/25) + 0.02 + noise

Usage:
    python3 scripts/loess_animation.py

Output:
    Lecture_Folder/images/week_08_S121_loess_animation.gif
    Lecture_Folder/images/week_08_S122_loess_multi_span_animation.gif
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import io
import os

# ── Settings ──────────────────────────────────────────────────────────
SEED = 42
N_POINTS = 150
FPS = 8
DPI = 120
HOLD_FRAMES_END = 12

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          '..', 'Lecture_Folder', 'images')

# ── Generate data (matches Week 9 R code exactly) ────────────────────
np.random.seed(SEED)
day = np.random.uniform(-155, 0, N_POINTS)
margin = 0.05 * np.sin(day / 25) + 0.02 + np.random.normal(0, 0.025, N_POINTS)


def tricube_weights(x_data, focal, span_frac, n):
    """Compute tri-cube kernel weights for LOESS."""
    dists = np.abs(x_data - focal)
    h = np.sort(dists)[int(np.ceil(span_frac * n)) - 1]
    if h == 0:
        h = 1e-10
    u = dists / h
    w = np.where(u < 1, (1 - u**3)**3, 0.0)
    w_max = w.max()
    if w_max > 0:
        w = w / w_max
    return w, h


def local_regression(x_data, y_data, weights, focal, h):
    """Fit weighted linear regression and return predicted line."""
    # Weighted least squares: y = a + b*x
    W = np.diag(weights)
    X = np.column_stack([np.ones_like(x_data), x_data])
    try:
        beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ y_data)
    except np.linalg.LinAlgError:
        return None, None

    # Prediction range (local neighborhood)
    x_range = np.linspace(focal - h * 0.5, focal + h * 0.5, 50)
    x_range = x_range[(x_range >= x_data.min()) & (x_range <= x_data.max())]
    y_pred = beta[0] + beta[1] * x_range
    return x_range, y_pred


def fig_to_image(fig, dpi=DPI):
    """Convert matplotlib figure to PIL Image."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, facecolor='white',
                bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).copy()


# ── Focal point positions for animation ──────────────────────────────
focal_points = np.linspace(-155, 0, 40)


# ======================================================================
#  GIF 1: Single-panel LOESS animation (span = 0.3)
# ======================================================================
def make_single_frame(x0, span=0.3):
    fig, ax = plt.subplots(figsize=(5.5, 5))
    fig.patch.set_facecolor('white')

    w, h = tricube_weights(day, x0, span, N_POINTS)

    # Point sizes and alpha based on weight
    sizes = 8 + 80 * w
    alphas = 0.15 + 0.85 * w
    colors = np.array([[0, 0, 0, a] for a in alphas])

    ax.scatter(day, margin, s=sizes, c=colors, edgecolors='none', zorder=2)

    # Local regression line
    x_line, y_line = local_regression(day, margin, w, x0, h)
    if x_line is not None:
        ax.plot(x_line, y_line, color='blue', linewidth=3, zorder=3)

    ax.set_xlabel('day', fontsize=13)
    ax.set_ylabel('margin', fontsize=13)
    ax.set_title(f'x0 = {x0:6.0f}', fontsize=14, fontweight='bold',
                 family='monospace')
    ax.set_xlim(-165, 10)
    ax.set_ylim(-0.06, 0.13)
    ax.tick_params(labelsize=11)

    # Weight legend
    legend_weights = [0.0, 0.25, 0.50, 0.75, 1.00]
    legend_x = 0.88
    ax.text(legend_x + 0.06, 0.97, 'weight', transform=ax.transAxes,
            fontsize=10, fontweight='bold', ha='center', va='top')
    for i, lw in enumerate(legend_weights):
        y_pos = 0.90 - i * 0.065
        sz = 8 + 80 * lw
        alpha = 0.15 + 0.85 * lw
        ax.scatter([legend_x], [y_pos], s=sz, c=[[0, 0, 0, alpha]],
                   transform=ax.transAxes, edgecolors='none', clip_on=False)
        ax.text(legend_x + 0.12, y_pos, f'{lw:.2f}', transform=ax.transAxes,
                fontsize=9, va='center')

    plt.tight_layout()
    return fig_to_image(fig)


print("Generating GIF 1: Single-panel LOESS animation...")
frames1 = []
for i, x0 in enumerate(focal_points):
    frames1.append(make_single_frame(x0))
    if (i + 1) % 10 == 0:
        print(f"  Frame {i + 1}/{len(focal_points)}")

# Hold final frame
for _ in range(HOLD_FRAMES_END):
    frames1.append(frames1[-1].copy())

path1 = os.path.join(OUTPUT_DIR, 'week_08_S121_loess_animation.gif')
frames1[0].save(path1, save_all=True, append_images=frames1[1:],
                duration=int(1000 / FPS), loop=0, optimize=True)
print(f"Saved: {path1} ({os.path.getsize(path1) / 1024:.0f} KB)")


# ======================================================================
#  GIF 2: Multi-span LOESS animation (2x2 facets)
# ======================================================================
SPANS = [0.1, 0.15, 0.25, 0.66]


def make_multi_frame(x0):
    fig, axes = plt.subplots(2, 2, figsize=(7, 6.5),
                             sharex=True, sharey=True)
    fig.patch.set_facecolor('white')

    for ax, span in zip(axes.flat, SPANS):
        w, h = tricube_weights(day, x0, span, N_POINTS)

        sizes = 8 + 80 * w
        alphas = 0.15 + 0.85 * w
        colors = np.array([[0, 0, 0, a] for a in alphas])

        ax.scatter(day, margin, s=sizes, c=colors, edgecolors='none', zorder=2)

        x_line, y_line = local_regression(day, margin, w, x0, h)
        if x_line is not None:
            ax.plot(x_line, y_line, color='blue', linewidth=2.5, zorder=3)

        # Span label in a header bar (mimicking ggplot2 facet strip)
        ax.set_title(f'{span}', fontsize=12, fontweight='bold',
                     backgroundcolor='#EBEBEB', pad=8)
        ax.set_xlim(-165, 10)
        ax.set_ylim(-0.06, 0.13)
        ax.tick_params(labelsize=9)

    # Shared axis labels
    fig.text(0.5, 0.02, 'day', ha='center', fontsize=13)
    fig.text(0.02, 0.5, 'margin', va='center', rotation='vertical', fontsize=13)

    # Main title
    fig.suptitle(f'x0 = {x0:6.0f}', fontsize=14, fontweight='bold',
                 family='monospace', y=0.99)

    # Weight legend (right side)
    legend_weights = [0.0, 0.25, 0.50, 0.75, 1.00]
    legend_ax = fig.add_axes([0.92, 0.35, 0.08, 0.35])
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.axis('off')

    legend_ax.text(0.5, 0.98, 'weight', fontsize=10, fontweight='bold',
                   ha='center', va='top')
    for i, lw in enumerate(legend_weights):
        y_pos = 0.85 - i * 0.18
        sz = 8 + 80 * lw
        alpha = 0.15 + 0.85 * lw
        legend_ax.scatter([0.2], [y_pos], s=sz, c=[[0, 0, 0, alpha]],
                          edgecolors='none', clip_on=False)
        legend_ax.text(0.55, y_pos, f'{lw:.2f}', fontsize=9, va='center')

    plt.subplots_adjust(left=0.08, right=0.90, top=0.92, bottom=0.08,
                        hspace=0.25, wspace=0.12)
    return fig_to_image(fig)


print("\nGenerating GIF 2: Multi-span LOESS animation...")
frames2 = []
for i, x0 in enumerate(focal_points):
    frames2.append(make_multi_frame(x0))
    if (i + 1) % 10 == 0:
        print(f"  Frame {i + 1}/{len(focal_points)}")

for _ in range(HOLD_FRAMES_END):
    frames2.append(frames2[-1].copy())

path2 = os.path.join(OUTPUT_DIR, 'week_08_S122_loess_multi_span_animation.gif')
frames2[0].save(path2, save_all=True, append_images=frames2[1:],
                duration=int(1000 / FPS), loop=0, optimize=True)
print(f"Saved: {path2} ({os.path.getsize(path2) / 1024:.0f} KB)")

print("\nDone!")
