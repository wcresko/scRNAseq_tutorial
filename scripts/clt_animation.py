#!/usr/bin/env python3
"""
Generate an animated GIF illustrating the Central Limit Theorem.

Shows repeated sampling from a skewed (exponential) population:
- Left panel: the population distribution with current sample highlighted
- Right panel: histogram of sample means building up over time
- Running estimate of the grand mean and SEM displayed and updated

The animation demonstrates how the distribution of sample means
becomes normal regardless of the population shape.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image
import io
import os

# ── Settings ──────────────────────────────────────────────────────────
np.random.seed(42)
N_POP = 100_000
SAMPLE_SIZE = 25
N_SAMPLES = 200          # total samples to draw
FRAMES_PER_SAMPLE = 1    # 1 frame per sample for smooth animation
HOLD_FRAMES_END = 30     # hold final frame
FPS = 12
DPI = 150
FIG_W, FIG_H = 10, 5

# Colors (matching course palette)
POP_COLOR = '#B0B0B0'
SAMPLE_HIGHLIGHT = '#E74C3C'
MEAN_HIST_COLOR = '#3498DB'
MEAN_LINE_COLOR = '#E74C3C'
SEM_COLOR = '#2ECC71'
TRUE_MEAN_COLOR = '#F39C12'

# ── Generate population ──────────────────────────────────────────────
population = np.random.exponential(scale=2.0, size=N_POP)
true_mean = np.mean(population)

# ── Precompute all samples and running stats ─────────────────────────
all_means = []
running_grand_mean = []
running_sem = []

for i in range(N_SAMPLES):
    samp = np.random.choice(population, size=SAMPLE_SIZE, replace=False)
    m = np.mean(samp)
    all_means.append(m)
    running_grand_mean.append(np.mean(all_means))
    if len(all_means) > 1:
        running_sem.append(np.std(all_means, ddof=1) / np.sqrt(len(all_means)))
    else:
        running_sem.append(0)

# Determine histogram bin edges for consistency across frames
hist_min = min(all_means) - 0.2
hist_max = max(all_means) + 0.2
bins = np.linspace(hist_min, hist_max, 35)

# ── Frame indices: show samples at 1,2,3,...,10, then every 2 up to 30,
#    then every 5 up to 100, then every 10 to 200
frame_indices = list(range(1, 11))
frame_indices += list(range(12, 31, 2))
frame_indices += list(range(35, 101, 5))
frame_indices += list(range(110, N_SAMPLES + 1, 10))
# Always include the final frame
if N_SAMPLES not in frame_indices:
    frame_indices.append(N_SAMPLES)

# ── Generate frames ──────────────────────────────────────────────────
frames = []

def make_frame(n_drawn):
    """Create a single frame showing n_drawn samples taken so far."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(FIG_W, FIG_H),
                                    gridspec_kw={'width_ratios': [1, 1.3]})
    fig.patch.set_facecolor('white')

    # ── Left panel: Population distribution ──────────────────────────
    ax1.hist(population, bins=60, density=True, color=POP_COLOR,
             edgecolor='white', linewidth=0.3, alpha=0.85)
    ax1.axvline(true_mean, color=TRUE_MEAN_COLOR, linewidth=2,
                linestyle='--', label=f'True μ = {true_mean:.2f}')

    # Highlight the most recent sample as red ticks on x-axis
    if n_drawn > 0:
        latest_sample = np.random.RandomState(42 + n_drawn).choice(
            population, size=SAMPLE_SIZE, replace=False)
        ax1.plot(latest_sample, np.zeros_like(latest_sample) - 0.01,
                 '|', color=SAMPLE_HIGHLIGHT, markersize=12, markeredgewidth=2,
                 alpha=0.8)
        latest_mean = all_means[n_drawn - 1]
        ax1.axvline(latest_mean, color=SAMPLE_HIGHLIGHT, linewidth=2,
                    alpha=0.7, label=f'Sample mean = {latest_mean:.2f}')

    ax1.set_xlabel('Value', fontsize=12)
    ax1.set_ylabel('Density', fontsize=12)
    ax1.set_title('Population (Exponential)', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=9, loc='upper right', framealpha=0.9)
    ax1.set_xlim(-0.5, 15)
    ax1.set_ylim(-0.03, 0.55)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # ── Right panel: Growing histogram of sample means ───────────────
    if n_drawn > 0:
        current_means = all_means[:n_drawn]
        ax2.hist(current_means, bins=bins, density=False,
                 color=MEAN_HIST_COLOR, edgecolor='white', linewidth=0.5,
                 alpha=0.85)

        gm = running_grand_mean[n_drawn - 1]
        sem = running_sem[n_drawn - 1]

        # Grand mean line
        ax2.axvline(gm, color=MEAN_LINE_COLOR, linewidth=2.5,
                    label=f'Estimate: x̄ = {gm:.3f}')
        # True mean for reference
        ax2.axvline(true_mean, color=TRUE_MEAN_COLOR, linewidth=2,
                    linestyle='--', alpha=0.7, label=f'True μ = {true_mean:.2f}')

        # SEM bracket
        if n_drawn > 1:
            ymax = ax2.get_ylim()[1]
            bracket_y = ymax * 0.88
            ax2.annotate('', xy=(gm - sem, bracket_y), xytext=(gm + sem, bracket_y),
                         arrowprops=dict(arrowstyle='<->', color=SEM_COLOR, lw=2.5))
            ax2.text(gm, bracket_y + ymax * 0.04,
                     f'SEM = {sem:.3f}',
                     ha='center', fontsize=11, fontweight='bold',
                     color=SEM_COLOR,
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                               edgecolor=SEM_COLOR, alpha=0.9))

        ax2.legend(fontsize=9, loc='upper right', framealpha=0.9)

    ax2.set_xlabel('Sample Mean', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_title('Distribution of Sample Means', fontsize=14, fontweight='bold')
    ax2.set_xlim(hist_min, hist_max)

    # Set consistent y-axis max
    max_count = max(np.histogram(all_means, bins=bins)[0]) * 1.25
    ax2.set_ylim(0, max_count)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # ── Counter badge ────────────────────────────────────────────────
    ax2.text(0.02, 0.95, f'Samples drawn: {n_drawn}/{N_SAMPLES}\nn = {SAMPLE_SIZE} per sample',
             transform=ax2.transAxes, fontsize=11, fontweight='bold',
             verticalalignment='top',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#ECF0F1',
                       edgecolor='#95A5A6', alpha=0.95))

    plt.tight_layout(pad=1.5)

    # Convert to PIL Image
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=DPI, facecolor='white',
                bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).copy()


print("Generating frames...")
for idx, n in enumerate(frame_indices):
    frames.append(make_frame(n))
    if (idx + 1) % 10 == 0:
        print(f"  Frame {idx + 1}/{len(frame_indices)} (sample {n})")

# Hold the final frame
for _ in range(HOLD_FRAMES_END):
    frames.append(frames[-1].copy())

print(f"Total frames: {len(frames)}")

# ── Save GIF ─────────────────────────────────────────────────────────
output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          '..', 'Lecture_Folder', 'images')
output_path = os.path.join(output_dir, 'week_02_clt_sampling_animation.gif')

frames[0].save(
    output_path,
    save_all=True,
    append_images=frames[1:],
    duration=int(1000 / FPS),
    loop=0,
    optimize=True,
)

print(f"Saved: {output_path}")
print(f"File size: {os.path.getsize(output_path) / 1024 / 1024:.1f} MB")
