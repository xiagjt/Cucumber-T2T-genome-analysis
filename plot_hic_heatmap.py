#!/usr/bin/env python3
"""
Usage:
    python plot_hic_heatmap.py -i input.cool -o output.png
    python plot_hic_heatmap.py -i input.cool -o output.png --chroms Chr1,Chr2,Chr3
    python plot_hic_heatmap.py -i input.cool -o output.png --mode oe --cmap coolwarm
"""

import argparse
import cooler
from cooler import balance_cooler
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot Hi-C/Pore-C contact heatmaps per chromosome')
    parser.add_argument('-i', '--input', required=True,
                        help='Input cool file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output image file (png/pdf/svg)')
    parser.add_argument('--chroms', default=None,
                        help='Chromosomes to plot, comma separated (default: all Chr* chroms)')
    parser.add_argument('--prefix', default='Chr',
                        help='Chromosome name prefix (default: Chr)')
    parser.add_argument('--ncols', type=int, default=4,
                        help='Number of columns in figure (default: 4)')
    parser.add_argument('--mode', default='raw', choices=['raw', 'oe'],
                        help='Matrix mode: raw=log1p(balanced), oe=log2(O/E) (default: raw)')
    parser.add_argument('--cmap', default=None,
                        help='Colormap (default: Reds for raw, coolwarm for oe)')
    parser.add_argument('--ignore-diags', type=int, default=3,
                        help='Number of diagonals to mask near main diagonal (default: 3)')
    parser.add_argument('--percentile', type=float, default=99.5,
                        help='Percentile for vmax (default: 99.5)')
    parser.add_argument('--panel-size', type=float, default=4.0,
                        help='Size of each panel in inches (default: 4.0)')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Output DPI (default: 300)')
    parser.add_argument('--title', default=None,
                        help='Figure title (default: input filename)')
    parser.add_argument('--nticks', type=int, default=5,
                        help='Number of axis ticks per panel (default: 5)')
    return parser.parse_args()

def ensure_balanced(c):
    if 'weight' in c.bins().columns:
        print("  Balance weights found, skipping balancing.")
        return
    print("  No balance weights found, running balance_cooler ...")
    balance_cooler(c, store=True, ignore_diags=2, min_nnz=10, min_count=0)
    print("  Balancing done.")

def mask_diagonals(mat, ignore_diags):
    if ignore_diags == 0:
        return mat
    mat = mat.copy()
    n = mat.shape[0]
    for d in range(ignore_diags + 1):
        idx = np.arange(n - d)
        mat[idx, idx + d] = np.nan
        mat[idx + d, idx] = np.nan
    return mat

def compute_oe(mat):
    n = mat.shape[0]
    expected = np.zeros_like(mat)
    for d in range(n):
        diag = np.diag(mat, d)
        mean_val = np.nanmean(diag)
        if mean_val > 0:
            idx = np.arange(n - d)
            expected[idx, idx + d] = mean_val
            expected[idx + d, idx] = mean_val
    with np.errstate(divide='ignore', invalid='ignore'):
        oe = np.where(expected > 0, np.log2(mat / expected), np.nan)
    return oe

def plot_heatmap(args):
    c = cooler.Cooler(args.input)
    res = c.info['bin-size']

    ensure_balanced(c)

    if args.chroms:
        chroms = args.chroms.split(',')
    else:
        chroms = [ch for ch in c.chromnames if ch.startswith(args.prefix)]
        chroms = sorted(chroms, key=lambda x: int(x.replace(args.prefix, '')))

    # colormap 选择
    if args.cmap:
        cmap = plt.get_cmap(args.cmap)
    elif args.mode == 'oe':
        cmap = plt.get_cmap('coolwarm')
    else:
        cmap = plt.get_cmap('Reds')

    print(f"Input:        {args.input}")
    print(f"Output:       {args.output}")
    print(f"Mode:         {args.mode}")
    print(f"Colormap:     {cmap.name}")
    print(f"Ignore diags: {args.ignore_diags}")
    print(f"Resolution:   {res/1000:.0f}kb")
    print(f"Chromosomes ({len(chroms)}): {', '.join(chroms)}")

    ncols = min(args.ncols, len(chroms))
    nrows = -(-len(chroms) // ncols)

    panel = args.panel_size
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * panel, nrows * (panel + 0.5)))
    if len(chroms) == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes.reshape(1, -1)
    axes_flat = axes.flatten()

    for i, chrom in enumerate(chroms):
        ax = axes_flat[i]
        try:
            mat = c.matrix(balance=True).fetch(chrom)
        except Exception as e:
            print(f"  Warning: skipping {chrom}: {e}")
            ax.set_visible(False)
            continue

        if args.mode == 'oe':
            mat = compute_oe(mat)

        mat_plot = mask_diagonals(mat, args.ignore_diags)

        if args.mode == 'oe':
            valid = mat_plot[~np.isnan(mat_plot)]
            vmax = np.percentile(np.abs(valid), args.percentile) if len(valid) > 0 else 1.0
            vmin = -vmax
            label = 'log2(O/E)'
        else:
            mat_plot = np.log1p(np.nan_to_num(mat_plot, nan=0))
            valid = mat_plot[mat_plot > 0]
            vmax = np.percentile(valid, args.percentile) if len(valid) > 0 else 1.0
            vmin = 0
            label = 'log1p(contacts)'

        n = mat.shape[0]
        im = ax.imshow(mat_plot, cmap=cmap, aspect='equal',
                       vmin=vmin, vmax=vmax, interpolation='none')

        cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        cbar.set_label(label, fontsize=7)
        cbar.ax.tick_params(labelsize=6)

        tick_bins = np.linspace(0, n-1, args.nticks).astype(int)
        tick_labels = [f"{b*res/1e6:.0f}Mb" for b in tick_bins]
        ax.set_xticks(tick_bins)
        ax.set_xticklabels(tick_labels, rotation=45, fontsize=7, ha='right')
        ax.set_yticks(tick_bins)
        ax.set_yticklabels(tick_labels, fontsize=7)

        chrom_len = n * res / 1e6
        ax.set_title(f"{chrom}  ({chrom_len:.1f} Mb)",
                     fontsize=11, fontweight='bold', pad=5)
        print(f"  {chrom}: {n} bins, vmax={vmax:.4f}")

    for j in range(len(chroms), len(axes_flat)):
        axes_flat[j].set_visible(False)

    title = args.title if args.title else args.input
    fig.suptitle(title, fontsize=13, fontweight='bold', y=1.01)

    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
    print(f"\nSaved: {args.output}")

if __name__ == '__main__':
    args = parse_args()
    plot_heatmap(args)
