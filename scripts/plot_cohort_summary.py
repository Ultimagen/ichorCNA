#!/usr/bin/env python3
"""
Cohort-level summary across N ichorCNA samples.

Two stacked panels:
  top    - tumor fraction per sample (sorted descending). Bar color is
           a quick visual cue: gray = below noise floor (~3%), red = above.
  bottom - stacked bar of total Mb covered by each call category, on the
           same x-axis order so you can visually correlate "high TF" vs
           "lots of CNV mass".

Usage:
    plot_cohort_summary.py --root /data/Runs/ichorCNA_omics \\
        --out /data/Runs/ichorCNA_omics/cohort_summary.png

The script auto-discovers samples by looking for <root>/<sample>/<sample>.seg.txt
and the matching .params.txt next to it.
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from plot_cnv_summary import (CALL_COLORS, CALL_ORDER, parse_tumor_fraction,
                              summarize_calls)

# Below this TF, ichorCNA's purity/ploidy estimates are unreliable (its
# documented noise floor is ~3%).
NOISE_FLOOR = 0.03


def collect(root: Path) -> pd.DataFrame:
    rows = []
    for d in sorted(root.iterdir()):
        if not d.is_dir():
            continue
        seg = d / f"{d.name}.seg.txt"
        params = d / f"{d.name}.params.txt"
        if not (seg.exists() and params.exists()):
            continue
        tf, ploidy = parse_tumor_fraction(params)
        summary = summarize_calls(seg)
        row = {"sample": d.name, "tumor_fraction": tf, "ploidy": ploidy}
        for call in CALL_ORDER:
            row[f"mb_{call}"] = summary.loc[call, "total_mb"]
            row[f"n_{call}"] = summary.loc[call, "count"]
        rows.append(row)
    return pd.DataFrame(rows)


def make_figure(df: pd.DataFrame, out: Path) -> None:
    df = df.sort_values("tumor_fraction", ascending=False).reset_index(drop=True)

    # Short labels — full sample IDs are 30+ chars and unreadable on x-axis
    short = df["sample"].str.split("-").str[1]
    fig, (ax_tf, ax_mb) = plt.subplots(
        2, 1, figsize=(max(10, 0.8 * len(df) + 4), 8),
        sharex=True, constrained_layout=True
    )

    # Top: tumor fraction
    bar_colors = [
        CALL_COLORS["AMP"] if tf is not None and tf >= NOISE_FLOOR else "#bdbdbd"
        for tf in df["tumor_fraction"]
    ]
    bars = ax_tf.bar(short, df["tumor_fraction"], color=bar_colors,
                     edgecolor="#333", linewidth=0.5)
    ax_tf.axhline(NOISE_FLOOR, ls="--", lw=1, color="#666",
                  label=f"Noise floor ({NOISE_FLOOR:.0%})")
    ax_tf.set_ylabel("Tumor Fraction")
    ax_tf.set_title("Estimated tumor fraction per sample (sorted)")
    ax_tf.set_axisbelow(True)
    ax_tf.grid(axis="y", alpha=0.3, linewidth=0.6)
    ax_tf.spines[["top", "right"]].set_visible(False)
    ax_tf.legend(loc="upper right", frameon=False, fontsize=9)
    for b, tf in zip(bars, df["tumor_fraction"]):
        if tf is not None:
            ax_tf.text(b.get_x() + b.get_width() / 2,
                       b.get_height() + ax_tf.get_ylim()[1] * 0.01,
                       f"{tf:.3f}", ha="center", va="bottom", fontsize=8)

    # Bottom: stacked Mb by call category
    bottoms = [0.0] * len(df)
    for call in CALL_ORDER:
        vals = df[f"mb_{call}"].astype(float).tolist()
        if sum(vals) == 0:
            continue  # skip categories that no sample has
        ax_mb.bar(short, vals, bottom=bottoms,
                  color=CALL_COLORS[call], label=call,
                  edgecolor="#333", linewidth=0.3)
        bottoms = [b + v for b, v in zip(bottoms, vals)]
    ax_mb.set_ylabel("Total genome covered (Mb)")
    ax_mb.set_title("CNV burden by call category")
    ax_mb.set_axisbelow(True)
    ax_mb.grid(axis="y", alpha=0.3, linewidth=0.6)
    ax_mb.spines[["top", "right"]].set_visible(False)
    ax_mb.legend(loc="upper right", frameon=False, fontsize=9, ncol=2)
    plt.setp(ax_mb.get_xticklabels(), rotation=30, ha="right")
    ax_mb.set_xlabel("Sample (short ID)")

    fig.suptitle(f"ichorCNA cohort summary — N={len(df)} samples",
                 fontsize=13, fontweight="bold")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", required=True, type=Path,
                   help="Directory containing one subdir per sample, each with "
                        "<sample>.seg.txt + <sample>.params.txt")
    p.add_argument("--out", required=True, type=Path,
                   help="Output figure path (.png or .pdf)")
    p.add_argument("--csv", type=Path, default=None,
                   help="Optional path to write the cohort summary table as CSV")
    args = p.parse_args()

    df = collect(args.root)
    if df.empty:
        sys.exit(f"No samples found under {args.root}")
    print(df.to_string(index=False), file=sys.stderr)
    make_figure(df, args.out)
    print(f"wrote {args.out}", file=sys.stderr)
    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"wrote {args.csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
