#!/usr/bin/env python3
"""
Build a CNV summary figure for one ichorCNA sample.

Inputs:  the per-sample <id>.seg.txt and <id>.params.txt produced by
         runIchorCNA.R.
Output:  a PNG (or PDF) with two panels —
           left  : count of segments by call category
           right : total Mb of genome covered by each category
         Title shows the estimated tumor fraction + ploidy.

Why two panels: a single bar with both count and length forces a dual y-axis
which is misleading. Side-by-side keeps each scale honest, and the bar colors
match across panels so the eye can connect them.

Usage:
    plot_cnv_summary.py --seg sample.seg.txt --params sample.params.txt \\
        --out sample.cnv_summary.png

Multi-sample (one figure per sample, useful when scripting a cohort):
    for d in /data/Runs/ichorCNA_*/ichorCNA/*/; do
        s=$(basename "$d")
        plot_cnv_summary.py --seg "$d/$s.seg.txt" \\
            --params "$d/$s.params.txt" \\
            --out "$d/$s.cnv_summary.png"
    done
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Stable order + colors so plots are comparable across samples.
# Category list matches ichorCNA call labels:
#   HOMD = homozygous deletion (CN=0)
#   HETD = heterozygous deletion (CN=1)
#   NEUT = neutral (CN=2)
#   GAIN = single-copy gain (CN=3)
#   AMP  = amplification (CN=4)
#   HLAMP = high-level amplification (CN>=5)
CALL_ORDER = ["HOMD", "HETD", "NEUT", "GAIN", "AMP", "HLAMP"]
CALL_COLORS = {
    "HOMD":  "#08306b",  # deepest blue — strongest loss
    "HETD":  "#4292c6",  # blue — loss
    "NEUT":  "#bdbdbd",  # gray — neutral
    "GAIN":  "#fcae91",  # light red — gain
    "AMP":   "#cb181d",  # red — amp
    "HLAMP": "#67000d",  # deepest red — high-level amp
}


def parse_tumor_fraction(params_path: Path) -> tuple[float | None, float | None]:
    """Pull the headline tumor fraction + ploidy from <id>.params.txt.

    The file's first non-empty line is the per-sample header:
        Sample\tTumor Fraction\tPloidy
        603747-...\t0.05393\t2.99
    """
    tf, ploidy = None, None
    with open(params_path) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    for i, ln in enumerate(lines):
        if ln.startswith("Sample\t") and i + 1 < len(lines):
            row = lines[i + 1].split("\t")
            if len(row) >= 3:
                try:
                    tf = float(row[1])
                    ploidy = float(row[2])
                except ValueError:
                    pass
            break
    return tf, ploidy


def summarize_calls(seg_path: Path, call_col: str = "Corrected_Call") -> pd.DataFrame:
    """Return a DataFrame indexed by CALL_ORDER with columns count, total_mb."""
    seg = pd.read_csv(seg_path, sep="\t")
    if call_col not in seg.columns:
        # Fallback for older ichorCNA outputs that don't have Corrected_Call
        call_col = "call"
    seg["length_mb"] = (seg["end"] - seg["start"]).clip(lower=0) / 1e6

    g = seg.groupby(call_col).agg(count=(call_col, "size"),
                                  total_mb=("length_mb", "sum"))
    # Reindex to canonical order, filling 0 for missing categories so every
    # plot has the same x axis — easier to compare across samples.
    g = g.reindex(CALL_ORDER, fill_value=0)
    return g


def make_figure(summary: pd.DataFrame, sample_id: str, tumor_frac: float | None,
                ploidy: float | None, out_path: Path) -> None:
    fig, (ax_count, ax_len) = plt.subplots(
        1, 2, figsize=(11, 4.2), constrained_layout=True
    )

    colors = [CALL_COLORS[c] for c in summary.index]

    # Left: segment count
    bars1 = ax_count.bar(summary.index, summary["count"], color=colors,
                         edgecolor="#333", linewidth=0.5)
    ax_count.set_ylabel("Number of segments")
    ax_count.set_title("Segment count by call")
    ax_count.set_axisbelow(True)
    ax_count.grid(axis="y", alpha=0.3, linewidth=0.6)
    ax_count.spines[["top", "right"]].set_visible(False)
    for b, v in zip(bars1, summary["count"]):
        if v > 0:
            ax_count.text(b.get_x() + b.get_width() / 2, v, f"{int(v)}",
                          ha="center", va="bottom", fontsize=9)

    # Right: total genome covered
    bars2 = ax_len.bar(summary.index, summary["total_mb"], color=colors,
                       edgecolor="#333", linewidth=0.5)
    ax_len.set_ylabel("Total genome covered (Mb)")
    ax_len.set_title("Genome coverage by call")
    ax_len.set_axisbelow(True)
    ax_len.grid(axis="y", alpha=0.3, linewidth=0.6)
    ax_len.spines[["top", "right"]].set_visible(False)
    for b, v in zip(bars2, summary["total_mb"]):
        if v > 0:
            ax_len.text(b.get_x() + b.get_width() / 2, v, f"{v:,.0f}",
                        ha="center", va="bottom", fontsize=9)

    # Title — center over both panels with key biology in the suptitle
    bits = [sample_id]
    if tumor_frac is not None:
        bits.append(f"Tumor Fraction = {tumor_frac:.4f}")
    fig.suptitle(" | ".join(bits), fontsize=12, fontweight="bold")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seg", required=True, type=Path,
                   help="Path to ichorCNA <id>.seg.txt")
    p.add_argument("--params", required=True, type=Path,
                   help="Path to ichorCNA <id>.params.txt")
    p.add_argument("--out", required=True, type=Path,
                   help="Output figure path (.png or .pdf)")
    p.add_argument("--sample-id", default=None,
                   help="Override sample id shown in the title "
                        "(default: derived from --seg filename)")
    p.add_argument("--call-col", default="Corrected_Call",
                   choices=["Corrected_Call", "call"],
                   help="Which call column to summarize (default: Corrected_Call)")
    args = p.parse_args()

    sample_id = args.sample_id or args.seg.name.replace(".seg.txt", "")
    tumor_frac, ploidy = parse_tumor_fraction(args.params)
    summary = summarize_calls(args.seg, args.call_col)

    print(f"sample: {sample_id}", file=sys.stderr)
    print(f"tumor fraction: {tumor_frac}", file=sys.stderr)
    print(f"ploidy: {ploidy}", file=sys.stderr)
    print(summary.to_string(), file=sys.stderr)

    make_figure(summary, sample_id, tumor_frac, ploidy, args.out)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
