#!/usr/bin/env python3
"""
Build a cohort ichorCNA report — Markdown + PDF — from a directory of
per-sample results.

Expected layout under --root:
    <root>/
        cohort_summary.csv        (from plot_cohort_summary.py)
        cohort_summary.png        (from plot_cohort_summary.py)
        <sample>/
            <sample>.seg.txt
            <sample>.params.txt
            <sample>.cnv_summary.png   (from plot_cnv_summary.py)

Outputs:
    <root>/COHORT_SUMMARY.md   self-contained markdown (base64-embedded PNGs)
    <root>/COHORT_SUMMARY.pdf  A4 landscape, GitHub-styled

Run plot_cnv_summary.py and plot_cohort_summary.py first (or pass --plot to do
both inline before building the report).

Usage:
    build_cohort_report.py --root /data/Runs/ichorCNA_omics
    build_cohort_report.py --root /data/Runs/ichorCNA_omics --plot
    build_cohort_report.py --root /data/Runs/ichorCNA_omics --skip-pdf
"""
import argparse
import base64
import csv
import re
import subprocess
import sys
import textwrap
from pathlib import Path

NOISE_FLOOR = 0.03  # ichorCNA's documented noise floor for tumor fraction

GITHUB_CSS = """\
@page { size: A4 landscape; margin: 12mm; }
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; max-width: 100%; margin: 0; padding: 0; color: #24292e; line-height: 1.5; }
h1 { font-size: 1.8em; border-bottom: 1px solid #eaecef; padding-bottom: 0.3em; }
h2 { font-size: 1.4em; border-bottom: 1px solid #eaecef; padding-bottom: 0.3em; margin-top: 1.5em; page-break-after: avoid; }
h3 { font-size: 1.15em; margin-top: 1em; page-break-after: avoid; page-break-before: auto; }
table { border-collapse: collapse; margin: 1em 0; font-size: 90%; page-break-inside: avoid; }
th, td { border: 1px solid #dfe2e5; padding: 5px 10px; }
th { background: #f6f8fa; font-weight: 600; }
tr:nth-child(2n) { background: #f6f8fa; }
img { border: 1px solid #eaecef; border-radius: 4px; margin: 0.5em 0; page-break-inside: avoid; display: block; }
code { background: #f6f8fa; padding: 0.2em 0.4em; border-radius: 3px; font-size: 85%; font-family: SFMono-Regular, Consolas, monospace; }
pre { background: #f6f8fa; padding: 12px; border-radius: 6px; overflow: auto; font-size: 85%; }
pre code { background: none; padding: 0; }
blockquote { border-left: 4px solid #dfe2e5; padding: 0 1em; color: #6a737d; margin: 0; }
ul, ol { padding-left: 2em; }
hr { border: none; border-top: 1px solid #eaecef; }
"""


def b64_data_uri(p: Path) -> str:
    return f"data:image/png;base64,{base64.b64encode(p.read_bytes()).decode('ascii')}"


def run_plotting(root: Path, scripts_dir: Path) -> None:
    """Generate per-sample + cohort figures from raw seg/params files."""
    plot_cnv = scripts_dir / "plot_cnv_summary.py"
    plot_cohort = scripts_dir / "plot_cohort_summary.py"
    if not plot_cnv.exists() or not plot_cohort.exists():
        sys.exit(f"--plot requires {plot_cnv.name} and {plot_cohort.name} "
                 f"in {scripts_dir}")

    for d in sorted(root.iterdir()):
        if not d.is_dir():
            continue
        seg = d / f"{d.name}.seg.txt"
        params = d / f"{d.name}.params.txt"
        out = d / f"{d.name}.cnv_summary.png"
        if seg.exists() and params.exists():
            subprocess.run(
                ["python3", str(plot_cnv), "--seg", str(seg),
                 "--params", str(params), "--out", str(out)],
                check=True,
            )

    subprocess.run(
        ["python3", str(plot_cohort),
         "--root", str(root),
         "--out", str(root / "cohort_summary.png"),
         "--csv", str(root / "cohort_summary.csv")],
        check=True,
    )


def build_markdown(root: Path, title: str, ga_run_id: str | None) -> Path:
    """Assemble COHORT_SUMMARY.md with base64-embedded figures."""
    csv_path = root / "cohort_summary.csv"
    cohort_png = root / "cohort_summary.png"
    out = root / "COHORT_SUMMARY.md"

    if not csv_path.exists() or not cohort_png.exists():
        sys.exit(f"missing {csv_path.name} or {cohort_png.name} in {root} "
                 f"(run plot_cohort_summary.py first, or pass --plot)")

    rows = list(csv.DictReader(csv_path.open()))
    rows.sort(key=lambda r: float(r["tumor_fraction"]), reverse=True)

    md: list[str] = []
    md.append(f"# {title}")
    md.append("")
    md.append(f"- **Cohort:** {len(rows)} cfDNA samples")
    md.append("- **Tool:** ichorCNA tumor-fraction / CNA caller, hg38 / 500 kb bins")
    md.append("- **PoN:** bundled hg38 HD-ULP `HD_ULP_PoN_hg38_500kb_median_normAutosome_median.rds`")
    if ga_run_id:
        md.append(f"- **Pipeline:** `IchorCNA` WDL on AWS HealthOmics, "
                  f"GitHub Actions run "
                  f"[{ga_run_id}](https://github.com/Ultimagen/terra_pipeline/actions/runs/{ga_run_id})")
    md.append("")
    md.append("---")
    md.append("")

    # Cohort overview figure
    md.append("## Cohort overview")
    md.append("")
    md.append(f"![cohort summary]({b64_data_uri(cohort_png)})")
    md.append("")

    # Per-sample summary table
    md.append("## Per-sample summary table")
    md.append("")
    md.append("Sorted by tumor fraction descending. **Bold TF** = above ichorCNA's "
              f"~{int(NOISE_FLOOR * 100)}% noise floor. Triploid-anchor flag = ploidy ≈ 3 with "
              "zero NEUT bins (all genome called GAIN/AMP/HLAMP) — strong indicator the "
              "HMM picked the wrong ploidy.")
    md.append("")
    md.append("| Sample (short) | Tumor Fraction | Ploidy | NEUT (Mb) | "
              "GAIN (Mb) | AMP (Mb) | HLAMP (Mb) | HETD (Mb) | Triploid-anchor flag |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        short = r["sample"].split("-")[1]
        tf = float(r["tumor_fraction"])
        pl = float(r["ploidy"])
        neut_mb = float(r["mb_NEUT"])
        tf_str = f"**{tf:.4f}**" if tf >= NOISE_FLOOR else f"{tf:.4f}"
        flag = "⚠️ yes" if (pl > 2.5 and neut_mb == 0) else ""
        md.append(
            f"| {short} | {tf_str} | {pl:.2f} | "
            f"{neut_mb:,.0f} | {float(r['mb_GAIN']):,.0f} | "
            f"{float(r['mb_AMP']):,.0f} | {float(r['mb_HLAMP']):,.0f} | "
            f"{float(r['mb_HETD']):,.0f} | {flag} |"
        )
    md.append("")

    # Auto-derive a short interpretation block from the data
    n_above_floor = sum(1 for r in rows if float(r["tumor_fraction"]) >= NOISE_FLOOR)
    n_zero = sum(1 for r in rows if float(r["tumor_fraction"]) == 0.0)
    n_marginal = n_above_floor - sum(
        1 for r in rows if float(r["tumor_fraction"]) >= 0.10
    )
    n_high = sum(1 for r in rows if float(r["tumor_fraction"]) >= 0.10)
    n_triploid_anchor = sum(
        1 for r in rows
        if float(r["ploidy"]) > 2.5 and float(r["mb_NEUT"]) == 0
    )

    md.append("## Interpretation")
    md.append("")
    md.append(f"- {n_high} sample(s) above 10% TF (confident positives).")
    md.append(f"- {n_marginal} sample(s) in the marginal {int(NOISE_FLOOR*100)}–10% TF band — calls unreliable.")
    md.append(f"- {n_zero} sample(s) returned TF = 0 — consistent with normal / non-tumor profiles.")
    if n_triploid_anchor:
        md.append(f"- ⚠️ {n_triploid_anchor} sample(s) hit the **triploid-anchor failure mode** "
                  "(ploidy ≈ 3 with zero NEUT bins). Consider re-running with "
                  "`--ploidy \"c(2)\"` to compare against the diploid solution.")
    md.append("")
    md.append("---")
    md.append("")

    # Per-sample figures (TF-descending)
    md.append("## Per-sample figures")
    md.append("")
    md.append("Each panel: left = segment count by call category; right = total "
              "genome covered. Title shows estimated tumor fraction. Categories "
              "ordered HOMD → HLAMP (loss → gain).")
    md.append("")
    for r in rows:
        sample = r["sample"]
        short = sample.split("-")[1]
        tf = float(r["tumor_fraction"])
        pl = float(r["ploidy"])
        fig = root / sample / f"{sample}.cnv_summary.png"
        if not fig.exists():
            print(f"WARN: missing per-sample figure for {sample}", file=sys.stderr)
            continue
        md.append(f"### {short} — TF = {tf:.4f}, ploidy = {pl:.2f}")
        md.append("")
        md.append(f"Full sample id: `{sample}`")
        md.append("")
        md.append(f"![{short}]({b64_data_uri(fig)})")
        md.append("")

    out.write_text("\n".join(md))
    print(f"wrote {out} ({out.stat().st_size/1e6:.2f} MB)", file=sys.stderr)
    return out


def build_pdf(md_path: Path, title: str) -> Path:
    """Convert the markdown to a GitHub-styled landscape A4 PDF.

    Two-stage: pandoc gfm → HTML, then wkhtmltopdf HTML → PDF. The base64
    data URIs in the .md are extracted to temp PNG files first because
    wkhtmltopdf doesn't handle large inline images reliably.
    """
    root = md_path.parent
    pdf_path = root / "COHORT_SUMMARY.pdf"
    tmpdir = root / "_pdf_imgs"
    tmpdir.mkdir(exist_ok=True)

    # 1) data: URIs → file paths so wkhtmltopdf can handle them
    src = md_path.read_text()

    def back_to_file(m: re.Match) -> str:
        alt = m.group(1)
        raw = base64.b64decode(m.group(2))
        h = abs(hash(raw))
        p = tmpdir / f"img_{h}.png"
        if not p.exists():
            p.write_bytes(raw)
        return f"![{alt}]({p})"

    md_for_pdf = re.sub(r"!\[([^\]]*)\]\(data:image/png;base64,([^)]+)\)",
                        back_to_file, src)
    md_pdf_path = root / "COHORT_SUMMARY_pdf.md"
    md_pdf_path.write_text(md_for_pdf)

    # 2) md → HTML via pandoc
    html_path = root / "COHORT_SUMMARY_pdf.html"
    subprocess.run(
        ["pandoc", str(md_pdf_path),
         "--from", "gfm", "--to", "html5", "--standalone",
         "--metadata", f"title={title}",
         "--css=-",
         "-o", str(html_path)],
        input=GITHUB_CSS, text=True, check=True,
    )

    # 3) Force img elements to absolute file:// + explicit width so wkhtmltopdf
    #    doesn't render at intrinsic pixel size and overflow the page
    html = html_path.read_text()

    def fix_img(m: re.Match) -> str:
        src = m.group(1)
        if src.startswith(("http", "data:", "/", "file:")):
            new_src = src
        else:
            new_src = f"file://{Path(src).resolve()}"
        return (f'<img src="{new_src}" '
                'style="width:260mm;max-width:100%;height:auto;'
                'display:block;margin:0.5em auto;'
                'border:1px solid #eaecef;border-radius:4px;">')

    html = re.sub(r'<img[^>]*src="([^"]+)"[^>]*>', fix_img, html)
    html_path.write_text(html)

    # 4) HTML → PDF via wkhtmltopdf. The "ContentNotFoundError" warning is a
    #    known false alarm from wkhtmltopdf and does not indicate failure.
    subprocess.run(
        ["wkhtmltopdf", "--enable-local-file-access", "--quiet",
         "--orientation", "Landscape", "--page-size", "A4",
         "--margin-top", "12mm", "--margin-bottom", "12mm",
         "--margin-left", "12mm", "--margin-right", "12mm",
         str(html_path), str(pdf_path)],
        check=False,  # wkhtmltopdf often returns 1 on harmless warnings
    )

    if not pdf_path.exists() or pdf_path.stat().st_size == 0:
        sys.exit(f"PDF generation failed; check {html_path}")

    # Clean up intermediates
    md_pdf_path.unlink(missing_ok=True)
    html_path.unlink(missing_ok=True)
    for img in tmpdir.glob("*.png"):
        img.unlink()
    tmpdir.rmdir()

    print(f"wrote {pdf_path} ({pdf_path.stat().st_size/1e6:.2f} MB)",
          file=sys.stderr)
    return pdf_path


def main():
    p = argparse.ArgumentParser(
        description=textwrap.dedent(__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--root", required=True, type=Path,
                   help="Cohort root directory (one subdir per sample)")
    p.add_argument("--title", default="ichorCNA cohort summary",
                   help="Top-level title for the report")
    p.add_argument("--ga-run-id", default=None,
                   help="Optional GitHub Actions run id (linked in the report header)")
    p.add_argument("--plot", action="store_true",
                   help="Run plot_cnv_summary.py + plot_cohort_summary.py first "
                        "(otherwise expects pre-generated figures)")
    p.add_argument("--skip-pdf", action="store_true",
                   help="Only build the markdown, skip wkhtmltopdf step")
    p.add_argument("--scripts-dir", type=Path,
                   default=Path(__file__).parent,
                   help="Directory containing plot_cnv_summary.py + "
                        "plot_cohort_summary.py (defaults to this script's directory)")
    args = p.parse_args()

    if not args.root.is_dir():
        sys.exit(f"--root not a directory: {args.root}")

    if args.plot:
        run_plotting(args.root, args.scripts_dir)

    md_path = build_markdown(args.root, args.title, args.ga_run_id)

    if not args.skip_pdf:
        build_pdf(md_path, args.title)


if __name__ == "__main__":
    main()
