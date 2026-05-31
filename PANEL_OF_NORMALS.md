# Panel of Normals (PoN) — guide

How the bundled ichorCNA PoN is structured, how the tool consumes it at runtime, and what to consider before replacing it with your own.

---

## File location

```
inst/extdata/HD_ULP_PoN_hg38_500kb_median_normAutosome_median.rds   661 KB
inst/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds     320 KB
```

`run_pipeline.sh` auto-selects one of these based on `bin_size` (500 kb → 500 kb PoN, 1 Mb → 1 Mb PoN, anything else → no PoN). The same auto-selection is wired into the WDL.

There are also two older hg19 / non-prefixed files in the same directory (`HD_ULP_PoN_500kb_…rds`, `HD_ULP_PoN_1Mb_…rds`) that the pipeline does not use.

The naming convention encodes how the file was built:
- `HD_ULP_` — Healthy Donor, Ultra-Low-Pass WGS
- `hg38_500kb_` — genome build and bin size
- `_median_` — bins aggregated across donors using a per-bin median (robust to outlier donors)
- `_normAutosome_` — each donor was normalized so its autosome (chr1–22) log-ratio mean is zero
- `_median.rds` — final per-bin statistic written into the `Median` column

---

## What the file contains

The `.rds` is an R serialized **`GRanges`** object.

| Property | Value |
|---|---|
| Number of normal samples | **27 healthy donors** (`HD1`, `HD2`, …, `HD30` with gaps at HD14, HD15, HD25, HD26) |
| Sample type | **Plasma cfDNA** (the column names `HDn.ctDNA.QC.FC<flowcell>.reads` make this explicit) |
| Coverage | Ultra-low-pass WGS (~0.1×) — the published Adalsteinsson et al. 2017 reference cohort |
| Bin size | exactly 500,000 bp (every bin) |
| Chromosomes | chr1–chr22 + chrX (no chrY, no decoys, no alt contigs) |
| Total bins | ~6,165 across the genome |
| Genome build | hg38 / UCSC ("chr" prefix) |
| Variant information | **None** — this is a depth-only / read-count-derived PoN |

### Per-bin metadata columns (`mcols`)

| Column | What it is |
|---|---|
| `reads` | Average read count for the bin (legacy/unused by ichorCNA's normalization path) |
| `gc` | GC fraction of the bin |
| `map` | Mappability score of the bin |
| `valid` | Boolean — bin passes basic quality filters |
| `ideal` | Boolean — bin passes stricter filters (high map, normal GC, etc.) |
| `cor.gc` | GC-corrected log-ratio for the median signal |
| `cor.map` | Map- and GC-corrected log-ratio |
| `HD1.ctDNA.QC.FC...reads` … `HD30.ctDNA.QC.FC...reads` | **27 columns**, one per donor — each holds the donor's per-bin GC + map-corrected log-ratio (already normalized to autosome median) |
| `Median` | The per-bin median across all 27 donors. **This is the "expected normal" track that ichorCNA subtracts from your tumor sample.** |

So **the PoN is depth-only**: per-bin log-ratios derived from read counts after GC and mappability correction. There are no SNP positions, no allele frequencies, and no variant calls. Allele balance / loss-of-heterozygosity (LOH) detection requires a BAF-aware tool such as TitanCNA.

### Whole genome or specific intervals?

**Whole genome, tiled at 500 kb fixed bins** (chr1–22 + chrX, no Y). Centromeres and telomeres are *not* explicitly excluded in the file — bins simply have lower `ideal=TRUE` rates there. There is no exome / target-panel layer.

---

## How ichorCNA uses the PoN at runtime

When you pass `--normalPanel /path/to/pon.rds` to `runIchorCNA.R`:

1. **Bin matching.** ichorCNA loads the PoN's `GRanges` and matches bins to the GC/map-corrected log-ratios it just computed for your sample. Bin coordinates and bin size must match exactly — that is why the PoN must be built at the same bin size as the GC/map WIGs and your sample's WIG.
2. **Median subtraction.** For each bin, ichorCNA subtracts the PoN's `Median` column from the sample's GC/map-corrected log-ratio. This removes systematic biases that recur across normal libraries — sequencing chemistry artifacts, repeated GC/map quirks the GC correction did not catch, fragment-size biases in cfDNA, etc.
3. **Bin filtering.** Bins where `ideal == FALSE` in the PoN are dropped from downstream HMM segmentation.
4. **Segmentation / tumor-fraction estimation.** The PoN-corrected log-ratios are passed to ichorCNA's HMM, which produces segment calls (HOMD/HETD/NEUT/GAIN/AMP/HLAMP), tumor-fraction estimate, and ploidy.

Without a PoN (`--normalPanel NULL`), step 2 is skipped and the sample is normalized using only its own GC/map correction. This is fine for very high-depth samples where library-specific biases are small relative to the signal, but at ULP-WGS depths (~0.1–1×) the PoN typically adds substantial signal-to-noise improvement.

---

## How to build your own PoN

The repo ships `scripts/createPanelOfNormals.R`, which takes a list of WIG files (one per normal sample) and produces an `.rds` file in the same format as the bundled PoN.

```bash
# 1. Stream each normal CRAM from S3 and produce a fixedStep WIG at the
#    desired bin size. These should be cfDNA libraries from healthy
#    donors processed identically to your tumor samples.
ls /path/to/normal_wigs/*.wig > normals.filelist

# 2. Build the PoN
Rscript /home/ubuntu/ichorCNA/scripts/createPanelOfNormals.R \
    --gcWig  inst/extdata/gc_hg38_500kb.wig \
    --mapWig inst/extdata/map_hg38_500kb.wig \
    --centromere inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
    --filelist normals.filelist \
    --chrs 'c(1:22,"X")' \
    --genomeStyle UCSC \
    --chrNormalize 'c(1:22)' \
    --method median \
    --outfile /path/to/my_PoN_hg38_500kb.rds
```

The resulting `.rds` is a drop-in replacement — pass it via `--normalPanel` to `runIchorCNA.R` (or `IchorCNA.normal_panel_override` in the WDL).

> **WIG format requirement.** The WIGs you feed to `createPanelOfNormals.R` must be **fixedStep with zero-fill** for empty bins, exactly as `scripts/bam_to_wig.py --fai <ref.fai>` produces them. ichorCNA's `wigToGRanges` parser crashes on `variableStep` WIGs (`Error in (breaks[i] + 1):(breaks[i + 1] - 1) : NA/NaN argument`) — using `--fai` is what enforces the correct format.

---

## Recommendations before replacing the PoN

1. **Match the protocol of your tumor samples.** The PoN's job is to remove *systematic* biases — sequencing chemistry, library prep, capture vs WGS, fragment size, GC profile from the sequencer's PCR chemistry. If your normals are on a different platform, you will *inject* new biases instead of removing them. Ultima sequencing's GC/coverage profile is different enough from Illumina that the bundled HD-ULP PoN may itself be a poor match for Ultima cfDNA samples — this could partially explain the triploid-anchor failures we saw in earlier runs.

2. **Match the sample type.** Plasma cfDNA tumor samples need plasma cfDNA normals — *not* buffy coat, not solid-tissue germline, not cell-line DNA. cfDNA has a characteristic ~167 bp fragment-size profile that produces distinctive GC/coverage patterns; tissue-mismatched normals will not cancel them out.

3. **Cohort size.** The published HD cohort has 27 samples. Aim for **≥10 high-quality matched normals**; the ichorCNA authors recommend ~15–40. Below ~10 the per-bin median becomes noisy and the PoN can introduce artifacts. Returns diminish above ~30.

4. **Use truly normal donors.** Avoid germline samples with known large-scale CNVs (Down syndrome, sex-chromosome aneuploidies). Do not include tumor-adjacent normals — copy-number aberrations in even one or two donors can leak into the median when the cohort is small.

5. **QC each candidate normal first.** Run each candidate through ichorCNA *with no PoN* and inspect the genome-wide plot. Drop any that show segments — they have CNVs that would contaminate the median.

6. **Match coverage to runtime.** The published PoN is ULP-WGS (~0.1×). For the 1× / 10× downsampled tumor samples we have been processing, build your PoN from normals **downsampled to the same target coverage**. Variance scales with depth, so coverage mismatch produces artifacts.

7. **Build separate PoNs per bin size.** The bundled hg38 PoN exists only at 500 kb and 1 Mb. For 50 kb / 10 kb bins you must build your own at that bin size — bins must match exactly.

8. **chrY is excluded.** If you care about chrY, build a male-only sub-PoN. Sex-chromosome handling generally is mixed-sex in the bundled PoN; it works for autosomes + chrX log-ratio interpretation only.

9. **Exclude problem regions.** The PoN inherits whatever bins the GC/map WIGs include. Bins with `ideal == FALSE` are filtered downstream, but if you have a known problem-region BED for your assay, pass it via `--exons.bed` so those bins are excluded up front.

---

## Quick sanity checks for an existing PoN

```r
# In R
pon <- readRDS("HD_ULP_PoN_hg38_500kb_median_normAutosome_median.rds")
class(pon)             # GRanges
length(pon)            # ~6165
table(seqnames(pon))   # chr1..chr22, chrX
table(width(pon))      # all 500000
colnames(mcols(pon))   # reads, gc, map, valid, ideal, cor.gc, cor.map,
                       # HD1.ctDNA.QC.FC..., ..., Median
sum(mcols(pon)$ideal)  # number of usable bins (~5000ish)
summary(mcols(pon)$Median)  # should be tightly centered on 0
```

If `summary(Median)` shows large deviations from zero or long tails, one or more donor samples are pulling the median — investigate before using the PoN.
