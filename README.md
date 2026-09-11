# qPCR ddCt Pipeline

Automated qPCR ddCt analysis from QuantStudio exports: SD-based technical-replicate
QC (optionally skippable), exploratory imputation of failed groups, ddCt/RQ/log2FC
calculation, statistical testing, and plotting.

## Workflow

1. **Read plates.** Reads one or more QuantStudio qPCR export `.xlsx` files (one file
   per plate). Each file must contain a `Results` sheet with the standard per-well
   columns (`Well`, `Well Position`, `Sample`, `Target`, `Cq`, ...).

2. **QC / outlier removal.** Groups technical replicates by `(plate, sample, target)`
   and applies an SD-based outlier removal rule:
   - compute SD of the replicate group
   - if SD <= threshold -> pass, keep all replicates
   - if SD > threshold -> remove the replicate farthest from the group mean,
     recompute SD, repeat until below threshold
   - if only 2 replicates remain and SD is still > threshold -> the whole
     `(sample, target)` group **FAILS QC** (but can still be salvaged later, see
     step 3)

   The SD threshold is set via `--sd-threshold` (default `0.3`).

   **To skip this step entirely**, pass `--skip-outlier-removal`. Every valid
   (non-NaN) replicate is then kept and simply averaged, `--sd-threshold` is
   ignored, and a group only fails if *none* of its replicates amplified at all
   (all `Cq` are `NaN`). Everything downstream (ddCt/RQ/log2FC, statistics, plots)
   runs exactly as it would otherwise — this just changes how the per-group Cq
   mean/SD is computed.

   Biological condition + replicate id are resolved per `Sample`, in priority order:
   1. an explicit `--sample-map` CSV
   2. DA3's own `Biogroup` field (`Replicate Group Result` sheet), if it was
      actually populated during plate setup
   3. a regex on the `Sample` text (`--condition-regex`, default: a trailing
      `" <number>"`)

   Sample names do **not** need to encode the replicate number themselves as long
   as DA3's Biogroup field was set, or you supply `--sample-map`.

3. **Exploratory imputation (optional).** With `--exploratory`, failed
   `(sample, target)` groups can be salvaged instead of dropped, using one of
   several strategies (`--fail-strategy`): `plate_mean`, `anchor_mean`,
   `best_subset`, or `drop`. Most useful for housekeeping genes, where losing a
   sample entirely because of one bad technical replicate is wasteful.

4. **ddCt method.** Computed per plate, using a user-specified anchor `SAMPLE`
   (one specific well/biological replicate used as the calibrator, e.g.
   `"Empty 1"`) and housekeeping gene (default `GAPDH`):

   ```
   dCt    = Cq(target) - Cq(housekeeping)
   ddCt   = dCt(sample) - dCt(anchor)
   RQ     = 2^-ddCt
   log2FC = -ddCt   (== log2(RQ))
   ```

   Omit `--housekeeping` entirely for experiments with no relative-quantification
   reference gene (e.g. a cloning/junction-validation qPCR) — the pipeline then
   runs QC + raw Cq plots only and skips ddCt/RQ/log2FC and statistics, rather than
   forcing a fake normalization that would just produce all-NaN results.

5. **Statistics.** Group-wise tests on log2FC per target vs. a control condition
   (all biological replicates of the anchor's condition, not just the single
   anchor well): per-group Shapiro-Wilk normality + Levene's test for equal
   variances decide parametric vs. non-parametric, then:
   - **parametric:** one-way ANOVA + Dunnett's post-hoc vs. control
   - **non-parametric:** Kruskal-Wallis + Dunn's post-hoc vs. control

6. **Outputs.** Writes result tables (technical-replicate level and
   sample/target-level with Cq, dCt, ddCt, RQ, log2FC) plus a QC log, and
   generates plots:
   - `plots/qc/` — per-plate Cq-per-technical-replicate plots with removed
     points marked (X), faceted by target
   - `plots/cq/` — Cq per sample per target, raw and anchor-referenced
   - `plots/log2fc/` — log2FC per condition per target with significance
     brackets vs. the control condition

## Usage

Standard run, with SD-based outlier removal and exploratory imputation for
failed groups:

```bash
python qpcr_pipeline.py \
    --input plate1.xlsx plate2.xlsx plate3.xlsx \
    --plate-labels Plate1 Plate2 Plate3 \
    --anchor "Empty 1" "E 1" "E 1" \
    --housekeeping GAPDH \
    --sd-threshold 0.3 \
    --exploratory --fail-strategy plate_mean \
    --outdir ./qpcr_results
```

Same run, but skipping SD-based outlier removal entirely (average every valid
replicate as-is; `--sd-threshold` is ignored) and going straight to ddCt:

```bash
python qpcr_pipeline.py \
    --input plate1.xlsx plate2.xlsx plate3.xlsx \
    --plate-labels Plate1 Plate2 Plate3 \
    --anchor "Empty 1" "E 1" "E 1" \
    --housekeeping GAPDH \
    --skip-outlier-removal \
    --outdir ./qpcr_results_no_qc
```

Run `python qpcr_pipeline.py --help` for the full argument list.

## CLI arguments

| Argument | Default | Description |
|---|---|---|
| `--input` (required) | — | One or more QuantStudio export `.xlsx` files, one per plate. |
| `--plate-labels` | file stem | Custom plate labels, matched by order to `--input`. |
| `--anchor` | — | Anchor sample name (exact `Sample` text, e.g. `"Empty 1"`), one per `--input` file. Required if `--housekeeping` is set; optional otherwise (highlights the sample as a reference line in the "with anchor" Cq plot, no ddCt math). |
| `--control-condition` | resolved from first plate's anchor | Baseline condition group for statistical comparisons. |
| `--condition-regex` | `^(?P<condition>.*?)\s+(?P<rep>\d+)$` | Regex splitting a `Sample` name into condition + replicate id, used when no `--sample-map` entry or DA3 Biogroup is available. |
| `--sample-map` | — | CSV with columns `Sample,Condition,BioRep` (optionally `Plate`) to explicitly declare condition/replicate per sample. Takes priority over DA3's Biogroup and over `--condition-regex`. |
| `--housekeeping` | — | Housekeeping gene/Target for the ddCt method (e.g. `GAPDH`). Omit for experiments with no reference gene — runs QC + raw Cq plots only. |
| `--sd-threshold` | `0.3` | Max allowed technical-replicate Cq SD before outlier removal kicks in. Ignored if `--skip-outlier-removal` is set. |
| `--skip-outlier-removal` | off | Bypass SD-based outlier removal entirely: average every valid replicate as-is; `--sd-threshold` is ignored. A group only fails if none of its replicates amplified. |
| `--exploratory` | off | Enable imputation of FAILED groups instead of dropping them. |
| `--fail-strategy` | `plate_mean` | Imputation strategy when `--exploratory` is set: `drop`, `plate_mean`, `anchor_mean`, or `best_subset`. |
| `--alpha` | `0.05` | Significance threshold for normality/variance/omnibus tests. |
| `--dpi` | `150` | Plot resolution. |
| `--outdir` (required) | — | Output folder (created if missing). |

## Output structure

```
<outdir>/
├── run.log
├── tables/
│   ├── technical_replicates_QC.csv
│   ├── sample_target_Cq_summary.csv
│   ├── ddCt_RQ_log2FC_results.csv        (only if --housekeeping given)
│   ├── statistics_posthoc.csv            (only if --housekeeping given)
│   └── statistics_report.txt             (only if --housekeeping given)
└── plots/
    ├── qc/
    ├── cq/
    └── log2fc/                           (only if --housekeeping given)
```

## Dependencies

```bash
pip install pandas numpy openpyxl scipy matplotlib seaborn scikit-posthocs statsmodels
```

`scipy.stats.dunnett` requires `scipy >= 1.11`.