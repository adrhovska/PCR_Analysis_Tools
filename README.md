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

4. **ddCt method.** Computed per plate, referenced to a user-specified anchor
   `SAMPLE` (one specific well/biological replicate used as the calibrator, e.g.
   `"Empty 1"` or `"Church wt"`). **Giving `--anchor` is what triggers this step**
   (and step 5) — omit it entirely to get QC + raw Cq outputs only. Whether
   `--housekeeping` is also given picks one of two modes:

   - **With a housekeeping gene** (e.g. `GAPDH`) — the standard two-step method:

     ```
     dCt    = Cq(target) - Cq(housekeeping)
     ddCt   = dCt(sample) - dCt(anchor)
     ```

   - **Without a housekeeping gene** (`--housekeeping` omitted) — ddCt is taken
     directly off each target's own Cq, with no reference-gene normalization step:

     ```
     dCt    = Cq(target)                 (no housekeeping subtraction)
     ddCt   = dCt(sample) - dCt(anchor)
     ```

     Use this mode when there's no housekeeping/endogenous-control assay on the
     plate by design — e.g. template input was already equalized another way
     (Qubit-based dilution to a fixed ng amount), so a reference-gene
     normalization would just add noise rather than remove it. This is the
     right mode for something like a Golden Gate cutting-efficiency plate with
     an uncut negative control (e.g. `"Church wt"`) as the anchor and no
     housekeeping gene on the plate.

   Either way:

   ```
   RQ     = 2^-ddCt
   log2FC = -ddCt   (== log2(RQ))
   ```

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
     brackets vs. the control condition, both as one file per target
     (`log2FC_<target>.png`) and as a single combined figure with all
     targets side by side (`log2FC_all_targets.png`) for an at-a-glance
     comparison across primer pairs
   - `plots/rq/` — the same comparison, drawn instead as a bar chart on the
     linear RQ scale: mean bar (`2^-mean(ddCt)`) + whiskers (± 1 SEM of
     ddCt) + individual biological-replicate points + a dashed RQ = 1 line,
     with the same significance brackets — one file per target
     (`RQ_<target>.png`) plus a combined figure (`RQ_all_targets.png`).
     Pass `--rq-highlight-target-match` to color a condition's bar
     differently when the target/gene name appears in the condition name
     (the "targets this gene" vs. "does not target this gene" convention
     used for dCas9/CRISPRoff-style gene-silencing experiments — leave it
     off for assays like digestion-efficiency plates where condition names
     don't encode a target gene).

   In `plots/qc/` and `plots/cq/`, samples are ordered left-to-right by their
   **physical position on the plate** (first well occupied, row letter then
   column number) rather than alphabetically or by condition — this matches
   how the plate was actually laid out, which makes visual interpretation
   easier.

   Plot titles are kept short (`"QC — {plate}"`, `"{plate} — Cq"`,
   `"{target} vs {control}"`) so they normally fit at the default figure
   size. If a plate label or target name is unusually long and would still
   get clipped, the figure widens itself just enough to fit the title
   rather than truncating it — you don't need to do anything for this, but
   a shorter `--plate-labels` value keeps the images from growing
   unnecessarily wide.

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

No housekeeping/endogenous-control assay at all (e.g. a Golden Gate
cutting-efficiency plate where template was Qubit-normalized): give `--anchor`
and simply omit `--housekeeping` — ddCt is computed directly off each target's
own Cq vs. the anchor sample:

```bash
python qpcr_pipeline.py \
    --input plate1.xlsx \
    --anchor "Church wt" \
    --outdir ./qpcr_results_no_housekeeping
```

Same idea, but with RQ bar plots colored by whether the condition targets the
plotted gene (matching a dCas9/CRISPRoff-style expression-silencing plot,
e.g. Alina's Day-14 qPCR):

```bash
python qpcr_pipeline.py \
    --input Plate1.xlsx Plate2.xlsx Plate3.xlsx \
    --plate-labels Plate1 Plate2 Plate3 \
    --anchor "Empty 1" "Empty 1" "Empty 1" \
    --control-condition "Empty" \
    --housekeeping GAPDH \
    --sd-threshold 0.3 \
    --exploratory --fail-strategy best_subset \
    --rq-highlight-target-match \
    --outdir ./qpcr_results_rq
```

Run `python qpcr_pipeline.py --help` for the full argument list.

## CLI arguments

| Argument | Default | Description |
|---|---|---|
| `--input` (required) | — | One or more QuantStudio export `.xlsx` files, one per plate. |
| `--plate-labels` | file stem | Custom plate labels, matched by order to `--input`. |
| `--anchor` | — | Anchor/calibrator sample name (exact `Sample` text, e.g. `"Empty 1"` or `"Church wt"`), one per `--input` file. **Giving this is what triggers ddCt/RQ/log2FC + statistics** — with `--housekeeping` also set, the standard two-step method is used; without it, ddCt is computed directly off the target's own Cq vs. the anchor. Omit `--anchor` entirely for QC + raw Cq outputs only. |
| `--control-condition` | resolved from first plate's anchor | Baseline condition group for statistical comparisons. |
| `--condition-regex` | `^(?P<condition>.*?)\s+(?P<rep>\d+)$` | Regex splitting a `Sample` name into condition + replicate id, used when no `--sample-map` entry or DA3 Biogroup is available. |
| `--sample-map` | — | CSV with columns `Sample,Condition,BioRep` (optionally `Plate`) to explicitly declare condition/replicate per sample. Takes priority over DA3's Biogroup and over `--condition-regex`. |
| `--housekeeping` | — | Housekeeping gene/Target for the standard two-step ddCt method (e.g. `GAPDH`). Omit for experiments with no reference gene: if `--anchor` is still given, ddCt is computed directly against the anchor's own Cq per target (no housekeeping subtraction); if `--anchor` is also omitted, runs QC + raw Cq plots only. |
| `--sd-threshold` | `0.3` | Max allowed technical-replicate Cq SD before outlier removal kicks in. Ignored if `--skip-outlier-removal` is set. |
| `--skip-outlier-removal` | off | Bypass SD-based outlier removal entirely: average every valid replicate as-is; `--sd-threshold` is ignored. A group only fails if none of its replicates amplified. |
| `--exploratory` | off | Enable imputation of FAILED groups instead of dropping them. **`--fail-strategy` has no effect at all unless this flag is also passed** — without `--exploratory`, FAILED groups are simply left failed/dropped regardless of what `--fail-strategy` says (including its default). |
| `--fail-strategy` | `plate_mean` | Imputation strategy, only applied when `--exploratory` is also set: `drop`, `plate_mean`, `anchor_mean`, or `best_subset`. |
| `--alpha` | `0.05` | Significance threshold for normality/variance/omnibus tests. |
| `--rq-highlight-target-match` | off | In `plots/rq/`, color a condition's bar/points differently when the plotted target/gene name is a case-insensitive substring of the condition name (e.g. `CRISPRoff-HER2` vs. `CRISPRoff-RGS9` on the `HER2` plot). Leave off for assays where condition names don't encode a target gene. |
| `--dpi` | `150` | Plot resolution. |
| `--outdir` (required) | — | Output folder (created if missing). |

## Output structure

```
<outdir>/
├── run.log
├── tables/
│   ├── technical_replicates_QC.csv
│   ├── sample_target_Cq_summary.csv
│   ├── ddCt_RQ_log2FC_results.csv        (only if --anchor given)
│   ├── statistics_posthoc.csv            (only if --anchor given)
│   └── statistics_report.txt             (only if --anchor given)
└── plots/
    ├── qc/
    ├── cq/
    ├── log2fc/                           (only if --anchor given)
    │   ├── log2FC_<target>.png           (one per target)
    │   └── log2FC_all_targets.png        (all targets side by side)
    └── rq/                               (only if --anchor given)
        ├── RQ_<target>.png               (one per target)
        └── RQ_all_targets.png            (all targets side by side)
```

## Dependencies

```bash
pip install pandas numpy openpyxl scipy matplotlib seaborn scikit-posthocs statsmodels
```

`scipy.stats.dunnett` requires `scipy >= 1.11`.