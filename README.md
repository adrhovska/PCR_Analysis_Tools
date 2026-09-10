# qPCR SD-outlier-removal & ddCt automation pipeline

Automates the workflow: QuantStudio (DA3) "Results" export followed by SD-based
outlier removal on technical replicates
(optionally) ddCt method and RQ generation with normality-gated statistics (ANOVA+Dunnett or Kruskal-Wallis+Dunn's)
Accompanied by relevant plots

The ddCt/statistics stage is optional, as the same script also works as a
standalone QC + raw-Cq tool for experiments that don't have a
housekeeping/reference gene (e.g. a cloning or junction-validation qPCR).
See "Two modes" below.

## Install

```bash
pip install pandas numpy openpyxl scipy matplotlib seaborn scikit-posthocs statsmodels
```
(`scipy.stats.dunnett` requires scipy ≥ 1.11 — `pip install -U scipy` if needed.)

## Two modes

### 1. Relative quantification (ddCt) — pass `--housekeeping`

Use this when you have a reference/housekeeping gene and want RQ, log2FC,
and statistics vs. a control condition.

```bash
python qpcr_pipeline.py \
  --input Plate1.xlsx Plate2.xlsx Plate3.xlsx \
  --plate-labels Plate1 Plate2 Plate3 \
  --anchor "Empty 1" "E 1" "E 1" \
  --control-condition Empty \
  --housekeeping GAPDH \
  --sd-threshold 0.3 \
  --exploratory --fail-strategy plate_mean \
  --outdir ./results
```

`--anchor` is **required** in this mode — one value per `--input` file, in
the same order, matching the exact `Sample` text in that plate's Results
sheet (e.g. `"E 1"`). It's the single well used as the ddCt calibrator for
that plate.

`--control-condition` is the *condition* (not a single well) used as the
baseline group in the statistics — e.g. all `Empty` replicates, not just
the one anchor well. If omitted, it's inferred from the first plate's
`--anchor` value via the condition-resolution logic below.

### 2. QC + raw Cq only — omit `--housekeeping`

Use this for experiments with no housekeeping gene at all.

```bash
python qpcr_pipeline.py \
  --input MyPlate.xlsx \
  --plate-labels PlateA \
  --sd-threshold 0.3 \
  --outdir ./results
```

No `--anchor` needed.

`--anchor` can still be given *without* `--housekeeping`, if you just want
that one sample highlighted with a dashed reference line on the Cq plot —
no ddCt math, purely visual. 

Passing `--housekeeping` without `--anchor` is rejected with a clear error
(ddCt can't be calibrated without an anchor well).

## What you get in `--outdir`

```
run.log                                  full parameter + step log
tables/
  technical_replicates_QC.csv            every well: kept/removed + reason
  sample_target_Cq_summary.csv           per (plate, sample, target): Cq mean/SD, QC status, imputation
  ddCt_RQ_log2FC_results.csv             [only if --housekeeping given] Cq, dCt, ddCt, RQ, log2FC
  statistics_posthoc.csv                 [only if --housekeeping given] post-hoc p-values/stars vs. control
  statistics_report.txt                  [only if --housekeeping given] normality, omnibus, post-hoc summary
plots/
  qc/qc_replicates_<plate>.png           Cq per technical replicate, outliers marked with X
  cq/cq_by_sample_raw_<plate>.png        Cq mean±SD per sample, no anchor reference
  cq/cq_by_sample_with_anchor_<plate>.png  [only if --anchor given] same, with anchor well highlighted
  log2fc/log2FC_<target>.png             [only if --housekeeping given] log2FC per condition with brackets
```

## How Sample names get grouped into biological conditions

You don't need your wells named with a trailing replicate number. Each
Sample is resolved to a (Condition, BioRep) pair using, in priority order:

1. **`--sample-map some.csv`** — an explicit override table with columns
   `Sample,Condition,BioRep` (and an optional `Plate` column; leave it
   blank to apply across all plates). Use this when neither of the options
   below applies — e.g. Sample names carry no replicate info at all and
   the DA3 "Biological Group" field was never set.
   ```csv
   Sample,Condition,BioRep
   HER2_KD_A,CRISPRoff-HER2,1
   HER2_KD_B,CRISPRoff-HER2,2
   Ctrl_A,Empty,1
   ```
2. **DA3's own "Biogroup" field** — QuantStudio Design & Analysis exports a
   `Biogroup` column in the "Replicate Group Result" sheet, populated only
   if you set a "Biological Group" during plate setup in the DA3 app. The
   script reads this automatically per plate; if it's populated, it's used
   with no extra flags needed. If it's blank (the common case — this field
   is easy to skip when setting up a run), the script falls through to (3)
   and says so in `run.log`.
3. **`--condition-regex`** — a regex applied to the Sample string, default
   `^(?P<condition>.*?)\s+(?P<rep>\d+)$` (i.e. "CRISPRoff-HER2 3" ->
   condition `CRISPRoff-HER2`, replicate `3`). Override this if your
   naming scheme differs (e.g. an underscore or letter suffix instead of a
   trailing number) — just make sure your regex has a `(?P<condition>...)`
   group and, optionally, a `(?P<rep>...)` group.

`run.log` always reports which of the three sources was actually used for
each plate, so you can confirm it did what you expect.

## Outlier rule (--sd-threshold)

For each (plate, sample, target) group of technical replicates:
1. Compute SD. If ≤ threshold → pass.
2. If not, remove whichever replicate is farthest from the group mean,
   recompute SD, repeat.
3. If SD is still above threshold with only 2 replicates left → the whole
   group FAILS (not silently averaged).

## Exploratory mode (--exploratory / --fail-strategy)

By default (no `--exploratory`), FAILED groups are left as missing (`NaN`
Cq) and drop out of any downstream ddCt/statistics — this is the
rigorous/conservative mode, and the right choice whenever "no valid
replicates" is itself a real, informative result (e.g. a target genuinely
didn't amplify) rather than a technical artifact to paper over.

With `--exploratory`, FAILED groups are instead imputed using one of:

- **`drop`** — same as rigorous mode (explicit no-op, for completeness).
- **`plate_mean`** — mean Cq of that target across the plate's other
  passing samples. Only appropriate for a **housekeeping/reference gene**,
  where "roughly the same across samples" is the expected biology — using
  it on a diagnostic target would erase the very difference you're trying
  to measure for the failed sample.
- **`anchor_mean`** — Cq of the anchor sample for that target. Same caveat
  as `plate_mean`: appropriate for a reference value, not for a target
  whose whole point is to vary between samples.
- **`best_subset`** — force-pass using whichever 2 of that *same sample's*
  original replicates have the smallest pairwise SD, even though it's
  still above threshold. This is the right choice for noisy-but-real
  diagnostic targets (the kind of failure where the sample did amplify,
  just with more technical scatter than your threshold allows), because it
  only uses that sample's own data — it never borrows a value from a
  different sample or condition. It correctly refuses to invent a value
  when fewer than 2 replicates actually amplified (a true negative stays
  missing, as it should).

Every imputed value is flagged in `Imputed_target` / `Imputed_note_target`
(or `Imputed` / `Imputed_note` in `sample_target_Cq_summary.csv`) so you
can always trace which numbers are real vs. salvaged, and how.

**Rule of thumb:** `plate_mean`/`anchor_mean` for housekeeping genes that
failed QC; `best_subset` for a diagnostic/target gene that failed QC due to
noise, not absence; `drop` (or just skip `--exploratory`) whenever a
missing value should stay missing.

## Statistics (only runs if --housekeeping is given)

Per target, per condition group:
- Shapiro-Wilk normality (only computable for n≥3; n<2 groups are flagged
  and forced into the non-parametric path since within-group variance
  can't be estimated at all)
- Levene's test for equal variances across groups with n≥2
- If all groups normal AND variances equal → one-way ANOVA + Dunnett's
  post-hoc vs. the control condition
- Otherwise → Kruskal-Wallis + Dunn's post-hoc (Holm-adjusted) vs. the
  control condition

Any group with fewer than 2 replicates is explicitly called out in the
report/plot notes rather than silently dropped, since it means the
p-values involving that group should be treated with caution.

## Troubleshooting

- **`FileNotFoundError` on your `.xlsx`**: the path in `--input` is
  relative to wherever your terminal's current directory is, not to the
  script's location. Use a full path, and quote it if the filename has
  spaces (e.g. `--input "/Users/you/Desktop/My File.xlsx"`). On macOS, you
  can get the exact path via Finder → right-click the file → hold
  **Option** → "Copy '...' as Pathname".
- **`--anchor is required when --housekeeping is set`**: ddCt can't be
  calibrated without an anchor well — either add `--anchor`, or drop
  `--housekeeping` if this experiment doesn't actually have a relative-
  quantification design (see "Two modes" above).
- **Every ddCt/log2FC value is `NaN`**: almost always means the
  `--housekeeping` name doesn't match any `Target` in the Results sheet —
  double check spelling/case, or reconsider whether this experiment has a
  housekeeping gene at all.
