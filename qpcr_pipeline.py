#!/usr/bin/env python3
"""
qPCR analysis tool which employs the updated SD-based replicate removal method.
Allows for exploratory imputation of failed (sample, target) groups using several strategies.
Includes statistical testing (parametric or non-parametric) and plotting of results.

WORKFLOW
1. Reads one or more QuantStudio qPCR export files (one file per plate). Each
   file must contain a "Results" sheet with the standard per-well columns
   (Well, Well Position, Sample, Target, Cq, ...).
2. Groups technical replicates by (plate, sample, target) and applies an
   SD-based outlier removal rule:
     - compute SD of the replicate group
     - if SD <= set threshold -> pass, keep all replicates
     - if SD > threshold  -> remove the replicate farthest from the group
       mean, recompute SD, repeat until below threshold
     - if only 2 replicates remain and SD is still > threshold -> the
       whole (sample, target) group FAILS QC (but allow exploration later)
   The SD threshold is a user parsed CLI argument (--sd-threshold).
3. In "exploratory" mode (--exploratory), failed (sample, target) groups
   can be salvaged instead of dropped, using one of several strategies
   (--fail-strategy): plate_mean, anchor_mean, best_subset, or drop.
   This is most useful for housekeeping genes, where losing a sample
   entirely because of one bad technical replicate is wasteful but
   exploratory imputation is reasonable/desired.
4. Computes the ddCt method per plate, using a user-specified anchor
   SAMPLE (one specific well/biological replicate used as the calibrator,
   e.g. "Empty 1") and housekeeping gene (default GAPDH):
       dCt        = Cq(target) - Cq(housekeeping)
       ddCt       = dCt(sample) - dCt(anchor)
       RQ         = 2^-ddCt
       log2FC     = -ddCt   (== log2(RQ))
5. Runs group-wise statistics on log2FC per target vs. a control condition
   (all biological replicates of the anchor's condition, not just the
   single anchor well): per-group Shapiro-Wilk normality + Levene's test
   for equal variances decide parametric vs. non-parametric, then:
     - parametric:      one-way ANOVA + Dunnett's post-hoc vs. control
     (can be changes but ANOVA used in previous analysis and this allows comparison)
     - non-parametric:  Kruskal-Wallis + Dunn's post-hoc vs. control
6. Writes result tables (technical-replicate level and
   sample/target-level with Cq, dCt, ddCt, RQ, log2FC) plus a QC log,
   and generates plots:
     - qc/            per-plate Cq-per-technical-replicate plots with
                       removed points marked (X), faceted by target
     - cq/             Cq per sample per target, raw and anchor-referenced
     - log2fc/         log2FC per condition per target with significance
                       brackets vs. the control condition

USAGE (PARSING CLI)
    python qpcr_pipeline.py \\
        --input plate1.xlsx plate2.xlsx plate3.xlsx \\
        --plate-labels Plate1 Plate2 Plate3 \\
        --anchor "Empty 1" "E 1" "E 1" \\
        --housekeeping GAPDH \\
        --sd-threshold 0.3 \\
        --exploratory --fail-strategy plate_mean \\
        --outdir ./qpcr_results

Run `python qpcr_pipeline.py --help` for the full argument list.

DEPENDENCIES
    pip install pandas numpy openpyxl scipy matplotlib seaborn scikit-posthocs statsmodels
(scipy.stats.dunnett requires scipy >= 1.11)
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import stats as sstats


sns.set_theme(style="whitegrid", context="talk")

LOG = logging.getLogger("qpcr_pipeline")

STAR_THRESHOLDS = [(1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (5e-2, "*")] 
def stars_for(p: Optional[float]) -> str:
    if p is None or np.isnan(p):
        return "n/a"
    for cutoff, label in STAR_THRESHOLDS:
        if p < cutoff:
            return label
    return "ns"


# 1. Reading QuantStudio "Results" sheets (has to be formatted beforehand)

REQUIRED_COLS = ["Well Position", "Sample", "Target", "Cq", "Omit", "Amp Status"]


def find_results_header_row(ws) -> int:
    for i, row in enumerate(ws.iter_rows(values_only=True), start=1):
        if row and row[0] == "Well":
            return i
    raise ValueError("Could not find the 'Well' header row in the Results sheet.")


def read_plate_results(path: Path, plate_label: str) -> pd.DataFrame:
    import openpyxl
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wb = openpyxl.load_workbook(path, data_only=True)
    if "Results" not in wb.sheetnames:
        raise ValueError(f"{path.name}: no 'Results' sheet found (sheets: {wb.sheetnames})")

    ws = wb["Results"]
    header_row = find_results_header_row(ws)
    headers = [c.value for c in ws[header_row]]

    records = []
    for row in ws.iter_rows(min_row=header_row + 1, values_only=True):
        if row[0] is None:
            continue
        records.append(row)

    df = pd.DataFrame(records, columns=headers)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: Results sheet is missing expected columns: {missing}")

    df["Cq"] = pd.to_numeric(df["Cq"], errors="coerce")
    df["Sample"] = df["Sample"].astype(str).str.strip()
    df["Target"] = df["Target"].astype(str).str.strip()
    df["Well Position"] = df["Well Position"].astype(str).str.strip()
    df["Omit"] = df["Omit"].astype(str).str.upper().eq("TRUE")
    df["Plate"] = plate_label
    df["SourceFile"] = path.name

    keep_cols = ["Plate", "SourceFile", "Well Position", "Sample", "Target",
                 "Cq", "Omit", "Amp Status"]
    return df[keep_cols].reset_index(drop=True)


# 2. Condition / biological-replicate parsing #has to e named accordingly

TRAILING_NUMBER_RE = re.compile(r"^(?P<condition>.*?)\s+(?P<rep>\d+)$")

def parse_condition(sample: str) -> tuple[str, str]:
    """Split a Sample label like 'CRISPRoff-HER2 3' into ('CRISPRoff-HER2', '3').

    If the Sample has no trailing replicate number, the whole string is
    treated as the condition and the replicate id defaults to '1'.
    """
    m = TRAILING_NUMBER_RE.match(sample)
    if m:
        return m.group("condition"), m.group("rep")
    return sample, "1"


# 3. SD-based outlier removal

@dataclass
class ReplicateResolution:
    plate: str
    sample: str
    target: str
    condition: str
    bio_rep: str
    kept_wells: list = field(default_factory=list)
    kept_cqs: list = field(default_factory=list)
    removed_wells: list = field(default_factory=list)
    removed_cqs: list = field(default_factory=list)
    removed_reasons: list = field(default_factory=list)
    n_start: int = 0
    final_sd: Optional[float] = None
    status: str = "PASS"          # PASS | PASS_AFTER_REMOVAL | FAILED
    fail_reason: Optional[str] = None
    cq_mean: Optional[float] = None
    imputed: bool = False
    imputed_note: Optional[str] = None


def resolve_replicate_group(plate, sample, target, condition, bio_rep,
                             wells, cqs, sd_threshold) -> ReplicateResolution:
    res = ReplicateResolution(plate=plate, sample=sample, target=target,
                               condition=condition, bio_rep=bio_rep)

    remaining = list(zip(wells, cqs))
    n_raw = len(remaining)

    # Drop replicates that never amplified (Cq is NaN) before SD logic starts, count as removed
    valid = [(w, c) for w, c in remaining if not pd.isna(c)]
    no_amp = [(w, c) for w, c in remaining if pd.isna(c)]
    for w, c in no_amp:
        res.removed_wells.append(w)
        res.removed_cqs.append(c)
        res.removed_reasons.append("no Cq / no amplification")

    remaining = valid
    res.n_start = n_raw

    while True:
        if len(remaining) < 2:
            res.status = "FAILED"
            res.fail_reason = "fewer than 2 valid replicates available"
            res.final_sd = np.nan
            break

        cq_vals = np.array([c for _, c in remaining], dtype=float)
        sd = float(np.std(cq_vals, ddof=1))

        if sd <= sd_threshold:
            res.status = "PASS" if not res.removed_wells else "PASS_AFTER_REMOVAL"
            res.final_sd = sd
            break

        if len(remaining) == 2:
            res.status = "FAILED"
            res.fail_reason = (
                f"SD={sd:.3f} still above threshold ({sd_threshold}) "
                f"with only 2 replicates left"
            )
            res.final_sd = sd
            break

        # remove the replicate farthest from the current mean
        mean = float(np.mean(cq_vals))
        idx = int(np.argmax(np.abs(cq_vals - mean)))
        bad_well, bad_cq = remaining.pop(idx)
        res.removed_wells.append(bad_well)
        res.removed_cqs.append(bad_cq)
        res.removed_reasons.append(
            f"largest deviation from mean ({abs(bad_cq - mean):.3f} Cq) while SD={sd:.3f}"
        )

    res.kept_wells = [w for w, _ in remaining]
    res.kept_cqs = [c for _, c in remaining]
    if remaining:
        res.cq_mean = float(np.mean([c for _, c in remaining]))
    return res


def run_qc(raw: pd.DataFrame, sd_threshold: float) -> tuple[pd.DataFrame, pd.DataFrame]:
_summary : one row per (plate, sample, target) with QC outcome + Cq mean
    rep_rows = []
    group_rows = []

    grouped = raw.groupby(["Plate", "Sample", "Target"], sort=False)
    for (plate, sample, target), sub in grouped:
        condition, bio_rep = parse_condition(sample)
        res = resolve_replicate_group(
            plate, sample, target, condition, bio_rep,
            sub["Well Position"].tolist(), sub["Cq"].tolist(),
            sd_threshold,
        )

        for w, c in zip(res.kept_wells, res.kept_cqs):
            rep_rows.append(dict(Plate=plate, Sample=sample, Condition=condition,
                                  BioRep=bio_rep, Target=target, Well=w, Cq=c,
                                  Removed=False, Reason=None))
        for w, c, r in zip(res.removed_wells, res.removed_cqs, res.removed_reasons):
            rep_rows.append(dict(Plate=plate, Sample=sample, Condition=condition,
                                  BioRep=bio_rep, Target=target, Well=w, Cq=c,
                                  Removed=True, Reason=r))

        group_rows.append(dict(
            Plate=plate, Sample=sample, Condition=condition, BioRep=bio_rep,
            Target=target, n_start=res.n_start, n_kept=len(res.kept_wells),
            n_removed=len(res.removed_wells), QC_status=res.status,
            fail_reason=res.fail_reason, Cq_mean=res.cq_mean,
            Cq_SD=res.final_sd, Imputed=False, Imputed_note=None,
        ))

    replicate_log = pd.DataFrame(rep_rows)
    group_summary = pd.DataFrame(group_rows)
    return replicate_log, group_summary

# 4. Exploratory imputation of FAILED groups

def apply_fail_strategy(group_summary: pd.DataFrame, raw: pd.DataFrame,
                         exploratory: bool, fail_strategy: str,
                         anchor_by_plate: dict) -> pd.DataFrame:
    gs = group_summary.copy()

    if not exploratory:
        # rigorous mode: failed groups simply stay NaN downstream
        return gs

    for idx, row in gs[gs["QC_status"] == "FAILED"].iterrows():
        plate, sample, target = row["Plate"], row["Sample"], row["Target"]

        if fail_strategy == "drop":
            continue

        elif fail_strategy == "plate_mean":
            others = gs[(gs["Plate"] == plate) & (gs["Target"] == target) &
                        (gs["Sample"] != sample) & (gs["QC_status"] != "FAILED")]
            if others.empty or others["Cq_mean"].isna().all():
                LOG.warning("plate_mean imputation impossible for %s / %s / %s "
                            "(no passing samples on this plate for this target)",
                            plate, sample, target)
                continue
            val = float(others["Cq_mean"].mean())
            gs.loc[idx, "Cq_mean"] = val
            gs.loc[idx, "Imputed"] = True
            gs.loc[idx, "Imputed_note"] = (
                f"plate_mean: mean Cq of {len(others)} other passing "
                f"'{target}' samples on {plate} = {val:.3f}"
            )

        elif fail_strategy == "anchor_mean":
            anchor_sample = anchor_by_plate.get(plate)
            anchor_row = gs[(gs["Plate"] == plate) & (gs["Sample"] == anchor_sample) &
                             (gs["Target"] == target)]
            if anchor_row.empty or pd.isna(anchor_row["Cq_mean"].iloc[0]):
                LOG.warning("anchor_mean imputation impossible for %s / %s / %s "
                            "(anchor sample has no valid Cq for this target)",
                            plate, sample, target)
                continue
            val = float(anchor_row["Cq_mean"].iloc[0])
            gs.loc[idx, "Cq_mean"] = val
            gs.loc[idx, "Imputed"] = True
            gs.loc[idx, "Imputed_note"] = f"anchor_mean: Cq of anchor '{anchor_sample}' = {val:.3f}"

        elif fail_strategy == "best_subset":
            sub = raw[(raw["Plate"] == plate) & (raw["Sample"] == sample) &
                      (raw["Target"] == target)]
            cqs = sub["Cq"].dropna().tolist()
            wells = sub.loc[sub["Cq"].notna(), "Well Position"].tolist()
            if len(cqs) < 2:
                LOG.warning("best_subset imputation impossible for %s / %s / %s "
                            "(fewer than 2 valid Cq values at all)", plate, sample, target)
                continue
            # choose the 2-replicate pair with the smallest difference
            best_pair, best_sd = None, np.inf
            for i in range(len(cqs)):
                for j in range(i + 1, len(cqs)):
                    pair_sd = float(np.std([cqs[i], cqs[j]], ddof=1))
                    if pair_sd < best_sd:
                        best_sd, best_pair = pair_sd, (cqs[i], cqs[j])
            val = float(np.mean(best_pair))
            gs.loc[idx, "Cq_mean"] = val
            gs.loc[idx, "Cq_SD"] = best_sd
            gs.loc[idx, "Imputed"] = True
            gs.loc[idx, "Imputed_note"] = (
                f"best_subset: forced-pass on best 2 of {len(cqs)} replicates, "
                f"SD={best_sd:.3f} (still above threshold)"
            )
        else:
            raise ValueError(f"Unknown fail-strategy: {fail_strategy}")

    return gs


# 5. ddCt method

def compute_ddct(group_summary: pd.DataFrame, housekeeping: str,
                  anchor_by_plate: dict) -> pd.DataFrame:
    rows = []

    for plate, plate_df in group_summary.groupby("Plate"):
        anchor_sample = anchor_by_plate.get(plate)
        if anchor_sample is None:
            raise ValueError(f"No anchor sample specified for plate '{plate}'.")

        hk = plate_df[plate_df["Target"] == housekeeping].set_index("Sample")["Cq_mean"]
        targets = sorted(t for t in plate_df["Target"].unique() if t != housekeeping)

        anchor_dct = {}
        for target in targets:
            tgt_series = plate_df[plate_df["Target"] == target].set_index("Sample")["Cq_mean"]
            if anchor_sample not in tgt_series.index or anchor_sample not in hk.index:
                anchor_dct[target] = np.nan
                continue
            anchor_dct[target] = tgt_series[anchor_sample] - hk[anchor_sample]

        for target in targets:
            tgt_df = plate_df[plate_df["Target"] == target]
            hk_df = plate_df[plate_df["Target"] == housekeeping]

            for _, row in tgt_df.iterrows():
                sample = row["Sample"]
                hk_row = hk_df[hk_df["Sample"] == sample]
                cq_hk = hk_row["Cq_mean"].iloc[0] if not hk_row.empty else np.nan
                cq_tgt = row["Cq_mean"]

                dct = cq_tgt - cq_hk if pd.notna(cq_tgt) and pd.notna(cq_hk) else np.nan
                a_dct = anchor_dct.get(target, np.nan)
                ddct = dct - a_dct if pd.notna(dct) and pd.notna(a_dct) else np.nan
                rq = 2 ** (-ddct) if pd.notna(ddct) else np.nan
                log2fc = -ddct if pd.notna(ddct) else np.nan

                rows.append(dict(
                    Plate=plate, Sample=sample, Condition=row["Condition"],
                    BioRep=row["BioRep"], Target=target,
                    Anchor_sample=anchor_sample, Housekeeping_gene=housekeeping,
                    Cq_target_mean=cq_tgt, Cq_target_SD=row["Cq_SD"],
                    Cq_housekeeping_mean=cq_hk,
                    QC_status_target=row["QC_status"], Imputed_target=row["Imputed"],
                    Imputed_note_target=row["Imputed_note"],
                    dCt=dct, anchor_dCt=a_dct, ddCt=ddct, RQ=rq, log2FC=log2fc,
                    is_anchor=(sample == anchor_sample),
                ))

    return pd.DataFrame(rows)


# 6. Statistics: normality -> parametric vs non-parametric -> vs control

@dataclass
class TargetStatsResult:
    target: str
    control_condition: str
    test_family: str                  
    omnibus_test: Optional[str] = None
    omnibus_p: Optional[float] = None
    normality: dict = field(default_factory=dict) 
    levene_p: Optional[float] = None
    posthoc: pd.DataFrame = None      
    note: Optional[str] = None


def run_stats_for_target(df_target: pd.DataFrame, control_condition: str,
                          alpha: float) -> TargetStatsResult:
    target_name = df_target["Target"].iloc[0]
    groups = {
        cond: sub["log2FC"].dropna().to_numpy()
        for cond, sub in df_target.groupby("Condition")
    }
    groups = {c: v for c, v in groups.items() if len(v) >= 1}

    if control_condition not in groups:
        return TargetStatsResult(target=target_name, control_condition=control_condition,
                                  test_family="insufficient_data",
                                  note="control condition has no usable log2FC values")

    if len(groups) < 2:
        return TargetStatsResult(target=target_name, control_condition=control_condition,
                                  test_family="insufficient_data",
                                  note="fewer than 2 conditions present for this target")

    # Conditions with <2 replicates can be compared with a
    # rank-based test, but variance-based methods (ANOVA/Dunnett) are not
    # valid for them -> force non-parametric 
    low_n = sorted(c for c, v in groups.items() if len(v) < 2)

    normality = {}
    all_normal = True
    for cond, vals in groups.items():
        if len(vals) >= 3:
            stat, p = sstats.shapiro(vals)
            normality[cond] = (float(stat), float(p))
            if p < alpha:
                all_normal = False
        else:
            normality[cond] = (np.nan, np.nan)
            all_normal = False  # can't confirm normality with n<3 -> be conservative

    groups_ge2 = {c: v for c, v in groups.items() if len(v) >= 2}
    levene_p = None
    if len(groups_ge2) >= 2:
        try:
            _, levene_p = sstats.levene(*groups_ge2.values())
        except Exception:
            levene_p = np.nan

    parametric = (all_normal and not low_n and
                  (levene_p is None or np.isnan(levene_p) or levene_p >= alpha))

    note_parts = []
    if low_n:
        note_parts.append(
            "condition(s) with <2 replicates (no within-group variance estimate; "
            f"forced non-parametric, inference for these is weak): {', '.join(low_n)}"
        )

    result = TargetStatsResult(
        target=target_name, control_condition=control_condition,
        test_family="parametric" if parametric else "nonparametric",
        normality=normality, levene_p=levene_p,
        note="; ".join(note_parts) or None,
    )

    other_conditions = [c for c in groups if c != control_condition]
    if not other_conditions:
        result.test_family = "insufficient_data"
        result.note = "no non-control conditions to compare"
        return result

    if parametric:
        try:
            f_stat, p_omni = sstats.f_oneway(*groups.values())
            result.omnibus_test = "one-way ANOVA"
            result.omnibus_p = float(p_omni)
        except Exception as exc:
            result.note = (result.note + " | " if result.note else "") + f"ANOVA failed: {exc}"

        try:
            dun = sstats.dunnett(
                *[groups[c] for c in other_conditions],
                control=groups[control_condition],
            )
            posthoc_rows = [
                dict(Condition=c, p_adj=float(p), stars=stars_for(p))
                for c, p in zip(other_conditions, dun.pvalue)
            ]
            result.posthoc = pd.DataFrame(posthoc_rows)
        except Exception as exc:
            result.note = (result.note + " | " if result.note else "") + f"Dunnett failed: {exc}"

    else:
        try:
            h_stat, p_omni = sstats.kruskal(*groups.values())
            result.omnibus_test = "Kruskal-Wallis"
            result.omnibus_p = float(p_omni)
        except Exception as exc:
            result.note = (result.note + " | " if result.note else "") + f"Kruskal-Wallis failed: {exc}"

        if sp is None:
            result.note = (result.note + " | " if result.note else "") + \
                "scikit-posthocs not installed; skipping Dunn's test"
        else:
            long_df = pd.concat(
                [pd.DataFrame({"Condition": c, "log2FC": v}) for c, v in groups.items()],
                ignore_index=True,
            )
            try:
                dunn = sp.posthoc_dunn(long_df, val_col="log2FC", group_col="Condition",
                                        p_adjust="holm")
                posthoc_rows = [
                    dict(Condition=c, p_adj=float(dunn.loc[control_condition, c]),
                         stars=stars_for(dunn.loc[control_condition, c]))
                    for c in other_conditions
                ]
                result.posthoc = pd.DataFrame(posthoc_rows)
            except Exception as exc:
                result.note = (result.note + " | " if result.note else "") + f"Dunn's test failed: {exc}"

    return result


def run_all_stats(ddct_df: pd.DataFrame, control_condition: str, alpha: float):
    results = {}
    for target, sub in ddct_df.groupby("Target"):
        results[target] = run_stats_for_target(sub, control_condition, alpha)
    return results

# 7. Plots

def plot_qc_replicates(replicate_log: pd.DataFrame, outdir: Path, dpi: int):
    outdir.mkdir(parents=True, exist_ok=True)
    for plate, pdf in replicate_log.groupby("Plate"):
        targets = sorted(pdf["Target"].unique())
        samples = sorted(pdf["Sample"].unique(),
                          key=lambda s: parse_condition(s))
        fig, axes = plt.subplots(len(targets), 1, figsize=(max(8, 0.5 * len(samples)), 4 * len(targets)),
                                  sharex=True)
        if len(targets) == 1:
            axes = [axes]

        for ax, target in zip(axes, targets):
            tdf = pdf[pdf["Target"] == target]
            for sample in samples:
                sdf = tdf[tdf["Sample"] == sample].sort_values("Well")
                x = [sample] * len(sdf)
                kept = sdf[~sdf["Removed"]]
                removed = sdf[sdf["Removed"] & sdf["Cq"].notna()]
                ax.scatter([sample] * len(kept), kept["Cq"], color="steelblue", zorder=3)
                ax.scatter([sample] * len(removed), removed["Cq"], color="crimson",
                           marker="x", s=90, zorder=4, label="Removed (outlier)")
            ax.set_ylabel(f"{target}\nCq")
            ax.grid(True, alpha=0.3)

        handles, labels = axes[0].get_legend_handles_labels()
        plt.setp(axes[-1].get_xticklabels(), rotation=60, ha="right")
        fig.suptitle(f"QC: Cq per technical replicate — {plate}", y=0.995)
        if handles:
            fig.legend(handles[:1], labels[:1], loc="lower center",
                       bbox_to_anchor=(0.5, 0.0), ncol=1)
        fig.tight_layout(rect=[0, 0.05, 1, 0.96])
        fig.savefig(outdir / f"qc_replicates_{plate}.png", dpi=dpi)
        plt.close(fig)


def plot_cq_by_sample(group_summary: pd.DataFrame, anchor_by_plate: dict,
                       outdir: Path, dpi: int):
    outdir.mkdir(parents=True, exist_ok=True)
    for plate, pdf in group_summary.groupby("Plate"):
        anchor_sample = anchor_by_plate.get(plate)
        targets = sorted(pdf["Target"].unique())
        samples = sorted(pdf["Sample"].unique(), key=lambda s: parse_condition(s))

        for suffix, with_anchor in [("raw", False), ("with_anchor", True)]:
            fig, axes = plt.subplots(len(targets), 1,
                                      figsize=(max(8, 0.5 * len(samples)), 4 * len(targets)),
                                      sharex=True)
            if len(targets) == 1:
                axes = [axes]
            for ax, target in zip(axes, targets):
                tdf = pdf[pdf["Target"] == target].set_index("Sample").reindex(samples)
                colors = ["crimson" if (with_anchor and s == anchor_sample) else "steelblue"
                          for s in samples]
                ax.errorbar(samples, tdf["Cq_mean"], yerr=tdf["Cq_SD"], fmt="o",
                             ecolor="gray", capsize=3, color="steelblue", zorder=3)
                ax.scatter(samples, tdf["Cq_mean"], c=colors, zorder=4)
                if with_anchor and anchor_sample in tdf.index and pd.notna(tdf.loc[anchor_sample, "Cq_mean"]):
                    ax.axhline(tdf.loc[anchor_sample, "Cq_mean"], color="crimson",
                                linestyle="--", alpha=0.6,
                                label=f"anchor ({anchor_sample})")
                    ax.legend(fontsize=9)
                ax.set_ylabel(f"{target}\nCq (mean \u00b1 SD)")
                ax.grid(True, alpha=0.3)
            plt.setp(axes[-1].get_xticklabels(), rotation=60, ha="right")
            title = f"Cq per sample — {plate}" + (" (anchor highlighted)" if with_anchor else " (raw)")
            fig.suptitle(title)
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            fig.savefig(outdir / f"cq_by_sample_{suffix}_{plate}.png", dpi=dpi)
            plt.close(fig)


def _draw_significance_brackets(ax, x_positions: dict, control: str,
                                 posthoc: pd.DataFrame, y_top: float, y_step: float):
    if posthoc is None or posthoc.empty:
        return
    level = 0
    for _, row in posthoc.sort_values("Condition").iterrows():
        cond = row["Condition"]
        if cond not in x_positions or control not in x_positions:
            continue
        x1, x2 = x_positions[control], x_positions[cond]
        y = y_top + level * y_step
        ax.plot([x1, x1, x2, x2], [y, y + y_step * 0.15, y + y_step * 0.15, y],
                color="black", linewidth=1.2)
        ax.text((x1 + x2) / 2, y + y_step * 0.18, row["stars"],
                ha="center", va="bottom", fontsize=13)
        level += 1


def plot_log2fc(ddct_df: pd.DataFrame, stats_results: dict, control_condition: str,
                 outdir: Path, dpi: int):
    outdir.mkdir(parents=True, exist_ok=True)
    for target, sub in ddct_df.groupby("Target"):
        sub = sub.dropna(subset=["log2FC"])
        if sub.empty:
            continue
        order = sorted(sub["Condition"].unique(),
                        key=lambda c: (c != control_condition, c))
        x_positions = {c: i for i, c in enumerate(order)}

        fig, ax = plt.subplots(figsize=(max(6, 1.3 * len(order)), 6))
        sns.stripplot(data=sub, x="Condition", y="log2FC", order=order,
                       ax=ax, size=8, jitter=0.15, alpha=0.85)
        means = sub.groupby("Condition")["log2FC"].mean().reindex(order)
        ax.scatter(range(len(order)), means, color="black", marker="_", s=800, zorder=5)
        ax.axhline(0, color="gray", linestyle=":", linewidth=1)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=60, ha="right")
        ax.set_ylabel("log2FC")
        ax.set_title(f"{target} — log2FC vs. {control_condition}")

        res = stats_results.get(target)
        if res is not None and res.posthoc is not None and not res.posthoc.empty:
            data_min, data_max = sub["log2FC"].min(), sub["log2FC"].max()
            data_range = max(data_max - data_min, 0.5)
            y_top = data_max + 0.15 * data_range
            y_step = 0.14 * data_range
            n_brackets = len(res.posthoc)
            _draw_significance_brackets(ax, x_positions, control_condition,
                                         res.posthoc, y_top, y_step)
            ax.set_ylim(data_min - 0.15 * data_range,
                        y_top + (n_brackets + 0.8) * y_step)

            fam = res.test_family
            omni = f"{res.omnibus_test}: p={res.omnibus_p:.4g}" if res.omnibus_p is not None else ""
            ax.text(0.01, 0.02, f"{fam}\n{omni}", transform=ax.transAxes,
                     ha="left", va="bottom", fontsize=9, color="dimgray")

        fig.tight_layout()
        safe_target = re.sub(r"[^\w\-.]", "_", target)
        fig.savefig(outdir / f"log2FC_{safe_target}.png", dpi=dpi)
        plt.close(fig)


# 8. Report

def write_stats_report(stats_results: dict, outpath: Path):
    lines = []
    for target, res in stats_results.items():
        lines.append(f"## Target: {target}")
        lines.append(f"Control condition: {res.control_condition}")
        lines.append(f"Test family: {res.test_family}")
        if res.note:
            lines.append(f"Note: {res.note}")
        if res.normality:
            lines.append("Shapiro-Wilk normality (per condition, n>=3 required):")
            for cond, (stat, p) in res.normality.items():
                p_str = f"{p:.4g}" if not np.isnan(p) else "n/a (n<3)"
                lines.append(f"  - {cond}: p={p_str}")
        if res.levene_p is not None:
            lp = f"{res.levene_p:.4g}" if not (res.levene_p is None or np.isnan(res.levene_p)) else "n/a"
            lines.append(f"Levene's test for equal variances: p={lp}")
        if res.omnibus_test:
            lines.append(f"Omnibus test: {res.omnibus_test}, p={res.omnibus_p:.4g}")
        if res.posthoc is not None and not res.posthoc.empty:
            lines.append("Post-hoc vs. control:")
            for _, row in res.posthoc.iterrows():
                lines.append(f"  - {row['Condition']}: p_adj={row['p_adj']:.4g} ({row['stars']})")
        lines.append("")
    outpath.write_text("\n".join(lines), encoding="utf-8")

# 9. CLI parsing and main entry point

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Automated qPCR ddCt analysis: QC/outlier removal, ddCt, "
                    "RQ/log2FC, statistics, and plots.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--input", nargs="+", required=True, type=Path,
                    help="One or more QuantStudio export .xlsx files, one per plate.")
    p.add_argument("--plate-labels", nargs="+", default=None,
                    help="Custom plate labels, matched by order to --input "
                         "(default: file stem).")
    p.add_argument("--anchor", nargs="+", required=True,
                    help="Anchor SAMPLE name (exact 'Sample' text, e.g. 'Empty 1'), "
                         "one per --input file in the same order.")
    p.add_argument("--control-condition", default=None,
                    help="Baseline condition group used for statistical comparisons "
                         "(default: condition parsed from the first plate's anchor "
                         "sample, e.g. 'E 1' -> 'E').")
    p.add_argument("--housekeeping", default="GAPDH",
                    help="Housekeeping gene / Target name (default: GAPDH).")
    p.add_argument("--sd-threshold", type=float, default=0.3,
                    help="Max allowed technical-replicate Cq SD before outlier "
                         "removal kicks in (default: 0.3).")
    p.add_argument("--exploratory", action="store_true", default=False, # turned off by default
                    help="Enable imputation of FAILED (sample, target) groups "
                         "instead of dropping them. Off by default (rigorous mode).")
    p.add_argument("--fail-strategy", choices=["drop", "plate_mean", "anchor_mean", "best_subset"],
                    default="plate_mean",
                    help="How to salvage FAILED groups when --exploratory is set "
                         "(default: plate_mean). Ignored if --exploratory is not set.")
    p.add_argument("--alpha", type=float, default=0.05,
                    help="Significance threshold for normality/variance/omnibus tests "
                         "(default: 0.05).")
    p.add_argument("--dpi", type=int, default=150, help="Plot resolution (default: 150).")
    p.add_argument("--outdir", required=True, type=Path,
                    help="Output folder (created if it doesn't exist).")
    return p


def main(argv=None):
    args = build_arg_parser().parse_args(argv)

    if len(args.anchor) != len(args.input):
        sys.exit(f"--anchor must have exactly one value per --input file "
                  f"({len(args.input)} inputs, {len(args.anchor)} anchors given).")

    plate_labels = args.plate_labels or [p.stem for p in args.input]
    if len(plate_labels) != len(args.input):
        sys.exit("--plate-labels must have exactly one value per --input file.")

    args.outdir.mkdir(parents=True, exist_ok=True)
    log_path = args.outdir / "run.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_path, mode="w"), logging.StreamHandler(sys.stdout)],
    )
    LOG.info("qpcr_pipeline starting")
    LOG.info("Inputs: %s", [str(p) for p in args.input])
    LOG.info("Plate labels: %s", plate_labels)
    LOG.info("Anchors: %s", args.anchor)
    LOG.info("SD threshold: %s | exploratory=%s | fail_strategy=%s",
              args.sd_threshold, args.exploratory, args.fail_strategy)

    anchor_by_plate = dict(zip(plate_labels, args.anchor))

    # 1. read & combine
    raw_frames = []
    for path, label in zip(args.input, plate_labels):
        LOG.info("Reading %s as plate '%s'", path, label)
        raw_frames.append(read_plate_results(path, label))
    raw = pd.concat(raw_frames, ignore_index=True)

    # 2. QC / outlier resolution 
    LOG.info("Running SD-based outlier resolution (threshold=%s)", args.sd_threshold)
    replicate_log, group_summary = run_qc(raw, args.sd_threshold)

    n_failed = (group_summary["QC_status"] == "FAILED").sum()
    LOG.info("QC complete: %d/%d (sample, target) groups FAILED",
              n_failed, len(group_summary))

    # 3. exploratory imputation 
    group_summary = apply_fail_strategy(group_summary, raw, args.exploratory,
                                         args.fail_strategy, anchor_by_plate)
    if args.exploratory:
        n_imputed = group_summary["Imputed"].sum()
        LOG.info("Exploratory mode ('%s'): %d groups imputed", args.fail_strategy, n_imputed)

    # 4. ddCt 
    LOG.info("Computing ddCt / RQ / log2FC (housekeeping=%s)", args.housekeeping)
    ddct_df = compute_ddct(group_summary, args.housekeeping, anchor_by_plate)

    control_condition = args.control_condition
    if control_condition is None:
        first_anchor = args.anchor[0]
        control_condition, _ = parse_condition(first_anchor)
        LOG.info("--control-condition not given; inferred '%s' from first anchor '%s'",
                  control_condition, first_anchor)

    # 5. statistics 
    LOG.info("Running statistics vs. control condition '%s'", control_condition)
    stats_results = run_all_stats(ddct_df, control_condition, args.alpha)

    # 6. write tables 
    tables_dir = args.outdir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    replicate_log.to_csv(tables_dir / "technical_replicates_QC.csv", index=False)
    group_summary.to_csv(tables_dir / "sample_target_Cq_summary.csv", index=False)
    ddct_df.to_csv(tables_dir / "ddCt_RQ_log2FC_results.csv", index=False)

    posthoc_frames = []
    for target, res in stats_results.items():
        if res.posthoc is not None and not res.posthoc.empty:
            f = res.posthoc.copy()
            f.insert(0, "Target", target)
            f.insert(1, "Control_condition", res.control_condition)
            f.insert(2, "Test_family", res.test_family)
            f.insert(3, "Omnibus_test", res.omnibus_test)
            f.insert(4, "Omnibus_p", res.omnibus_p)
            posthoc_frames.append(f)
    if posthoc_frames:
        pd.concat(posthoc_frames, ignore_index=True).to_csv(
            tables_dir / "statistics_posthoc.csv", index=False)
    write_stats_report(stats_results, tables_dir / "statistics_report.txt")

    # 7. plots 
    LOG.info("Generating plots")
    plot_qc_replicates(replicate_log, args.outdir / "plots" / "qc", args.dpi)
    plot_cq_by_sample(group_summary, anchor_by_plate, args.outdir / "plots" / "cq", args.dpi)
    plot_log2fc(ddct_df, stats_results, control_condition, args.outdir / "plots" / "log2fc", args.dpi)

    LOG.info("Done. Results written to %s", args.outdir.resolve())
    print(f"\nDone. See:\n  {tables_dir}\n  {args.outdir / 'plots'}\n  {log_path}")


if __name__ == "__main__":
    main()
