#!/usr/bin/env python3
"""Build the Savant LCA stable analysis master v1.

The builder is deliberately additive: it preserves frozen v2.11 fields, fills
only the 48 SSC Module 4 CSS gaps, and maps the existing Korean clinical CSS
without recalculation.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


SSC = "SSC"
KOREAN = "SNUBH/Korean"
SNUBH = "SNU_Bundang_Hospital"


def _numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _metric(rows: list[dict], metric: str, value: int | float | str, detail: str = "") -> None:
    rows.append({"metric": metric, "value": value, "detail": detail})


def combine_css(
    subjects: pd.DataFrame,
    m13: pd.DataFrame,
    m4: pd.DataFrame,
    ssc_pheno: pd.DataFrame,
    korean_pheno: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Preserve M1-3 CSS, fill only M4 gaps, and map Korean source CSS."""
    out = subjects.copy()
    out["individual"] = out["individual"].astype("string")
    ssc_ids = set(out.loc[out["dataset"] == SSC, "individual"].dropna())

    m13_map = (
        m13.assign(individual=m13["individual"].astype("string"), css=_numeric(m13["total_css"]))
        .dropna(subset=["individual", "css"])
        .drop_duplicates("individual")
        .set_index("individual")["css"]
    )
    m4_map = (
        m4.assign(individual=m4["individual"].astype("string"), css=_numeric(m4["Total_CSS"]))
        .dropna(subset=["individual", "css"])
        .drop_duplicates("individual")
        .set_index("individual")["css"]
    )
    m13_current = set(m13_map.index) & ssc_ids
    m4_current = set(m4_map.index) & ssc_ids
    overlap = m13_current & m4_current
    if overlap:
        raise ValueError(f"SSC M1-3/Module4 overlap detected: {len(overlap)}")
    missing_ids = ssc_ids - (m13_current | m4_current)
    if missing_ids:
        raise ValueError(f"SSC CSS union leaves {len(missing_ids)} current subjects missing")

    old_ssc = _numeric(out.loc[out["dataset"] == SSC, "ados_css"])
    old_ids = out.loc[out["dataset"] == SSC, "individual"]
    expected_old = old_ids.map(m13_map)
    changed = old_ssc.notna() & ~np.isclose(old_ssc, expected_old, equal_nan=False)
    if changed.any():
        raise ValueError(f"existing SSC M1-3 CSS would change for {int(changed.sum())} subjects")

    out["ados_css_final"] = np.nan
    out["ados_css_source"] = pd.NA
    is_ssc = out["dataset"] == SSC
    out.loc[is_ssc, "ados_css_final"] = out.loc[is_ssc, "individual"].map(m13_map)
    out.loc[is_ssc & out["ados_css_final"].notna(), "ados_css_source"] = "SSC_M1to3"
    need_m4 = is_ssc & out["ados_css_final"].isna()
    out.loc[need_m4, "ados_css_final"] = out.loc[need_m4, "individual"].map(m4_map)
    out.loc[need_m4 & out["ados_css_final"].notna(), "ados_css_source"] = "SSC_Module4"

    is_korean = out["dataset"] == KOREAN
    korean_map = (
        korean_pheno.assign(
            Sample_ID=korean_pheno["sample_id"].astype("string"),
            source_css=_numeric(korean_pheno["ADOS_Total"]),
        )
        .drop_duplicates("Sample_ID")
        .set_index("Sample_ID")["source_css"]
    )
    korean_source_css = out.loc[is_korean, "Sample_ID"].map(korean_map)
    korean_intermediate = _numeric(out.loc[is_korean, "ados_total"])
    intermediate_mismatch = korean_intermediate.notna() & ~np.isclose(
        korean_intermediate, korean_source_css, equal_nan=False
    )
    out.loc[is_korean, "ados_css_final"] = korean_source_css
    out.loc[is_korean & out["ados_css_final"].notna(), "ados_css_source"] = "Korean_clinical_source"

    pheno_map = (
        ssc_pheno.assign(
            individual=ssc_pheno["sample_id"].astype("string"),
            pheno_css=_numeric(ssc_pheno["ADOS_total"]),
        )
        .drop_duplicates("individual")
        .set_index("individual")["pheno_css"]
    )
    ssc_final = out.loc[is_ssc, ["individual", "ados_css_final"]].copy()
    ssc_final["pheno_css"] = ssc_final["individual"].map(pheno_map)
    pheno_match = np.isclose(ssc_final["ados_css_final"], ssc_final["pheno_css"], equal_nan=False)
    if not pheno_match.all():
        raise ValueError(f"SSC union differs from phenotype ADOS_total for {int((~pheno_match).sum())} subjects")
    if out["ados_css_final"].isna().any():
        raise ValueError(f"final CSS missing for {int(out['ados_css_final'].isna().sum())} subjects")

    rows: list[dict] = []
    _metric(rows, "ssc_current_n", len(ssc_ids))
    _metric(rows, "ssc_m1to3_n", len(m13_current))
    _metric(rows, "ssc_module4_n", len(m4_current))
    _metric(rows, "ssc_source_overlap", len(overlap))
    _metric(rows, "ssc_union_gap", len(missing_ids))
    _metric(rows, "ssc_existing_value_changes", int(changed.sum()))
    _metric(rows, "ssc_pheno_exact", int(pheno_match.sum()))
    _metric(rows, "korean_existing_css_n", int(out.loc[is_korean, "ados_css_final"].notna().sum()))
    _metric(rows, "korean_intermediate_mismatches", int(intermediate_mismatch.sum()))
    _metric(rows, "final_css_missing", int(out["ados_css_final"].isna().sum()))
    return out, pd.DataFrame(rows)


def normalize_korean_clinical(raw: pd.DataFrame) -> pd.DataFrame:
    """Normalize the two-part SNP-D identifier and source comparison score."""
    family_col, subject_col = raw.columns[:2]

    def fmt(value: object, width: int) -> str | None:
        if pd.isna(value):
            return None
        text = str(value).strip()
        if text.endswith(".0"):
            text = text[:-2]
        return text.zfill(width)

    family = raw[family_col].map(lambda x: fmt(x, 3))
    subject = raw[subject_col].map(lambda x: fmt(x, 2))
    normalized = pd.DataFrame(
        {
            "Sample_ID": [
                f"SNP-D-{fam}-{sub}" if fam is not None and sub is not None else pd.NA
                for fam, sub in zip(family, subject)
            ],
            "source_css": _numeric(raw["비교점수"]),
        }
    )
    normalized = normalized.dropna(subset=["Sample_ID"]).drop_duplicates("Sample_ID")
    return normalized


def validate_korean_clinical(subjects: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    korean = subjects.loc[subjects["dataset"] == KOREAN, ["Sample_ID", "ados_css_final"]].copy()
    joined = korean.merge(clinical, on="Sample_ID", how="left", validate="one_to_one")
    current = _numeric(joined["ados_css_final"])
    source = _numeric(joined["source_css"])
    matched = source.notna()
    exact = matched & np.isclose(current, source, equal_nan=False)
    if not matched.all() or not exact.all():
        raise ValueError(
            "Korean clinical CSS validation failed: "
            f"matched={int(matched.sum())}/{len(joined)}, exact={int(exact.sum())}/{len(joined)}"
        )
    rows: list[dict] = []
    _metric(rows, "korean_current_n", len(joined))
    _metric(rows, "korean_source_css_matched", int(matched.sum()))
    _metric(rows, "korean_source_css_exact", int(exact.sum()))
    return pd.DataFrame(rows)


def resolve_institution(
    subjects: pd.DataFrame, supertable: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Resolve Korean institution by VCF ID, then by Sample_ID fallback."""
    out = subjects.copy()
    cols = ["VCF_ID_2025", "Sample_ID", "Collection_Institution"]
    source = supertable[cols].copy()

    def unique_map(key: str) -> pd.Series:
        subset = source.dropna(subset=[key, "Collection_Institution"])[[key, "Collection_Institution"]]
        conflicts = subset.groupby(key)["Collection_Institution"].nunique()
        if (conflicts > 1).any():
            raise ValueError(f"institution conflict for {key}")
        return subset.drop_duplicates(key).set_index(key)["Collection_Institution"]

    by_vcf = unique_map("VCF_ID_2025")
    by_sample = unique_map("Sample_ID")
    out["institution"] = "SSC"
    is_korean = out["dataset"] == KOREAN
    out.loc[is_korean, "institution"] = out.loc[is_korean, "vcf_iid"].map(by_vcf)
    unresolved = is_korean & out["institution"].isna()
    out.loc[unresolved, "institution"] = out.loc[unresolved, "Sample_ID"].map(by_sample)
    unresolved_n = int((is_korean & out["institution"].isna()).sum())
    if unresolved_n:
        raise ValueError(f"institution unresolved for {unresolved_n} Korean subjects")
    non_snubh = is_korean & (out["institution"] != SNUBH)
    if non_snubh.any():
        raise ValueError(f"current master contains {int(non_snubh.sum())} non-SNUBH Korean subjects")
    out["hospital_included"] = (~is_korean) | (out["institution"] == SNUBH)
    out["exclusion_reason"] = pd.NA
    out.loc[is_korean & ~out["hospital_included"], "exclusion_reason"] = "IRB_scope_excluded"

    rows: list[dict] = []
    _metric(rows, "korean_institution_resolved", int(is_korean.sum()) - unresolved_n)
    _metric(rows, "korean_institution_unresolved", unresolved_n)
    _metric(rows, "korean_snubh_included", int((is_korean & out["hospital_included"]).sum()))
    _metric(rows, "current_non_snubh_included", int(non_snubh.sum()))
    return out, pd.DataFrame(rows)


def add_eligibility(master: pd.DataFrame) -> pd.DataFrame:
    out = master.copy()
    out["lca_eligible"] = out["class_label"].notna()
    out["pgs_eligible"] = out["updated_pgs_mapped"].fillna(False).astype(bool)
    common = ["age", "sex", "dataset", "ados_css_final", "vineland"]
    srs = ["age", "sex", "dataset", "srs", "vineland"]
    pcs = common + [f"PC{i}" for i in range(1, 6)]
    out["eligible_primary_base"] = out[common].notna().all(axis=1)
    out["eligible_srs_base"] = out[srs].notna().all(axis=1)
    out["eligible_primary_pc_base"] = out[pcs].notna().all(axis=1)
    return out


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--subject-table", type=Path, required=True)
    parser.add_argument("--ssc-m13", type=Path, required=True)
    parser.add_argument("--ssc-m4", type=Path, required=True)
    parser.add_argument("--ssc-pheno", type=Path, required=True)
    parser.add_argument("--korean-pheno", type=Path, required=True)
    parser.add_argument("--korean-clinical", type=Path, required=True)
    parser.add_argument("--supertable", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    inputs = {
        "frozen_subject_table": args.subject_table,
        "ssc_css_m1to3": args.ssc_m13,
        "ssc_css_module4": args.ssc_m4,
        "ssc_phenotype": args.ssc_pheno,
        "korean_phenotype": args.korean_pheno,
        "korean_clinical_css": args.korean_clinical,
        "korean_supertable": args.supertable,
    }

    subjects = pd.read_csv(args.subject_table)
    m13 = pd.read_excel(args.ssc_m13)
    m4 = pd.read_excel(args.ssc_m4)
    ssc_pheno = pd.read_csv(args.ssc_pheno, sep="\t", low_memory=False)
    korean_pheno = pd.read_csv(args.korean_pheno, sep="\t", low_memory=False)
    supertable = pd.read_excel(args.supertable)
    korean_raw = pd.read_excel(args.korean_clinical, sheet_name="Sheet1", header=1)
    korean_clinical = normalize_korean_clinical(korean_raw)

    master, css_report = combine_css(subjects, m13, m4, ssc_pheno, korean_pheno)
    korean_report = validate_korean_clinical(master, korean_clinical)
    master, hospital_report = resolve_institution(master, supertable)
    master = add_eligibility(master)

    validation = pd.concat([css_report, korean_report, hospital_report], ignore_index=True)
    master.to_csv(args.output_dir / "analysis_master_subject_level_v1.tsv", sep="\t", index=False)
    validation.to_csv(args.output_dir / "css_validation_report_v1.tsv", sep="\t", index=False)

    exclusions = master.loc[
        ~master["hospital_included"],
        ["individual", "dataset", "institution", "exclusion_reason"],
    ]
    exclusions.to_csv(args.output_dir / "sample_exclusion_ledger_v1.tsv", sep="\t", index=False)

    coverage_rows = []
    for dataset, group in list(master.groupby("dataset", observed=True)) + [("Combined", master)]:
        for variable in ["ados_css_final", "srs", "vineland", "age"]:
            coverage_rows.append(
                {
                    "dataset": dataset,
                    "variable": variable,
                    "n_total": len(group),
                    "n_nonmissing": int(group[variable].notna().sum()),
                    "n_missing": int(group[variable].isna().sum()),
                }
            )
    pd.DataFrame(coverage_rows).to_csv(
        args.output_dir / "variable_coverage_v1.tsv", sep="\t", index=False
    )

    registry = []
    for role, path in inputs.items():
        registry.append(
            {
                "role": role,
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "change_policy": "read_only",
            }
        )
    pd.DataFrame(registry).to_csv(
        args.output_dir / "input_source_manifest_v1.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
