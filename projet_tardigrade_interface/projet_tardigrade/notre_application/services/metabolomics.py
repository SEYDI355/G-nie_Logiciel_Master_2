"""Service d'analyse statistique métabolomique et de mapping vers les IDs KEGG."""

from __future__ import annotations

import math
import os
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from django.db import transaction
from scipy import stats

from ..models import Dataset, MetaboliteMeasurement

MIN_REPLICATS = 3
SEUIL_P = 0.05
SEUIL_LFC = 1.0

KEGG_ID_RE = re.compile(r"^(?:cpd:)?(C\d{5})$", re.IGNORECASE)

# Dictionnaire local de synonymes normalisés -> identifiant KEGG (fonctionne hors-ligne).
BUILTIN_KEGG_SYNONYMS: dict[str, str] = {
    "2,3-dpg": "C01159",
    "trehalose": "C01083",
    "alpha,alpha-trehalose": "C01083",
    "trehalose-6-phosphate": "C00689",
    "trehalose 6-phosphate": "C00689",
    "udp-glucose": "C00029",
    "glucose": "C00031",
    "d-glucose": "C00031",
    "sucrose": "C00089",
    "glucose-6-phosphate": "C00092",
    "d-glucose 6-phosphate": "C00092",
    "glucose-1-phosphate": "C00103",
    "fructose-6-phosphate": "C00085",
    "d-fructose 6-phosphate": "C00085",
    "maltose": "C00208",
    "glycerol": "C00116",
    "pyruvate": "C00022",
    "lactate": "C00186",
    "phosphoenolpyruvate": "C00074",
    "3-phosphoglycerate": "C00197",
    "3-phospho-d-glycerate": "C00197",
    "fructose-1,6-bisphosphate": "C00354",
    "d-fructose 1,6-bisphosphate": "C00354",
    "dihydroxyacetone phosphate": "C00111",
    "glycerone phosphate": "C00111",
    "citrate": "C00158",
    "succinate": "C00042",
    "fumarate": "C00122",
    "malate": "C00149",
    "alpha-ketoglutarate": "C00026",
    "2-oxoglutarate": "C00026",
    "oxaloacetate": "C00036",
    "ribose-5-phosphate": "C00117",
    "d-ribose 5-phosphate": "C00117",
    "sedoheptulose-7-phosphate": "C05382",
    "sedoheptulose 7-phosphate": "C05382",
    "6-phosphogluconate": "C00345",
    "6-phospho-d-gluconate": "C00345",
    "erythrose-4-phosphate": "C00279",
    "d-erythrose 4-phosphate": "C00279",
    "alanine": "C00041",
    "glycine": "C00037",
    "serine": "C00065",
    "proline": "C00148",
    "glutamate": "C00025",
    "glutamine": "C00064",
    "aspartate": "C00049",
    "leucine": "C00123",
    "valine": "C00183",
    "lysine": "C00047",
    "tryptophan": "C00078",
    "kynurenine": "C00328",
    "atp": "C00002",
    "adp": "C00008",
    "amp": "C00020",
    "gtp": "C00044",
    "utp": "C00075",
    "adenosine": "C00212",
    "inosine": "C00294",
    "uridine": "C00299",
    "glutathione": "C00051",
    "gsh": "C00051",
    "glutathione disulfide": "C00127",
    "gssg": "C00127",
    "ascorbate": "C00072",
    "nad+": "C00003",
    "nadh": "C00004",
    "nadp+": "C00006",
    "nadph": "C00005",
    "cysteine": "C00097",
    "choline": "C00114",
    "phosphocholine": "C00588",
    "choline phosphate": "C00588",
    "palmitate": "C00249",
    "oleate": "C00712",
    "glycerol-3-phosphate": "C00093",
    "sn-glycerol 3-phosphate": "C00093",
}

CONTROL_CANDIDATES = ("control", "ctrl", "hydraté", "hydrate", "hydrated", "témoin", "temoin")


def _clean_float(val: object) -> float | None:
    """Convertit un float numpy/pandas en float Python valide JSON (None si NaN/Inf)."""
    if val is None:
        return None
    try:
        f = float(val)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def to_numeric_series(series: pd.Series) -> pd.Series:
    """Convertit une colonne en nombres : 'N.D.' -> NaN, '9,8' -> 9.8."""
    text = series.astype(str).str.strip().str.replace(",", ".", regex=False)
    return pd.to_numeric(text, errors="coerce")


def normalise_name(name: object) -> str:
    return re.sub(r"\s+", " ", str(name).strip().lower())


def name_variants(name: object) -> list[str]:
    n = normalise_name(name)
    variants = [n, f"l-{n}", f"d-{n}"]
    if n.endswith("ate"):
        variants.append(n[:-3] + "ic acid")
    return variants


def load_optional_kegg_cache(cache_paths: Iterable[Path] | None = None) -> dict[str, list[str]]:
    """Charge un éventuel cache local kegg_compounds.tsv s'il est présent sur le disque."""
    if cache_paths is None:
        base = Path(__file__).resolve().parents[4]
        cache_paths = [
            Path(__file__).resolve().parent.parent / "data" / "kegg_compounds.tsv",
            base / "kegg_compounds.tsv",
        ]

    synonyms: dict[str, list[str]] = {}
    for cid_key, cid in BUILTIN_KEGG_SYNONYMS.items():
        synonyms[cid_key] = [cid]

    for path in cache_paths:
        if path.exists() and path.stat().st_size > 100_000:
            with open(path, encoding="utf-8") as fh:
                for line in fh:
                    if "\t" not in line:
                        continue
                    cid_raw, names = line.split("\t", 1)
                    cid = cid_raw.replace("cpd:", "").strip()
                    for alias in names.split(";"):
                        key = normalise_name(alias)
                        if not key:
                            continue
                        ids = synonyms.setdefault(key, [])
                        if cid not in ids:
                            ids.append(cid)
            break

    return synonyms


def map_to_kegg_id(name: object, synonyms: dict[str, list[str]], explicit_id: object = None) -> tuple[str, str]:
    """Associe un nom ou identifiant de métabolite vers un ID KEGG Cxxxxx et son statut."""
    for candidate in (explicit_id, name):
        if candidate is None or (isinstance(candidate, float) and math.isnan(candidate)):
            continue
        match = KEGG_ID_RE.match(str(candidate).strip())
        if match:
            return match.group(1).upper(), "direct"

    clean_key = normalise_name(name)
    if not clean_key:
        return "", "non trouvé"
    if clean_key in BUILTIN_KEGG_SYNONYMS:
        return BUILTIN_KEGG_SYNONYMS[clean_key], "manuel"
    if " + " in str(name):
        return "", "plusieurs composés"

    for i, variant in enumerate(name_variants(name)):
        ids = synonyms.get(variant)
        if ids:
            if len(ids) > 1:
                return ";".join(ids), "ambigu"
            return ids[0], ("exact" if i == 0 else "approché")

    return "", "non trouvé"


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    p = pd.Series(pvalues, dtype=float)
    valid = p.dropna().sort_values()
    n = len(valid)
    if n == 0:
        return p
    adjusted = valid * n / np.arange(1, n + 1)
    adjusted = adjusted[::-1].cummin()[::-1].clip(upper=1)
    return adjusted.reindex(p.index)


def _read_tabular(path: str) -> pd.DataFrame:
    """Lit un CSV/TSV en détectant automatiquement le séparateur."""
    return pd.read_csv(path, sep=None, engine="python")


def read_dataset_tables(matrix_path: str, metadata_path: str) -> tuple[pd.DataFrame, dict[str, list[str]], dict[str, str], dict[str, str]]:
    """Lit la matrice et les métadonnées d'un Dataset.

    Retourne (valeurs, conditions, pathways_par_metabolite, kegg_ids_explicites).
    """
    df_meta = _read_tabular(metadata_path)
    meta_cols = {c.strip().lower(): c for c in df_meta.columns}
    if "sample" not in meta_cols or "condition" not in meta_cols:
        raise ValueError("Le fichier de métadonnées doit contenir les colonnes 'sample' et 'condition'.")

    sample_col = meta_cols["sample"]
    cond_col = meta_cols["condition"]

    conditions: dict[str, list[str]] = {}
    for _, row in df_meta.iterrows():
        sample = str(row[sample_col]).strip()
        cond = str(row[cond_col]).strip()
        if sample and cond and sample.lower() != "nan" and cond.lower() != "nan":
            conditions.setdefault(cond, []).append(sample)

    if len(conditions) < 2:
        raise ValueError("Au moins deux conditions distinctes sont requises dans les métadonnées.")

    df_mat = _read_tabular(matrix_path)
    mat_cols = {c.strip().lower(): c for c in df_mat.columns}
    met_col = (
        mat_cols.get("metabolite")
        or mat_cols.get("compound")
        or mat_cols.get("name")
        or mat_cols.get("compound_id")
        or mat_cols.get("kegg_id")
        or df_mat.columns[0]
    )
    pw_col = mat_cols.get("pathway")
    kegg_col = mat_cols.get("kegg_id") or mat_cols.get("compound_id")
    if kegg_col == met_col:
        kegg_col = None

    df_mat[met_col] = df_mat[met_col].astype(str).str.strip()
    df_mat = df_mat[~df_mat[met_col].isin(["", "nan", "None"])].drop_duplicates(met_col)

    all_samples = [s for cols in conditions.values() for s in cols]
    missing_samples = [s for s in all_samples if s not in df_mat.columns]
    if missing_samples:
        raise ValueError(f"Échantillons absents de la matrice : {', '.join(missing_samples[:5])}")

    # Filtrer chaque condition aux colonnes réellement présentes
    conditions = {
        cond: [s for s in cols if s in df_mat.columns]
        for cond, cols in conditions.items()
    }

    indexed = df_mat.set_index(met_col)
    pathways = (
        indexed[pw_col].fillna("").astype(str).str.strip().to_dict()
        if pw_col else {}
    )
    explicit_ids = (
        indexed[kegg_col].fillna("").astype(str).str.strip().to_dict()
        if kegg_col else {}
    )
    values = indexed[all_samples].apply(to_numeric_series)
    return values, conditions, pathways, explicit_ids


def pick_control_condition(conditions: dict[str, list[str]]) -> str:
    """Choisit la condition de référence (Control / hydraté en priorité)."""
    by_lower = {c.lower(): c for c in conditions}
    for candidate in CONTROL_CANDIDATES:
        if candidate in by_lower:
            return by_lower[candidate]
    return next(iter(conditions))


def compute_dataset_measurements(
    dataset: Dataset,
    min_replicates: int = MIN_REPLICATS,
) -> list[MetaboliteMeasurement]:
    """Calcule les statistiques + mapping KEGG et retourne les objets non encore sauvegardés."""
    values, conditions, pathways, explicit_ids = read_dataset_tables(
        dataset.matrix_file.path,
        dataset.metadata_file.path,
    )
    control_name = pick_control_condition(conditions)
    control_cols = conditions[control_name]
    synonyms = load_optional_kegg_cache()

    log2_vals = np.log2(values.where(values > 0))
    ctrl_ok = values[control_cols].notna().sum(axis=1) >= min_replicates

    rows: list[MetaboliteMeasurement] = []
    for cond_name, cond_cols in conditions.items():
        if cond_name == control_name:
            continue

        keep = ctrl_ok & (values[cond_cols].notna().sum(axis=1) >= min_replicates)
        if not keep.any():
            continue

        raw_sub = values[keep]
        lg_sub = log2_vals[keep]

        ctrl_means = raw_sub[control_cols].mean(axis=1)
        cond_means = raw_sub[cond_cols].mean(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            log2fc_series = np.log2(cond_means / ctrl_means)

        test = stats.ttest_ind(
            lg_sub[cond_cols].to_numpy(),
            lg_sub[control_cols].to_numpy(),
            axis=1,
            equal_var=False,
            nan_policy="omit",
        )
        pvals = pd.Series(np.asarray(test.pvalue, dtype=float), index=lg_sub.index)
        padj = benjamini_hochberg(pvals)

        for met_name in lg_sub.index:
            cid, status = map_to_kegg_id(met_name, synonyms, explicit_ids.get(met_name))
            lfc = _clean_float(log2fc_series.loc[met_name])
            pv = _clean_float(pvals.loc[met_name])
            bh = _clean_float(padj.loc[met_name])

            if pv is not None and lfc is not None and pv < SEUIL_P and abs(lfc) > SEUIL_LFC:
                sig = "hausse" if lfc > 0 else "baisse"
            else:
                sig = "ns"

            rows.append(
                MetaboliteMeasurement(
                    dataset=dataset,
                    metabolite_name=str(met_name),
                    pathway_label=pathways.get(met_name, ""),
                    compound_id=cid,
                    mapping_status=status,
                    condition=cond_name,
                    control_condition=control_name,
                    value=_clean_float(cond_means.loc[met_name]),
                    control_value=_clean_float(ctrl_means.loc[met_name]),
                    log2fc=lfc,
                    pvalue=pv,
                    padj_bh=bh,
                    significance=sig,
                )
            )

    return rows


@transaction.atomic
def analyze_and_store_dataset(dataset: Dataset, min_replicates: int = MIN_REPLICATS) -> list[MetaboliteMeasurement]:
    """Calcule et enregistre en base les mesures d'un Dataset."""
    measurements = compute_dataset_measurements(dataset, min_replicates=min_replicates)
    dataset.measurements.all().delete()
    if measurements:
        MetaboliteMeasurement.objects.bulk_create(measurements)
    return measurements
