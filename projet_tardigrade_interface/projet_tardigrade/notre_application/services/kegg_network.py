"""Service d'alignement des mesures d'un Dataset sur le réseau Cytoscape KEGG."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from ..models import Dataset
from .metabolomics import analyze_and_store_dataset

CYTOSCAPE_JSON_PATH = (
    Path(__file__).resolve().parent.parent
    / "static"
    / "notre_application"
    / "data"
    / "pmeo_cytoscape.json"
)


@lru_cache(maxsize=1)
def get_kegg_metabolite_ids() -> frozenset[str]:
    """Charge une fois les IDs de métabolites présents dans pmeo_cytoscape.json."""
    if not CYTOSCAPE_JSON_PATH.exists():
        return frozenset()
    with open(CYTOSCAPE_JSON_PATH, encoding="utf-8") as fh:
        payload = json.load(fh)
    nodes = payload.get("bipartite", {}).get("nodes", [])
    return frozenset(
        str(n["data"]["id"])
        for n in nodes
        if n.get("data", {}).get("node_type") == "metabolite" and n.get("data", {}).get("id")
    )


def build_dataset_network_payload(dataset: Dataset, condition: str | None = None) -> dict:
    """Retourne les mesures KEGG prêtes à être injectées dans Cytoscape.js."""
    qs = dataset.measurements.all()
    if not qs.exists():
        analyze_and_store_dataset(dataset)
        qs = dataset.measurements.all()

    conditions = list(qs.values_list("condition", flat=True).distinct())
    selected_condition = condition if condition in conditions else (conditions[0] if conditions else None)
    if selected_condition:
        qs = qs.filter(condition=selected_condition)

    graph_ids = get_kegg_metabolite_ids()
    measurements: dict[str, dict] = {}
    unmapped_examples: list[str] = []
    mapped_count = 0
    in_graph_count = 0
    significant_count = 0
    total_rows = 0
    control_condition = ""

    for m in qs:
        total_rows += 1
        control_condition = control_condition or m.control_condition
        if m.significance in ("hausse", "baisse"):
            significant_count += 1

        # Ignorer les IDs non résolus ou ambigus (contenant ';') pour la coloration du graphe
        if not m.compound_id or ";" in m.compound_id:
            if len(unmapped_examples) < 10:
                unmapped_examples.append(m.metabolite_name)
            continue

        mapped_count += 1
        in_graph = m.compound_id in graph_ids
        if in_graph:
            in_graph_count += 1

        # Si plusieurs lignes pointent vers le même ID KEGG, garder la plus significative (|log2fc| max)
        existing = measurements.get(m.compound_id)
        if existing is not None:
            old_fc = abs(existing["log2fc"] or 0.0)
            new_fc = abs(m.log2fc or 0.0)
            if new_fc <= old_fc:
                continue

        measurements[m.compound_id] = {
            "metabolite_name": m.metabolite_name,
            "pathway_label": m.pathway_label,
            "value": m.value,
            "control_value": m.control_value,
            "log2fc": m.log2fc,
            "pvalue": m.pvalue,
            "padj_bh": m.padj_bh,
            "significance": m.significance,
            "condition": m.condition,
            "control_condition": m.control_condition,
            "mapping_status": m.mapping_status,
            "in_kegg_graph": in_graph,
        }

    return {
        "dataset_id": dataset.pk,
        "dataset_name": dataset.name,
        "summary": {
            "total_rows": total_rows,
            "mapped_count": mapped_count,
            "in_graph_count": in_graph_count,
            "significant_count": significant_count,
            "conditions": conditions,
            "selected_condition": selected_condition,
            "control_condition": control_condition,
            "unmapped_examples": unmapped_examples,
        },
        "measurements": measurements,
    }
