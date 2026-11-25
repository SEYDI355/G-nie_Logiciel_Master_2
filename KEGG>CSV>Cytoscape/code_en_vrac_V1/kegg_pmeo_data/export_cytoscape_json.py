import os
import json
import pandas as pd

GRAPH_DIR = os.path.join("kegg_pmeo_data", "graph_tables")
APP_DIR = "cytoscape_app"
os.makedirs(APP_DIR, exist_ok=True)

df_met = pd.read_csv(os.path.join(GRAPH_DIR, "nodes_metabolites.csv"))
df_rxn = pd.read_csv(os.path.join(GRAPH_DIR, "nodes_reactions.csv"))
df_edges = pd.read_csv(os.path.join(GRAPH_DIR, "edges.csv"))

nodes = []
edges = []
pathway_catalog = {}  # pathway_id -> titre lisible

# 1. Convertir les métabolites en nœuds Cytoscape.js
for _, row in df_met.iterrows():
    pws = [p for p in str(row.get("pathways", "")).split(";") if p and p != "nan"]
    full_name = str(row.get("name", row["node_id"]))
    short_label = full_name.split(";")[0][:24]

    node_obj = {
        "group": "nodes",
        "data": {
            "id": str(row["node_id"]),
            "label": short_label,
            "full_name": full_name,
            "node_type": "metabolite",
            "is_currency": bool(row.get("is_currency", False)),
            "evidence_status": str(row.get("evidence_status", "inferred_kegg")),
            "pathways": pws,
        },
    }
    # Conserver les coordonnées (X, Y) issues du KGML pour le layout 'preset'
    if pd.notna(row.get("x_kegg")) and pd.notna(row.get("y_kegg")):
        node_obj["position"] = {
            "x": float(row["x_kegg"]),
            "y": float(row["y_kegg"]),
        }
    nodes.append(node_obj)

# 2. Convertir les réactions en nœuds Cytoscape.js
for _, row in df_rxn.iterrows():
    pws = [p for p in str(row.get("pathways", "")).split(";") if p and p != "nan"]
    titles = [t for t in str(row.get("pathway_titles", "")).split(";") if t and t != "nan"]
    for pw in pws:
        if pw not in pathway_catalog:
            pathway_catalog[pw] = titles[0] if titles else pw

    genes_raw = str(row.get("genes_pmeo", ""))
    genes = "" if genes_raw == "nan" else genes_raw

    node_obj = {
        "group": "nodes",
        "data": {
            "id": str(row["node_id"]),
            "label": str(row["node_id"]),
            "node_type": "reaction",
            "genes_pmeo": genes,
            "reversible": bool(row.get("reversible", False)),
            "evidence_source": str(row.get("evidence_source", "KEGG_reference")),
            "pathways": pws,
        },
    }
    if pd.notna(row.get("x_kegg")) and pd.notna(row.get("y_kegg")):
        node_obj["position"] = {
            "x": float(row["x_kegg"]),
            "y": float(row["y_kegg"]),
        }
    nodes.append(node_obj)

# 3. Convertir les arêtes bipartites
for idx, row in df_edges.iterrows():
    pws = [p for p in str(row.get("pathways", "")).split(";") if p and p != "nan"]
    edges.append({
        "group": "edges",
        "data": {
            "id": f"e_{idx}",
            "source": str(row["source"]),
            "target": str(row["target"]),
            "role": str(row["role"]),
            "reversible": bool(row.get("reversible", False)),
            "pathways": pws,
        },
    })

payload = {
    "pathways": [
        {"id": k, "title": v}
        for k, v in sorted(pathway_catalog.items())
    ],
    "elements": {
        "nodes": nodes,
        "edges": edges,
    },
}

out_json = os.path.join(APP_DIR, "pmeo_graph.json")
with open(out_json, "w", encoding="utf-8") as f:
    json.dump(payload, f, ensure_ascii=False, indent=2)

print(f"Export réussi : {out_json}")
print(f" - {len(nodes)} nœuds, {len(edges)} arêtes, {len(pathway_catalog)} voies.")