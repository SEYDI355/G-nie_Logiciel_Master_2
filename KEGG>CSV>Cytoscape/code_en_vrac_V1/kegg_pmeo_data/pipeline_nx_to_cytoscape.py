import os
import glob
import json
import pandas as pd
import networkx as nx
from networkx.readwrite import json_graph
from Bio.KEGG.KGML import KGML_parser

# --- CONFIGURATION ---
DATA_DIR = "kegg_pmeo_data"
KGML_DIR = os.path.join(DATA_DIR, "kgml")
QUANT_CSV = os.path.join(DATA_DIR, "metabolites_quant.csv")  # Optionnel (vos mesures GC-MS / dosages)
APP_DIR = "cytoscape_app"
os.makedirs(APP_DIR, exist_ok=True)

# Cofacteurs / currency metabolites (à filtrer pour les calculs topologiques et la vue monopartite)
CURRENCY_COMPOUNDS = {
    "C00001", "C00002", "C00003", "C00004", "C00005", "C00006",
    "C00007", "C00008", "C00009", "C00010", "C00011", "C00013",
    "C00014", "C00020", "C00080",
}


def load_lookup_tables(data_dir):
    """Charge les tables KEGG extraites (noms chimiques, gènes -> KO/EC, titres de voies)."""
    cpd_names = {}
    cpd_path = os.path.join(data_dir, "metabolites_uniques.csv")
    if os.path.exists(cpd_path):
        df = pd.read_csv(cpd_path)
        for _, r in df.iterrows():
            cid = str(r["compound_id"]).replace("cpd:", "").strip()
            cpd_names[cid] = str(r["compound_name"])

    gene_to_ko = {}
    ko_path = os.path.join(data_dir, "gene_to_ko.csv")
    if os.path.exists(ko_path):
        df = pd.read_csv(ko_path)
        for g, grp in df.groupby("gene_id"):
            gene_to_ko[str(g).strip()] = sorted({str(x).replace("ko:", "").strip() for x in grp["ko_id"]})

    gene_to_ec = {}
    ec_path = os.path.join(data_dir, "gene_to_ec.csv")
    if os.path.exists(ec_path):
        df = pd.read_csv(ec_path)
        for g, grp in df.groupby("gene_id"):
            gene_to_ec[str(g).strip()] = sorted({str(x).replace("ec:", "").strip() for x in grp["ec_id"]})

    pw_titles = {}
    pw_path = os.path.join(data_dir, "pathways.csv")
    if os.path.exists(pw_path):
        df = pd.read_csv(pw_path)
        for _, r in df.iterrows():
            pid = str(r["pathway_id"]).replace("path:", "").strip()
            pw_titles[pid] = str(r["pathway_name"])

    # Optionnel : mesures quantitatives de métabolites (GC-MS / dosages)
    quant_data = {}
    if os.path.exists(QUANT_CSV):
        df_q = pd.read_csv(QUANT_CSV)
        for _, r in df_q.iterrows():
            cid = str(r["compound_id"]).replace("cpd:", "").strip()
            quant_data[cid] = {
                "value": None if pd.isna(r.get("value")) else float(r["value"]),
                "log2fc": None if pd.isna(r.get("log2fc")) else float(r["log2fc"]),
                "pvalue": None if pd.isna(r.get("pvalue")) else float(r["pvalue"]),
                "unit": "" if pd.isna(r.get("unit")) else str(r["unit"]),
            }

    return cpd_names, gene_to_ko, gene_to_ec, pw_titles, quant_data


def build_bipartite_nx(kgml_dir, cpd_names, gene_to_ko, gene_to_ec, pw_titles, quant_data):
    """Construit l'objet NetworkX biparti (Métabolites <-> Réactions) avec tous les attributs."""
    G = nx.DiGraph()
    pathway_catalog = {}

    kgml_files = sorted(glob.glob(os.path.join(kgml_dir, "*.kgml")))
    print(f"Lecture de {len(kgml_files)} fichiers KGML...")

    for filepath in kgml_files:
        pw_id = os.path.splitext(os.path.basename(filepath))[0]
        if pw_id.endswith("01100"):  # Exclure la carte globale redondante
            continue

        try:
            with open(filepath, "r", encoding="utf-8") as f:
                pathway = KGML_parser.read(f)
        except Exception as e:
            print(f"  [Ignoré] {pw_id}: {e}")
            continue

        if not pathway.reactions:
            continue

        pw_title = pw_titles.get(pw_id, pathway.title or pw_id)
        pathway_catalog[pw_id] = pw_title

        entries_by_id = {}
        for eid, entry in pathway.entries.items():
            g = entry.graphics[0] if entry.graphics else None
            entries_by_id[eid] = {
                "name": entry.name,
                "type": entry.type,
                "label": g.name if g else "",
                "x": float(g.x) if g and g.x is not None else None,
                "y": float(g.y) if g and g.y is not None else None,
                "bgcolor": getattr(g, "bgcolor", ""),
            }

        for rxn in pathway.reactions:
            enz = entries_by_id.get(rxn.id, {})
            rxn_ids = [r.replace("rn:", "").strip() for r in rxn.name.split() if r.startswith("rn:")] or [rxn.name.strip()]
            is_rev = (rxn.type == "reversible")
            genes = [g.strip() for g in enz.get("name", "").split() if g.strip() and ":" in g]
            is_pmeo = (enz.get("type") == "gene") or (enz.get("bgcolor") == "#BFFFBF")

            kos = sorted({k for g in genes for k in gene_to_ko.get(g, [])})
            ecs = sorted({e for g in genes for e in gene_to_ec.get(g, [])})

            substrates, products = [], []
            for sub in rxn.substrates:
                se = entries_by_id.get(sub.id, {})
                for cid in se.get("name", "").split():
                    c = cid.replace("cpd:", "").strip()
                    if c:
                        substrates.append((c, se))

            for prod in rxn.products:
                pe = entries_by_id.get(prod.id, {})
                for cid in pe.get("name", "").split():
                    c = cid.replace("cpd:", "").strip()
                    if c:
                        products.append((c, pe))

            # 1. Ajouter / mettre à jour les nœuds Métabolites
            for cid, ce in substrates + products:
                if cid not in G:
                    full_name = cpd_names.get(cid, ce.get("label", cid))
                    q = quant_data.get(cid, {})
                    G.add_node(
                        cid,
                        node_type="metabolite",
                        label=full_name.split(";")[0][:24],
                        full_name=full_name,
                        is_currency=(cid in CURRENCY_COMPOUNDS),
                        evidence_status="measured" if cid in quant_data else "inferred_kegg",
                        quant_value=q.get("value"),
                        log2fc=q.get("log2fc"),
                        pvalue=q.get("pvalue"),
                        unit=q.get("unit", ""),
                        pathways=set(),
                        pos_by_pathway={},
                    )
                G.nodes[cid]["pathways"].add(pw_id)
                if ce.get("x") is not None and ce.get("y") is not None:
                    G.nodes[cid]["pos_by_pathway"][pw_id] = {"x": ce["x"], "y": ce["y"]}

            # 2. Ajouter / mettre à jour les nœuds Réactions et les arêtes bipartites
            for rid in rxn_ids:
                if rid not in G:
                    G.add_node(
                        rid,
                        node_type="reaction",
                        label=rid,
                        genes_pmeo=set(),
                        ko_ids=set(),
                        ec_ids=set(),
                        reversible=is_rev,
                        evidence_source="KEGG_pmeo" if is_pmeo else "KEGG_reference",
                        pathways=set(),
                        pos_by_pathway={},
                    )
                node_r = G.nodes[rid]
                node_r["genes_pmeo"].update(genes)
                node_r["ko_ids"].update(kos)
                node_r["ec_ids"].update(ecs)
                node_r["pathways"].add(pw_id)
                if is_pmeo:
                    node_r["evidence_source"] = "KEGG_pmeo"
                if is_rev:
                    node_r["reversible"] = True
                if enz.get("x") is not None and enz.get("y") is not None:
                    node_r["pos_by_pathway"][pw_id] = {"x": enz["x"], "y": enz["y"]}

                # Arêtes Substrat -> Réaction
                for cid, _ in substrates:
                    if G.has_edge(cid, rid):
                        G[cid][rid]["pathways"].add(pw_id)
                        G[cid][rid]["reversible"] = G[cid][rid]["reversible"] or is_rev
                    else:
                        G.add_edge(cid, rid, role="substrate", reversible=is_rev, pathways={pw_id})

                # Arêtes Réaction -> Produit
                for cid, _ in products:
                    if G.has_edge(rid, cid):
                        G[rid][cid]["pathways"].add(pw_id)
                        G[rid][cid]["reversible"] = G[rid][cid]["reversible"] or is_rev
                    else:
                        G.add_edge(rid, cid, role="product", reversible=is_rev, pathways={pw_id})

    # Supprimer d'éventuels nœuds isolés
    G.remove_nodes_from(list(nx.isolates(G)))
    return G, pathway_catalog


def compute_topology_metrics(G):
    """Calcule les métriques NetworkX en excluant les cofacteurs pour éviter les faux hubs."""
    non_currency_nodes = [
        n for n, d in G.nodes(data=True)
        if not (d.get("node_type") == "metabolite" and d.get("is_currency"))
    ]
    G_clean = G.subgraph(non_currency_nodes)

    deg = dict(G_clean.degree())
    in_deg = dict(G_clean.in_degree())
    out_deg = dict(G_clean.out_degree())

    # Centralité d'intermédiarité (échantillonnée si graphe très grand pour rester rapide)
    k_sample = min(250, G_clean.number_of_nodes()) if G_clean.number_of_nodes() > 250 else None
    betw = nx.betweenness_centrality(G_clean, k=k_sample, normalized=True, seed=42)

    for n in G.nodes():
        G.nodes[n]["degree"] = int(deg.get(n, G.degree(n)))
        G.nodes[n]["in_degree"] = int(in_deg.get(n, G.in_degree(n)))
        G.nodes[n]["out_degree"] = int(out_deg.get(n, G.out_degree(n)))
        G.nodes[n]["betweenness"] = round(float(betw.get(n, 0.0)), 6)


def build_monopartite_projection(G_bip):
    """Projette le graphe biparti en graphe monoparti (Métabolites -> Métabolites) sans cofacteurs."""
    G_mono = nx.DiGraph()

    # Copier tous les nœuds métabolites non-cofacteurs
    for n, d in G_bip.nodes(data=True):
        if d.get("node_type") == "metabolite" and not d.get("is_currency"):
            G_mono.add_node(n, **{k: (set(v) if isinstance(v, set) else v) for k, v in d.items()})

    # Relier chaque substrat principal à chaque produit principal via le nœud réaction
    for r, rd in G_bip.nodes(data=True):
        if rd.get("node_type") != "reaction":
            continue
        subs = [u for u in G_bip.predecessors(r) if u in G_mono]
        prods = [v for v in G_bip.successors(r) if v in G_mono]
        is_rev = bool(rd.get("reversible", False))
        pws = set(rd.get("pathways", set()))
        ev = rd.get("evidence_source", "KEGG_reference")
        genes = set(rd.get("genes_pmeo", set()))

        for s in subs:
            for p in prods:
                if s == p:
                    continue
                if G_mono.has_edge(s, p):
                    e = G_mono[s][p]
                    e["reactions"].add(r)
                    e["genes_pmeo"].update(genes)
                    e["pathways"].update(pws)
                    e["reversible"] = e["reversible"] or is_rev
                    if ev == "KEGG_pmeo":
                        e["evidence_source"] = "KEGG_pmeo"
                else:
                    G_mono.add_edge(
                        s, p,
                        reactions={r},
                        genes_pmeo=set(genes),
                        pathways=set(pws),
                        reversible=is_rev,
                        evidence_source=ev,
                    )

    G_mono.remove_nodes_from(list(nx.isolates(G_mono)))
    return G_mono


def sanitize_graph_for_json(G):
    """Convertit les sets en listes triées pour que nx.cytoscape_data() soit 100 % sérialisable en JSON."""
    G_copy = G.copy()
    for n, d in G_copy.nodes(data=True):
        for k, v in list(d.items()):
            if isinstance(v, set):
                d[k] = sorted(v)
        if d.get("node_type") == "reaction":
            d["gene_count"] = len(d.get("genes_pmeo", []))
    for u, v, d in G_copy.edges(data=True):
        for k, val in list(d.items()):
            if isinstance(val, set):
                d[k] = sorted(val)
    return G_copy


def to_cytoscape_elements(G, prefix):
    """Utilise la fonction native nx.cytoscape_data(G) et nettoie les IDs d'arêtes."""
    cy = json_graph.cytoscape_data(sanitize_graph_for_json(G))
    for i, edge in enumerate(cy["elements"]["edges"]):
        edge["data"]["id"] = f"{prefix}_e_{i}"
    return cy["elements"]


# --- EXÉCUTION DU PIPELINE ---
if __name__ == "__main__":
    cpd_names, gene_to_ko, gene_to_ec, pw_titles, quant_data = load_lookup_tables(DATA_DIR)

    # 1. Construction de l'objet NetworkX biparti (source de vérité)
    G_bip = build_bipartite_nx(KGML_DIR, cpd_names, gene_to_ko, gene_to_ec, pw_titles, quant_data)
    G_bip, pathway_catalog = G_bip

    # 2. Enrichissement topologique par NetworkX
    compute_topology_metrics(G_bip)

    # 3. Projection monopartite par NetworkX (Métabolites uniquement, sans cofacteurs)
    G_mono = build_monopartite_projection(G_bip)

    print(f"NetworkX Biparti   : {G_bip.number_of_nodes()} nœuds, {G_bip.number_of_edges()} arêtes")
    print(f"NetworkX Monoparti : {G_mono.number_of_nodes()} métabolites, {G_mono.number_of_edges()} arêtes")

    # 4. Export JSON natif Cytoscape.js via NetworkX
    payload = {
        "meta": {
            "organism": "Paramacrobiotus metropolitanus (pmeo)",
            "bipartite_nodes": G_bip.number_of_nodes(),
            "bipartite_edges": G_bip.number_of_edges(),
            "monopartite_nodes": G_mono.number_of_nodes(),
            "monopartite_edges": G_mono.number_of_edges(),
        },
        "pathways": [{"id": k, "title": v} for k, v in sorted(pathway_catalog.items())],
        "bipartite": to_cytoscape_elements(G_bip, "bip"),
        "monopartite": to_cytoscape_elements(G_mono, "mono"),
    }

    out_json = os.path.join(APP_DIR, "pmeo_cytoscape.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False)

    print(f"JSON Cytoscape.js exporté avec succès : {out_json}")