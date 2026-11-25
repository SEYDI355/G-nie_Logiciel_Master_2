"""
analyse_metabo.py
Lecture, analyse statistique et mapping KEGG des données métabolomiques
(LC-MS / IC-MS) du tardigrade.

Règles (alignées sur la présentation) :
- métabolite supprimé si le contrôle a moins de 3 réplicats mesurés
- comparaison ignorée pour un métabolite si la condition a moins de 3 réplicats mesurés
- log2FC = log2(moyenne condition / moyenne contrôle)
- test de Welch sur les valeurs log2, p-value brute
  (p-value corrigée Benjamini-Hochberg fournie à titre indicatif)
- significatif si p < 0.05 et |log2FC| > 1 ; rouge = hausse, bleu = baisse

Utilisation seule (sans Django) :
    py analyse_metabo.py LC-MS.xlsx IC-MS.xlsx
    py analyse_metabo.py LC-MS.xlsx IC-MS.xlsx --reseau reseau.json

Utilisation dans Django :
    from .analyse_metabo import analyser, en_json
    resultat = analyser([request.FILES["lcms"], request.FILES["icms"]])
    return JsonResponse(en_json(resultat))
"""

import argparse
import json
import os
import re

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import requests
from scipy import stats

# ===================== PARAMÈTRES =====================
CONTROLE = "Control"   # nom de la condition de référence dans l'Excel
MIN_REPLICATS = 3
SEUIL_P = 0.05
SEUIL_LFC = 1.0
CACHE_KEGG = "kegg_compounds.tsv"   # liste KEGG enregistrée localement

COULEURS = {"hausse": "#d62728", "baisse": "#1f77b4", "ns": "#bbbbbb"}

# Corrections manuelles : nom de l'Excel (en minuscules) -> identifiant KEGG.
# Vérifier chaque identifiant sur https://www.kegg.jp avant de l'ajouter.
CORRESPONDANCES = {
    "2,3-dpg": "C01159",   # 2,3-Bisphospho-D-glycerate
}
# ======================================================


# ---------- 1. Lecture du fichier Excel ----------

def en_nombre(serie):
    """Convertit une colonne en nombres : 'N.D.' devient NaN, '9,8' devient 9.8."""
    texte = serie.astype(str).str.strip().str.replace(",", ".", regex=False)
    return pd.to_numeric(texte, errors="coerce")


def lire_fichier(fichier):
    """Lit un Excel au format : titre / ligne 'Compound' / ligne 'File name' / données.

    Accepte un chemin ou un fichier uploadé (Django).
    Renvoie (valeurs, conditions, source) :
      - valeurs : DataFrame (lignes = métabolites, colonnes = échantillons)
      - conditions : {"Control": ["ara_S1", ...], "Bleomycin": [...], ...}
      - source : nom du fichier
    """
    source = os.path.basename(getattr(fichier, "name", str(fichier)))
    brut = pd.read_excel(fichier, header=None, dtype=object)

    premiere_col = brut.iloc[:, 0].astype(str).str.strip()
    lignes_compound = brut.index[premiere_col == "Compound"]
    if len(lignes_compound) == 0:
        raise ValueError(f"{source} : ligne 'Compound' introuvable, format non reconnu")

    ligne = lignes_compound[0]
    groupes = brut.iloc[ligne, 1:]
    debut = ligne + 1

    if str(brut.iloc[debut, 0]).strip() == "File name":
        echantillons = brut.iloc[debut, 1:]
        debut += 1
    else:
     echantillons = groupes

    colonnes_utiles = groupes.notna()
    groupes = groupes[colonnes_utiles].astype(str).str.strip()
    echantillons = echantillons[colonnes_utiles].astype(str).str.strip()

    donnees = brut.iloc[debut:, [0] + list(groupes.index)].copy()
    donnees.columns = ["metabolite"] + list(echantillons)
    donnees["metabolite"] = donnees["metabolite"].astype(str).str.strip()
    donnees = donnees[~donnees["metabolite"].isin(["", "nan", "None"])]
    donnees = donnees.drop_duplicates("metabolite").set_index("metabolite")

    valeurs = donnees.apply(en_nombre)

    # "Control #1" -> "Control"
    conditions = {}
    for groupe, echantillon in zip(groupes, echantillons):
        condition = re.sub(r"\s*#\s*\d+\s*$", "", groupe)
        conditions.setdefault(condition, []).append(echantillon)

    if CONTROLE not in conditions:
        raise ValueError(f"{source} : condition '{CONTROLE}' absente "
           f"(conditions trouvées : {list(conditions)})")
    return valeurs, conditions, source


# ---------- 2. Statistiques ----------

def benjamini_hochberg(pvalues):
    p = pd.Series(pvalues, dtype=float)
    valides = p.dropna().sort_values()
    n = len(valides)
    if n == 0:
        return p
    ajustees = valides * n / np.arange(1, n + 1)
    ajustees = ajustees[::-1].cummin()[::-1].clip(upper=1)
    return ajustees.reindex(p.index)


def analyser_tableau(valeurs, conditions, source):
    """Compare chaque condition au contrôle. Renvoie (résultats, résumé)."""
    controle = conditions[CONTROLE]
    log2 = np.log2(valeurs.where(valeurs > 0))
    n_valides = valeurs.notna().sum(axis=1)

    controle_ok = valeurs[controle].notna().sum(axis=1) >= MIN_REPLICATS
    tables, resume = [], []

    for condition, colonnes in conditions.items():
        if condition == CONTROLE:
            continue

        garde = controle_ok & (valeurs[colonnes].notna().sum(axis=1) >= MIN_REPLICATS)
        resume.append({"source": source, "comparaison": f"{CONTROLE} vs {condition}",
         "total": len(valeurs), "analysables": int(garde.sum())})
        if not garde.any():
            continue

        brut = valeurs[garde]
        lg = log2[garde]

        log2fc = np.log2(brut[colonnes].mean(axis=1) / brut[controle].mean(axis=1))
        test = stats.ttest_ind(lg[colonnes].to_numpy(), lg[controle].to_numpy(),
   axis=1, equal_var=False, nan_policy="omit")
        pval = pd.Series(np.asarray(test.pvalue, dtype=float), index=lg.index)

        table = pd.DataFrame({
     "metabolite": lg.index,
     "source": source,
     "condition": condition,
     "log2fc": log2fc.values,
     "pval": pval.values,
     "padj_bh": benjamini_hochberg(pval).values,
            "n_valides": n_valides[garde].values,
        })

        significatif = (table["pval"] < SEUIL_P) & (table["log2fc"].abs() > SEUIL_LFC)
        table["sens"] = np.where(significatif & (table["log2fc"] > 0), "hausse",
          np.where(significatif & (table["log2fc"] < 0), "baisse", "ns"))
        tables.append(table)

    resultats = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    return resultats, resume


# ---------- 3. Mapping vers les identifiants KEGG ----------

def normaliser(nom):
    return re.sub(r"\s+", " ", str(nom).strip().lower())


def variantes(nom):
    n = normaliser(nom)
    liste = [n, "l-" + n, "d-" + n]
    if n.endswith("ate"):
        liste.append(n[:-3] + "ic acid")
    return liste


def charger_kegg():
    """Dictionnaire {synonyme normalisé: [identifiants KEGG]} (mis en cache)."""
    # On ne réutilise le cache que s'il semble complet (> 100 Ko)
    if os.path.exists(CACHE_KEGG) and os.path.getsize(CACHE_KEGG) > 100_000:
        with open(CACHE_KEGG, encoding="utf-8") as f:
            texte = f.read()
    else:
        reponse = requests.get("https://rest.kegg.jp/list/compound", timeout=60)
        reponse.raise_for_status()
        texte = reponse.text
        with open(CACHE_KEGG, "w", encoding="utf-8") as f:
            f.write(texte)

    synonymes = {}
    for ligne in texte.splitlines():
        if "\t" not in ligne:
            continue
        cid, noms = ligne.split("\t", 1)
        cid = cid.replace("cpd:", "").strip()
        for nom in noms.split(";"):
            cle = normaliser(nom)
            if not cle:
                continue
            ids = synonymes.setdefault(cle, [])
            if cid not in ids:
                ids.append(cid)

    print(f"KEGG : {len(synonymes)} synonymes chargés")
    return synonymes


def mapper_kegg(nom, synonymes):
    """Renvoie (identifiant KEGG ou None, statut)."""
    cle = normaliser(nom)
    if cle in CORRESPONDANCES:
        return CORRESPONDANCES[cle], "manuel"
    if " + " in str(nom):
        return None, "plusieurs composés"

    for i, variante in enumerate(variantes(nom)):
        ids = synonymes.get(variante)
        if ids:
            if len(ids) > 1:
                return ";".join(ids), "ambigu"
            return ids[0], ("exact" if i == 0 else "approché")
    return None, "non trouvé"



# ---------- 4. Visualisation ----------

def volcano(resultats, condition):
    d = resultats[resultats["condition"] == condition]
    fig = go.Figure()

    for sens, libelle in [("ns", "Non significatif"), ("baisse", "Baisse"), ("hausse", "Hausse")]:
        s = d[d["sens"] == sens]
        fig.add_trace(go.Scatter(
            x=s["log2fc"],
     y=-np.log10(s["pval"]),
     mode="markers",
     name=libelle,
     marker=dict(color=COULEURS[sens], size=8, opacity=0.85),
     text=s["metabolite"],
     customdata=s["kegg_id"].fillna("non associé"),
     hovertemplate="<b>%{text}</b><br>KEGG : %{customdata}"
     "<br>log2FC = %{x:.2f}<br>-log10(p) = %{y:.2f}<extra></extra>",
 ))

    fig.add_hline(y=-np.log10(SEUIL_P), line_dash="dash", line_color="grey")
    fig.add_vline(x=SEUIL_LFC, line_dash="dash", line_color="grey")
    fig.add_vline(x=-SEUIL_LFC, line_dash="dash", line_color="grey")
    fig.update_layout(
        title=f"{CONTROLE} vs {condition}",
 xaxis_title="log2 fold change",
 yaxis_title="-log10(p-value, Welch)",
 template="simple_white",
    )
    return fig


def pour_reseau(resultats):
    """{kegg_id: {"nom": ..., condition: {log2fc, pval, sens, couleur}}} pour colorer le graphe."""
    sortie = {}
    for _, ligne in resultats.dropna(subset=["kegg_id"]).iterrows():
        if ";" in ligne["kegg_id"]:
            continue  # identifiant ambigu : pas de coloration
        noeud = sortie.setdefault(ligne["kegg_id"], {"nom": ligne["metabolite"]})
        noeud[ligne["condition"]] = {
     "log2fc": round(float(ligne["log2fc"]), 3),
     "pval": None if pd.isna(ligne["pval"]) else float(ligne["pval"]),
     "sens": ligne["sens"],
     "couleur": COULEURS[ligne["sens"]],
        }
    return sortie


# ---------- 5. Fonction principale (appelée par Django) ----------

def analyser(fichiers):
    """fichiers : liste de chemins ou de fichiers uploadés (LC-MS, IC-MS)."""
    synonymes = charger_kegg()
    tables, resume = [], []

    for fichier in fichiers:
        valeurs, conditions, source = lire_fichier(fichier)
        table, r = analyser_tableau(valeurs, conditions, source)
        tables.append(table)
        resume.extend(r)

    resultats = pd.concat(tables, ignore_index=True)

    # Métabolite présent dans les deux fichiers : on garde la source la plus complète
    meilleure_source = (resultats.sort_values("n_valides", ascending=False)
          .drop_duplicates("metabolite")[["metabolite", "source"]])
    resultats = resultats.merge(meilleure_source, on=["metabolite", "source"])

    mapping = {m: mapper_kegg(m, synonymes) for m in resultats["metabolite"].unique()}
    resultats["kegg_id"] = resultats["metabolite"].map(lambda m: mapping[m][0])
    resultats["mapping"] = resultats["metabolite"].map(lambda m: mapping[m][1])

    figures = {c: volcano(resultats, c) for c in resultats["condition"].unique()}

    return {
        "resultats": resultats,
        "resume": pd.DataFrame(resume),
        "figures": figures,
        "reseau": pour_reseau(resultats),
    }


def en_json(resultat):
    """Convertit le résultat en dictionnaire sérialisable (JsonResponse)."""
    table = resultat["resultats"].astype(object).where(resultat["resultats"].notna(), None)
    return {
 "resultats": table.to_dict("records"),
        "resume": resultat["resume"].to_dict("records"),
        "figures": {c: json.loads(fig.to_json()) for c, fig in resultat["figures"].items()},
        "reseau": resultat["reseau"],
    }


# ---------- 6. Vérification du recouvrement avec le réseau ----------

def ids_du_reseau(chemin):
    """Accepte un JSON Cytoscape ({"nodes": [{"data": {"id"}}]}) ou NetworkX node_link."""
    with open(chemin, encoding="utf-8") as f:
        graphe = json.load(f)
    noeuds = graphe.get("nodes") or graphe.get("elements", {}).get("nodes", [])
    ids = set()
    for noeud in noeuds:
        identifiant = noeud.get("data", noeud).get("id")
        if identifiant:
            ids.add(str(identifiant).replace("cpd:", ""))
    return ids


# ---------- 7. Lancement en ligne de commande ----------

def main():
    parser = argparse.ArgumentParser(description="Analyse métabolomique du tardigrade")
    parser.add_argument("fichiers", nargs="+", help="fichiers Excel (LC-MS, IC-MS)")
    parser.add_argument("--reseau", help="JSON du réseau pour tester le recouvrement")
    args = parser.parse_args()

    resultat = analyser(args.fichiers)
    resultats = resultat["resultats"]

    print("\n=== Métabolites analysables ===")
    print(resultat["resume"].to_string(index=False))

    resultats.to_csv("resultats.csv", index=False, encoding="utf-8-sig")

    a_verifier = resultats.loc[~resultats["mapping"].isin(["exact", "manuel"]),
          ["metabolite", "kegg_id", "mapping"]].drop_duplicates()
    a_verifier.to_csv("a_verifier.csv", index=False, encoding="utf-8-sig")

    associes = resultats.loc[resultats["mapping"].isin(["exact", "manuel"]), "metabolite"].nunique()
    print(f"\nMapping KEGG : {associes} métabolites associés, "
   f"{len(a_verifier)} à vérifier (voir a_verifier.csv)")

    for condition, fig in resultat["figures"].items():
        nom = f"volcano_{condition.replace(' ', '_')}.html"
        fig.write_html(nom)
        n_sig = (resultats[(resultats["condition"] == condition)]["sens"] != "ns").sum()
        print(f"{nom} : {n_sig} métabolites significatifs")

    if args.reseau:
        ids = ids_du_reseau(args.reseau)
        communs = set(resultat["reseau"]) & ids
        print(f"\nRéseau : {len(communs)} métabolites mesurés présents "
       f"sur {len(ids)} nœuds du graphe")

    print("\nFichiers écrits : resultats.csv, a_verifier.csv, volcano_*.html")


if __name__ == "__main__":
    main()