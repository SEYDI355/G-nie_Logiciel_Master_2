import os
import re
import time
from urllib.parse import urlparse, parse_qs
import requests
import pandas as pd

# 1. Extraire automatiquement le code organisme depuis l'URL KEGG
URL_KEGG = "https://www.kegg.jp/kegg-bin/show_pathway?map=pmeo01100"
map_param = parse_qs(urlparse(URL_KEGG).query).get("map", [""])[0]
ORG = re.match(r"([a-zA-Z]+)", map_param).group(1)  # Donne 'pmeo'

BASE_API = "https://rest.kegg.jp"
OUT_DIR = f"kegg_{ORG}_data"
KGML_DIR = os.path.join(OUT_DIR, "kgml")
os.makedirs(KGML_DIR, exist_ok=True)


def kegg_get_text(endpoint):
    """Interroge l'API REST de KEGG avec une courte pause pour respecter le serveur."""
    url = f"{BASE_API}/{endpoint}"
    r = requests.get(url, timeout=30)
    time.sleep(0.35)  # Pause de courtoisie (limite KEGG)
    if r.status_code == 200 and r.text.strip():
        return r.text.strip()
    return ""


def parse_two_cols(text, col1, col2):
    """Convertit une réponse KEGG à 2 colonnes (tabulée) en DataFrame."""
    if not text:
        return pd.DataFrame(columns=[col1, col2])
    rows = [line.split("\t")[:2] for line in text.split("\n") if "\t" in line]
    return pd.DataFrame(rows, columns=[col1, col2])


print(f"Extraction des données KEGG pour l'organisme : {ORG}...\n")

# 2. Informations générales sur l'organisme
info = kegg_get_text(f"info/{ORG}")
with open(os.path.join(OUT_DIR, f"{ORG}_info.txt"), "w", encoding="utf-8") as f:
    f.write(info)
print("1/6 - Infos générales sauvegardées.")

# 3. Liste des voies métaboliques (pathways) et des gènes
df_pathways = parse_two_cols(kegg_get_text(f"list/pathway/{ORG}"), "pathway_id", "pathway_name")
df_pathways.to_csv(os.path.join(OUT_DIR, "pathways.csv"), index=False)

genes_raw = kegg_get_text(f"list/{ORG}")
genes_rows = [line.split("\t") for line in genes_raw.split("\n") if line]
df_genes = pd.DataFrame(genes_rows)
df_genes.to_csv(os.path.join(OUT_DIR, "genes.csv"), index=False, header=False)
print(f"2/6 - {len(df_pathways)} voies et {len(df_genes)} gènes récupérés.")

# 4. Liens Gènes -> Voies, Gènes -> KO (enzymes), Gènes -> EC
df_gene_pw = parse_two_cols(kegg_get_text(f"link/pathway/{ORG}"), "gene_id", "pathway_id")
df_gene_ko = parse_two_cols(kegg_get_text(f"link/ko/{ORG}"), "gene_id", "ko_id")
df_gene_ec = parse_two_cols(kegg_get_text(f"link/ec/{ORG}"), "gene_id", "ec_id")

df_gene_pw.to_csv(os.path.join(OUT_DIR, "gene_to_pathway.csv"), index=False)
df_gene_ko.to_csv(os.path.join(OUT_DIR, "gene_to_ko.csv"), index=False)
df_gene_ec.to_csv(os.path.join(OUT_DIR, "gene_to_ec.csv"), index=False)
print(f"3/6 - Annotations récupérées ({df_gene_ko['ko_id'].nunique()} KO uniques).")

# 5. Reconstruire le pont KO -> Réactions (Rxxxxx) -> Métabolites (Cxxxxx)
df_ko_rn = parse_two_cols(kegg_get_text("link/reaction/ko"), "ko_id", "reaction_id")
df_rn_cpd = parse_two_cols(kegg_get_text("link/compound/reaction"), "reaction_id", "compound_id")
df_cpd_names = parse_two_cols(kegg_get_text("list/compound"), "compound_id", "compound_name")

# Filtrer uniquement pour les KO présents chez pmeo
df_pmeo_rxn = df_gene_ko.merge(df_ko_rn, on="ko_id", how="inner")
df_pmeo_rxn.to_csv(os.path.join(OUT_DIR, "gene_ko_reaction.csv"), index=False)

# Extraire les métabolites uniques avec leurs noms chimiques
df_pmeo_cpd = (
    df_pmeo_rxn[["reaction_id"]].drop_duplicates()
    .merge(df_rn_cpd, on="reaction_id", how="inner")
    .merge(df_cpd_names, on="compound_id", how="left")
)
df_pmeo_cpd.to_csv(os.path.join(OUT_DIR, "reaction_to_compound.csv"), index=False)

df_unique_cpd = df_pmeo_cpd[["compound_id", "compound_name"]].drop_duplicates()
df_unique_cpd.to_csv(os.path.join(OUT_DIR, "metabolites_uniques.csv"), index=False)
print(f"4/6 - {df_pmeo_rxn['reaction_id'].nunique()} réactions et {len(df_unique_cpd)} métabolites uniques extraits.")

# 6. Télécharger les fichiers KGML (topologie Substrat -> Produit et coordonnées X, Y)
# Remarque : la carte globale 01100 est très lourde ; on télécharge ici toutes les voies de l'organisme.
print("5/6 - Téléchargement des fichiers KGML des voies...")
kgml_count = 0
for pw_id in df_pathways["pathway_id"]:
    clean_id = pw_id.replace("path:", "")
    kgml_text = kegg_get_text(f"get/{clean_id}/kgml")
    if kgml_text:
        with open(os.path.join(KGML_DIR, f"{clean_id}.kgml"), "w", encoding="utf-8") as f:
            f.write(kgml_text)
        kgml_count += 1

print(f"6/6 - Terminé ! {kgml_count} fichiers KGML enregistrés dans '{KGML_DIR}'.")
print(f"Toutes les tables CSV sont disponibles dans le dossier '{OUT_DIR}/'.")