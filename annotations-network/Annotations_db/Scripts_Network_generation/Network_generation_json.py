#authors: Thomas Jaylet - Université Paris Cité - France
#         Karine Audouze - Université Paris Cité - France

#contact: systox@paris-descartes.fr


#AOP-helpFinder is provided without any warranty. But if you have any probleme please feel free to contact us by mail.

#------- WHAT IS AOPHELPFINDER? -------------

#AOP-helpFinder is a tool developed to help AOP development (Jean-Charles Carvaillo: https://github.com/jecarvaill/aop-helpFinder)(Environ Health Perspect. 2019 Apr;127(4):47005).

#It is based on text mining and parsing process on scientific abstracts. AOP-helpFinder identify links between stressors and molecular initiating event, key events and adverse outcomes through abstracts from the PubMed database (https://pubmed.ncbi.nlm.nih.gov/).

#AOP-helpFinder was implemented under the H2020 Human Biomonintoring in Europe (HBM4EU) project, Work Package 13.1.
#HBM4EU has received funding from the European Union’s H2020 research and innovation programme under grant agreement No 733032.

#------- LICENCE ---------------------------

#This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can  use,  modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL

# http://cecill.info/licences/Licence_CeCILL_V2.1-en.txt

#######-----------Parameters-----------######
import pandas as pd
import os 
import argparse
import json



#######-----------Parser-----------######
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--aophf_folder', type=str, help='Le chemin vers le dossier de resultat AOPhf')
parser.add_argument('-m', '--mode', type=str, help='Mode de recherche (ee ou se)')
parser.add_argument('-o', '--output', type=str, help='Chemin pour ecrire l output')
parser.add_argument('-e', '--enrichment', type=str, help='Chemin vers le dossier d enrichissement')
parser.add_argument('-r', '--raw_db', type=str, help='Chemin vers le dossier contenant les donnees brutes des bdd d enrichissement')



args = parser.parse_args()

#Parse variables
folder_AOPhf = args.aophf_folder
mode = args.mode #ee ou se
output_network = args.output
enrichment_folder = args.enrichment
rawdb_folder = args.raw_db



#######-----------Functions-----------######
def network_creation(scoring_df, mode):
    nodes = []
    edges = []
    if mode =="se":
        for i, row in scoring_df.iterrows():
            # Créer des nœuds pour le Stressor et l'Event
            nodes.append({ 'data': { 'id': row['Stressor'].capitalize(), 'type': 'Stressor' } })
            nodes.append({ 'data': { 'id': row['Event'].capitalize(), 'type': 'Event' } })
            # Créer une arête entre le Stressor et l'Event
            edges.append({ 'data': { 'id': f"e{i}", 'source': row['Stressor'].capitalize(), 'target': row['Event'].capitalize(), 'confidence': row['Confidence'], 'link': row['Link'] } })
    elif mode =="ee":       
        for i, row in scoring_df.iterrows():
            # Créer des nœuds pour les paires d'events
            nodes.append({ 'data': { 'id': row['Event 1'].capitalize(), 'type': 'Event' } })
            nodes.append({ 'data': { 'id': row['Event 2'].capitalize(), 'type': 'Event' } })
            # Créer une arête les events
            edges.append({ 'data': { 'id': f"e{i}", 'source': row['Event 1'].capitalize(), 'target': row['Event 2'].capitalize(), 'confidence': row['Confidence'], 'link': row['Link'] } })
    return nodes, edges, i



def process_pairs(df, num_rows):
    """
    Process pairs of columns, keeping the top `num_rows` entries sorted by date or pubdate.
    
    :param df: DataFrame to process
    :param num_rows: Number of top rows to keep for each pair
    :return: Processed DataFrame with top `num_rows` entries for each pair
    """
    # Rename 'date' to 'pubdate' if 'date' column exists
    if 'date' in df.columns:
        df = df.rename(columns={'date': 'pubdate'})
    
    # Determine the date column to use
    if 'pubdate' in df.columns:
        date_column = 'pubdate'
    else:
        raise ValueError("Neither 'date' nor 'pubdate' column found in DataFrame.")
    
    # Determine the pair of columns to use
    if 'stressor' in df.columns and 'event' in df.columns:
        pair_columns = ('stressor', 'event')
    elif 'event_1' in df.columns and 'event_2' in df.columns:
        pair_columns = ('event_1', 'event_2')
    else:
        raise ValueError("Neither 'stressor' and 'event' nor 'event_1' and 'event_2' columns found in DataFrame.")
    
    # Group by the determined pair of columns
    grouped = df.groupby(list(pair_columns))
    
    # Initialize a list to store the processed DataFrames
    frames = []
    
    # Process each group
    for name, group in grouped:
        # Sort by the determined date column
        sorted_group = group.sort_values(by=date_column, ascending=False)
        # Keep the top num_rows elements
        top_group = sorted_group.head(num_rows)
        frames.append(top_group)
    
    # Combine all the processed groups into a single DataFrame
    result_df = pd.concat(frames)
    
    return result_df


def AOPhf_edges_info(df_xlsx_aophf):
    df_xlsx_aophf['title'] = df_xlsx_aophf['title'].str.replace('"', '')
    if mode=="se":
        df_xlsx_aophf['stressor'] = df_xlsx_aophf['stressor'].str.replace('_', ' ')
        df_xlsx_aophf.rename(columns={'stressor': 'source_node', 'event': 'target_node'}, inplace=True)
    if mode=="ee":
        df_xlsx_aophf.rename(columns={'event_1': 'source_node', 'event_2': 'target_node'}, inplace=True)
    df_xlsx_aophf['source_node'] = df_xlsx_aophf['source_node'].str.capitalize()
    df_xlsx_aophf['target_node'] = df_xlsx_aophf['target_node'].str.capitalize()

    return df_xlsx_aophf.to_dict(orient='records')


def aopwiki_enrichment(aopwiki_enrichment_df, nodes, edges, i):
    j=i
    for i, row in aopwiki_enrichment_df.iterrows():
        j=j+i
        # Créer des nœuds pour les paires d'events
        nodes.append({ 'data': { 'id': row['Event_AOP-helpFinder'].capitalize(), 'type': 'Event' } })
        nodes.append({ 'data': { 'id': row['N°AOP'], 'type': 'AOP' } })
        # Créer une arête les events
        edges.append({ 'data': { 'id': f"e{j}", 'source': row['Event_AOP-helpFinder'].capitalize(), 'target': row['N°AOP'], 'AOPWiki_event': row['Event_AOP-wiki'], 'num_event_AOPWiki': row['N°Event'] } })
    return nodes, edges, j


def AOPWiki_nodes_info(aopwiki_enrichment_df, rawtsv_AOPWiki_df, rawjson_AOPWiki_df):
    unique_aop_ids = set(aopwiki_enrichment_df['N°AOP'].astype(str))

    rawtsv_AOPWiki_df['AOP_ID'] = rawtsv_AOPWiki_df['AOP'].apply(lambda x: str(x.split(':')[1]))
    rawtsv_AOPWiki_df['Event_ID'] = rawtsv_AOPWiki_df['Event'].apply(lambda x: x.split(':')[1])


    # Filtrer pour garder uniquement les lignes avec un AOP_ID présent dans unique_aop_ids
    rawtsv_AOPWiki_filtered_df = rawtsv_AOPWiki_df[rawtsv_AOPWiki_df['AOP_ID'].isin(unique_aop_ids)]

    
    # Transformer en dictionnaire {AOP_ID: [events]} puis en liste de dictionnaires
    events_dict = rawtsv_AOPWiki_filtered_df.groupby('AOP_ID').apply(lambda x: x[['Event_ID', 'Type', 'Description']].to_dict('records')).to_dict()
    aop_events_list = [{'AOP_ID': k, 'Events': v} for k, v in events_dict.items()]

    filtered_rawjson_AOPWiki_df = rawjson_AOPWiki_df[rawjson_AOPWiki_df['id'].astype(str).isin(unique_aop_ids)]
    # Convertir en {AOP_ID: Title} puis en liste de dictionnaires
    aop_titles_dict = pd.Series(filtered_rawjson_AOPWiki_df.title.values,index=filtered_rawjson_AOPWiki_df.id.astype(str)).to_dict()
    aop_titles_list = [{'AOP_ID': k, 'Title': v} for k, v in aop_titles_dict.items()]

    return aop_events_list, aop_titles_list


def HPA_nodes_info(HPA_enrichment_df):
    # Sélectionner les colonnes nécessaires et convertir en liste de dictionnaires
    HPA_enrichment_df['Event_AOP-helpFinder'] = HPA_enrichment_df['Event_AOP-helpFinder'].str.capitalize()
    HPA_enrichment_df['RNA tissue specific'] = HPA_enrichment_df['RNA tissue specific'].str.capitalize()
    HPA_enrichment_df['Gene'] = HPA_enrichment_df['Gene'].str.strip()
    hpa_tissue_info = HPA_enrichment_df[['Event_AOP-helpFinder', 'RNA tissue specific', 'Gene']].to_dict(orient='records')
    
    # Retourner le résultat sous forme de dictionnaire
    return hpa_tissue_info

def uniprot_enrichment(Uniprot_enrichment_df, nodes, edges, j):
    k=j
    for i, row in Uniprot_enrichment_df.iterrows():
        k=k+i
        # Créer des nœuds pour les paires d'events
        nodes.append({ 'data': { 'id': row['Event_AOP-helpFinder'].capitalize(), 'type': 'Event' } })
        nodes.append({ 'data': { 'id': row['Biological process'], 'type': 'Uniprot' } })
        # Créer une arête les events
        edges.append({ 'data': { 'id': f"e{k}", 'source': row['Event_AOP-helpFinder'].capitalize(), 'target': row['Biological process'], 'Gene_symbol': row['Gene'] } })
    return nodes, edges, k


def KEGG_enrichment(KEGG_enrichment_df, nodes, edges, k):
    l=k
    for i, row in KEGG_enrichment_df.iterrows():
        l=l+i
        # Créer des nœuds pour les paires d'events
        nodes.append({ 'data': { 'id': row['Event_AOP-helpFinder'].capitalize(), 'type': 'Event' } })
        nodes.append({ 'data': { 'id': row['Pathway'], 'type': 'KEGG' } })
        # Créer une arête les events
        edges.append({ 'data': { 'id': f"e{l}", 'source': row['Event_AOP-helpFinder'].capitalize(), 'target': row['Pathway'], 'Gene_symbol': row['Gene_associated'],'KEGG_ID': row['KEGG_Id'] } })
    return nodes, edges, l



def DisG_enrichment(DisG_enrichment_df, nodes, edges, l):
    m=l
    for i, row in DisG_enrichment_df.iterrows():
        m=m+i
        # Créer des nœuds pour les paires d'events
        nodes.append({ 'data': { 'id': row['Event_AOP-helpFinder'].capitalize(), 'type': 'Event' } })
        nodes.append({ 'data': { 'id': row['diseaseName'], 'type': 'DisGeNET' } })
        # Créer une arête les events
        edges.append({ 'data': { 'id': f"e{m}", 'source': row['Event_AOP-helpFinder'].capitalize(), 'target': row['diseaseName'], 'Gene_symbol': row['geneSymbol'],'Disease_ID': row['diseaseId'] } })
    return nodes, edges, m


def DIS_enrichment(DIS_enrichment_df, nodes, edges, m):
    n=m
    for i, row in DIS_enrichment_df.iterrows():
        n=n+i
        # Créer des nœuds pour les paires d'events
        nodes.append({ 'data': { 'id': row['Event_AOP-helpFinder'].capitalize(), 'type': 'Event' } })
        nodes.append({ 'data': { 'id': row['diseaseName'], 'type': 'DISEASES' } })
        # Créer une arête les events
        edges.append({ 'data': { 'id': f"e{n}", 'source': row['Event_AOP-helpFinder'].capitalize(), 'target': row['diseaseName'], 'Gene_symbol': row['geneSymbol'],'Disease_ID': row['diseaseId'] } })
    return nodes, edges, n


def Reactome_enrichment(Reactome_enrichment_df, nodes, edges, n):
    o=n
    for i, row in Reactome_enrichment_df.iterrows():
        o=o+i
        # Créer des nœuds pour les paires d'events
        nodes.append({ 'data': { 'id': row['Event_AOP-helpFinder'].capitalize(), 'type': 'Event' } })
        nodes.append({ 'data': { 'id': row['Pathway'], 'type': 'Reactome' } })
        # Créer une arête les events
        edges.append({ 'data': { 'id': f"e{o}", 'source': row['Event_AOP-helpFinder'].capitalize(), 'target': row['Pathway'], 'Gene_symbol': row['Gene_associated'],'Reactome_ID': row['Reactome_Id'] } })
    return nodes, edges, o

def WikiPathways_enrichment(WikiPathways_enrichment_df, nodes, edges, o):
    p=o
    for i, row in WikiPathways_enrichment_df.iterrows():
        p=p+i
        # Créer des nœuds pour les paires d'events
        nodes.append({ 'data': { 'id': row['Event_AOP-helpFinder'].capitalize(), 'type': 'Event' } })
        nodes.append({ 'data': { 'id': row['Pathway'], 'type': 'WikiPathways' } })
        # Créer une arête les events
        edges.append({ 'data': { 'id': f"e{p}", 'source': row['Event_AOP-helpFinder'].capitalize(), 'target': row['Pathway'], 'Gene_symbol': row['Gene_associated'],'WikiPathways_ID': row['WikiPathways_Id'] } })
    return nodes, edges, p

#######-----------Open and prep files-----------######
#Scoring from AOPhf
if mode=="se":
    tsv_scoring = os.path.join(folder_AOPhf, "raw_scoring", "scoring_stressor-event.tsv")
if mode=="ee":
    tsv_scoring = os.path.join(folder_AOPhf, "raw_scoring", "scoring_event-event.tsv")

scoring_df = pd.read_table(tsv_scoring)

# Mapping des valeurs de la colonne 'Confidence'
confidence_mapping = {
    'Very High': 5,
    'High': 4,
    'Moderate': 3,
    'Quite low': 2,
    'Low': 1
}
scoring_df['Confidence'] = scoring_df['Confidence'].map(confidence_mapping)
scoring_df['Link'] = scoring_df['Link'].astype(int)


#Output from AOPhf in xlsx
AOPhf_xlsx_name = "AOPhF.xlsx"
xlsx_resAOPhf = os.path.join(folder_AOPhf, AOPhf_xlsx_name)
df_xlsx_aophf = pd.read_excel(xlsx_resAOPhf)
df_xlsx_aophf = process_pairs(df_xlsx_aophf, 30) #Garder les 30 premiers liens pour chaque pairs



# AOP_Wiki enrichment
AOPWiki_enrichment_name = 'annotations_AOPWiki_ntwrk.tsv' #Quentin Capdet file
AOPWiki_rawdb_name = 'raw_AOPwiki.tsv' #raw AOPwiki db in tsv
AOPWiki_raw_json_name = 'AOPwiki_json.json' #raw AOPwiki json to retrieve AOP name
tsv_AOPWiki_enrichment = os.path.join(enrichment_folder, AOPWiki_enrichment_name)
is_AOPWiki_enrichment = os.path.isfile(tsv_AOPWiki_enrichment)


if is_AOPWiki_enrichment:
    aopwiki_enrichment_df = pd.read_table(tsv_AOPWiki_enrichment)
    aopwiki_enrichment_df.index = aopwiki_enrichment_df.index + 1
    tsv_raw_AOPWiki = os.path.join(rawdb_folder, AOPWiki_rawdb_name)
    rawtsv_AOPWiki_df = pd.read_csv(tsv_raw_AOPWiki, delimiter='\t', header=None, names=['AOP', 'Event', 'Type', 'Description'])
    json_raw_AOPWiki = os.path.join(rawdb_folder, AOPWiki_raw_json_name)
    rawjson_AOPWiki_df = pd.read_json(json_raw_AOPWiki)

else :
     aopwiki_enrichment_df = None
     awtsv_AOPWiki_df = None
     rawjson_AOPWiki_df = None

# HPA enrichment
HPA_enrichment_name = 'annotations_HPA_tissue_info.tsv'
tsv_HPA_enrichment = os.path.join(enrichment_folder, HPA_enrichment_name)
is_HPA_enrichment = os.path.isfile(tsv_HPA_enrichment)

if is_HPA_enrichment:
    HPA_enrichment_df = pd.read_table(tsv_HPA_enrichment)
    HPA_enrichment_df.index = HPA_enrichment_df.index + 1
    HPA_enrichment_df["RNA tissue specific"] = HPA_enrichment_df["RNA tissue specific"].str.strip()
else :
    HPA_enrichment_df = None


# Uniprot enrichment (from HPA)
Uniprot_enrichment_name = 'annotations_Uniprot_bio_process.tsv'
tsv_Uniprot_enrichment = os.path.join(enrichment_folder, Uniprot_enrichment_name)
is_Uniprot_enrichment = os.path.isfile(tsv_Uniprot_enrichment)

if is_Uniprot_enrichment:
    Uniprot_enrichment_df = pd.read_table(tsv_Uniprot_enrichment)
    Uniprot_enrichment_df.index = Uniprot_enrichment_df.index + 1
    Uniprot_enrichment_df = Uniprot_enrichment_df[Uniprot_enrichment_df['Biological process'] != 'Not specified']
    Uniprot_enrichment_df["Biological process"] = Uniprot_enrichment_df["Biological process"].str.strip() + " "

else :
    Uniprot_enrichment_df = None


# KEGG enrichment
KEGG_enrichment_name = 'annotations_KEGG.tsv'
tsv_KEGG_enrichment = os.path.join(enrichment_folder, KEGG_enrichment_name)
is_KEGG_enrichment = os.path.isfile(tsv_KEGG_enrichment)

if is_KEGG_enrichment:
    KEGG_enrichment_df = pd.read_table(tsv_KEGG_enrichment)
    KEGG_enrichment_df.index = KEGG_enrichment_df.index + 1
    KEGG_enrichment_df['Pathway'] = KEGG_enrichment_df['Pathway'].str.replace('_', ' ') + "  "

else :
    KEGG_enrichment_df = None


# Reactome enrichment
Reactome_enrichment_name = 'annotations_Reactome.tsv'
tsv_Reactome_enrichment = os.path.join(enrichment_folder, Reactome_enrichment_name)
is_Reactome_enrichment = os.path.isfile(tsv_Reactome_enrichment)

if is_Reactome_enrichment:
    Reactome_enrichment_df = pd.read_table(tsv_Reactome_enrichment)
    Reactome_enrichment_df.index = Reactome_enrichment_df.index + 1
    Reactome_enrichment_df['Pathway'] = Reactome_enrichment_df['Pathway'].str.replace('_', ' ') + "   "

else :
    Reactome_enrichment_df = None

# WikiPathways enrichment
WikiPathways_enrichment_name = 'annotations_WikiPathways.tsv'
tsv_WikiPathways_enrichment = os.path.join(enrichment_folder, WikiPathways_enrichment_name)
is_WikiPathways_enrichment = os.path.isfile(tsv_WikiPathways_enrichment)

if is_WikiPathways_enrichment:
    WikiPathways_enrichment_df = pd.read_table(tsv_WikiPathways_enrichment)
    WikiPathways_enrichment_df.index = WikiPathways_enrichment_df.index + 1
    WikiPathways_enrichment_df['Pathway'] = WikiPathways_enrichment_df['Pathway'].str.replace('_', ' ') + "    "

else :
    WikiPathways_enrichment_df = None

# DisGeNET enrichment
DisG_enrichment_name = 'annotations_DisGeNET_2015.tsv'
tsv_DisG_enrichment = os.path.join(enrichment_folder, DisG_enrichment_name)
is_DisG_enrichment = os.path.isfile(tsv_DisG_enrichment)

if is_DisG_enrichment:
    DisG_enrichment_df = pd.read_table(tsv_DisG_enrichment)
    DisG_enrichment_df.index = DisG_enrichment_df.index + 1
    DisG_enrichment_df['diseaseName'] = DisG_enrichment_df['diseaseName'].str.strip() + "     "
else :
    DisG_enrichment_df = None

# Disease enrichment
DIS_enrichment_name = 'annotations_DISEASES_jensenlab.tsv'
tsv_DIS_enrichment = os.path.join(enrichment_folder, DIS_enrichment_name)
is_DIS_enrichment = os.path.isfile(tsv_DIS_enrichment)

if is_DIS_enrichment:
    DIS_enrichment_df = pd.read_table(tsv_DIS_enrichment)
    DIS_enrichment_df.index = DIS_enrichment_df.index + 1
    DIS_enrichment_df['diseaseName'] = DisG_enrichment_df['diseaseName'].str.strip() + "      "
else :
    DIS_enrichment_df = None




#######-----------Main-----------######

# Calling functions - creating the network
nodes, edges, i = network_creation(scoring_df,mode)
cpt_edge = i

if is_AOPWiki_enrichment: 
    nodes, edges, cpt_edge = aopwiki_enrichment(aopwiki_enrichment_df, nodes, edges, cpt_edge)

if is_Uniprot_enrichment: 
    nodes, edges, cpt_edge = uniprot_enrichment(Uniprot_enrichment_df, nodes, edges, cpt_edge)

if is_KEGG_enrichment:
    nodes, edges, cpt_edge = KEGG_enrichment(KEGG_enrichment_df, nodes, edges, cpt_edge)

if is_DisG_enrichment:
    nodes, edges, cpt_edge = DisG_enrichment(DisG_enrichment_df, nodes, edges, cpt_edge)

if is_DIS_enrichment:
    nodes, edges, cpt_edge = DIS_enrichment(DIS_enrichment_df, nodes, edges, cpt_edge)

if is_Reactome_enrichment:
    nodes, edges, cpt_edge = Reactome_enrichment(Reactome_enrichment_df, nodes, edges, cpt_edge)

if is_WikiPathways_enrichment:
    nodes, edges, cpt_edge = WikiPathways_enrichment(WikiPathways_enrichment_df, nodes, edges, cpt_edge)


elements = nodes + edges

json_dict_AOPhf_edges_info = AOPhf_edges_info(df_xlsx_aophf)


aop_events_list = []
aop_titles_list = []
if is_AOPWiki_enrichment:
    aop_events_list, aop_titles_list = AOPWiki_nodes_info(aopwiki_enrichment_df, rawtsv_AOPWiki_df, rawjson_AOPWiki_df)

hpa_tissue_info_list = []
if is_HPA_enrichment:
    hpa_tissue_info_list = HPA_nodes_info(HPA_enrichment_df)




# Writing results
combined_data = {
    "Network_elements": elements,  # Assurez-vous que `elements` est une structure de données Python sérialisable en JSON
    "Link_information": json_dict_AOPhf_edges_info,
    "AOP_events_info" : aop_events_list,
    "AOP_name_info" : aop_titles_list,
    "HPA_tissue_info" : hpa_tissue_info_list
}


# Écrire les données combinées dans un fichier JSON
with open(output_network, 'w') as combined_file:
    json.dump(combined_data, combined_file, indent=4)




print("Network generated")



#python Network_generation_json.py --aophf_folder ../Case_study/Event_Event --mode ee --output ../../web/aophf-v3/from_tassinca/123/data.json --enrichment ../Enrichment --raw_db ../Folder_raw_db/updated_raw_db
#python Network_generation_json.py --aophf_folder ../Case_study/Stressor_Event --mode se --output ../../web/aophf-v3/from_tassinca/456/data.json --enrichment ../Enrichment --raw_db ../Folder_raw_db/updated_raw_db

