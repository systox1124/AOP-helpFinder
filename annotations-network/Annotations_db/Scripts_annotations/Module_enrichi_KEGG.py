#!/usr/bin/env python3

#authors: Quentin Capdet - Université Paris Cité - France
#         Thomas Jaylet - Université Paris Cité - France
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

# Libraries

import sys
import os
import argparse
import csv
import pandas as pd
import re
import string
#import nltk
#nltk.download('punkt')
#nltk.download('stopwords')
import nltk.data
nltk.data.path.append('./data/nltk_data/')
import nltk.corpus
import nltk.stem
from nltk.tokenize import sent_tokenize
from nltk import word_tokenize

# Functions

def remove_duplicates_hF_se(file_name):
    """
    INPUT : "file scoring"
	Stressor	Event	Article	Link	Pvalue	pvar	Confidence
	Cigarette smoke	Lung cancer	83707	4015	0.000e+00	0.000e+00	Very High
	CCigarette smoke	Oxidative Stress	83707	2323	0.000e+00	0.000e+00	Very High
	Radon	Lung cancer	9591	1297	0.000e+00	0.000e+00	Very High
	Radon	Mutation	9591	112	1.000e+00	NA	Low
	OUTPUT :
	Liste = ["Lung cancer", "Oxidative Stress", "Mutation"]
	Avec stopwords :
	Liste_sw =  ["Lung cancer", "Oxidative Stress", "Mutation"]
	
    Open tsv file (stressor-event) and remove event present more than once (duplicates) and stopwords.
    """
    liste_hf = []
    with open(file_name, 'r') as filin:
        csvreader = csv.reader(filin, delimiter='\t')
        header = next(csvreader)
        if header[0] == "Stressor":
        	rows = [row for row in csvreader]
        	columns = [clean_column(row[1]) for row in rows]
        	columns_raw = [row[1].lower() for row in rows]
        list_with_sw = list(set(columns))
        #print(list_with_sw)
        for event in list_with_sw:
        	event_without_sw = remove_stopwords(event)
        	#print(event_without_sw, '!!!!!!!!!!!')
        	liste_hf.extend(event_without_sw)
        #print(liste_hf)
        
        #create dico_raw
        dico_raw = {}
        for event1, event2 in zip(columns, columns_raw):
        	dico_raw[event1] = event2

        return liste_hf, list_with_sw, dico_raw

def remove_duplicates_hF_ee(file_name):
	"""
	Open tsv file (event-event) and remove event present more than once (duplicates) and stopwords.
	"""
	liste_hf = []
	with open(file_name, 'r') as filin:
		csvreader = csv.reader(filin, delimiter='\t')
		header = next(csvreader)
		if header[0] == "Event 1":
			rows = [row for row in csvreader]
			columns1 = [clean_column(row[0]) for row in rows]
			filin.seek(0)
			next(csvreader)
			columns2 = [clean_column(row[1]) for row in rows]
		list_with_sw = list(set(columns1 + columns2))
		for event in list_with_sw:
			event_without_sw = remove_stopwords(event)
			liste_hf.extend(event_without_sw)
		
		#create dico_raw
		columns1_raw = [row[0].lower() for row in rows]
		columns2_raw = [row[1].lower() for row in rows]
		dico_raw = {}
		for col1, col1_raw in zip(columns1, columns1_raw):
			dico_raw[col1] = col1_raw
		for col2, col2_raw in zip(columns2, columns2_raw):
			dico_raw[col2] = col2_raw
		
		return liste_hf, list_with_sw, dico_raw

def clean_column(column):
    """
    Clean a column value, removing everything after '|'.
    """
    if '|' in column:
        return column.split('|')[0].strip()
    return column.strip()

def process_db(file_name): # (+ liste des noms des éléments à récupérer) / renommer col1, col2 par nom en-tête / file de HPA dispo dans ce dossier et datant de mars 2024.
    """
    Open tsv file (database : HPA), select columns chosen and create a list of genes and other features.
    """
    pro_db = []
    genes_db = []
    with open(file_name, 'r') as filin:
    	csvreader = csv.reader(filin, delimiter='\t')
    	#next(csvreader)
    	for ligne in csvreader:
    		col0 = [ligne[0].split('_')[0]] if ligne[0] else ["Not specified"] # KEGG_ID
    		col1 = ligne[0].split('_')[1:] if ligne[0] else ["Not specified"] # Pathway
    		col1 = ['_'.join(col1)]
    		col3 = ligne[2:] if ligne[2] else ["Not specified"] # Gene_associated
    		col3b = ligne[2:] if ligne[2] else ["Not specified"] # Gene_associated
    		
    		columns = col0 + col1 + col3
    		#print(columns)
    		cleaned_columns = [item.strip() for item in columns if item != '']
    		pro_db.append(cleaned_columns)
    		gens = col3b
    		#print(gens)
    		gens_clean = [item.strip() for item in gens if item != '']
    		genes_db.append(gens_clean)
    	#print(genes_db)
    	#genes_db = [','.join(inner_list) for inner_list in genes_db]
    	#print(genes_db)
	
    return pro_db, genes_db

def remove_stopwords(list_events):
    """
    remove stopwords ('a', 'the', etc..) and add space between terms like 'injury/death'.
    """
    sents = sent_tokenize(list_events)
    sents = [sent for sent in sents if not (bool(re.search('\d', sent) and 'body weight' in sent))]
    list_events = [word_tokenize(sent) for sent in sents]
    punctuation = list(string.punctuation)
    stop = nltk.corpus.stopwords.words('english') + punctuation + ['\'s']
    for i in range(len(list_events)):
        # Check for terms like 'injury/death' and separate them with spaces
        new_list = []
        for word in list_events[i]:
            if '/' in word:
                # Split terms like 'injury/death' into separate words
                new_list.extend(re.split(r'/', word))
            else:
                new_list.append(word)
        list_events[i] = [word for word in new_list if word not in stop]
     
    return list_events

def find_words_btw_hF_db(liste_hf, genes_db, list_with_sw, pro_db):
	"""
	find words in common between hF and database events.
	"""

	HPA_df = pd.read_csv(gene_synonym, sep='\t')
	
	HPA_df['Gene'] = HPA_df['Gene'].fillna('Not specified').str.split(',')
	HPA_df['Gene synonym'] = HPA_df['Gene synonym'].fillna('Not specified').str.split(',')
	HPA_df['Gene description'] = HPA_df['Gene description'].fillna('Not specified').str.split(',')
	
	gene_syn = [list(map(str.strip, a + b)) for a, b in zip(HPA_df['Gene'], HPA_df['Gene synonym'])]
	
	gene_name = [list(map(str.strip, a + b)) for a, b in zip(HPA_df['Gene'], HPA_df['Gene description'])]
	
	#Gene and gene_synonym identification
	gene_dico = {}
	liste_gene = []
	for i, pro_list in enumerate(liste_hf):
		for sublist_gene in gene_syn:
			for gene in sublist_gene:
				#print(gene)
				for word in pro_list:
					if gene.lower() == word.lower():
						liste_gene.append(gene)
						gene_symbol = [g for list_g in genes_db for g in list_g if g in sublist_gene]
						if len(gene_symbol) > 0:
							gene_dico[list_with_sw[i]] = gene_symbol[0]
						else:
							pass
	
	#print(liste_gene)
	
	#Gene description identification
	# Créer un dico des noms de gènes avec les id correspondants
	gene_name_dict = {sublist_pro[1].lower(): sublist_pro[0] for sublist_pro in gene_name}
	# Générer et compiler une expression régulière unique pour correspondre à tous les noms de gènes
	gene_pattern = r'\b(?:' + '|'.join(map(re.escape, gene_name_dict.keys())) + r')\b'
	gene_regex = re.compile(gene_pattern)

	gene_name_dico = {}
	for event in list_with_sw:
		# Rechercher tous les noms de gènes dans l'événement
		matches = gene_regex.findall(event.lower())
		for match in matches:
		    gene_name_dico[event] = gene_name_dict[match]

	#print(gene_name_dico)

	
	gene_event = gene_dico.copy()
	gene_event.update(gene_name_dico)
	#print(gene_event)
	
	return gene_event

def creation_table_genes(pro_db, output_filename, gene_event) :
	"""
	Create the table with events of hF, gene and features of db
	"""
	#print(pro_db)
	remove_name = ["DAMAGE", "Not specified"] #NOT RESEARCHED IN HPA
	data = []
	for genes, event in gene_event.items():
		#print(event, '!!!!!!!!!!!', genes)
		if event not in remove_name:
			for slist in pro_db:
				for gene in slist[2:]:
					g = gene.split(' ')
					for w in g:
						if re.match(r'\b' + re.escape(event.lower()) + r'\b(?![-])', w.lower()):
							#print(event, w, slist[2:], slist[0])
							gene_info_cleaned = [re.sub(r'[\d:.]', '', item) for item in slist[4:]]
							data.append({
										"Event_AOP-helpFinder":genes,
										"KEGG_Id": slist[0],
										"Pathway": ','.join(slist[1:2]), 
										"Gene_associated": event #','.join(slist[2:]) #event
									})

	table_complete = pd.DataFrame(data)
	table_complete.drop_duplicates(keep = 'first', inplace=True)
	for i, event in enumerate(table_complete.iloc[:,0]):
		for val, val_raw in dico_raw.items():
			if event == val:
				table_complete.iloc[i, 0] = dico_raw[val]
	table_complete.to_csv(output_filename, sep='\t', index=False)
	
	return output_filename, table_complete

# Main Program

parser = argparse.ArgumentParser(description="Module d'enrichissement (base disponible : AOPWiki)")
parser.add_argument("-hf", "--AOPHF", help="Résultats trouvés par AOP-helpFinder pour les cas stressor-event : liste d'events")
parser.add_argument("-db", "--database", help="Base de données pour l'enrichissement des events obtenu par AOPHF")
parser.add_argument("-o", "--output", help="Fichier de sortie qui sera un tableau tsv")
parser.add_argument("-mod", "--MODE", help="Sélection du mode selon la nature des résultats d'AOP_hF obtenu avec des stressor-event (se) ou des event-event (ee)")
parser.add_argument("-syn", "--synonyme", help="Liste des synonymes des gènes")
args = parser.parse_args()

select_mode = args.MODE #se ou ee
list_events = args.AOPHF #AOP_hF "scoring_stressor-event.tsv"
database = args.database #database "proteinatlas.tsv"
output_filename = args.output #output "{name}.tsv"
gene_synonym = args.synonyme #liste des synonymes des gènes

# Function

if select_mode == "se":
	liste_hf, list_with_sw, dico_raw = remove_duplicates_hF_se(list_events)
	
elif select_mode == "ee":
	liste_hf, list_with_sw, dico_raw = remove_duplicates_hF_ee(list_events)

else:
	print("ERROR ! Option du mode incorrecte lire l'instruction '-help' de Module_enrichi_HPA.py") 

pro_db, genes_db = process_db(database)
#print(pro_db)
gene_event = find_words_btw_hF_db(liste_hf, genes_db, list_with_sw, pro_db)  
#print(gene_event)
output_file, table_complete = creation_table_genes(pro_db, output_filename, gene_event)

