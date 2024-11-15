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

#import sys
#import os
import argparse
import csv
import pandas as pd
import re
import string #sw
#import nltk  #sw & stem
#nltk.download('punkt')
#nltk.download('stopwords')
import nltk.data
nltk.data.path.append('./data/nltk_data/')
import nltk.corpus
import nltk.stem
from nltk.tokenize import sent_tokenize
from nltk import word_tokenize
#import subprocess

# Functions

def remove_duplicates_hF_se(file_name):
    """
    Open tsv file and remove event present more than once (duplicates).
    """
    with open(file_name, 'r') as filin:
        csvreader = csv.reader(filin, delimiter='\t')
        header = next(csvreader)
        if header[0] == "Stressor":
        	rows = [row for row in csvreader]
        	columns = [clean_column(row[1].lower()) for row in rows]
        	columns_raw = [row[1].lower() for row in rows]
        liste = list(set(columns))
        
        #create dico_raw
        dico_raw = {}
        for event1, event2 in zip(columns, columns_raw):
        	dico_raw[event1] = event2
        return liste, dico_raw

def remove_duplicates_hF_ee(file_name):

	with open(file_name, 'r') as filin:
		csvreader = csv.reader(filin, delimiter='\t')
		header = next(csvreader)
		if header[0] == "Event 1":
			rows = [row for row in csvreader]
			columns1 = [clean_column(row[0].lower()) for row in rows]
			filin.seek(0)
			next(csvreader)
			columns2 = [clean_column(row[1].lower()) for row in rows]
		liste = list(set(columns1 + columns2))
		
		#create dico_raw
		columns1_raw = [row[0].lower() for row in rows]
		columns2_raw = [row[1].lower() for row in rows]
		dico_raw = {}
		for col1, col1_raw in zip(columns1, columns1_raw):
			dico_raw[col1] = col1_raw
		for col2, col2_raw in zip(columns2, columns2_raw):
			dico_raw[col2] = col2_raw
		
		return liste, dico_raw

def clean_column(column):
    """
    Clean a column value, removing everything after '|'.
    """
    if '|' in column:
        return column.split('|')[0].strip()
    return column.strip()

def process_file(file_name):
    """
    Open tsv file, remove duplicates from column 2 and make a word process (stopwords + stemming).
    """
    res = []
    seen_values = set()
    with open(file_name, 'r') as filin:
        for ligne in filin.readlines():
        	columns = ligne.strip().split('\t')
        	if columns[1] not in seen_values:
        		seen_values.add(columns[1])
		        sw = remove_stopwords(ligne)
		        for words_list in sw:
		            res.append(stem_process(sorted(words_list)))
    
    return res

def process_list(liste):
    """
    Take a list and make a word process (stopwords + stemming).
    """
    res = []
    for ligne in liste:
            sw = remove_stopwords(ligne.lower())
            for words_list in sw:
                res.append(stem_process(sorted(words_list)))
    
    return res

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

def stem_process(words):
    """Method to stem each words in a list.
    Return:
        words (list): a list which contains stemmed words.
    """
    # snowball stemmers developed by Martin Poter
    sno = nltk.stem.SnowballStemmer('english')
    # search if a particular word is in words
    for i in range(len(words)):
        words[i] = sno.stem(words[i])
    
    return words

def check_words_presence(word, resdb):
    """
    Determine the links between AOP-hF and AOP-wiki events if present.
    """
    set1 = set(word)
    result_list = []
    found = False
    for sublist in resdb:
        set2 = set(sublist)
        if set1.issubset(set2):
	    # Filter words of the form 'event:numéro' from set2
            filtered_set2 = {word for word in set2 if word.startswith('event:')}
            result_list.append(filtered_set2)
            #print(set2, filtered_set2)
            found = True
    
    return found, result_list

def find_words_btw_hF_db(reshF, resdb, uniq_hF):
	"""
	find words in common between hF and database events (after remove stemming and stopwords)
	"""
	found_all_words = True
	nbr_event = []
	all_filtered_sets = []
	remove_line = []
	
	for i, word in enumerate(reshF):
	    #print(word)
	    found, result_list = check_words_presence(word, resdb)

	    if not found:
		    #print(f"L'un des mots de la liste n'est présents dans aucune sous-liste : '{i}, {word}'")
		    found_all_words = False
		    remove_line.append(i)
	    else:
	    	nbr_event.append(len(result_list))
	    	all_filtered_sets.extend(result_list)
	    	#print(nbr_event)

	#if found_all_words:
		#print(f"Tous les mots de la liste sont présents dans au moins une sous-liste.")

	# Create a list of event numbers
	num_events = [item.pop().capitalize() for item in all_filtered_sets]
	#print(num_events)
	
	# Delete the specified row(s)
	new_events_hF = [word for i, word in enumerate(uniq_hF) if i not in remove_line]
	#print(remove_line)
	#print(new_events_hF, len(new_events_hF))
	
	occ_event = [event for event, nbr in zip(new_events_hF, nbr_event) for _ in range(nbr)] 
	#print(occ_event)
	num_events_df = pd.DataFrame({'Colonne1': num_events, 'Colonne2': occ_event})
	#tsvfile = "num_events.tsv"
	#col1 = num_events_df.iloc[:,0]
	#col1.to_csv(tsvfile, sep='\t', index=False)

	return num_events_df, new_events_hF

def update_event_with_synonym(file_syn, liste_events) :
	"""
	Updating the event list with synonyms based on associated words with a list of synonyms created from a tsv file.
	"""
	# Build a list of synonym
	list_syn = []
	corresp = {}
	liste_events = [event.lower() for event in liste_events]
	with open(file_syn, 'r') as filin:
		csvreader = csv.reader(filin, delimiter='\t')
		for row in csvreader:
			row_lower = [item.lower() for item in row]
			list_syn.append(row_lower)
	
	Liste_maj = set()
	# Check if all elements of event list are present in synonym list or not
	for element in liste_events:
		for syn in list_syn:
		    for terme in syn:
		        if terme in element and all(val in element.split() for val in terme.split()):
		            # Replace synonyms in the element, except the one that is already present
		            for nouveau_syn in syn:
		                if nouveau_syn != terme and nouveau_syn not in liste_events:
		                    nouvel_element = element.replace(terme, nouveau_syn)
		                    Liste_maj.add(nouvel_element)
		                    #print(f"La sous-liste '{syn}' a été trouvée dans la liste complète '{element}'. Mise à jour en '{nouvel_element}'.")
		                    corresp[nouvel_element] = element
		            else:
		                # If all the synonyms are already in the full list, do not update
		                Liste_maj.add(element)
		else:
		    Liste_maj.add(element)
	
	return list(Liste_maj), corresp

def write_to_file(file_path, data):
    """
    Write a tsv file.
    """
    with open(file_path, 'w') as file:
        for item in data:
            file.write(f"{item}\n")

def search_pattern_in_db(database, num_events_df):
	"""
	Make a list of event number (pattern) from a tsv file and retrieve the lines in the database (events_wiki) corresponding to each pattern 
	"""
	data_events_wiki = []
	search_patterns = num_events_df.iloc[:, 0].tolist()
	
	results = {pattern: {'aops': set(), 'event_type': set()} for pattern in search_patterns}
	aop_event_mapping = {}
	
	with open(database, 'r', newline='', encoding='utf-8') as tsv_file:
		tsv_reader = csv.reader(tsv_file, delimiter='\t')
		for row in tsv_reader:
			new_columns = [row[3], row[1], row[0], row[2]]
			line = '\t'.join(new_columns)
			
			for pattern in search_patterns:
				if re.search(r'\b{}\b'.format(re.escape(pattern)), line):
					before_event_match = re.search(r'(.+?)\s*Event:', line)
					if before_event_match:
						before_event = before_event_match.group(1).strip()
						if pattern not in aop_event_mapping:
							aop_event_mapping[pattern] = before_event
					
					aop_match = re.search(r'Aop:(\d+)', line)
					if aop_match:
						aop = aop_match.group(1)
						results[pattern]['aops'].add(aop)
					
					event_type_match = re.search(r'(KeyEvent|MolecularInitiatingEvent|AdverseOutcome)', line)
					if event_type_match:
						event_type = event_type_match.group(1)
						results[pattern]['event_type'].add(event_type)
		
		for pattern in search_patterns:
			aops_str = ', '.join(sorted(results[pattern]['aops']))
			event_type_str = ', '.join(sorted(results[pattern]['event_type']))
			
			if aop_event_mapping[pattern] and aops_str and event_type_str:
				event_number = pattern.replace("Event:", "")
				data_events_wiki.append({
            		"Event_AOP-wiki": aop_event_mapping[pattern],
            		"N°Event": event_number,
            		"N°AOP": aops_str,
            		"Type_event": event_type_str
            	})
	events_wiki_df = pd.DataFrame(data_events_wiki)
	return  events_wiki_df

def replace_syn_by_right_word(num_events_df, correspondance):
	num_events_df.iloc[:,1] = num_events_df.iloc[:,1].replace(correspondance)
	return num_events_df

def creation_table(events_wiki_df, num_events_df, output_file):
	"""
	Creation of the enrichment table
	"""
	col2 = num_events_df.iloc[:,1]
	col2 = col2.rename('Event_AOP-helpFinder')
	Table_resultat = pd.concat([col2, events_wiki_df], axis=1)
	Table_resultat = Table_resultat.drop_duplicates(subset="N°Event", keep='first')
	
	for i, event in enumerate(Table_resultat.iloc[:,0]):
		for val, val_raw in dico_raw.items():
			if event == val:
				Table_resultat.iloc[i, 0] = dico_raw[val]
		
	Table_resultat.to_csv(output_file, sep='\t', index=False)
	
	return output_file

def table_for_interactiv_network(output_file, database, interact_ntwrk):

	# Charger le tableau initial
	df = pd.read_csv(output_file, delimiter='\t')

	# Diviser les valeurs de la colonne 'N°AOP' qui sont des listes séparées par des virgules
	df['N°AOP'] = df['N°AOP'].apply(lambda x: [int(num) for num in str(x).split(', ')] if pd.notna(x) else [])

	# Utiliser explode pour répéter les lignes avec plusieurs N°AOP et reset_index
	df_expanded = df.explode('N°AOP')
	df_expanded = df_expanded.reset_index(drop=True)

	for index, row in df_expanded.iterrows():
		# Récupérer les informations de la 3ème et 4ème colonne
		event_info = row.iloc[2]
		aop_info = row.iloc[3]

		# Formater les informations 
		formatted_event = f"Event:{event_info}"
		formatted_aop = f"Aop:{aop_info}"

		# Faire une recherche dans un fichier avec les termes formatés
		search_terms = [formatted_event, formatted_aop]

		# fichier de recherche
		with open(database, "r") as file:
			lines = file.readlines()
			for line in lines:
				if all(term in line for term in search_terms):
				    new_value = line.split('\t')[2].strip()
				    df_expanded.at[index, df_expanded.columns[4]] = new_value

	# Sauvegarder le nouveau tableau dans un fichier tsv
	df_expanded.to_csv(interact_ntwrk, sep='\t', index=False)


# Main Program (pseudocode)

parser = argparse.ArgumentParser(description="Module d'enrichissement (base disponible : AOPWiki)")
parser.add_argument("-hf", "--AOPHF", help="Résultats trouvés par AOP-helpFinder pour les cas stressor-event : liste d'events")
parser.add_argument("-db", "--database", help="Base de données pour l'enrichissement des events obtenu par AOPHF")
parser.add_argument("-syn", "--synonym", help="Liste de synonymes")
parser.add_argument("-o", "--output", help="Fichier de sortie qui sera un tableau tsv")
parser.add_argument("-net", "--NETWORKTAB", help="Tableau tsv dans un format adapté au réseau interactif")
parser.add_argument("-mod", "--MODE", help="Sélection du mode selon la nature des résultats d'AOP_hF obtenu avec des stressor-event (se) ou des event-event (ee)")
args = parser.parse_args()

select_mode = args.MODE #se ou ee
list_events = args.AOPHF #AOP_hF "scoring_stressor-event.tsv"
database = args.database #database "raw_AOPwiki.tsv / ..."
file_syn = args.synonym #synonym "AOPWiki_synonym.tsv"
output_file = args.output #output "{name}.tsv"
interact_ntwrk = args.NETWORKTAB #NETWORKTAB "Table_enrichi_network.tsv"

if select_mode == "se":
	uniq_hF, dico_raw = remove_duplicates_hF_se(list_events)
	
elif select_mode == "ee":
	uniq_hF, dico_raw = remove_duplicates_hF_ee(list_events)

else:
	print("ERROR ! Option du mode incorrecte lire l'instruction help '-h' de Module_enrichi.py") 
	
#uniq_hF = remove_duplicates_hF_se(list_events)
stem_hF = process_list(uniq_hF)
stem_db = process_file(database)
num_events_df, new_events_hF = find_words_btw_hF_db(stem_hF, stem_db, uniq_hF)

if file_syn:
	Liste_maj, corresp = update_event_with_synonym(file_syn, uniq_hF)
	stem_hF = process_list(Liste_maj)
	num_events_df, new_events_hF = find_words_btw_hF_db(stem_hF, stem_db, Liste_maj)
	num_events_df = replace_syn_by_right_word(num_events_df, corresp)

events_wiki_df = search_pattern_in_db(database, num_events_df)
output_file = creation_table(events_wiki_df, num_events_df, output_file)

if interact_ntwrk:
	table_for_interactiv_network(output_file, database, interact_ntwrk)








