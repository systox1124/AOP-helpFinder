# -*- coding: utf-8 -*-
#authors: Florence Jornod - INSERM UMRS 1124
#         Thomas Jaylet - Université de Paris - France
#         Karine Audouze - Université de Paris - France

#contact: systox@paris-descartes.fr


#AOP-helpFinder is provided without any warranty. But if you have any probleme please feel free to contact us by mail.

#------- WHAT IS AOPHELPFINDER? -------------

#AOP-helpFinder is a tool developed to help AOP development (Jean-Charles Carvaillo: https://github.com/jecarvaill/aop-helpFinder)(Environ Health Perspect. 2019 Apr;127(4):47005).

#It is based on text mining and parsing process on scientific abstracts. AOP-helpFinder identify links between stressors and molecular initiating event, key events and adverse outcomes through abstracts from the PubMed database (https://pubmed.ncbi.nlm.nih.gov/).

#AOP-helpFinder was implemented under the H2020 Human Biomonintoring in Europe (HBM4EU) project, Work Package 13.1.
#HBM4EU has received funding from the European Union’s H2020 research and innovation programme under grant agreement No 733032.

#------- LICENCE ---------------------------

#This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can  use,  modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL

#http://cecill.info/licences/Licence_CeCILL_V2.1-en.txt

#-------------------------------------------
#-------------------------------------------
import csv
import tm_module as tm


#-------------------------------------------------------------------------------
def create_collection_event(events_file, events_type):
    """
    Launch the cleaning according to the parameters provided in param.py
    """
    if events_type == "homemade":
        list_events = get_collection_event_homemade(events_file)
    if events_type == "homemade_synonym":
        list_events = get_collection_event_homemade_synonym(events_file)
    list_events=cleaning_event_list(list_events)
    return list_events


#-------------------------------------------------------------------------------
def get_collection_event_homemade(event_filename):
    """
    Method to take from a file (event_filename) a list of event
    The file must be a .txt with one event per line
    ARGUMENT: event_filename -> the name of the file with the information about the event
    RETURN: list_events -> a dictionnary with all the events.
    """
    list_events = []
    with open(event_filename,mode="r") as infile:
        reader = csv.reader(infile, delimiter = "\t")
        for rows in reader:
            if rows:
                dict_event = {}
                event = rows[0].strip()
                dict_event["name"] = [event]
                dict_event["stem_synonym"] = tm.clean_abstract(event, "stem")
                list_events.append(dict_event)
    return list_events


def get_collection_event_homemade_synonym(event_filename):
    """Method to take from a file (event_filename) a list of event
    The file must be a .txt with each different event per line.
    Synonyms must be written on the same line and separated by a tabulation
    ARGUMENT: event_filename -> the name of the file with the information about the event
    RETURN: list_events -> a dictionnary with all the events.
    """
    list_events = []
    with open(event_filename,mode="r") as infile:
        reader = csv.reader(infile, delimiter = "\t")
        for rows in reader:
            if rows:
                dict_event = {}
                for event in rows:
                    event = event.strip()
                    if 'name' not in dict_event:
                        dict_event['name'] = [event]
                        dict_event["stem_synonym"] = tm.clean_abstract(event,"stem")
                    else :
                        dict_event['name'].append(event)
                        dict_event["stem_synonym"].append(tm.clean_abstract(event,"stem")[0])
                list_events.append(dict_event)
    return list_events


################################################################################
def cleaning_event_list(dico_event):
    """
    Clean event list to remove duplicates (based on stemnames)
    e.g. Increased or Increase -> Increas (same stame)
    """
    dico_event_cleaned = []
    for i in range(len(dico_event)):
        if dico_event[i]["name"] != "Delete event" :
            for j in range(i+1, len(dico_event)):
                if len(set(dico_event[i]["stem_synonym"]) & set(dico_event[j]["stem_synonym"])) > 0:
                    new_stem= list(set(dico_event[i]["stem_synonym"]) | set(dico_event[j]["stem_synonym"]))
                    new_name= list(set(dico_event[i]["name"]) | set(dico_event[j]["name"]))
                    dico_event[i]["name"] = new_name
                    dico_event[i]["stem_synonym"] = new_stem
                    dico_event[j]["name"] = "Delete event"
            dico_event_cleaned.append(dico_event[i])
    return(dico_event_cleaned)

#-------------------------------------------------------------------------------
#Function test regardless of AOPhF
"""
test = create_collection_event("event.txt", "homemade", lemma = False)
print("\nBefore cleaning: \n",test)
test = cleaning_event_list(test)
print("\nAfter cleaning: \n",test)

test = create_collection_event("event.txt", "homemade_synonym", lemma = False)
print("\nBefore cleaning: \n",test)
test = cleaning_event_list(test)
print("\nAfter cleaning: \n",test)
"""
