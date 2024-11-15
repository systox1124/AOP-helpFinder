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
import networkx as nx
import tm_module as tm
from nltk.tokenize import sent_tokenize
from nltk import word_tokenize



#-------------------------------------------------------------------------------
def compute_scores(abstract, events, ratio_intro, search_method):
    """
    Launch the scoring according to the parameters provided
    """
    if search_method =="event-event":
        scoring = compute_scores_event_event(abstract, events, ratio_intro)
        return scoring
    if search_method == "stressor-event":
        return compute_scores_stressor_event(abstract, events, ratio_intro)



#-------------------------------------------------------------------------------
def compute_scores_stressor_event(abstract, events, ratio_intro):
    ''' METHOD FOR STRESSOR-EVENT: compute the score of presence of the event in the abstract, return also the different position in the abstract.
	ARGUMENTS:
    abstract -> a string with the abstract
    events -> the list of the events
    ratio_intro -> The percentage of first sentences of the abstract that we want to ignore (20% adapted to reduce false positive rates)
	RETURN: dico -> a dictionnary with the name if the event as key and the scores (multiscre ; localisation ; sentence found)
    '''
    abstract_stem = tm.clean_abstract(abstract["abstract"], "stem")
    len_abstract_stem = len(abstract_stem)
    abstract_raw_without_neg = tm.remove_negative_sentences(abstract["abstract"])
    abstract_raw_in_sentences = sent_tokenize(abstract["abstract"])
    dico={}

    for event in events:
        score = []
        words_found = []
        for synonym in event["stem_synonym"]:
            if synonym!='':
                score, localisation, sentence, words_found = find_event(abstract_stem, abstract_raw_without_neg, synonym)
            if(score):
                # If there is a score, we then check the position in the abstract
                localisation, ratio_localisation = real_localisation(abstract_raw_in_sentences, sentence)
                score, localisation, sentence, words_found, ratio_localisation = check_position(score, localisation, sentence, words_found, ratio_localisation, ratio_intro)
                if(score):
                    name = " | ".join(event["name"])
                    if name in dico :
                        for i in range(len(sentence)):
                            if sentence[i] not in dico[name]["sentence"]:
                                dico[name]["multiscore"].append(score[i])
                                dico[name]["localisation"].append(localisation[i])
                                dico[name]["sentence"].append(sentence[i])
                                dico[name]["words_found"].append(words_found[i]) #d
                    else :
                        dico[name]={}
                        dico[name]["multiscore"]=score
                        dico[name]["localisation"]=localisation
                        dico[name]["sentence"]=sentence
                        dico[name]["words_to_find"]=len(synonym.split()) #d
                        dico[name]["words_found"]= words_found #d
    return(dico)

################################################################################
def compute_scores_event_event(abstract, events, ratio_intro):
    ''' METHOD FOR EVENT-EVENT: compute the score of presence of the event in the abstract, return also the different position in the abstract.
    ARGUMENTS:
    abstract -> a string with the abstract
    events -> the list of the events
    ratio_intro -> The percentage of first sentences of the abstract that we want to ignore (20% adapted to reduce false positive rates)
    RETURN: dico -> a dictionnary with the name if the event as key and the scores (multiscre ; localisation ; sentence found)
    '''
    abstract_stem = tm.clean_abstract(abstract["abstract"], "stem")
    len_abstract_stem = len(abstract_stem)
    abstract_raw_without_neg = tm.remove_negative_sentences(abstract["abstract"])
    abstract_raw_in_sentences = sent_tokenize(abstract["abstract"])
    dico = {}
    for i in range(len(events)):
        score_i = [] ; words_found_i = []
        for synonym_i in events[i]["stem_synonym"]:
            if synonym_i !='':
                score_i, localisation, sentence, words_found_i = find_event(abstract_stem, abstract_raw_without_neg, synonym_i)
                if(score_i):
                    localisation, ratio_localisation = real_localisation(abstract_raw_in_sentences, sentence)
                    score_i, localisation, sentence, words_found_i, ratio_localisation = check_position(score_i, localisation, sentence, words_found_i, ratio_localisation, ratio_intro)
                    if(score_i):
                        for j in range(i+1, len(events)):
                            for synonym_j in events[j]["stem_synonym"]:
                                sentence_event_i = []
                                for sentence_i in sentence :
                                    sentence_event_i.append(" ".join(tm.clean_abstract(sentence_i, "stem")))
                                event_split = list(set(synonym_j.split()))
                                score_j = [] ; words_found_j = [] ; index_j = []
                                for w in range(len(sentence_event_i)):
                                    best, words = presence_score(sentence_event_i, event_split,w)
                                    if (abs(1-best/len(event_split)) <= 1 and best!=0):
                                        score_j.append(best/len(event_split))
                                        words_found_j.append(words)
                                        index_j.append(w)
                                score_i_update = [score_i[x] for x in index_j]
                                cooc_multiscore = [score_i_update, score_j]
                                words_found_i_update = [words_found_i[x] for x in index_j]
                                #cooc_wfound = [words_found_i_update, words_found_j]
                                cooc_sentence = [sentence[x] for x in index_j]
                                cooc_loca = [localisation[x] for x in index_j]
                                if (score_i_update and score_j):
                                    name = "{} with {}".format(" | ".join(events[i]["name"]), " | ".join(events[j]["name"]))
                                    if name in dico:
                                        for s in range(len(cooc_sentence)):
                                            if cooc_sentence[s] not in dico[name]["sentence"]:
                                                dico[name]["multiscore"][0].append(cooc_multiscore[0][s])
                                                dico[name]["multiscore"][1].append(cooc_multiscore[1][s])
                                                #dico[name]["cooc_wfound"][0].append(cooc_multiscore[0][s]) #d
                                                #dico[name]["cooc_wfound"][1].append(cooc_multiscore[1][s]) #d
                                                dico[name]["localisation"].append(cooc_loca[s])
                                                dico[name]["sentence"].append(cooc_sentence[s])
                                    else :
                                        dico[name]={}
                                        dico[name]["multiscore"]=cooc_multiscore
                                        dico[name]["localisation"]=cooc_loca
                                        dico[name]["sentence"]=cooc_sentence
                                        dico[name]["words_to_find"]=[[len(synonym_i.split())], [len(synonym_j.split())]] #d
                                        #dico[name]["words_found"]= cooc_wfound #d
    #print(dico)
    return dico


#-------------------------------------------------------------------------------
def find_event(sentences_stem, sentences_raw, event):
    ''' METHOD: find a event in a list of sentences and compute the score
    ARGUMENTS: sentences -> a string with the sentences (for example an abstract)
    event -> a string with one event
    RETURN: score, localisation -> the score associated to the event and the abstract, and the localisation of the event (index of the list of sentences) , sentence and number of words found in the sentence
    '''
    event_split = list(set(event.split()))
    if not event_split:
        return
    score = []
    words_found = []
    localisation = []
    phrase = []
    ind = 0
    for i in range(len(sentences_stem)):
        best, cpt = presence_score(sentences_stem,event_split,i)
        if (abs(1-best/len(event_split)) <= 1 and best!=0):
            score.append(best/len(event_split))
            localisation.append(i)
            words_found.append(cpt)
            phrase.append(sentences_raw[i])
        ind += 1
    return score, localisation, phrase, words_found


################################################################################
def presence_score(sentence, event,i):
	''' METHOD: check if a event is in a sentence in a text (sentences)
	    ARGUMENTS: sentences -> a string with the sentences (for example an abstract)
		       event -> a list of the words from one event
		       i -> the index of the sentence
	    RETURN: best -> the best score found ; cpt -> number of words found
    '''
	best = 0
	cpt = 0
	index = {}
	words = sentence[i].split()
	for e in event:
		if e in words:
			cpt += 1
			index[e] = [pos+1 for pos, value in enumerate(words) if value == e]
	proportion = cpt/len(event)
	if 0.75 <= proportion and len(event) != 1:
		values = index.values()
		values = [sorted(value) for value in values]
		values = sorted(values, key = lambda k: k[-1])
		G = nx.Graph()
		for i in range(len(values) -1):
			for j in range(len(values[i])):
				for l in range(len(values[i + 1])):
					G.add_edge(values[i][j], values[i+1][l], weight = abs(values[i+1][l]-values[i][j]))
		shortest_paths = {}
		for start in values[0]:
			for end in values[-1]:
				shortest_path = nx.dijkstra_path(G, start, end)
				tot_weight = nx.dijkstra_path_length(G, start, end)
				shortest_paths[tot_weight +1] = shortest_path
		best = min(shortest_paths, key = shortest_paths.get)
	elif proportion == 1.0 and len(event) ==1 :
		best = 1.0
	return best, cpt

################################################################################
def real_localisation(abstract_raw_without_tm, sentences_found):
    """
    Find the real localization of the sentence found in the complete abstract (without preprocessing)
    Returns the position of the sentence in the abstract as an index and a ratio between 0 and 1
    """
    len_abstract_raw = len(abstract_raw_without_tm)
    ratio_localisation = []
    index_localisation = []
    for sentence_found in sentences_found :
        real_localisation = - 1
        for i in range(len_abstract_raw):
            if ' '.join(word_tokenize(abstract_raw_without_tm[i]))== sentence_found :
                real_localistation = i
        ratio = real_localistation/len_abstract_raw
        index_localisation.append(real_localistation)
        ratio_localisation.append(ratio)
    return index_localisation, ratio_localisation

################################################################################
def check_position(score, localisation, sentence, words_found, r_localisation, limit):
    """
    Remove the score if the found sentence is below the limit
    """
    cpt = 0
    for j in range(len(r_localisation)):
        i = j - cpt
        if r_localisation[i] <= limit :
            del score[i]
            del localisation[i]
            del sentence[i]
            del r_localisation[i]
            del words_found[i]
            cpt += 1
    return score, localisation, sentence, words_found, r_localisation
