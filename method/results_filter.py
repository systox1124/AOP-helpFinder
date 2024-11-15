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
import tm_module as tm
import compute_scores as cscores
import spacy
nlp = spacy.load("en_core_web_sm")

import nltk.data
nltk.data.path.append('./data/nltk_data/')

import nltk.corpus
import nltk.stem
from nltk.tokenize import sent_tokenize
from nltk import word_tokenize
import re
import string
import param

#-------------------------------------------------------------------------------
def refinement_filter(dict_abstracts, list_events, search_mode):
    """
    Launch the filter according to the parameters provided in param.py
    """
    print("\nRefinement_filter\n")
    if param.context_choice is True:
        context_choice(dict_abstracts, search_mode)
    if param.tm_lemma is True:
        lemma_filter(dict_abstracts, list_events, search_mode)

#-------------------------------------------------------------------------------
def context_choice(dict_abstracts, search_mode):
    """
    -Filter on AOP-HelpFinder results: eliminates sentences containing context words to improve accuracy
    -Input: dictionary containing the results; list of events to lemmatize and type of search (event-event or stressor-event)
    -Output : Filtered dictionary (deletion of links not found by lemmatization)
    """
    context_word = [ "investigated", "investigate", "examined", "examine", "studied",
                    "tested", "assayed", "measured", "evaluated", "objective", "objectives",
                    "method", "methods", "background", "aimed", "context", "assesses", "explored",
                    "aim","assays", "evaluate", "explore", "goal", "analyzed", "monitored",
                    "abbreviations", "performed", "examining", "conducted", "compiles",
                    "purpose", "recruiter", "assessed", "assess", "explores"]
    contexts = [context.lower() for context in context_word]

    for element in dict_abstracts.copy():
        if 'score' in element:
            for key, val in element["score"].copy().items():
                sents = val["sentence"]
                abstract = [word_tokenize(sent.lower()) for sent in sents]
                index = []
                for i in range(len(abstract)):
                    flag = False
                    for word in abstract[i]:
                        if word in contexts:
                            flag = True
                    if flag is False:
                        index.append(i)
                val["sentence"] = [val["sentence"][x] for x in index]
                val["localisation"] = [val["localisation"][x] for x in index]
                if search_mode == "event-event":
                    score_i = [val["multiscore"][0][x] for x in index]
                    score_j = [val["multiscore"][1][x] for x in index]
                    val["multiscore"] = [score_i, score_j]
                if search_mode == "stressor-event" :
                    val["multiscore"] = [val["multiscore"][x] for x in index]
                if not val["sentence"]:
                    del element["score"][key]
            if not element["score"]:
                del element["score"]

#-------------------------------------------------------------------------------
def lemma_filter(dict_abstracts, list_events, search_mode):
    """
    -Filter on AOP-HelpFinder results: Lemmatization of sentences found by AOP-helpFinder to improve accuracy compared to stemming
    -Input: dictionary containing the results and type of search (event-event or stressor-event)
    -Output : Filtered dictionary
    """
    for event in list_events:
        event["lemmaname"] = list(set([tm.clean_abstract(x, "lemma")[0] for x in event["name"]]))
    for element in dict_abstracts.copy():
        if 'score' in element:
            for key, val in element["score"].copy().items():
                sentence = val["sentence"]
                sentence_event_i = []
                for sentence_i in sentence:
                    sentence_event_i.append(" ".join(tm.clean_abstract(sentence_i, "lemma")))
                #----
                if search_mode == "stressor-event":
                    for event in list_events:
                        event_i = key.split("; ")
                        if event_i == event["name"]:
                            event_lemma = event["lemmaname"]
                    score_update, words_found_update, ind = lemma_rescoring(sentence_event_i, event_lemma)
                    val["multiscore"] = [score_update[x] for x in ind]
                    val["localisation"] = [val["localisation"][x] for x in ind]
                    val["sentence"] = [val["sentence"][x] for x in ind]
                    #val["words_found"] = [words_found_update[x] for x in ind]
                #----
                if search_mode == "event-event":
                    evs = key.split(" with ")
                    event_i = evs[0].split("; ")
                    event_j = evs[1].split("; ")
                    for event in list_events:
                        if event_i == event["name"]:
                            event_i_lemma = event["lemmaname"]
                        if event_j == event["name"]:
                            event_j_lemma = event["lemmaname"]
                    score_i_update, words_found_i_update, ind_i = lemma_rescoring(sentence_event_i, event_i_lemma)
                    score_j_update, words_found_j_update, ind_j = lemma_rescoring(sentence_event_i, event_j_lemma)
                    commun_index = list(set(ind_i) & set(ind_j))
                    score_i_update = [score_i_update[x] for x in commun_index]
                    score_j_update = [score_j_update[x] for x in commun_index]
                    val["multiscore"] = [score_i_update, score_j_update]
                    val["localisation"] = [val["localisation"][x] for x in commun_index]
                    val["sentence"] = [val["sentence"][x] for x in commun_index]
                    #words_found_i_update = [words_found_i_update[x] for x in commun_index]
                    #words_found_j_update = [words_found_j_update[x] for x in commun_index]
                    #val["words_found"] = [words_found_i_update, words_found_j_update]
                #----
                if not val["sentence"]:
                    del element["score"][key]
            if not element["score"]:
                del element["score"]



#-------------------------------------------------------------------------------

def lemma_rescoring(sentence_event, event_lemma):
    """
    Input: Sentence and Event (list)
    Output: new score and number found after recalculation. Index of sentences where a score was found
    """
    score_update = [] ; words_found_update = []
    for i in range(len(sentence_event)):
        score = [] ; words_found = []
        for synonym in event_lemma:
            synonym_split = list(set(synonym.split()))
            best, words = cscores.presence_score(sentence_event, synonym_split, i)
            if (abs(1-best/len(synonym_split)) <= 1 and best!=0):
                score.append(best/len(synonym_split))
                words_found.append(words)
        if score :
            best_score = score[score.index(max(score))]
            best_words = words_found[words_found.index(max(words_found))]
        else :
            best_score = None ; best_words = None ; sentences_update = None
        score_update.append(best_score)
        words_found_update.append(best_words)
    ind = [i for i, x in enumerate(score_update) if x is not None]
    return score_update, words_found_update, ind
