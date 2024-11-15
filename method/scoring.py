# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
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

from Bio import Entrez
import os
import scipy.stats as stats
import pandas as pd
import numpy as np
import scipy.stats as stats
from pyexcel_ods3 import save_data
from pyexcel_ods3 import get_data
from collections import OrderedDict
import argparse
import time
Entrez.email = "example@example.org"
#nb_pubmed_article =  33000000
alpha = 0.05
Retries = 4


query = '("1750"[Date - Modification] : "3000"[Date - Modification])'
handle = Entrez.esearch(db="pubmed", term=query)
record = Entrez.read(handle)
handle.close()
nb_pubmed_article = int(record["Count"])
print("Article in PubMed:", int(record["Count"]))


parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="path of the input")
parser.add_argument("--mode", type=str, help="stressor-event or event-event")
args = parser.parse_args()

input = args.input
mode = args.mode


"""
- "stressor-event" : Calculates the strength of a link between a stressor and an event.
- "event-event" : Calculates the strength of a link between a pair of events.
In both cases, a Fisher exact test is applied to calculate the score.
This score makes it possible to identify whether a link is significant or not,
i.e. whether it is well documented in the literature or not.
/!\ If the p-value is not significant, it does not necessarily mean that the link does not exist,
but simply that it is very little documented/studied. It will therefore be necessary to check the results manually
"""

# Fisher's exact test is used to calculate (assuming equivalence between the 2 groups)
# the exact probability of obtaining, between the 2 groups, a difference in the
# distribution of the 2 modalities greater than or equal to that observed (p unilateral)


#################################### STRESSOR-EVENT ##############################################
if mode == "stressor-event":
    print("Scoring stressor-event in progress")
    """
    Calculates the strength of a link between a stressor and an event.
    It calculates the probability of obtaining, between two groups (PubMed vs Stressor) a difference
    in the distribution of the 2 modalities (event; non-event) greater than or equal to that observed
    """
    res_ods = []
    res_ods.append(["Stressor","Event","Article","Link","Pvalue","pvar","Confidence"])

    input_read=pd.read_csv(input,sep='\t')
    dico_article_stressor = {}
    for index, row in input_read.iterrows():
        stressor=row["stressor"].replace('_',' ')
        if stressor not in dico_article_stressor:
            dico_article_stressor[stressor] = {}
            dico_article_stressor[stressor] = {"article" : row["nb_traited_AOPhf"]}
        if "link" not in dico_article_stressor[stressor]:
            dico_article_stressor[stressor]["link"] = []
        dico_article_stressor[stressor]["link"].append({ row["event"] : int(row["link"])})


    #print(dico_article_stressor)
    # Calculation of the score: fisher's exact test for each event-stressor pair.
    # Comparison with PubMed distribution
    for stressor in dico_article_stressor:
        nb_stressor_article = dico_article_stressor[stressor]["article"]
        query = stressor
        for event in dico_article_stressor[stressor]["link"]:
            ### PubMed request
            ev_name, link = list(event.items())[0]
            query = ev_name.replace("|", "OR")
            #handle = Entrez.esearch(db="pubmed", term=query.replace(" ","+"))
            for x in range(Retries):
                try:
                    handle = Entrez.esearch(db="pubmed", term=query)
                    record = Entrez.read(handle)
                    handle.close()
                except Exception as e:
                    print("Error, retry request...")
                    time.sleep(5)
                else:
                    break

            ### Creation of contingence table
            nb_event_stressor_AOPhf = int(link)
            nb_non_event_stressor = int(nb_stressor_article) - int(nb_event_stressor_AOPhf)
            if (int(nb_event_stressor_AOPhf) > int(nb_stressor_article)):
                nb_non_event_stressor = 0
            nb_event_pubmed = int(record["Count"]) - nb_event_stressor_AOPhf #Tous les articles de Pubmed traitant de l'event sans compter ceux qui ont été retrouvé par AOPhf pour le stresseur d'interet (independance des donnees)
            nb_non_event_pubmed = int(int(nb_pubmed_article)-int(nb_stressor_article)) - int(nb_event_pubmed) #Tout le reste de pubmed (pubmed - stressor) qui ne traite pas de l'event d'interet (-event)
            if nb_event_pubmed < 0:
                nb_event_pubmed = 0
            ### Calculation Fisher 1
            dico_stressor = {}
            df = pd.DataFrame({'PubMed':[nb_non_event_pubmed, nb_event_pubmed], stressor:[nb_non_event_stressor, nb_event_stressor_AOPhf]}, index=pd.Index(['Non {}'.format(ev_name), '{}'.format(ev_name)]))
            oddsratio, pvalue = stats.fisher_exact(df.to_numpy(), alternative="greater")
            dico_stressor["fisher_table"] = df
            dico_stressor["pvalue"] = format(pvalue, ".3e")
            print(dico_stressor["fisher_table"])
            dico_article_stressor[stressor]["first_test"] = dico_stressor

            ### Calculation Fisher 2
            dico_stressor_2 = {}
            if (pvalue > alpha):
                dico_stressor_2["fisher_table"] = "NA"
                dico_stressor_2["pvalue"] = "NA"
                dico_article_stressor[stressor]["p_var"] = "NA"
            else :
                df2 = pd.DataFrame({'PubMed':[nb_non_event_pubmed, nb_event_pubmed], "stressor":[nb_non_event_stressor+1, nb_event_stressor_AOPhf-1]})
                oddsratio2, pvalue2 = stats.fisher_exact(df2.to_numpy(), alternative="greater")
                dico_stressor_2["fisher_table"] = df2
                dico_stressor_2["pvalue"] = pvalue2
                dico_article_stressor[stressor]["p_var"] = format((pvalue2 - pvalue), ".3e")
            dico_article_stressor[stressor]["second_test"] = dico_stressor_2

            ### Save results
            res_liste = []
            res_liste.append(stressor)
            res_liste.append(ev_name)
            res_liste.append(int(dico_article_stressor[stressor]["article"]))
            res_liste.append(link)
            res_liste.append(dico_article_stressor[stressor]["first_test"]["pvalue"])
            res_liste.append(dico_article_stressor[stressor]["p_var"])

            if dico_article_stressor[stressor]["p_var"] == "NA":
                if link > 999:
                    res_liste.append("Very High") 
                elif link > 99:
                    res_liste.append("High")
                else:
                    res_liste.append("Low")
            elif float(dico_article_stressor[stressor]["p_var"]) >= 5e-02:
                if link > 999:
                    res_liste.append("Very High") 
                elif link > 99:
                    res_liste.append("High")
                else:
                    res_liste.append("Quite low")
            elif 1e-03 <= float(dico_article_stressor[stressor]["p_var"]) < 5e-02:
                if link > 999:
                    res_liste.append("Very High") 
                elif link > 99:
                    res_liste.append("High")
                else:
                    res_liste.append("Moderate")
            elif 1e-05 <= float(dico_article_stressor[stressor]["p_var"]) < 1e-03:
                if link > 999:
                    res_liste.append("Very High") 
                else:
                    res_liste.append("High")
            elif float(dico_article_stressor[stressor]["p_var"]) < 1e-05:
                res_liste.append("Very High")

            """
            if (dico_article_stressor[stressor]["p_var"] == "NA"):
                res_liste.append("Low")
            elif float(dico_article_stressor[stressor]["p_var"]) >= 5e-02:
                res_liste.append("Quite low")
            elif 1e-03 <= float(dico_article_stressor[stressor]["p_var"]) < 5e-02:
                res_liste.append("Moderate")
            elif 1e-05 <= float(dico_article_stressor[stressor]["p_var"]) < 1e-03:
                res_liste.append("High")
            elif float(dico_article_stressor[stressor]["p_var"]) < 1e-05:
                res_liste.append("Very High")
            """

            print(res_liste)
            res_ods.append(res_liste)


    ### Writing results
    data = OrderedDict()
    data.update({"Stressor-Event": res_ods})
    save_data("scoring_stressor-event.ods", data)

    with open("scoring_stressor-event.tsv", "w") as output_csv:
        for element in res_ods:
            output_csv.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(element[0],element[1],element[2],element[3],element[4],element[5],element[6]))


################################### EVENT-EVENT #############################################

if mode == "event-event":
    """
    Calculates the strength of a link between a pair of events.
    Calculates the probability of obtaining, between the 2 groups (non event1; event1),
    a difference in distribution of the 2 modalities (non event2; event2) greater than or equal to that observed
    """
    # Open the stats file to retrieve the different Event-Event information of AOP-helpfinder ; Fisher test
    print("\nScoring event-event in progress\n")

    res_ods = []
    res_ods.append(["Event 1","Event 2","Link","Pvalue","pvar", "Confidence"])

    input_read=pd.read_csv(input,sep='\t')
    dico_article_event = {}
    for index, row in input_read.iterrows():
        dico_event_pair = {}
        Event_1 = row["event_1"]
        Event_2 = row["event_2"]
        ### PubMed query
        query_ev1 = Event_1.replace("|" , "OR")
        for x in range(Retries):
            try:
                handle_ev1 = Entrez.esearch(db="pubmed", term=query_ev1)
                record_ev1 = Entrez.read(handle_ev1)
                handle_ev1.close()
            except Exception as e:
                print("Error, retry request event 1...")
                time.sleep(5)
            else :
                break
        query_ev2 = Event_2.replace("|" , "OR")
        for x in range(Retries):
            try:
                handle_ev2 = Entrez.esearch(db="pubmed", term=query_ev2)
                record_ev2 = Entrez.read(handle_ev2)
                handle_ev2.close()
            except Exception as e:
                print("Error, retry request event 2...")
                time.sleep(5)
            else :
                break

        ### Creation of contigence table
        nb_ev1 = int(record_ev1["Count"])
        nb_ev2 = int(record_ev2["Count"])

        nb_ev1_ev2 = int(row["link"])
        nb_ev1_noev2 = int(nb_ev1) - int(nb_ev1_ev2)

        if (int(nb_ev1_ev2) > int(nb_ev1)):
            nb_ev1_noev2 = 0
        nb_pubmed_noev1_ev2 = int(int(nb_pubmed_article)-int(nb_ev1)) #Tous les articles
        nb_pubmed_noev1_noev2 = int(nb_pubmed_noev1_ev2) - int(nb_ev2)
        nb_pubmed_noev1_ev2 = int(nb_ev2) - int(nb_ev1_ev2) #Tous les articles qui parlent de l'ev2 sur PubMed moins ceux retrouvés par AOPhf pour pas prendre en compte les articles de l'ev1
        nb_pubmed_noev1_noev2 = int((nb_pubmed_article) - int(nb_ev1)) - int(nb_pubmed_noev1_ev2) #Tous les articles qui ne parlent pas de l'ev2 sur (pubmed - ev1)
        if nb_pubmed_noev1_ev2 < 0 :
            nb_pubmed_noev1_ev2 = 0

        Ev1_Ev2 = "{} AND {}".format(Event_1, Event_2)
        Ev1_noEv2 = "{} NOT {}".format(Event_1, Event_2)
        noEv1_Ev2 = "NOT {} AND {}".format(Event_1, Event_2)
        noEv1_noEv2 = "NOT {} NOT {}".format(Event_1, Event_2)

        dico_event_pair[Ev1_Ev2] = nb_ev1_ev2
        dico_event_pair[Ev1_noEv2] = nb_ev1_noev2
        dico_event_pair[noEv1_Ev2] = nb_pubmed_noev1_ev2
        dico_event_pair[noEv1_noEv2] = nb_pubmed_noev1_noev2

        pair = "{} ~~ {}".format(Event_1, Event_2)
        dico_article_event[pair] = {}

        ### Calculation Fisher 1
        df = pd.DataFrame({"Non {}".format(Event_1) : [nb_pubmed_noev1_noev2, nb_pubmed_noev1_ev2], Event_1 : [nb_ev1_noev2, nb_ev1_ev2]}, index = pd.Index(["Non {}".format(Event_2), Event_2]))
        oddsratio, pvalue = stats.fisher_exact(df.to_numpy(), alternative="greater")
        dico_event = {}
        dico_event["fisher_table"] = df
        dico_event["pvalue"] = format(pvalue, ".3e")
        dico_article_event[pair]["first_test"] = dico_event


        ### Calculation Fisher 2
        dico_event_2 = {}
        if (pvalue > alpha):
            dico_event_2["fisher_table"] = "NA"
            dico_event_2["pvalue"] = "NA"
            dico_article_event[pair]["p_var"] = "NA"

        else :
            df2 = pd.DataFrame({"Non {}".format(Event_1) : [nb_pubmed_noev1_noev2, nb_pubmed_noev1_ev2], Event_1 : [nb_ev1_noev2+1, nb_ev1_ev2-1]}, index = pd.Index(["Non {}".format(Event_2), Event_2]))
            oddsratio2, pvalue2 = stats.fisher_exact(df2.to_numpy(), alternative="greater")
            dico_event_2["fisher_table"] = df2
            dico_event_2["pvalue"] = format(pvalue2, ".3e")
            dico_article_event[pair]["p_var"] = format((pvalue2 - pvalue), ".3e")

        dico_article_event[pair]["second_test"] = dico_event_2


        ###Save results
        res_liste = []
        res_liste.append(Event_1)
        res_liste.append(Event_2)
        res_liste.append(nb_ev1_ev2)
        res_liste.append(str(dico_article_event[pair]["first_test"]["pvalue"]))
        res_liste.append(str(dico_article_event[pair]["p_var"]))
        
        """
        if dico_article_event[pair]["p_var"] == "NA":
            res_liste.append("Low")
        elif float(dico_article_event[pair]["p_var"]) >= 5e-02:
            res_liste.append("Quite low")
        elif 1e-03 <= float(dico_article_event[pair]["p_var"]) < 5e-02:
            res_liste.append("Moderate")
        elif 1e-05 <= float(dico_article_event[pair]["p_var"]) < 1e-03:
            res_liste.append("High")
        elif float(dico_article_event[pair]["p_var"]) < 1e-05:
            res_liste.append("Very High")
        """
         
        if dico_article_event[pair]["p_var"] == "NA":
            if nb_ev1_ev2 > 999:
                res_liste.append("Very High") 
            elif nb_ev1_ev2 > 99:
                res_liste.append("High")
            else:
                res_liste.append("Low")
        elif float(dico_article_event[pair]["p_var"]) >= 5e-02:
            if nb_ev1_ev2 > 999:
                res_liste.append("Very High") 
            elif nb_ev1_ev2 > 99:
                res_liste.append("High")
            else:
                res_liste.append("Quite low")
        elif 1e-03 <= float(dico_article_event[pair]["p_var"]) < 5e-02:
            if nb_ev1_ev2 > 999:
                res_liste.append("Very High") 
            elif nb_ev1_ev2 > 99:
                res_liste.append("High")
            else:
                res_liste.append("Moderate")
        elif 1e-05 <= float(dico_article_event[pair]["p_var"]) < 1e-03:
            if nb_ev1_ev2 > 999:
                res_liste.append("Very High") 
            else:
                res_liste.append("High")
        elif float(dico_article_event[pair]["p_var"]) < 1e-05:
            res_liste.append("Very High")

        res_ods.append(res_liste)


    ### Writing results
    data = OrderedDict()
    data.update({"Event-Event": res_ods})
    save_data("scoring_event-event.ods", data)

    with open("scoring_event-event.tsv", "w") as output_csv:
        for element in res_ods:
            output_csv.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(element[0],element[1],element[2],element[3],element[4],element[5]))

