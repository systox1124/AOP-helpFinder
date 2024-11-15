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

import argparse
import pandas as pd
import warnings
import os

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument("--mode", type=str, help="stressor-event or event-event")
parser.add_argument("--uid", type=str, help="path of the file with all uid")
parser.add_argument("--stats", type=str, help="path of the file with effective number of abstracts read by AOPhf")
parser.add_argument("--output", type=str, help="path of the output file (table)")
parser.add_argument("--stressor", type=str, help="name of the stressor")
args = parser.parse_args()

uid_file = args.uid
stats_file = args.stats
output_file = args.output
stressor_name = os.path.basename(output_file)
# stressor_name = args.stressor
mode = args.mode

###### uid
with open(uid_file) as id_file:
    n_uid = 0
    for uid in id_file:
        if(stressor_name == uid.split(" ")[0]):
            n_uid=uid.split(" ")[1]

###### traited abstract
with open(stats_file) as s_file:
    total_abs = 0
    for abs in s_file:
        total_abs += int(abs)

###### tsv output
o_df = pd.read_csv(output_file, sep='\t', header=0)
o_df = o_df.rename(columns=lambda x: x.strip())
print(o_df)
# count link total
nb_link = len(o_df["pmid"])
# count abstract total
nb_abst = len(o_df["pmid"].unique())
if mode == "stressor-event":
    #ev_links = o_df.groupby('event')['pmid'].count()
    ev_links = o_df.event.value_counts()
if mode =="event-event":
    ev_links = o_df.groupby(['event_1', 'event_2'])["pmid"].count()

###### Writing results
#stressor-event
if mode =="stressor-event":
    # Table stats
    resume_file = "resume_{}.txt".format(stressor_name)
    with open(resume_file, "w") as resume:
        resume.write("name of compound : {}.\n".format(stressor_name))
        resume.write("We found {} links between {} and your list of events. It correspond to {} abstracts and {} events found.\n".format(nb_link, stressor_name, nb_abst, len(ev_links)))
        resume.write("Please find here a resume of what we found:\n")
        df_for_stats = [{'stressor': stressor_name, 'nb_uid_extracted':str(n_uid)+"("+str(total_abs)+")", 'nb_abstract':nb_abst, 'nb_link':nb_link}]
        #df_for_stats = [{'stressor': stressor_name, 'nb_uid_extracted':n_uid, 'nb_traited_AOPhf':total_abs, 'nb_abstract':nb_abst, 'nb_link':nb_link}]

        df_for_stats= pd.DataFrame(data=df_for_stats)
        df_for_stats.to_csv('stats_stressor-event.tsv',sep='\t',index=False)
        # df which will be used for scoring
        df_for_scoring = pd.DataFrame(columns=['stressor','nb_traited_AOPhf','event','link'])
        for ev, link in ev_links.iteritems():
        #for ev, link in ev_links.items():
            new_row = {'stressor':stressor_name, 'nb_traited_AOPhf':total_abs, 'event':ev, 'link':link}
            df_for_scoring = df_for_scoring.append(new_row, ignore_index=True)
            resume.write("{}\t{}\n".format(ev,link))
        df_for_scoring.to_csv('prep_scoring_stressor-event.tsv',sep='\t',index=False) 
#event-event
if mode =="event-event":
    #Table stats
    df_for_stats = pd.DataFrame(columns=['event_1','event_2','nb_uid_extracted','nb_traited_AOPhf', 'nb_abstract', 'nb_link'])
    df_for_scoring = pd.DataFrame(columns=['event_1','event_2','link'])
    for ev, link in ev_links.iteritems():
        new_row_stats = [{'event_1': ev[0], 'event_2' : ev[1], 'nb_uid_extracted':n_uid, 'nb_traited_AOPhf':total_abs, 'nb_abstract':nb_abst, 'nb_link':nb_link}]
        new_row_scoring = [{'event_1': ev[0], 'event_2' : ev[1], 'link':link}]
        df_for_stats = df_for_stats.append(new_row_stats, ignore_index=True)
        df_for_scoring = df_for_scoring.append(new_row_scoring, ignore_index=True)
    df_for_stats.to_csv('stats_event-event.tsv',sep='\t', index=False)
    df_for_scoring.to_csv('prep_scoring_event-event.tsv',sep='\t', index=False)
