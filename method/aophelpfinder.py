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
import sys
import re
import os
import get_abstracts as getabs
import compute_scores as cscores
import compute_output as cout
import create_collection_event as ev
import results_filter as rf
import tm_module as tm

import importlib.util
paramfile = sys.argv[1]
spec = importlib.util.spec_from_file_location("param", paramfile)
param = importlib.util.module_from_spec(spec)
sys.modules["param"] = param
spec.loader.exec_module(param)

#-------------------------------------------

chem_name = param.abstracts_file.split(".")[0].replace("-"," ")
try:
    dict_abstracts = getabs.get_abstracts_from_file(param.abstracts_file, param.abstracts_type)
except:
    dict_abstracts = {}
    stats_file =open(param.stats_file,"a")
    stats_file.write(str(len(dict_abstracts))+"\n")
    cout.output_writing(dict_abstracts, param.output_name, param.search_mode, param.outputtype)
    exit()

list_events = ev.create_collection_event(param.events_file, param.events_type)
print("\nListe d'events:\n",list_events)
print(len(dict_abstracts))

#-------------------------------------------

stats_file =open(param.stats_file,"a")
stats_file.write(str(len(dict_abstracts))+"\n")
#stats_file.write(chem_name + "\t" + str(len(dict_abstracts)))

#-------------------------------------------
for abstract in dict_abstracts:
    #print(abstract)
    score = cscores.compute_scores(abstract,list_events,param.intro, param.search_mode)
    if score:
        print("score")
        abstract['score']=score
#print(dict_abstracts)
#-------------------------------------------
if param.res_filter is True :
    rf.refinement_filter(dict_abstracts, list_events, param.search_mode)

#-------------------------------------------
cout.output_writing(dict_abstracts, param.output_name, param.search_mode, param.outputtype)
