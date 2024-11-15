# -*- coding: utf-8 -*-
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

# http://cecill.info/licences/Licence_CeCILL_V2.1-en.txt

import pandas as pd

#-------------------------------------------------------------------------------
def output_writing(dict_abstracts, output_name, research_mode, output_type):
    """
    Compute output according to the parameters provided in param.py
    """
    output_file = output_name + ".tsv"
    tsv_format(dict_abstracts, output_file, research_mode, output_type)
   # tsv_format_EWS(dict_abstracts, output_file, research_mode, output_type) # !!!!ONLY FOR stressor-event - created for EWS!!!!
    if "txt" in output_type:
        output_file = output_name + ".txt"
        txt_format(dict_abstracts, output_file, research_mode, output_type)
    return

#-------------------------------------------------------------------------------
def abst_info(abstract):
    """
    Retrieves the info from the pubmed abstract and puts NA if the info is not available
    """
    # Retrives abstract informations. "NA" if not
    dico_attribute = {
            'title' : "NA",
            'pmid' : "NA",
            'pubdate' : "NA",
            'abstract' : "NA"
            }
    if 'title' in abstract : dico_attribute['title'] = abstract['title']
    if 'pmid' in abstract : dico_attribute['pmid'] = abstract['pmid']
    if 'pubdate' in abstract : dico_attribute['pubdate'] = abstract['pubdate']
    if 'abstract' in abstract : dico_attribute['abstract'] = " ".join(abstract['abstract'].split('\n'))
    return dico_attribute


#Original version (work for event-event and stressor-event)
#-------------------------------------------------------------------------------
def tsv_format(dict_abstracts, output_file, research_mode, output_type):
    """
    Prepare output when the requested output type is "tsv" or "tsv_abstract"
    """
    #
    df = {
            'title' : [],
            'pmid' : [],
            'date' : [],
            'abstract' : [],
        }
    ### Research mode specificity
    if research_mode == "stressor-event":
        df["event"] = []
    if research_mode == "event-event":
        df["event_1"] = []
        df["event_2"] = []
    ### Data filling
    for abstract in dict_abstracts:
        if 'score' in abstract:
            abst_attribute = abst_info(abstract)
            for event, event_values in abstract['score'].items():
                df['title'].append(abst_attribute['title'])
                df['pmid'].append(abst_attribute['pmid'])
                df['date'].append(abst_attribute['pubdate'])
                df['abstract'].append(abst_attribute['abstract'])
                if research_mode == "stressor-event":
                    df["event"].append(event)
                if research_mode == "event-event":
                    co_event = event.split(" with ")
                    df["event_1"].append(co_event[0])
                    df["event_2"].append(co_event[1])
    ### Transform to pandas
    df_output = pd.DataFrame(df)
    ### Modality according to output type and research mode
    if research_mode == "stressor-event":
        if output_type == "tsv_abstracts" : df_output = df_output[["date", "title", "pmid", "event", "abstract"]]
        else : df_output = df_output[["date", "title", "pmid", "event"]]
    if research_mode == "event-event":
        if output_type == "tsv_abstracts" : df_output = df_output[["date", "title", "pmid", "event_1", "event_2", "abstract"]]
        else : df_output = df_output[["date", "title", "pmid", "event_1", "event_2"]]
    ### Writting
    print("\nOutput:\n")
    print(df_output)
    df_output.to_csv(output_file, sep = '\t', mode='a', index=False, header = True)
    return


#Original version (work for event-event and stressor-event)
#-------------------------------------------------------------------------------
def txt_format(dict_abstracts, output_file, research_mode, output_type):
    """
    Prepare output when the requested output type is "txt"
    """
    df = {
        'title' : [],
        'pmid' : [],
        'pubdate' : [],
        'abstract' : [],
        }
    ### Data filling
    with open(output_file, "w") as txt_output:
        for abstract in dict_abstracts:
            if 'score' in abstract:
                abst_attribute = abst_info(abstract)
                abst_title = abst_attribute['title']
                abst_date = abst_attribute['pubdate']
                abst_id = abst_attribute['pmid']
                abst = abst_attribute["abstract"]
                event_list = ""
                for event, event_values in abstract['score'].items():
                    ### Research mode specificity
                    if research_mode == "stressor-event" :
                        event_list = event_list + event + ", "
                    if research_mode == "event-event":
                        event_list = event_list + event + ", "
                info_abst = "{}\n{}\n{}\n{}\n{}\n\n --------------------------------------- \n\n".format(abst_date, abst_title, abst_id, event_list, abst)
                txt_output.write(info_abst)
    return



#Version for EWS + AOP-Wiki (only work for stressor-event)
#-------------------------------------------------------------------------------
def tsv_format_EWS(dict_abstracts, output_file, research_mode, output_type):
    """
    Prepare output when the requested output type is "tsv" or "tsv_abstract"
    """
    #
    df = {
            'title' : [],
            'pmid' : [],
            'pubdate' : [],
            'abstract' : [],
            'sentence' : [],
            'localisation' : [],
            'words_find' : [],
            'total_words' : [],
            'aophf_score' : []
        }
    df["event"] = []
    ### Data filling
    for abstract in dict_abstracts:
        if 'score' in abstract:
            abst_attribute = abst_info(abstract)
            for event, event_values in abstract['score'].items():
                nb_sentence = len(event_values['sentence'])
                for i in range(nb_sentence): #d
                    df['title'].append(abst_attribute['title'])
                    df['pmid'].append(abst_attribute['pmid'])
                    df['pubdate'].append(abst_attribute['pubdate'])
                    df['abstract'].append(abst_attribute['abstract'])
                    df['sentence'].append(event_values['sentence'][i])
                    df['localisation'].append(event_values['localisation'][i])
                    df['words_find'].append(event_values['words_found'][i])
                    df['total_words'].append(event_values['words_to_find'])
                    df['aophf_score'].append(event_values['multiscore'][i])
                    df["event"].append(event)
    ### Transform to pandas
    df_output = pd.DataFrame(df)
    ### Modality according to output type and research mode
    if output_type == "tsv_abstracts" : df_output = df_output[["pubdate", "title", "pmid", "event", "abstract","localisation", "sentence", "words_find", "total_words", "aophf_score"]]
    else : df_output = df_output[["pubdate", "title", "pmid", "event","localisation", "sentence", "words_find", "total_words", "aophf_score"]]
    ### Writting
    print("\nOutput:\n")
    print(df_output)
    # if(df_output.empty):
    #    noresultfile=open("no-results.txt","w")
    #    noresultfile.write("no results \n")
    # else:
    df_output.to_csv(output_file, sep = '\t', index=False, header = True)
        
    return
