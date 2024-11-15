#!/bin/bash


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
START_TIME=$(date +%s)
 
db_dir=$1
scripts_annotations_dir=$2
output_dir=$3
scoring=$4
mode=$5


if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
    echo "Répertoire $output_dir créé."
else
    echo "Répertoire $output_dir existe déjà."
fi

#liste de résultats d'AOP-hf

#se "stressor - event": 
#scoring_stressor-event.tsv - scoring_stressor-event2.tsv - scoring_stressor-event_lung_radon.tsv#

#ee "event - event":
#scoring_event-event.tsv
 
 
#AOP-Wiki
python3 "$scripts_annotations_dir/Module_enrichi_AOPWiki.py" -hf $scoring -db $db_dir/raw_AOPwiki.tsv -syn $db_dir/AOPWiki_synonym.tsv -o $output_dir/annotations_AOPWiki.tsv -net $output_dir/annotations_AOPWiki_ntwrk.tsv -mod $mode 

#HPA
python3 "$scripts_annotations_dir/Module_enrichi_HPA.py" -hf $scoring -db $db_dir/proteinatlas.tsv -o $output_dir/annotations_HPA_Uniprot.tsv -o_net1 $output_dir/annotations_HPA_tissue_info.tsv -o_net2 $output_dir/annotations_Uniprot_bio_process.tsv -mod $mode

#DISEASES
python3 "$scripts_annotations_dir/Module_enrichi_DISEAS.py" -hf $scoring -db $db_dir/human_disease_knowledge_filtered.tsv -syn $db_dir/proteinatlas.tsv -o $output_dir/annotations_DISEASES_jensenlab.tsv -mod $mode

#DisGeNET
python3 "$scripts_annotations_dir/Module_enrichi_DisGeNet.py" -hf $scoring -db $db_dir/DisGenNet-curated_gene_disease_associations.tsv -syn $db_dir/proteinatlas.tsv -o $output_dir/annotations_DisGeNET_2015.tsv -mod $mode

#KEGG
python3 "$scripts_annotations_dir/Module_enrichi_KEGG.py" -hf $scoring -db $db_dir/2024_04_12_kegg_hsa.gmt -syn $db_dir/proteinatlas.tsv -o $output_dir/annotations_KEGG.tsv -mod $mode

#Reactome
python3 "$scripts_annotations_dir/Module_enrichi_Reactome.py" -hf $scoring -db $db_dir/ReactomePathways.gmt -syn $db_dir/proteinatlas.tsv -o $output_dir/annotations_Reactome.tsv -mod $mode

#WikiPathways
python3 "$scripts_annotations_dir/Module_enrichi_WikiPathways.py" -hf $scoring -db $db_dir/converted_data_WikiPathways.gmt -syn $db_dir/proteinatlas.tsv -o $output_dir/annotations_WikiPathways.tsv -mod $mode

END_TIME=$(date +%s)

EXECUTION_TIME=$((END_TIME - START_TIME))

echo -e "\033[0;32mLe module d'annotation a mis $EXECUTION_TIME secondes à s'exécuter.\033[0m"
