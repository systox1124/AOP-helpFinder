# AOP-helpFinder 3.0

## Authors:   
Thomas Jaylet - Université Paris Cité, T3S, INSERM UMR-S 1124, Paris, France  
Florence Jornod - Université Paris Cité, T3S, INSERM UMR-S 1124, Paris, France 
Quentin Capdet - Université Paris Cité, T3S, INSERM UMR-S 1124, Paris, France
Thibaut Coustillet - Université Paris Cité, T3S, INSERM UMR-S 1124, Paris, France  
Karine Audouze - Université Paris Cité, T3S, INSERM UMR-S 1124, Paris, France  

contact: systox@parisdescartes.fr

AOP-helpFinder 3.0 is provided without any warranty. But if you have any problem please feel free to contact us by e-mail.


## HOW TO CITE 
If you use AOP-helpFinder V3.0, please cite :
BIORXIV/2025/648318. AOP-helpFinder 3.0: from text mining to network visualization of key event relationships, and knowledge integration from multi-sources. Thomas Jaylet, Florence Jornod, Quentin Capdet, Olivier Armant, and Karine Audouze*

*Corresponding author: Karine Audouze (Email: karine.audouze [at] u-paris.fr)

*for the method(s)*:  
* Carvaillo, J. C., Barouki, R., Coumoul, X., & Audouze, K. (2019). Linking Bisphenol S to Adverse Outcome Pathways Using a Combined Text Mining and Systems Biology Approach. Environmental health perspectives, 127(4), 47005. https://doi.org/10.1289/EHP4200  
&nbsp;
* Jaylet, T., Coustillet, T., Jornod, F., Margaritte-Jeannin, P., & Audouze, K. (2023). AOP-helpFinder 2.0: Integration of an event-event searches module. Environment international, 177, 108017. https://doi.org/10.1016/j.envint.2023.108017

*for the webserver*:
* Jornod, F., Jaylet, T., Blaha, L., Sarigiannis, D., Tamisier, L., & Audouze, K. (2022). AOP-helpFinder webserver: a tool for comprehensive analysis of the literature to support adverse outcome pathways development. Bioinformatics (Oxford, England), 38(4), 1173–1175. https://doi.org/10.1093/bioinformatics/btab750
* Jaylet, T., Coustillet, T., Jornod, F., Margaritte-Jeannin, P., & Audouze, K. (2023). AOP-helpFinder 2.0: Integration of an event-event searches module. Environment international, 177, 108017. https://doi.org/10.1016/j.envint.2023.108017


## WHAT IS AOP-helpFinder ? 

AOP-helpFinder is a tool developed to help AOP development (Jean-Charles Carvaillo: https://github.com/systox1124/aop-helpFinder) (Environ Health Perspect. 2019 Apr;127(4):47005).

It is based on text mining and parsing process on scientific abstracts. AOP-helpFinder allows the identification and extraction of associations between 'Prototypical stressor & event' and 'event & event' at various level of the biological organization (molecular initiating event, key event, and adverse outcome) through analysis of abstracts from the PubMed database (https://pubmed.ncbi.nlm.nih.gov/). 

AOP-helpFinder 2.0 allow the identification and extraction of associations between 'prototypical stressor - event' and 'event - event' associations at various level of the biological organization (Molecular Initiating Event (MIE), Key Event (KE), and Adverse Outcome (AO)).  
&nbsp;
The version 2.0 proposes:
1. A new module to identify event-event linkage, 
2. A classification system (Confidence Score) to prioritize the relationships in order to support the weight of evidence, 
3. New visualization options for easier interpretation of the results,
4. Interactive tables to be forwarded to the PubMed page of articles hosting links for better data curation.

AOP-helpFinder was implemented under the European HBM4EU, OBERON, RADONORM and PARC projects
The HBM4EU, OBERON and RadonNorm projects have received funding from the European Union's Horizon 2020 research and innovation programme, under grant agreement No 733032, 825712, and 900009. The RadoNorm project from the Euratom programme 2019-2020 under grant agreement No 733032, 825712, 101057014, 900009. The PARC projet has received funding from the European Union's Horizon Europe framework under the gratn agreement 101057014.


AOP-helpFinder is freely available at : http://aop-helpfinder-v3.u-paris-sciences.fr/



##  WHAT IS AOP-helpFinder 3.0 ? 
With the new version of AOP-helpFinder you can visualize you results as an annotated network. You will have the annotations from 8 databases : AOP-wiki, HPA, KEGG, Wikipathways, REACTOME, uniprot, DisGeNET, DISEASES.


## LICENCE 
This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can  use,  modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL

http://cecill.info/licences/Licence_CeCILL_V2.1-en.txt

AOP-helpFinder figures are governed by CC-by-nc license


## How To Use 

require for AOP-helpFinder
Package Python:
* nltk
* spacy
* pandas
* scipy
* numpy
* pyexcel_ods3
* pubmed_parser (http://github.com/titipata/pubmed_parser)

You need to complete the parameters with completed the param.py file  
You will find example file in the example_file folder

`$ python3 aophelpfinder.py param.py` 

`$ python3 scoring.py --input (the file from aophelpfinder) --mode (event-event or stressor-event)`

Some figures are also available with the different figures scripts.
