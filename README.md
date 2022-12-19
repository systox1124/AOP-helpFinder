authors: Florence Jornod - Université Paris Cité, T3S, INSERM UMR-S 1124, Paris, France
         Thomas Jaylet - Université Paris Cité, T3S, INSERM UMR-S 1124, Paris, France
	 Thibaut Coustillet - Université Paris Cité, T3S, INSERM UMR-S 1124, Paris, France
         Karine Audouze - Université Paris Cité, T3S, INSERM UMR-S 1124, Paris, France

contact: systox@parisdescartes.fr


AOP-helpFinder is provided without any warranty. But if you have any problem please feel free to contact us by e-mail.


# HOW TO CITE 
If you use AOP-helpFinder, please cite :

* Jornod et al. AOP-helpFinder webserver: a tool for comprehensive analysis of the literature to support adverse outcome pathways development. Bioinformatics. 2021 oct 30; doi: 10.1093/bioinformatics/btab750.
* XXXXXXXXXXXX. AOP-helpFinder 2.0 : ...

# WHAT IS AOP-helpFinder ? 

AOP-helpFinder is a tool developed to help AOP development (Jean-Charles Carvaillo: https://github.com/jecarvaill/aop-helpFinder)(Environ Health Perspect. 2019 Apr;127(4):47005).

It is based on text mining and parsing process on scientific abstracts. AOP-helpFinder allows the identification and extraction of associations between 'Prototypical stressor & event' and 'event & event' at various level of the biological organization (molecular initiating event, key event, and adverse outcome) through analysis of abstracts from the PubMed database. (https://pubmed.ncbi.nlm.nih.gov/). 


AOP-helpFinder was implemented under the H2020 Human Biomonintoring in Europe (HBM4EU) project, Work Package 13.1.
HBM4EU has received funding from the European Union’s H2020 research and innovation programme under grant agreement No 733032.

#  WHAT IS AOP-helpFinder web server 2.0 ? 

AOP-helpFinder 2.0 is a web server for identification and extraction of associations between 'prototypical stressor - event' and and 'event - event' associations at various level of the biological organization (molecular initiating event (MIE), key event (KE), and adverse outcome (AO)).
The web server AOP-helpFinder 2.0 proposes a new module to identify event-event linkage, and to prioritize them by a classification system, supporting the weight of evidence. For easier interpretation of the results, new visualization options are also available.

For more information please read :
* Jornod et al. AOP-helpFinder webserver: a tool for comprehensive analysis of the literature to support adverse outcome pathways development. Bioinformatics. 2021 oct 30; doi: 10.1093/bioinformatics/btab750.
* XXXXXXXXXXXX. AOP-helpFinder 2.0 : ...

# LICENCE 
This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can  use,  modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL

http://cecill.info/licences/Licence_CeCILL_V2.1-en.txt

AOP-helpFinder figures are governed by CC-by-nc license


# How To Use 

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

python3 aophelpfinder.py param.py 

python3 scoring.py --input (the file from aophelpfinder) --mode (event-event or stressor-event)

Some figures are also available with the different figures scripts

