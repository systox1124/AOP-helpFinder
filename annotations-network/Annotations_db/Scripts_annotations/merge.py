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


import pandas as pd
import argparse

#######-----------Parser-----------######
parser = argparse.ArgumentParser()
parser.add_argument('-sef', '--se_scoring', type=str, help='Path to the stressor_event file (TSV format)', required=True)
parser.add_argument('-eef', '--ee_scoring', type=str, help='Path to the event_event file (TSV format)', required=True)
parser.add_argument('-o', '--output', type=str, help='Path to write the merged output file (TSV format)', required=True)
args = parser.parse_args()

#######-----------Load Data-----------######
# Load the stressor_event and event_event TSV files
se_df = pd.read_csv(args.se_scoring, sep='\t')
ee_df = pd.read_csv(args.ee_scoring, sep='\t')

#######-----------Prepare Data-----------######
# Remove the 'Article' column from stressor_event file
se_df.drop(columns=['Article'], inplace=True)

# Rename columns in the event_event file to match the stressor_event file
ee_df.rename(columns={'Event 1': 'Stressor', 'Event 2': 'Event'}, inplace=True)

# Concatenate the two dataframes
merged_df = pd.concat([se_df, ee_df], ignore_index=True)

#######-----------Write Output-----------######
# Save the merged dataframe to the output TSV file
merged_df.to_csv(args.output, sep='\t', index=False)

print(f"Merged data has been written to {args.output}")
