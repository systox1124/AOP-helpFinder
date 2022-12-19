#!/usr/bin/python3
# -*- coding: utf-8 -*-

######################################
# DESCRIPTION
######################################

# Text Mining Module
# Coded by Jean-Charles Carvaillo, Florence Jornod, Thomas Jaylet, Karine Audouze

#contact: systox@parisdescartes.fr


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
"""
Text mining module
"""

######################################
# IMPORT
######################################

# set the path of nltk_data
import nltk.data
nltk.data.path.append('./data/nltk_data/')

import nltk.corpus
import nltk.stem
from nltk.tokenize import sent_tokenize
from nltk import word_tokenize
import re
import string


import spacy
"""
Spacy :
    - Reference : https://github.com/explosion/spaCy
    - Installation : https://spacy.io/usage
"""

######################################
# FUNCTION
######################################

def remove_negative_sentences(abstract):
    """
    delete the sentences containing a negation word
    """
    sents = sent_tokenize(abstract)
    sents = [sent for sent in sents if not (bool(re.search('\d', sent) and
             'body weight' in sent))]
    abstract = [word_tokenize(sent) for sent in sents]
    negations = ['never', 'neither', 'no', 'none', 'nor', 'not', 'ain',
                     'aren', 'couldn', 'didn', 'doesn', 'hadn', 'hasn', 'haven',
                     'isn', 'mightn', 'mustn', 'needn', 'shan', 'shouldn', 'wasn',
                     'weren', 'won', 'wouldn']
    abstract = [sent for sent in abstract if not
                any(negation in sent for negation in negations)]
    for i in range(len(abstract)):
        abstract[i] = ' '.join(abstract[i])
    return abstract


def remove_stopwords(abstract):
    """
    remove stopwords ('a', 'the', etc..)
    """
    sents = sent_tokenize(abstract)
    sents = [sent for sent in sents if not (bool(re.search('\d', sent) and
             'body weight' in sent))]

    abstract = [word_tokenize(sent) for sent in sents]

    punctuation = list(string.punctuation)
    stop = nltk.corpus.stopwords.words('english') + punctuation + ['\'s']
    for i in range(len(abstract)):
        abstract[i] = [word for word in abstract[i] if word not in stop]
    for i in range(len(abstract)):
        abstract[i] = ' '.join(abstract[i])
    # Return
    return abstract



def clean_abstract(abstract, process):
    """
    - Input: a text (ex abstract) in the form of a string AND the pre-processing method ("lemma" OR "stem")
    - Steps: cut the text into sentences, cut the sentences into words, delete the sentences containing a negation word, delete the stopwords, then perform the pre-processing (lemmatization or stemming)
    - Output: text (abstract) cleaned as a list containing the sentences
    """
    # 1. split abstract by sentences
    sents = sent_tokenize(abstract)

    sents = [sent for sent in sents if not (bool(re.search('\d', sent) and
             'body weight' in sent))]

    # 2. split sentences by words
    abstract = [word_tokenize(sent) for sent in sents]

    # 3. remove sentences which contain a negation

    negations = ['never', 'neither', 'no', 'none', 'nor', 'not', 'ain',
                 'aren', 'couldn', 'didn', 'doesn', 'hadn', 'hasn', 'haven',
                 'isn', 'mightn', 'mustn', 'needn', 'shan', 'shouldn', 'wasn',
                 'weren', 'won', 'wouldn']

    abstract = [sent for sent in abstract if not
                any(negation in sent for negation in negations)]

    #4. remove Stopwords
    punctuation = list(string.punctuation)
    stop = nltk.corpus.stopwords.words('english') + punctuation + ['\'s']

    for i in range(len(abstract)):
        abstract[i] = [word for word in abstract[i] if word not in stop]

    #5. Preprocessing
    for i in range(len(abstract)):
        if process == "stem" :
            abstract[i] = stem_process(abstract[i])
        elif process == "lemma" :
            abstract[i] = lemma_process(abstract[i])
        abstract[i] = ' '.join(abstract[i])

    # Return
    return abstract


def stem_process(words):
    """Method to stem each words in a list.
    Return:
        words (list): a list which contains stemmed words.
    """
    # snowball stemmers developed by Martin Poter
    sno = nltk.stem.SnowballStemmer('english')
    # search if a particular word is in words
    for i in range(len(words)):
        words[i] = sno.stem(words[i])
    return words


def lemma_process(words):
    """
    Method to lemm each words in a list.
    Return: words in canonical form
    """
    nlp = spacy.load("en_core_web_sm")
    words = ' '.join(words)
    doc = nlp(words.lower())
    lemma_list = []
    for token in doc:
        lemma_list.append(token.lemma_)
    return lemma_list




######################################
# Main()
######################################

if __name__ == '__main__':
    abstract = "BACKGROUND Bisphenol S (BPS) was introduced in the market as a potentially safer alternative to bisphenol A (BPA). However, there are limited studies on health effects of BPS and no epidemiologic studies on its relationship with male reproductive health outcomes, specifically semen quality.   OBJECTIVE To investigate predictors of urinary BPS concentrations and its association with semen parameters among men attending a fertility center.   METHODS This cross-sectional analysis included 158 men of couples seeking fertility treatment (2011-2017) contributing 338 paired semen and urine samples. At the time of sample collection, men completed a questionnaire on self-reported use of household products and food intake within the previous 24 h. Urinary concentrations of BPA, BPS and bisphenol F were quantified using isotope-dilution tandem mass spectrometry. Semen samples were analyzed following WHO guidelines. Multivariable mixed models were used to investigate predictors of urinary BPS concentrations and to evaluate associations between urinary BPS concentrations and semen parameters, using random intercept to account for correlation in outcomes across multiple observations per man and adjusting for abstinence time, specific gravity, age, body mass index (BMI), year of sample collection and BPA concentrations. Analyses were also stratified by BMI (≥25 vs <25 kg/m2).   RESULTS Median (IQR) urinary BPS concentration was 0.30 (0.20, 0.90) μg/L, and 76% of samples had detectable (>0.1 μg/L) concentrations. Self-reported fabric softener and paint/solvent use as well as intake of beef and cheese within 24 h before urine collection were positively associated with BPS concentrations. Men with higher BPS concentrations also had significantly higher BMI. Lower semen parameters were found among men with detectable BPS concentrations, compared to men with non-detectable BPS [2.66 vs. 2.91 mL for volume (p = 0.03), 30.7 vs. 38.3 mil/mL for concentration (p = 0.03), 76.8 vs. 90.0 mil for total count (p = 0.09), 43.7 vs. 47.0% for motility (p = 0.06), and 5.42 vs. 6.77% for morphologically normal sperm (p = 0.24)]. Some associations of BPS with lower semen parameters were only found among men with a BMI ≥ 25 kg/m2.   CONCLUSIONS We identified dietary and lifestyle factors associated with BPS exposure, suggesting potential avenues for reducing exposures. We also observed negative associations between BPS and semen parameters, especially among overweight and obesity men."
    abstract1 = clean_abstract(abstract, "lemma")
    print("\n", abstract1)
    abstract2 = clean_abstract(abstract, "stem")
    print("\n", abstract2)

    abstract2 = remove_negative_sentences(abstract)
    print("\n", abstract2)

    print("\n", sent_tokenize(abstract))
