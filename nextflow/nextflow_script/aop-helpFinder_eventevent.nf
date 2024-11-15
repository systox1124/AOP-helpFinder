#!/usr/bin/env nextflow
nextflow.enable.dsl=1 
params.event = "event-event.txt"
params.ncbiKey=""

// #authors: Florence Jornod - INSERM UMRS 1124
// #         Thomas Jaylet - Université de Paris - France
// #         Karine Audouze - Université de Paris - France

// #contact: systox@paris-descartes.fr


// #AOP-helpFinder is provided without any warranty. But if you have any probleme please feel free to contact us by mail.

// #------- WHAT IS AOPHELPFINDER? -------------

// #AOP-helpFinder is a tool developed to help AOP development (Jean-Charles Carvaillo: https://github.com/jecarvaill/aop-helpFinder)(Environ Health Perspect. 2019 Apr;127(4):47005).

// #It is based on text mining and parsing process on scientific abstracts. AOP-helpFinder identify links between stressors and molecular initiating event, key events and adverse outcomes through abstracts from the PubMed database (https://pubmed.ncbi.nlm.nih.gov/).

// #AOP-helpFinder was implemented under the H2020 Human Biomonintoring in Europe (HBM4EU) project, Work Package 13.1.
// #HBM4EU has received funding from the European Union’s H2020 research and innovation programme under grant agreement No 733032.

// #------- LICENCE ---------------------------

// #This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can  use,  modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL

// #http://cecill.info/licences/Licence_CeCILL_V2.1-en.txt

// #-------------


event = file(params.event)
RESULTSDIR = file(params.resultspath)
RESULTSDIR_AOPhFout = file("${RESULTSDIR}/raw_AOP-helpFinder/")
RESULTSDIR_scoring = file("${RESULTSDIR}/raw_scoring/")
RESULTSDIR_stats = file("${RESULTSDIR}/resume/")
paramUser = params.userParam

process cleanEventFile{
	input:
	file(event) from event.splitText()

	output:
	file("eventCleaned.txt") into eventCleaned

	script:
	"""
	sed '/^ *\$/d' ${event} > eventCleaned.txt
	"""
}

process extractPubmedUidEvent {
	label "extractPubmedUidEvent"
	maxForks 5
	
	input:
	val(event) from eventCleaned.splitText()
	
	output:
	file("event-uid.txt") into uids
	
	script:
	"""
	
	#update PATH 
	PATH=/PATH/TO/perl5/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/PATH/TO/edirect
	#update our PUBMED DataBase

	archive-pubmed -path /mnt/ssd

	EDIRECT_PUBMED_MASTER=/mnt/ssd

	export NCBI_API_KEY=${params.ncbiKey}
	esearch -db pubmed -query "${event.trim().replaceAll(/\t/," OR ")}" | efetch -format uid -db pubmed > event-uid.txt
	"""
}

process catUid{
    publishDir "${RESULTSDIR}/uid", mode: "copy"
    
    input:
    file(uidFile) from uids.collectFile()

    output:
    file('event-uid_uniq.txt') into uniqUid, uids2

    script:
    """
    cat ${uidFile} | sort | uniq > event-uid_uniq.txt
    """
}
process countUids{
	input:
	file(uidFile) from uids2

	output:
	file('stats_uid.txt') into statsUid

	script:
	"""
	nbLine=\$(wc -l ${uidFile})
	echo "eventevent \$nbLine" > stats_uid.txt
	"""
}

statsUid.collectFile(name: 'stats_uid.txt').into{statUidFile;statUidFile2}

process downloadAbstractXML{
	label 'DLXML'
	
	input:
	file(uid) from uniqUid.splitText(by:750)
	
	output:
	file("abstract.xml") into xml
	
	script:
	"""
	PATH=/PATH/TO/perl5/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/PATH/TO/edirect
	EDIRECT_PUBMED_MASTER=/mnt/ssd
	export NCBI_API_KEY=${params.ncbiKey}
	cat ${uid} | /home/fjornod/edirect/fetch-pubmed -path '/mnt/ssd/Archive'  | tr -dc '[:alnum:]\\n\\r<>/ !?"=.,-' | sed '3,\${/<?xml/d}' | sed '3,\${/<!DOCTYPE/d}'  > abstract.xml
	"""	
}

process generateParamsFile{
	label 'paramFile'

	input:
	file(abstractFile) from xml
	val(parameters) from paramUser
	file(eventFile) from event

	output:
	tuple file('abstract.xml'), file('params.py') into dataAOPhF

	shell:
	'''
	IFS=$';' read -r user_id events_type output_type intro cleaning Fig search_mode word_number user_id2 encours  <<< "!{parameters}"
	echo "abstracts_file = '!{abstractFile}'" > params.py
	echo "abstracts_type = 'xml'" >> params.py
	echo "events_file = '!{event}'" >> params.py
	echo "events_type = '$events_type'" >> params.py
	echo "output_name = 'output-eventevent'" >> params.py
	echo "outputtype = '${output_type}'" >> params.py
	echo "intro = ${intro}" >> params.py
	echo "stats_file='stats.txt'" >> params.py
	echo "tm_lemma = False" >> params.py
	echo "res_filter = False" >> params.py
	echo "context_choice = False" >> params.py
	echo "search_mode = '${search_mode}'" >> params.py
	echo "word_number = '${word_number}'" >> params.py
	'''
// /!\ TM_LEMMA, RES_FILTER & CONTEXT_CHOICE VIENNENT DE CLEANING QUI DOIT ETRE TRUE OU FALSE */
}

process aopHelpFinder{
 	label 'aophelpfinder'

	input:
	tuple file(abstractFile), file(param) from dataAOPhF
	file(eventFile) from event

	output:
	file("output-eventevent.tsv") into AOPhF_out optional true
    file("no-results.txt") into AOPhF_outnr optional true
	file("stats.txt") into  AOPhF_UidEffect

	script:
	"""
	python3 aophelpfinder.py ${param}
	"""
}

AOPhF_out.collectFile(keepHeader: true, skip: 1).into{AOPhF_outFile; AOPhF_outFile2}

AOPhF_outFile2.subscribe{
	f -> f.copyTo(RESULTSDIR_AOPhFout.resolve("${f.name}"))
}

AOPhF_UidEffect.collectFile(keepHeader: false).into{AOPhF_outUidEffectFile; AOPhF_outUidEffectFile2}



process PrepScoringandStatFile{
	input:
	file(result) from AOPhF_outFile
	file(UidEffect) from AOPhF_outUidEffectFile
	each file(statsUid) from statUidFile
	
	output:
	file("stats_event-event.tsv") into StatsResults
	file("prep_scoring_event-event.tsv") into ScoreResults
	script:
	"""
	
	python3 stats.py  --uid ${statsUid} --stats ${UidEffect} --output "${result}" --mode "event-event"
	"""
}


ScoreResults.collectFile(keepHeader: true, skip: 1).into{ToScoring; ToScoring2}


StatsResults.collectFile(keepHeader: true, skip: 1).into{StatsFile; StatsFile2}
StatsFile.subscribe{
	f -> f.copyTo(RESULTSDIR_stats.resolve(f.name))
}


process ScoringandStatFile{

	input:
	file(result) from ToScoring

	output:
	file("scoring_event-event.tsv") into ScoringResultsTSV
	file("scoring_event-event.ods") into ScoringResultsODS

	script:
	"""
	python3 scoring.py --input ${result} --mode event-event
	touch scoring_event-event.ods
	"""
}

ScoringResultsTSV.collectFile(keepHeader: true, skip: 1).into{ScoringResultsTSVall; ScoringResultsTSVall2}

ScoringResultsTSVall2.subscribe{
	f -> f.copyTo(RESULTSDIR_scoring.resolve(f.name))
}

