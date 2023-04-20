#!/usr/bin/env nextflow
nextflow.enable.dsl=1 
params.stressors = "stressor.txt"
params.event = "event.txt"
params.resultspath= 'results/'
params.userParam = "user_id;homemade;tsv;-1;False;stressor-event;encours"

stressors = file(params.stressors)
event = file(params.event)
RESULTSDIR = file(params.resultspath)
RESULTSDIR_AOPhFout = file("${RESULTSDIR}/raw_AOP-helpFinder/")
RESULTSDIR_scoring = file("${RESULTSDIR}/raw_scoring/")
RESULTSDIR_stats = file("${RESULTSDIR}/resume/")
paramUser = params.userParam

process cleanStressorFile{
	input:
	file(stressor) from stressors.splitText()

	output:
	file("stressorCleaned.txt") into stressorsCleaned

	script:
	"""
	sed '/^ *\$/d' ${stressor} > stressorCleaned.txt
	"""
}


process extractPubmedUidStressor {
	publishDir "${RESULTSDIR}/uid", mode: "copy"
	label "extractPubmedUidStressor"
	
	input:
	val(stressor) from stressorsCleaned.splitText()
	
	output:
	tuple val(stressor), file("*uid.txt") into uids, uids2
	
	script:
	"""
	
esearch -db pubmed -query "${stressor.trim()}" | efetch -chunk 10000 -format uid > ${stressor.trim().replaceAll(/[^a-zA-Z0-9\-\_]/,'_')}-uid.txt
	
	"""
}

process countUids{
	input:
	tuple val(stressor), file(uidFile) from uids2

	output:
	file('stats_uid.txt') into statsUid

	script:
	"""
	nbLine=\$(wc -l ${uidFile})
	echo "${stressor.trim().replaceAll(/[^a-zA-Z0-9\-\_]/,'_')} \$nbLine" > stats_uid.txt
	"""
}

statsUid.collectFile(name: 'stats_uid.txt').into{statUidFile;statUidFile2}
//statUidFile2.subscribe(f -> f.copyTo(RESULTSDIR.resolve(f.name)))


process downloadAbstractXML{
	label 'DLXML'
	
	input:
	tuple val(stressor), file(uid) from uids.splitText(by:750)
	
	output:
	tuple val(stressor), file("abstract.xml") into xml
	
	script:
	"""
	cat ${uid} | fetch-pubmed -path | tr -dc '[:alnum:]\\n\\r<>/ !?"=.,-' | sed '3,\${/<?xml/d}' | sed '3,\${/<!DOCTYPE/d}'  > abstract.xml
	"""	
}

process generateParamsFile{
	label 'paramFile'

	input:
	tuple val(stressor), file(abstractFile) from xml
	val(parameters) from paramUser
	file(eventFile) from event

	output:
	tuple val(stressor), file('abstract.xml'), file('params.py') into dataAOPhF

	shell:
	'''
	IFS=$';' read -r user_id events_type output_type intro cleaning Fig search_mode encours <<< "!{parameters}"
	echo "abstracts_file = '!{abstractFile}'" > params.py
	echo "abstracts_type = 'xml'" >> params.py
	echo "events_file = '!{event}'" >> params.py
	echo "events_type = '$events_type'" >> params.py
	echo "output_name = '!{stressor.trim().replaceAll(/[^a-zA-Z0-9\\-\\_]/,'_')}'" >> params.py
	echo "outputtype = '${output_type}'" >> params.py
	echo "intro = ${intro}" >> params.py
	echo "stats_file='stats.txt'" >> params.py
	echo "tm_lemma = ${cleaning}" >> params.py
	echo "res_filter = ${cleaning}" >> params.py
	echo "context_choice = ${cleaning}" >> params.py
	echo "search_mode = '${search_mode}'" >> params.py
	'''
// /!\ TM_LEMMA, RES_FILTER & CONTEXT_CHOICE VIENNENT DE CLEANING QUI DOIT ETRE TRUE OU FALSE */
}

process aopHelpFinder{
 	label 'aophelpfinder'

	input:
	tuple val(stressor), file(abstractFile), file(param) from dataAOPhF
	file(eventFile) from event

	output:
	tuple val("${stressor.trim().replaceAll(/[^a-zA-Z0-9\-\_]/,'_')}"), file("${stressor.trim().replaceAll(/[^a-zA-Z0-9\-\_]/,'_')}*") into AOPhF_out optional true
	file("no-results.txt") into AOPhF_outnr optional true
	tuple val("${stressor.trim().replaceAll(/[^a-zA-Z0-9\-\_]/,'_')}_stats"), file("stats.txt") into  AOPhF_UidEffect

	script:
	"""
	python3 aophelpfinder.py ${param}
	"""
}


AOPhF_out.collectFile(keepHeader: true, skip: 1).into{AOPhF_outFile; AOPhF_outFile2}

AOPhF_outFile2.subscribe{
	f -> f.copyTo(RESULTSDIR_AOPhFout.resolve("${f.name}.tsv"))
}


AOPhF_UidEffect.collectFile(keepHeader: false).into{AOPhF_outUidEffectFile; AOPhF_outUidEffectFile2}


process PrepScoringandStatFile{
	//publishDir "${RESULTSDIR}/stat", mode: "copy"
	publishDir "${RESULTSDIR_stats}", pattern:"resume*", mode: "copy"


	input:
	file(result) from AOPhF_outFile
	file(UidEffect) from AOPhF_outUidEffectFile
	each file(statsUid) from statUidFile

	output:
	file("stats_stressor-event.tsv") into StatsResults
	file("prep_scoring_stressor-event.tsv") into ScoreResults
	file("resume*.txt") into ResumeStressor
	script:
	"""
	
	python3 stats.py  --uid ${statsUid} --stats ${UidEffect} --output "${result.getSimpleName().trim()}" --mode "stressor-event"
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
	file("scoring_stressor-event.tsv") into ScoringResultsTSV
	file("scoring_stressor-event.ods") into ScoringResultsODS

	script:
	"""
	python3 scoring.py --input ${result} --mode stressor-event
	touch scoring_event-event.ods
	"""
}

ScoringResultsTSV.collectFile(keepHeader: true, skip: 1).into{ScoringResultsTSVall; ScoringResultsTSVall2}

ScoringResultsTSVall2.subscribe{
	f -> f.copyTo(RESULTSDIR_scoring.resolve(f.name))
}


