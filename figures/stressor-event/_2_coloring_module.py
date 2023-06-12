#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:34:29 2022

@author: <Thibaut Coustillet>
__lab__: SysTox Team (U1124 | T3S) FRANCE
__version__: 1.0.0
"""

####################

import csv
import glob
import numpy as np
import os
import pandas as pd
import re
import xlsxwriter
import sys
from pathlib import Path

####################

# Get the parent directory of the directory containing this script.
#pp_directory = Path(__file__).resolve().parents[1].as_posix()
pp_directory = sys.argv[5]

def cell_coloring(workbook, color_hex):
    """
    """
    
    F = workbook.add_format({
        "bold": True,
        "fg_color": color_hex,
        "border": 1, 
        "align": "center",
        "font_name": "Myriad Pro"})
    
    return F


####################

aophf_color = "#AEB1BA"

####################


def coloring_csv(file_name):
    """To color the output file of AOP-helpFinder (after preprocessing).
    
    Args:
    -----
        file_name: A .csv file with AOP-helpFinder data. = "all_stressors.csv"
    
    Returns:
    --------
        The same DataFrame (colored and in .xlsx format).
    
    Note:
    -----
        file_name is the output of 'add_stressors' function 
        (_1_preprocessing_module.py script).
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data = pd.read_csv(file_name, sep=";")
    data.fillna(" ", inplace=True)
    data[" "] = " "
    
    if "abstract" not in data.columns:
        data = data[["PMID", "stressor", "event", "pubdate", "title"]]
        abstract_p = None
    else:
        data = data[["PMID", "stressor", "event", "pubdate", "title", "abstract", " "]]
        abstract_p = True
    
    n_rows = len(data)
    
    event_id = data.columns.get_loc("event")
    pmid_id = data.columns.get_loc("PMID")
    stressor_id = data.columns.get_loc("stressor")
    title_id = data.columns.get_loc("title")
    date_id = data.columns.get_loc("pubdate")
    aophf_ind = [data.columns.to_list().index(c) for c in data.columns]
    
    output_name = "AOPhF.xlsx" 
    output_path = os.path.join(pp_directory, output_name)
    
    writer = pd.ExcelWriter(output_path, engine="xlsxwriter")
    data.to_excel(writer, sheet_name="AOP-helpFinder results", index=None)
    
    info_AOPhF = pd.DataFrame([])
    info_AOPhF.to_excel(writer, sheet_name="Useful Information", index=None)
    
    workbook = writer.book
    ws1 = writer.sheets["AOP-helpFinder results"]
    ws_info = writer.sheets["Useful Information"]
    
    cell_format1 = workbook.add_format({
        "font_name": "Myriad Pro", 
        "border": 1, 
        "border_color": "#D6D6D6"})
    
    alternating_color1 = workbook.add_format({
        "bg_color": "#F5F5F5", 
        "font_name": "Myriad Pro", 
        "border": 1, 
        "border_color": "#D6D6D6"})
    
    alternating_color2 = workbook.add_format({ 
        "font_name": "Myriad Pro", 
        "border": 1, 
        "border_color": "#D6D6D6"})
    
    alternating_color_c = workbook.add_format({
        "bg_color": "#F5F5F5", 
        "font_name": "Myriad Pro", 
        "border": 1, 
        "border_color": "#D6D6D6", 
        "align": "center"})
    
    alternating_color2_c = workbook.add_format({ 
        "font_name": "Myriad Pro", 
        "border": 1, 
        "border_color": "#D6D6D6", 
        "align": "center"})
    
    url_format = workbook.add_format({
        "font_color": "blue",
        "bold": 1,
        "underline": 1,
        "font_size": 12, 
        "font_name": "Myriad Pro"})
    
    ws1.freeze_panes(1, 0)
    
    ws1.set_column(event_id, event_id, 30, cell_format1)
    ws1.set_column(pmid_id, pmid_id, 14, cell_format1)
    ws1.set_column(stressor_id, stressor_id, 20, cell_format1)
    ws1.set_column(title_id, title_id, 45, cell_format1)
    
    for col_num, value in zip(aophf_ind, data.columns.values):
        ws1.write(0, col_num, value, cell_coloring(workbook, aophf_color))
    
    for row in range(n_rows+1):
        ws1.set_row(row, cell_format=(alternating_color1 if row%2==0 else alternating_color2))
    
    if n_rows < 65_530:
        for row in range(len(data)):
            if row % 2 == 0:
                pmid = str(data.loc[row]["PMID"])
                ws1.write_url(row+1, 0, 
                                    f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", 
                                    alternating_color_c,
                                    string=pmid)
            else:
                pmid = str(data.loc[row]["PMID"])
                ws1.write_url(row+1, 0, 
                                    f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", 
                                    alternating_color2_c,
                                    string=pmid)
    else:
        pass
    
    for row in range(len(data)):
        if row % 2 == 0:
            ws1.write(row+1, stressor_id, str(data.loc[row]["stressor"]), alternating_color_c)
            ws1.write(row+1, event_id, str(data.loc[row]["event"]), alternating_color_c)
            ws1.write(row+1, date_id, str(data.loc[row]["pubdate"]), alternating_color_c)
            ws1.write(row+1, title_id, str(data.loc[row]["title"]), alternating_color1)
        else:
            ws1.write(row+1, stressor_id, str(data.loc[row]["stressor"]), alternating_color2_c)
            ws1.write(row+1, event_id, str(data.loc[row]["event"]), alternating_color2_c)
            ws1.write(row+1, date_id, str(data.loc[row]["pubdate"]), alternating_color2_c)
            ws1.write(row+1, title_id, str(data.loc[row]["title"]), alternating_color2)
    
    myriad12 = workbook.add_format({"font_size": 12, "font_name": "Myriad Pro",})
    bold12 = workbook.add_format({"bold": True, "font_size": 12, "font_name": "Myriad Pro",})
    bold12c = workbook.add_format({"bold": True, "font_size": 12, "font_name": "Myriad Pro", "align": "center"})
    
    for row in [0, 1, 11, 12]:
        ws_info.set_row(row, 8)
    
    ws_info.set_column("A:B", 1)
    ws_info.set_column("S:T", 1)
    ws_info.set_column("C:R", 20)
    ws_info.set_column("F:F", 26)
    
    for c in [f"A{i}" for i in range(1,13)]:
        ws_info.write(c, "", cell_coloring(workbook, "#17203C"))
    for c in [f"{l}13" for l in "ABCDEFGHIJKLMNOPQRST"]:
        ws_info.write(c, "", cell_coloring(workbook, "#17203C"))
    for c in [f"T{i}" for i in range(1,13)]:
        ws_info.write(c, "", cell_coloring(workbook, "#17203C"))
    for c in [f"{l}1" for l in "ABCDEFGHIJKLMNOPQRST"]:
        ws_info.write(c, "", cell_coloring(workbook, "#17203C"))
    
    for c in [f"B{i}" for i in list(range(2,12))]:
        ws_info.write(c, "", cell_coloring(workbook, "#EEA83C"))
    for c in [f"{l}12" for l in "BCDEFGHIJKLMNOPQRS"]:
        ws_info.write(c, "", cell_coloring(workbook, "#EEA83C"))
    for c in [f"S{i}" for i in range(2,12)]:
        ws_info.write(c, "", cell_coloring(workbook, "#EEA83C"))
    for c in [f"{l}2" for l in "BCDEFGHIJKLMNOPQRS"]:
        ws_info.write(c, "", cell_coloring(workbook, "#EEA83C"))
    
    if abstract_p:
        provided_by_AOPhF = ["PMID", "stressor", "event", "title", "date", "abstract"]
    else:
        provided_by_AOPhF = ["PMID", "stressor", "event", "title", "date"]
        
    positions = [f"C{i}" for i in range(5, 5+len(provided_by_AOPhF))]
    positions2 = [f"D{i}" for i in range(5, 5+len(provided_by_AOPhF))]
    for i, j in zip(positions, provided_by_AOPhF):
        ws_info.write(i, j, cell_coloring(workbook, aophf_color))
    for ii in positions2:
        ws_info.write(ii, "provided by SysTox", cell_coloring(workbook, aophf_color))
    
    ws_info.write("F5", "How to cite the AOP-helpFinder webserver:", bold12)
    ws_info.write("H5", "Jornod F, Jaylet T, Blaha L, Sarigiannis D, Tamisier L, Audouze K. AOP-helpFinder webserver: a tool for comprehensive analysis of the literature to support adverse outcome pathways development. Bioinformatics. 2021 oct 30.", myriad12)
    ws_info.write_url("Q5", "https://doi.org/10.1093/bioinformatics/btab750", url_format)
    
    ws_info.write("F6", "How to cite the AOP-helpFinder original method:", bold12)
    ws_info.write("H6", "Carvaillo JC, Barouki R, Coumoul X, Audouze K. Linking Bisphenol S to Adverse Outcome Pathways Using a Combined Text Mining and Systems Biology Approach. Environ Health Perspect. 2019 Apr;127(4):47005.", myriad12)
    ws_info.write_url("Q6", "https://doi.org/10.1289/EHP4200", url_format)
    
    ws_info.write("F7", "How to cite the AOP-helpFinder v2.0:", bold12)
    ws_info.write("H7", "Jaylet T, Coustillet T, Jornod F, Margaritte-Jeannin P, Audouze K. AOP-helpFinder 2.0: Integration of an event-event searches module. Environ Int. 2023;177:108017.", myriad12)
    ws_info.write_url("Q7", "https://doi.org/10.1016/j.envint.2023.108017", url_format)
    
    ws_info.write("F9", "Our website:", bold12c)
    ws_info.write_url("G9", "https://systox.u-paris-sciences.fr/", url_format)
    
    ws_info.write("F10", "AOP-helpFinder web server:", bold12c)
    ws_info.write_url("G10", "http://aop-helpfinder.u-paris-sciences.fr", url_format)
    
    
    
    writer.save()
    return None


AOPhF_1s1 = pp_directory + "AOPhF_1s1s.xlsx"
writer_1s1s = pd.ExcelWriter(AOPhF_1s1, engine="xlsxwriter")


def coloring_1s(file_name, directory_name, writer):
    """To color output files of AOP-helpFinder (raw output files).
    
    Args:
    -----
        file_name: file with AOP-helpFinder data
        = WWWWW.tsv where WWWWW is the stressor name.
        It is hosted in the AOPhF results directory, usually "result/".
    
    Returns:
    --------
        The same file (but colored and in a .xlsx format.)
    """
    
    #file_name = os.path.join(pp_directory, file_name)
    
    col_names4 = ["pubdate", "title", "PMID", "event"]
    col_names5 = ["pubdate", "title", "PMID", "event", "abstract"]
    
    with open(file_name) as f:
        reader = csv.reader(f, delimiter="\t")
        first_row = next(reader)
    
    if "abstract" in first_row:
        data = pd.read_csv(file_name, sep="\t", names=col_names5)
        data.fillna(" ", inplace=True)
        data[" "] = " "
    else:
        data = pd.read_csv(file_name, sep="\t", names=col_names4)
        data.fillna(" ", inplace=True)
        
    data.drop(data.head(1).index, inplace=True)
    data.reset_index(drop=True, inplace=True)
    data["PMID"] = data["PMID"].astype(int)
    
    n_rows = len(data)
    stressor = re.search(f"{directory_name}(.*).tsv", file_name).group(1)
    data.insert(0, "stressor", stressor)
    
    if "abstract" not in first_row:
        data = data[["PMID", "stressor", "event", "pubdate", "title"]]
        abstract_p = None
    else:
        data = data[["PMID", "stressor", "event", "pubdate", "title", "abstract", " "]]
        abstract_p = True
    
    if len(stressor) > 31:
        stressor = stressor[:30]
    else:
        pass
    
    data.to_excel(writer, sheet_name=stressor, index=None)
    workbook = writer.book
    ws1 = writer.sheets[stressor]
    
    cell_format1 = workbook.add_format({
        "font_name": "Myriad Pro", 
        "border": 1, 
        "border_color": "#D6D6D6"})
    
    alternating_color1 = workbook.add_format({
        "bg_color": "#F5F5F5", 
        "font_name": "Myriad Pro", 
        "border": 1, 
        "border_color": "#D6D6D6"})
    
    alternating_color2 = workbook.add_format({
        "bg_color": "#FFFFFF", 
        "font_name": "Myriad Pro", 
        "border": 1, 
        "border_color": "#D6D6D6"})
    
    alternating_color_c = workbook.add_format({
        "bg_color": "#F5F5F5", 
        "font_name": "Myriad Pro", 
        "border": 1, 
        "border_color": "#D6D6D6", 
        "align": "center"})
    
    alternating_color2_c = workbook.add_format({ 
        "font_name": "Myriad Pro", 
        "border": 1, 
        "border_color": "#D6D6D6", 
        "align": "center"})
    
    ws1.freeze_panes(1, 0)
    
    for row in range(n_rows+1):
        ws1.set_row(row, cell_format=(alternating_color1 if row%2==0 else alternating_color2))
    
    event_id = data.columns.get_loc("event")
    pmid_id = data.columns.get_loc("PMID")
    title_id = data.columns.get_loc("title")
    stressor_id = data.columns.get_loc("stressor")
    date_id = data.columns.get_loc("pubdate")
    aophf_ind = [data.columns.to_list().index(c) for c in data.columns]
    
    for col_num, value in zip(aophf_ind, data.columns.values):
        ws1.write(0, col_num, value, cell_coloring(workbook, aophf_color))
    
    ws1.set_column(stressor_id, stressor_id, 20, cell_format1)
    ws1.set_column(event_id, event_id, 30, cell_format1)
    ws1.set_column(pmid_id, pmid_id, 14, cell_format1)
    ws1.set_column(title_id, title_id, 45, cell_format1)
    
    if n_rows < 65_530:
        for row in range(len(data)):
            if row % 2 == 0:
                pmid = str(data.loc[row]["PMID"])
                ws1.write_url(row+1, 0, 
                                    f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", 
                                    alternating_color_c,
                                    string=pmid)
            else:
                pmid = str(data.loc[row]["PMID"])
                ws1.write_url(row+1, 0, 
                                    f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", 
                                    alternating_color2_c,
                                    string=pmid)
    else:
        pass
    
    for row in range(len(data)):
        if row % 2 == 0:
            ws1.write(row+1, stressor_id, str(data.loc[row]["stressor"]), alternating_color_c)
            ws1.write(row+1, event_id, str(data.loc[row]["event"]), alternating_color_c)
            ws1.write(row+1, date_id, str(data.loc[row]["pubdate"]), alternating_color_c)
            ws1.write(row+1, title_id, str(data.loc[row]["title"]), alternating_color1)
        else:
            ws1.write(row+1, stressor_id, str(data.loc[row]["stressor"]), alternating_color2_c)
            ws1.write(row+1, event_id, str(data.loc[row]["event"]), alternating_color2_c)
            ws1.write(row+1, date_id, str(data.loc[row]["pubdate"]), alternating_color2_c)
            ws1.write(row+1, title_id, str(data.loc[row]["title"]), alternating_color2)
    
    return None


def coloring_alls(directory_name):
    """To color all files from AOP-helpFinder output (raw output files) 
    and add a 'Useful information' sheet.
    
    Args:
    -----
        directory_name: directory containing files from AOP-helpFinder output,
        usually "result/".
    
    Returns:
    --------
        A colored .xlsx file. 1 sheet = 1 stressor.
    """
    
    directory_name = os.path.join(pp_directory, directory_name)
    directory_content = directory_name + "*.tsv"
    
    all_s = []
    for s in glob.glob(directory_content):
        FILE  = s.split("/")[-1]
        if FILE.startswith("scoring"):
            pass
        else:
            all_s += [s]
    
    with open(all_s[0]) as f:
        reader = csv.reader(f, delimiter="\t")
        first_row = next(reader)
    
    if "abstract" in first_row:
        abstract_p = True
    else:
        abstract_p = None
    
    all_s = sorted(all_s, key=str.casefold)
    for s in all_s:
        coloring_1s(s, directory_name, writer_1s1s)
    
    info_AOPhF = pd.DataFrame([])
    info_AOPhF.to_excel(writer_1s1s, sheet_name="Useful Information", index=None)
    ws_info = writer_1s1s.sheets["Useful Information"]
    workbook = writer_1s1s.book
    
    myriad12 = workbook.add_format({"font_size": 12, "font_name": "Myriad Pro",})
    bold12 = workbook.add_format({"bold": True, "font_size": 12, "font_name": "Myriad Pro",})
    bold12c = workbook.add_format({"bold": True, "font_size": 12, "font_name": "Myriad Pro", "align": "center"})
    
    url_format = workbook.add_format({
        "font_color": "blue",
        "bold": 1,
        "underline": 1,
        "font_size": 12, 
        "font_name": "Myriad Pro"})
    
    for row in [0, 1, 11, 12]:
        ws_info.set_row(row, 8)
    
    ws_info.set_column("A:B", 1)
    ws_info.set_column("S:T", 1)
    ws_info.set_column("C:R", 20)
    ws_info.set_column("F:F", 26)
    
    for c in [f"A{i}" for i in range(1,13)]:
        ws_info.write(c, "", cell_coloring(workbook, "#17203C"))
    for c in [f"{l}13" for l in "ABCDEFGHIJKLMNOPQRST"]:
        ws_info.write(c, "", cell_coloring(workbook, "#17203C"))
    for c in [f"T{i}" for i in range(1,13)]:
        ws_info.write(c, "", cell_coloring(workbook, "#17203C"))
    for c in [f"{l}1" for l in "ABCDEFGHIJKLMNOPQRST"]:
        ws_info.write(c, "", cell_coloring(workbook, "#17203C"))
    
    for c in [f"B{i}" for i in list(range(2,12))]:
        ws_info.write(c, "", cell_coloring(workbook, "#EEA83C"))
    for c in [f"{l}12" for l in "BCDEFGHIJKLMNOPQRS"]:
        ws_info.write(c, "", cell_coloring(workbook, "#EEA83C"))
    for c in [f"S{i}" for i in range(2,12)]:
        ws_info.write(c, "", cell_coloring(workbook, "#EEA83C"))
    for c in [f"{l}2" for l in "BCDEFGHIJKLMNOPQRS"]:
        ws_info.write(c, "", cell_coloring(workbook, "#EEA83C"))
    
    if abstract_p:
        provided_by_AOPhF = ["PMID", "stressor", "event", "title", "date", "abstract"]
    else:
        provided_by_AOPhF = ["PMID", "stressor", "event", "title", "date"]
    
    positions = [f"C{i}" for i in range(5, 5+len(provided_by_AOPhF))]
    positions2 = [f"D{i}" for i in range(5, 5+len(provided_by_AOPhF))]
    for i, j in zip(positions, provided_by_AOPhF):
        ws_info.write(i, j, cell_coloring(workbook, aophf_color))
    for ii in positions2:
        ws_info.write(ii, "provided by SysTox", cell_coloring(workbook, aophf_color))
    
    ws_info.write("F5", "How to cite the AOP-helpFinder webserver:", bold12)
    ws_info.write("H5", "Jornod F, Jaylet T, Blaha L, Sarigiannis D, Tamisier L, Audouze K. AOP-helpFinder webserver: a tool for comprehensive analysis of the literature to support adverse outcome pathways development. Bioinformatics. 2021 oct 30.", myriad12)
    ws_info.write_url("Q5", "https://doi.org/10.1093/bioinformatics/btab750", url_format)
    
    ws_info.write("F6", "How to cite the AOP-helpFinder method:", bold12)
    ws_info.write("H6", "Carvaillo JC, Barouki R, Coumoul X, Audouze K. Linking Bisphenol S to Adverse Outcome Pathways Using a Combined Text Mining and Systems Biology Approach. Environ Health Perspect. 2019 Apr;127(4):47005.", myriad12)
    ws_info.write_url("Q6", "https://doi.org/10.1289/EHP4200", url_format)
    
    ws_info.write("F8", "Our website:", bold12c)
    ws_info.write_url("G8", "https://systox.u-paris-sciences.fr/", url_format)
    
    ws_info.write("F9", "AOP-helpFinder web server:", bold12c)
    ws_info.write_url("G9", "http://aop-helpfinder.u-paris-sciences.fr", url_format)
    
    writer_1s1s.save()
    return None


####################

low_color = "#fe9929"
quite_low_color = "#fec44f"
moderate_color = "#e6f598"
high_color = "#91cf60"
very_high_color = "#1a9850"

####################


def make_data_colors(file_name):
    """Create the color file necessary for the coloring of the stressor-event links file.
    
    Args:
    -----
        file_name: scoring file in AOP-helpFinder output directory
        = "scoring_stressor-event.csv"
    
    Returns:
    --------
        A DataFrame with stressors in columns and events in row.
        The location (i,j) holds the color linking the stressor i to the event j, 
        i.e. the strength of association corresponding to the confidence score.
    
    Note:
    -----
        It is the same file as the ouput of the 'make_file4heatmap' function 
        (_1_preprocessing_module.py script), except that the value of the i-th row
        and j-th column is the color of the future cell in HEX format 
        (instead of the number of links).
    """
    
    data = pd.read_csv(file_name, sep="\t")
    data_colors = data.pivot(index="Stressor", columns="Event", values="Confidence")
    data_colors.fillna(0, inplace=True)
    data_colors.index.name = "stressor \ event"
    data_colors.columns.name = ""
    
    data_colors.replace("Low", low_color, inplace=True)
    data_colors.replace("Quite low", quite_low_color, inplace=True)
    data_colors.replace("Moderate", moderate_color, inplace=True)
    data_colors.replace("High", high_color, inplace=True)
    data_colors.replace("Very High", very_high_color, inplace=True)
    data_colors.replace(0, "#f0f0f0", inplace=True)
    
    data_colors.sort_values(by="stressor \ event", inplace=True, key=lambda col: col.str.lower())
    
    return data_colors


def coloring_stressor_event_links(file4heatmap, color_file):
    """Create the colored file hosting stressor-event links.
    
    Args:
    -----
        file4heatmap: file to be colored 
        = output of 'make_file4heatmap' function, = "file4heatmap.csv"
        color_file: output of 'make_data_colors' function
        = the same file as the "file4heatmap" file, except that the value of 
        the i-th row and j-th column is the color of the future cell in HEX format 
        (instead of the number of links).
    
    Returns:
    --------
        The same DataFrame as input (colored and in .xlsx format).
    
    Note:
    -----
        The color indicates the confidence score.
    """
    
    file4heatmap = os.path.join(pp_directory, file4heatmap)
    color_file = os.path.join(pp_directory, color_file)
    data_scoring = pd.read_csv(file4heatmap, sep=";", index_col=0)
    data_colors = make_data_colors(color_file)
    
    output_excel_name = "stressor_event_links.xlsx"
    output_excel_path = os.path.join(pp_directory, output_excel_name)
    writer = pd.ExcelWriter(output_excel_path, engine="xlsxwriter")
    data_scoring.to_excel(writer, sheet_name="stressor-event links")
    workbook = writer.book
    worksheet = writer.sheets["stressor-event links"]
    
    n_cols = len(data_scoring.columns)
    n_rows = len(data_scoring)
    
    worksheet.set_column(f"{xlsxwriter.utility.xl_col_to_name(n_cols+1)}:{xlsxwriter.utility.xl_col_to_name(n_cols+1)}", 1)
    worksheet.set_row(n_rows+1, 10)
    
    worksheet.freeze_panes(1, 1)
    
    worksheet.set_column(f"A:{xlsxwriter.utility.xl_col_to_name(n_cols)}", 25)
    for row in range(n_rows):
        worksheet.set_row(row, 20)
    
    cell_format1 = workbook.add_format({"font_name": "Myriad Pro", 
                                        "bold": True, 
                                        "font_size": 14, 
                                        "align": "center", 
                                        "border": 1, 
                                        "border_color": "#878787"})
    
    worksheet.write(0, 0, data_scoring.index.name, cell_format1)
    
    for row in range(data_scoring.shape[0]):
        for col in range(data_scoring.shape[1]):
            color_numbers = workbook.add_format({"bg_color": f"{data_colors.iloc[row,col]}", 
                                             "font_name": "Myriad Pro", 
                                             "bold": True, 
                                             "font_size": 14, 
                                             "align": "center", 
                                             "border": 1, 
                                             "border_color": "#878787"})
            
            worksheet.write(row+1, col+1, int(data_scoring.iloc[row,col]), color_numbers)
    
    for row in range(data_scoring.shape[0]):
        format1_col0 = workbook.add_format({"font_name": "Myriad Pro", 
                             "bold": False, 
                             "font_size": 12, 
                             "align": "left", 
                             "border": 1, 
                             "border_color": "#878787"})
        
        format2_col0 = workbook.add_format({'bg_color': "#F5F5F5", 
                                            "font_name": "Myriad Pro", 
                                            "font_size": 12, 
                                            "align": "left", 
                                            "border": 1, 
                                            "border_color": "#878787"})
        
        if row%2 == 0:
            worksheet.write(row+1, 0, data_scoring.index[row], format1_col0)
        else:
            worksheet.write(row+1, 0, data_scoring.index[row], format2_col0)
    
    for col in range(data_scoring.shape[1]):
        format1_row0 = workbook.add_format({"font_name": "Myriad Pro", 
                                            "bold": False, 
                                            "font_size": 12, 
                                            "align": "center", 
                                            "border": 1, 
                                            "border_color": "#878787"})
        
        format2_row0 = workbook.add_format({"bg_color": "#F5F5F5", 
                                            "font_name": "Myriad Pro", 
                                            "font_size": 12, "align": 
                                            "center", "border": 1, 
                                            "border_color": "#878787"})
        
        if col%2 == 0:
            worksheet.write(0, col+1, data_scoring.columns[col], format1_row0)
        else:
            worksheet.write(0, col+1, data_scoring.columns[col], format2_row0)
    
    ttotal = 0
    
    worksheet.write(0, data_scoring.shape[1]+2, "Total", workbook.add_format({"font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787", "bold": True}))
    for row in data_scoring.index.values.tolist():
        total4stressors = sum(data_scoring.loc[row])
        ttotal += total4stressors
        worksheet.write(data_scoring.index.get_loc(row)+1, data_scoring.shape[1]+2, total4stressors, workbook.add_format({"font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    
    worksheet.write(data_scoring.shape[0]+2, 0, "Total", workbook.add_format({"font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787", "bold": True}))
    for col in data_scoring.columns.tolist():
        total4events = sum(data_scoring[col])
        worksheet.write(data_scoring.shape[0]+2, data_scoring.columns.get_loc(col)+1, total4events, workbook.add_format({"font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    
    worksheet.write(data_scoring.shape[0]+2, data_scoring.shape[1]+2, ttotal, workbook.add_format({"font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    
    worksheet.write(n_rows+5, 0, "Confidence score:", workbook.add_format({'bg_color': "#F5F5F5", "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    worksheet.write(n_rows+6, 0, "Low", workbook.add_format({'bg_color': low_color, "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    worksheet.write(n_rows+7, 0, "Quite Low", workbook.add_format({'bg_color': quite_low_color, "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    worksheet.write(n_rows+8, 0, "Moderate", workbook.add_format({'bg_color': moderate_color, "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    worksheet.write(n_rows+9, 0, "High", workbook.add_format({'bg_color': high_color, "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    worksheet.write(n_rows+10, 0, "Very High", workbook.add_format({'bg_color': very_high_color, "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    
    writer.save()
    return None
    