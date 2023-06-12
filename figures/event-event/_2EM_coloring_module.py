#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:34:29 2022

@author: <Thibaut Coustillet>
__lab__: SysTox Team (U1124 | T3S) FRANCE
__version__: 1.0.0
"""

####################

import numpy as np
import os
import pandas as pd
import re
import xlsxwriter
import sys
from pathlib import Path

import _1EM_preprocessing_module as ppm

####################

# Get the parent directory of the directory containing this script.
pp_directory = sys.argv[3]

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
    """To color the output file of AOP-helpFinder (event-event search).
    
    Args:
    -----
        file_name: A .csv file with AOP-helpFinder data = "event-event.tsv".
    
    Returns:
    --------
        The same DataFrame (colored and in .xlsx format).
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data = pd.read_csv(file_name, sep="\t")
    data.fillna(" ", inplace=True)
    data[" "] = " "
    
    if "abstract" not in data.columns:
        data = data[["pmid", "event_1", "event_2", "pubdate", "title"]]
        abstract_p = None
    else:
        data = data[["pmid", "event_1", "event_2", "pubdate", "title", "abstract", " "]]
        abstract_p = True
    
    n_rows = len(data)
    data.rename({"pmid": "PMID" }, axis=1, inplace=True)
    
    event1_id = data.columns.get_loc("event_1")
    pmid_id = data.columns.get_loc("PMID")
    event2_id = data.columns.get_loc("event_2")
    title_id = data.columns.get_loc("title")
    date_id = data.columns.get_loc("pubdate")
    aophf_ind = [data.columns.to_list().index(c) for c in data.columns]
    
    output_file_name = "AOPhF.xlsx" 
    output_path = os.path.join(pp_directory, output_file_name)
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
    
    alternating_color = workbook.add_format({
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
    
    ws1.set_column(event1_id, event1_id, 30, cell_format1)
    ws1.set_column(event2_id, event2_id, 30, cell_format1)
    ws1.set_column(pmid_id, pmid_id, 14, cell_format1)
    ws1.set_column(title_id, title_id, 45, cell_format1)
    
    for col_num, value in zip(aophf_ind, data.columns.values):
        ws1.write(0, col_num, value, cell_coloring(workbook, aophf_color))
        
    for row in range(n_rows+1):
        ws1.set_row(row, cell_format=(alternating_color if row%2==0 else alternating_color2))
    
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
            ws1.write(row+1, event1_id, str(data.loc[row]["event_1"]), alternating_color_c)
            ws1.write(row+1, event2_id, str(data.loc[row]["event_2"]), alternating_color_c)
            ws1.write(row+1, date_id, str(data.loc[row]["pubdate"]), alternating_color_c)
            ws1.write(row+1, title_id, str(data.loc[row]["title"]), alternating_color)
        else:
            ws1.write(row+1, event1_id, str(data.loc[row]["event_1"]), alternating_color2_c)
            ws1.write(row+1, event2_id, str(data.loc[row]["event_2"]), alternating_color2_c)
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
        provided_by_AOPhF = ["PMID", "event_1", "event_2", "title", "pubdate", "abstract"]
    else:
        provided_by_AOPhF = ["PMID", "event_1", "event_2", "title", "pubdate"]
    
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

####################

low_color = "#fe9929"
quite_low_color = "#fec44f"
moderate_color = "#e6f598"
high_color = "#91cf60"
very_high_color = "#1a9850"
    
####################


def make_data_colors(file_name):
    """Create the color file necessary for the coloring of the event-event links file.
    
    Args:
    -----
        file_name: scoring file in AOP-helpFinder output directory
        = "scoring_event-event.csv"
    
    Returns:
    --------
        A DataFrame with events in rows and columns (double-input table).
        The location (i,j) holds the color linking the event i to the event j, 
        i.e. the strength of association corresponding to the confidence score.
    
    Note:
    -----
        It is the same file as the ouput of the 'preprocess' function 
        (_1ev_preprocessing_module.py script), except that the value of the i-th row
        and j-th column is the color of the future cell in HEX format 
        (instead of the number of links).
    """
    
    data_colors1 = pd.read_csv(file_name, sep="\t")
    data_colors2 = data_colors1.copy()
    data_colors2.rename({"Event 1": "Event 2", "Event 2": "Event 1"}, axis="columns", inplace=True)
    data_colors3 = pd.concat([data_colors1, data_colors2], axis=0)
    data_colors3.reset_index(drop=True, inplace=True)
    
    df = data_colors3.pivot(index="Event 1", columns="Event 2", values="Confidence")
    df.fillna(0, inplace=True)
    data_colors = df.where(np.tril(np.ones(df.shape), -1).astype(bool))
    data_colors.fillna("x", inplace=True)
    data_colors.index.name = "event1 \ event2"
    data_colors.columns.name = ""
    
    data_colors.replace("Low", low_color, inplace=True)
    data_colors.replace("Quite low", quite_low_color, inplace=True)
    data_colors.replace("Moderate", moderate_color, inplace=True)
    data_colors.replace("High", high_color, inplace=True)
    data_colors.replace("Very High", very_high_color, inplace=True)
    data_colors.replace("x", "#bdbdbd", inplace=True)
    data_colors.replace(0, "#f0f0f0", inplace=True)
    
    return data_colors


def coloring_event_event_links(file_name):
    """Create the colored file hosting event-event links.
    
    Args:
    -----
        file_name: file to be colored 
        = output of 'preprocess' function, = "event_event_links.csv"
    
    Returns:
    --------
        A file (same as input) colored and in a .xslx format.
    
    Note:
    -----
        The color indicates the confidence score.
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data_scoring = ppm.preprocess(file_name)
    data_colors = make_data_colors(file_name)
    
    output_excel_name = "event_event_links.xlsx"
    output_excel_path = os.path.join(pp_directory, output_excel_name)
    writer = pd.ExcelWriter(output_excel_path, engine="xlsxwriter")
    data_scoring.to_excel(writer, sheet_name="event-event links")
    workbook = writer.book
    ws1 = writer.sheets["event-event links"]
    
    n_cols = len(data_scoring.columns)
    n_rows = len(data_scoring)
    
    ws1.freeze_panes(1, 1)
    
    ws1.set_column(f"A:{xlsxwriter.utility.xl_col_to_name(n_cols)}", 25)
    for row in range(n_rows):
        ws1.set_row(row, 20)
        
    cell_format1 = workbook.add_format({"font_name": "Myriad Pro", 
                                        "bold": True, 
                                        "font_size": 14, 
                                        "align": "center", 
                                        "border": 1, 
                                        "border_color": "#878787"})
    
    ws1.write(0, 0, data_scoring.index.name, cell_format1)
    
    for row in range(data_scoring.shape[0]):
        for col in range(data_scoring.shape[1]):
            color_x = workbook.add_format({"bg_color": f"{data_colors.iloc[row,col]}", 
                                           "font_name": "Myriad Pro", 
                                           "font_size": 12, 
                                           "align": "center"})
            
            color_numbers = workbook.add_format({"bg_color": f"{data_colors.iloc[row,col]}", 
                                             "font_name": "Myriad Pro", 
                                             "bold": True, 
                                             "font_size": 14, 
                                             "align": "center", 
                                             "border": 1, 
                                             "border_color": "#878787"})
            
            if row != col:
                if row < col:
                    ws1.write(row+1, col+1, data_scoring.iloc[row,col], color_x)
                else:
                    ws1.write(row+1, col+1, int(data_scoring.iloc[row,col]), color_numbers)
            else:
                ws1.write(row+1, col+1, data_scoring.iloc[row,col], workbook.add_format({"bg_color": "black"}))
            
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
            ws1.write(row+1, 0, data_scoring.index[row], format1_col0)
        else:
            ws1.write(row+1, 0, data_scoring.index[row], format2_col0)
              
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
            ws1.write(0, col+1, data_scoring.columns[col], format1_row0)
        else:
            ws1.write(0, col+1, data_scoring.columns[col], format2_row0)
    
    ws1.write(n_rows+3, 0, "Confidence score:", workbook.add_format({'bg_color': "#F5F5F5", "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    ws1.write(n_rows+5, 0, "Quite Low", workbook.add_format({'bg_color': quite_low_color, "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    ws1.write(n_rows+4, 0, "Low", workbook.add_format({'bg_color': low_color, "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    ws1.write(n_rows+6, 0, "Moderate", workbook.add_format({'bg_color': moderate_color, "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    ws1.write(n_rows+7, 0, "High", workbook.add_format({'bg_color': high_color, "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    ws1.write(n_rows+8, 0, "Very High", workbook.add_format({'bg_color': very_high_color, "font_name": "Myriad Pro", "font_size": 14, "align": "center", "border": 1, "border_color": "#878787"}))
    
    writer.save()
    return None
    