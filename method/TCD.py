# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
#authors: Thibaut Coustillet - Université Paris Cité - France
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

####################

import pandas as pd
import sys
import time
import xlsxwriter

####################

low_color = "#fe9929"
quite_low_color = "#fec44f"
moderate_color = "#e6f598"
high_color = "#91cf60"
very_high_color = "#1a9850"

####################

date_name = "pubdate"


def colors(row):
    """To match a color to each confidence level.
    """
    
    if row["confidence"] == "Low":
        val = low_color
    elif row["confidence"] == "Quite low":
        val = quite_low_color
    elif row["confidence"] == "Moderate":
        val = moderate_color
    elif row["confidence"] == "High":
        val = high_color
    elif row["confidence"] == "Very High":
        val = very_high_color
    else:
        val = "#636363"
    return val


def stars(row):
    """To match star(s) to each confidence level.
    """
    
    if row["confidence"] == "Low":
        val = "★☆☆☆☆"
    elif row["confidence"] == "Quite low":
        val = "★★☆☆☆" 
    elif row["confidence"] == "Moderate":
        val = "★★★☆☆"
    elif row["confidence"] == "High":
        val = "★★★★☆"
    elif row["confidence"] == "Very High":
        val = "★★★★★"
    else:
        val = "?????"
    return val


def add_scoring(tcd_file, data_scoring):
    """To merge data from AOP-helpFinder output : add the scoring data to 
       the main data (after pre-processing of the data).
    """
    
    data = tcd_file.copy()
    data = data.where(data != 1, data.columns.to_series(), axis=1)
    data["PMIDs"] = data[data.columns].values.tolist()
    data = pd.DataFrame(data["PMIDs"])
    data["PMIDs"] = data["PMIDs"].apply(lambda x: sorted(x, reverse=True))
    data['t_PMIDs'] = data['PMIDs'].apply(lambda row: [val for val in row if val != 0])
    data.insert(0, "Total", data["t_PMIDs"].str.len())
    del data["t_PMIDs"]
    
    data[[f"PMID {i+1}" for i in range(len(data["PMIDs"][0]))]] = pd.DataFrame(data["PMIDs"].tolist(), index=data.index)
    del data["PMIDs"]
    
    nunique = data.nunique()
    df_ws = data.drop(nunique[nunique == 1].index, axis=1)
    df_ws.replace(0, "", inplace=True)
    
    output_df = pd.concat([df_ws, data_scoring], axis=1)
    output_df.insert(0, 'confidence', output_df.pop('confidence'))
    output_df["confidence"] = output_df.apply(stars, axis=1)
    
    return output_df


def make_TCD(file_name, scoring_file):
    """To generate TCD-like files from AOPhF output.
    
    Args:
    -----
        file_name: output file from AOP-helpFinder 
        (usually called all_stressors.csv or output-eventevent.tsv).
        scoring_file: scoring file from AOP-helpFinder output.
    
    Returns:
    --------
        A DataFrame : One row per stressor-event or (event-event) link.
        As many columns as there are PMIDs that capture the link.
    """
    
    scoreF = pd.read_csv(scoring_file, sep="\t")

    if "Event 1" and "Event 2" in scoreF.columns:
        search_type = "EV-EV"
    else:
        search_type = "ST-EV"

    if search_type == "ST-EV":
        data0 = pd.read_csv(file_name, sep=";")

        if 'date' in data0.columns:
            data0.rename(columns={'date': date_name}, inplace=True)
        data0 = data0[["PMID", "stressor", "event", date_name, "title"]]
        
        scoring_cols = ["stressor", "event", "article", "link", "pvalue", "pvar", "confidence"]
        data_scoring = pd.read_csv(scoring_file, sep="\t", names=scoring_cols, header=0)
        data_scoring["stressor"] = data_scoring["stressor"].str.replace(" ", "_")
        
        data_scoringA = data_scoring.set_index(["stressor", "event"])
        data_scoringA = data_scoringA[["confidence"]]
        data_scoringA["colors"] = data_scoringA.apply(colors, axis=1)
        
        data_scoringB = data_scoring.set_index(["event", "stressor"])
        data_scoringB = data_scoringB[["confidence"]]
        data_scoringB["colors"] = data_scoringB.apply(colors, axis=1)
        
        data0A = pd.crosstab([data0["stressor"], data0["event"]], data0["PMID"])
        data0B = pd.crosstab([data0["event"], data0["stressor"]], data0["PMID"])
        
        dataA = add_scoring(data0A, data_scoringA)
        dataB = add_scoring(data0B, data_scoringB)
        
        dataA.sort_values(by=["stressor", "Total"], key=lambda col: col.str.lower() if col.name == "stressor" else col, 
                          ascending=[True, False], inplace=True)
        dataB.sort_values(by=["event", "Total"], key=lambda col: col.str.lower() if col.name == "event" else col, 
                          ascending=[True, False], inplace=True)
        
        return dataA, dataB
    
    elif search_type == "EV-EV":
        data00 = pd.read_csv(file_name, sep="\t")
        
        if 'date' in data00.columns:
            data00.rename(columns={'date': date_name}, inplace=True)

        data00 = data00[["pmid", "event_1", "event_2", date_name, "title"]]
        data01 = data00.copy()
        data01.rename(columns={"event_1": "event_2", "event_2": "event_1"}, inplace=True)
        data0 = pd.concat([data00, data01], ignore_index=True)
        
        scoring_cols = ["event_1", "event_2", "link", "pvalue", "pvar", "confidence"]
        data_scoring0 = pd.read_csv(scoring_file, sep="\t", names=scoring_cols, header=0)
        data_scoring1 = data_scoring0.copy()
        data_scoring1.rename(columns={"event_1": "event_2", "event_2": "event_1"}, inplace=True)
        data_scoring = pd.concat([data_scoring0, data_scoring1], ignore_index=True)
        data_scoring.set_index(["event_1", "event_2"], inplace=True)
        
        data_scoring = data_scoring[["confidence"]]
        data_scoring["colors"] = data_scoring.apply(colors, axis=1)
        
        data_tcd = pd.crosstab([data0["event_1"], data0["event_2"]], data0["pmid"])
        data = add_scoring(data_tcd, data_scoring)
        
        data.sort_values(by=['event_1', "Total"], key=lambda col: col.str.lower() if col.name == "event_1" else col, 
                         ascending=[True, False], inplace=True)
        
        return data


# Setting up the color functions to format the output table.
def ccoloring(workbook, color, font_size):
    f = workbook.add_format({
        "align": "center",
        "font_size":font_size,
        "bold": True,
        "border": 1,
        "fg_color": color,
        "font_name": "Myriad Pro"})
    return f


def ccoloring_cs(workbook, color_hex):
    F = workbook.add_format({
        "fg_color": color_hex,
        "border": 1, 
        "align": "right"})
    return F


def falternating_color(workbook, color, border_color, align, font_size):
    f = workbook.add_format({
        "align": align,
        "bold": False,
        "border": 1,
        "border_color": border_color,
        "bg_color": color,
        "font_color": "#000000",
        "font_name": "Myriad Pro", 
        "font_size":font_size})
    return f


def coloring(writer, df_name, sheet_name):
    """To color output DataFrame (TCD-like file).
    
    Args:
    -----
        writer: name of the output file.
        df_name: Dataframe to be colored.
        sheet_name: sheet name of the .xlsx output file 
    
    Returns:
    --------
        A .xlsx colored file.
    
    Notes:
    ------
    "sheet_name" must belong to {"stressor-event", "event-stressor", "event-event"}.
    For the TCD-like files related to stressor-event search, two sheets:
    One where events depend on stressors and the other where stressors depend on events.
    """
    
    df_name.drop(labels="colors", axis=1).to_excel(writer, sheet_name=sheet_name, index=True)
    workbook = writer.book
    ws1 = writer.sheets[sheet_name]
    
    n_rows = df_name.shape[0]
    n_cols = df_name.shape[1]
    
    CF1 = workbook.add_format({
            "font_name": "Myriad Pro", 
            "border": 1, 
            "border_color": "#D6D6D6"})
    
    CF2 = workbook.add_format({
            "bold": True,
            "font_size" : 12,
            "font_name": "Myriad Pro", 
            "border": 1, 
            "border_color": "#919191", 
            "align": "center", 
            "valign": "vcenter"})
    
    FALTERNATING1 = falternating_color(workbook, "#F5F5F5", "#D6D6D6", "left", 11)
    FALTERNATING2 = falternating_color(workbook, "#FFFFFF", "#D6D6D6", "left", 11)
    FALTERNATINGC = falternating_color(workbook, "#F5F5F5", "#D6D6D6", "center", 11)
    FALTERNATINGC2 = falternating_color(workbook, "#FFFFFF", "#D6D6D6", "center", 11)
    
    PINK1 = ccoloring(workbook, "#F8C5FC", 12)
    PINK2 = ccoloring(workbook, "#F5A2C7", 12)
    YELLOW = ccoloring(workbook, "#FFFB00", 12)
    PURPLE1 = ccoloring(workbook, "#BEC1FF", 11)
    PURPLE2 = ccoloring(workbook, "#AC95CB", 11)
    
    ws1.freeze_panes(1, 4)
    ws1.set_zoom(120)
    
    ws1.set_column(0, 0, 28, CF1)
    ws1.set_column(1, 1, 28, CF1)
    
    if sheet_name == "stressor-event":
        ws1.write(0, 0, "stressor", PINK1)
        ws1.write(0, 1, "event", PINK2)
    elif sheet_name == "event-stressor":
        ws1.write(0, 0, "event", PINK2)
        ws1.write(0, 1, "stressor", PINK1)
    elif sheet_name == "event-event":
        ws1.write(0, 0, "event 1", PINK1)
        ws1.write(0, 1, "event 2", PINK2)
    else:
        pass
    
    ws1.set_column(2, 2, 10.5, CF1)
    ws1.set_column(3, n_cols, 10)
    
    ws1.write(0, 2, "Confidence", YELLOW)
    ws1.write(0, 3, "Total", YELLOW)
    
    for row in range(n_rows):
        ws1.write(row+1, 0, df_name.index[row][0], CF2)
        if row % 2 != 0:
            ws1.write(row+1, 1, df_name.index[row][1], FALTERNATING1)
        else:
            ws1.write(row+1, 1, df_name.index[row][1], FALTERNATING2)
    
    for col in range(2, n_cols-1):
        if col % 2 == 0:
            ws1.write(0, col+2, df_name.columns[col], PURPLE1)
        else:
            ws1.write(0, col+2, df_name.columns[col], PURPLE2)
    
    for row in range(n_rows+1):
        ws1.set_row(row, cell_format=(FALTERNATING1 if row%2==0 else FALTERNATING2))
    
    if sum(df_name["Total"]) < 65_530:
        for row in range(n_rows):
            for col in range(2, n_cols-1):
                pmid = str(df_name.iloc[row, col])
                if pmid != "" and row % 2 != 0:
                    ws1.write_url(row+1, col+2, f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", FALTERNATINGC, string=pmid)
                elif pmid != "" and row % 2 == 0:
                    ws1.write_url(row+1, col+2, f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", FALTERNATINGC2, string=pmid)
                else:
                    pass
    else:
        pass
    
    for row in range(n_rows):
        if row % 2 != 0:
            ws1.write(row+1, 3, int(df_name["Total"][row]), FALTERNATINGC)
            ws1.write(row+1, 2, df_name["confidence"][row], ccoloring_cs(workbook, f"{df_name['colors'][row]}"))
        else:
            ws1.write(row+1, 3, int(df_name["Total"][row]), FALTERNATINGC2)
            ws1.write(row+1, 2, df_name["confidence"][row], ccoloring_cs(workbook, f"{df_name['colors'][row]}"))
    
    for row in range(n_rows):
        for col in range(2, n_cols-1):
            pmid = str(df_name.iloc[row, col])
            if pmid.isdigit() and row % 2 != 0:
                ws1.write(row+1, col+2, int(pmid), FALTERNATINGC)
            elif pmid.isdigit() and row % 2 == 0:
                ws1.write(row+1, col+2, int(pmid), FALTERNATINGC2)
            else:
                pass
    
    FCS = falternating_color(workbook, "#F5F5F5", "#878787", "center", 12)
    FLOW = falternating_color(workbook, low_color, "#878787", "center", 12)
    FQLOW = falternating_color(workbook, quite_low_color, "#878787", "center", 12)
    FMOD = falternating_color(workbook, moderate_color, "#878787", "center", 12)
    FHIG = falternating_color(workbook, high_color, "#878787", "center", 12)
    FVHI = falternating_color(workbook, very_high_color, "#878787", "center", 12)
    
    ws1.write(n_rows+3, 0, "Confidence score:", FCS)
    
    ws1.write(n_rows+8, 0, "★☆☆☆☆", FLOW)
    ws1.write(n_rows+7, 0, "★★☆☆☆", FQLOW)
    ws1.write(n_rows+6, 0, "★★★☆☆", FMOD)
    ws1.write(n_rows+5, 0, "★★★★☆", FHIG)
    ws1.write(n_rows+4, 0, "★★★★★", FVHI)
    
    ws1.write(n_rows+8, 1, "Low", FLOW)
    ws1.write(n_rows+7, 1, "Quite Low", FQLOW)
    ws1.write(n_rows+6, 1, "Moderate", FMOD)
    ws1.write(n_rows+5, 1, "High", FHIG)
    ws1.write(n_rows+4, 1, "Very High", FVHI)
    
    return None


def make_xlsx(file_s):
    """Main function to run the coloring script.
    """
    
    if isinstance(file_s, tuple):
        output_path = "TCD-stressor_event.xlsx"
        writer = pd.ExcelWriter(output_path, engine="xlsxwriter")
        stressor_event = file_s[0]
        event_stressor = file_s[1]
        coloring(writer, stressor_event, "stressor-event")
        coloring(writer, event_stressor, "event-stressor")
    
    else:
        output_path = "TCD-event_event.xlsx"
        writer = pd.ExcelWriter(output_path, engine="xlsxwriter")
        coloring(writer, file_s, "event-event")
    
    writer.save()
    return None


if __name__ == "__main__":
    all_stressors = sys.argv[1]
    scoring_file = sys.argv[2]
    start_time = time.time()
    
    TCD = make_TCD(all_stressors, scoring_file)
    make_xlsx(TCD)
    
    print(f"\u2023 TCD-like file generated in {(time.time() - start_time):.3f} seconds.")
    