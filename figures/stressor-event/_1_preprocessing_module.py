#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:00:28 2022

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
import sys
from pathlib import Path

####################

# Get the parent directory of the directory containing this script.
#pp_directory = Path(__file__).resolve().parents[1].as_posix()
pp_directory = sys.argv[2]

def add_stressors(directory_name):
    """Concatenate all files from AOP-helpFinder output directory 
    and add a column 'stressor' in position 1.
    
    Args:
    -----
        directory_name: name of the directory storing AOP-helpFinder 
        output files (in .tsv format), usually "results/".
    
    Returns:
    --------
        A .csv file with all events and all stressors (column 1 = stressor)
        named "all_stressors.csv".
    
    Note:
    -----
        It is not a problem if the <directory_name> does not contain only .tsv files.
        <<!!>> Do not forget the "/" after the directory name ! <<!!>>
     """
    

    df = pd.DataFrame([])

    directory_name = os.path.join(pp_directory, directory_name)
    files = directory_name + "*.tsv"
    data = []
    col_names4 = ["pubdate", "title", "pmid", "event"]
    col_names5 = ["pubdate", "title", "pmid", "event", "abstract"]
    
    for file in glob.glob(files):
        FILE = file.split("/")[-1]
        if FILE.startswith("scoring"):
            pass
        else:
            with open(file) as f:
                reader = csv.reader(f, delimiter="\t")
                first_row = next(reader)
        
            if "abstract" in first_row:
                df = pd.read_csv(file, sep="\t", names=col_names5)
            else:
                df = pd.read_csv(file, sep="\t", names=col_names4)
        
            df.rename({"pmid": "PMID"}, axis=1, inplace=True)
            pattern = directory_name + "(.*?).tsv"
            stressor = re.search(pattern, file) 
            print(pattern)
            print(stressor.group(1))
            df.drop(df.head(1).index, inplace=True)
            df.insert(0, "stressor", stressor.group(1))
            data.append(df)
    
    data = pd.concat(data, ignore_index=True)
    data.sort_values(by="stressor", inplace=True)
    output_name = "all_stressors.csv"
    output_path = os.path.join(pp_directory, output_name)
    data.to_csv(output_path, sep=";", index=False)
    
    return None


def make_file4heatmap(directory_name):
    """Setting up the input file for the heatmap.
    
    Args:
    -----
        directory_name: name of the directory storing AOP-helpFinder 
        output files (in .tsv format), usually "results/".
    
    Returns:
    --------
        A .csv file "file4heatmap.csv" with stressors in columns and events in rows.
        The location (i,j) holds the number of papers linking the stressor i to the event j.
    """
    
    directory_name = os.path.join(pp_directory, directory_name)
    all_files = directory_name + '*.tsv'
    col_names4 = ["date", "title", "PMID", "event"]
    col_names5 = ["date", "title", "PMID", "event", "abstract"]
    res_tmp = pd.DataFrame()
    
    for file in glob.glob(all_files):
        FILE = file.split("/")[-1]
        if FILE.startswith("scoring"):
            pass
        else:
            stressor_name = file.split("/")[-1].split(".")[0]
            with open(file) as f:
                reader = csv.reader(f, delimiter="\t")
                first_row = next(reader)
            
            if "abstract" in first_row:
                table = pd.read_csv(file, sep="\t", names=col_names5)
            else:
                table = pd.read_csv(file, sep="\t", names=col_names4)
        
            table.drop(table.head(1).index, inplace=True)
            table["event"].str.strip
            table = table.assign(stressor=stressor_name)
            res = table.groupby(["stressor", "event"]).count()
            res_tmp = pd.concat([res_tmp, res])
        
    df_stressor_event = res_tmp.reset_index()
    df_stressor_event = df_stressor_event[["stressor", "event", "PMID"]].pivot(index="stressor", columns="event")
    df_stressor_event = df_stressor_event.fillna(0)
    df_stressor_event.columns = df_stressor_event.columns.droplevel()
    df_stressor_event = df_stressor_event.astype(int)
    df_stressor_event.index.name = "stressor \ event"
    
    df_stressor_event.sort_values(by="stressor \ event", inplace=True, key=lambda col: col.str.lower())
    
    output_name = "file4heatmap.csv"
    output_path = os.path.join(pp_directory, output_name)
    df_stressor_event.to_csv(output_path, sep=";")
    
    return None
    
