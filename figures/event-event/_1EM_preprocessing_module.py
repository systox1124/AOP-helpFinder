#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 10:02:10 2022

@author: <Thibaut Coustillet>
__lab__: SysTox Team (U1124 | T3S) FRANCE
__version__: 1.0.0
"""

####################

import csv
import numpy as np
import os
import pandas as pd
import re
import sys
from pathlib import Path

####################

# Get the parent directory of the directory containing this script.
pp_directory = sys.argv[2]

def preprocess(file_name):
    """Formatting the file to count the event-event links.
    
    Args:
    -----
        file_name: A .csv file with event-event links and associated scores
        = "scoring_event-event.csv".
    
    Returns:
    --------
        A .csv file (double-input table) with events in rows ans columns.
        The location (i,j) holds the number of papers linking the event i to the event j.
        
    Note:
    -----
        The .csv output file is a lower triangular matrix with the 1st diagonal 
        equal to 0 since data are symmetrical.
    """
    
    data_scoring = pd.read_csv(file_name, sep="\t")
    data_scoring2 = data_scoring.copy()
    data_scoring2.rename({"Event 1":"Event 2", "Event 2":"Event 1"}, axis="columns", inplace=True)
    data_scoring3 = pd.concat([data_scoring, data_scoring2], axis=0)
    data_scoring3.reset_index(drop=True, inplace=True)
    
    df = data_scoring3.pivot(index="Event 1", columns="Event 2", values="Link")
    df.fillna(0, inplace=True)
    output_file = df.where(np.tril(np.ones(df.shape), -1).astype(bool))
    output_file.fillna(-1, inplace=True)
    output_file.index.name = "event1 \ event2"
    output_file = output_file.astype(int)
    output_file.replace(-1, "x", inplace=True)
    output_file.columns.name = ""
    
    output_name = "event_event_links.csv"
    output_path = os.path.join(pp_directory, output_name)
        
    return output_file
    