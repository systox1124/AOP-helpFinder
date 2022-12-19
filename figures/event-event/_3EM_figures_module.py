#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 09:38:45 2022

@author: <Thibaut Coustillet>
__lab__: SysTox Team (U1124 | T3S) FRANCE
__version__: 1.0.0
"""

####################

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys
from matplotlib import ticker
from pathlib import Path

####################


# Get the parent directory of the directory containing this script.
pp_directory = sys.argv[3]


def plot_years(file_name):
    """Plot the number of articles published per year for a given data set.
    
    Args:
    -----
        file_name: A .csv file with AOPhF data: "event-event.tsv".
    
    Returns:
    --------
        A barplot "distribution_years.pdf" showing the number of articles published per year.
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data = pd.read_csv(file_name, sep="\t")
    data["pubdate"] = data["pubdate"].astype(str)
    data["pubdate"] = data["pubdate"].str.replace("-(.*)", "", regex=True)
    data = data.drop_duplicates(subset=["pmid"])
    data.reset_index(drop=True, inplace=True)
    papers = data["pubdate"].value_counts().sort_index()
    n_dates = data["pubdate"].nunique()
    date_size = 13 - ((((n_dates//10)*10+5)/5)//2+1)
    
    my_dpi = 192
    fig, ax = plt.subplots(figsize=(14, 3), dpi=my_dpi)
    
    papers.plot.bar(ax=ax, color="#F7A617")
    ax.grid(axis="y", linestyle=":", linewidth=0.2, alpha=0.8, color="grey")
    ax.bar_label(ax.containers[0], fontsize=date_size)
    
    plt.xlabel("year of publication", fontsize=11, labelpad=12)
    plt.ylabel("number of papers", fontsize=11, labelpad=12)
    plt.xticks(fontsize=9, rotation=90, ha='center')
    plt.yticks(fontsize=9, rotation=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(False)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    
    plt.title(f"Distribution of articles by year of publication (n = {len(data)} articles)", 
              fontsize=12, pad=20)
    
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    # ax.xaxis.set_minor_formatter(ticker.IndexFormatter(papers.index.tolist()))
    # ax.tick_params(axis="x", which="minor", length=18)
    # ax.tick_params(axis="x", which="both", color="lightgrey")
    # ax.autoscale(enable=True, axis="x", tight=True)
    
    output_name = "distribution_years.pdf"
    output_path = os.path.join(pp_directory, output_name)
    
    fig.savefig(output_path, bbox_inches="tight", format="pdf")
    plt.close()
    
    return None


####################

low_color = "#fe9929"
quite_low_color = "#fec44f"
moderate_color = "#e6f598"
high_color = "#91cf60"
very_high_color = "#1a9850"

####################


def plot_event_event_links(file_name):
    """Plot the event-event links found according to the number of papers.
    
    Args:
    -----
        file_name: A .csv file with AOPhF scoring data: "scoring_event-event.csv".
    
    Returns:
    --------
        A barplot "distribution_scores.pdf" showing the most common event-event links.
    
    Note:
    -----
        If the number of event-event links is > 30, 
        then the figure displays only the first 30.
        The color indicates the confidence score.
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data_scoring = pd.read_csv(file_name, sep="\t")
    data_scoring.insert(2, "event-event link", (data_scoring["Event 1"].astype(str) + " + " + data_scoring["Event 2"].astype(str)))
    data_scoring.sort_values(by="Link", ascending=False, inplace=True)
    
    def make_confidence_colors(row):
        """To match a color to each confidence level.
        """
        
        val = ""
        if row["Confidence"] == "Low":
            val = low_color
        elif row["Confidence"] == "Quite low":
                val = quite_low_color
        elif row["Confidence"] == "Moderate":
            val = moderate_color
        elif row["Confidence"] == "High":
            val = high_color
        elif row["Confidence"] == "Very High":
            val = very_high_color
        else:
            val = "#636363"
        
        return val
    
    data_scoring["colors"] = data_scoring.apply(make_confidence_colors, axis=1)
    
    n_linksT = len(data_scoring)
    sum_links = np.sum(data_scoring["Link"])
    ind = np.arange(n_linksT)
    data_scoring4plot = data_scoring.copy()
    
    percentage = [0 for k in ind]
    tot_p = 0
    for s in ind:
        tot_p += data_scoring["Link"][s] / sum_links * 100
        percentage[s] = tot_p
    
    diff_percentage = [percentage[0]]
    diff_percentage += [j-i for i,j in zip(percentage[:-1], percentage[1:])]
    
    cut_off = 30
    if n_linksT > cut_off:
        data_scoring4plot = data_scoring4plot.head(cut_off)
        n_links = len(data_scoring4plot)
        threshold = float(f"{(np.sum(data_scoring4plot['Link'].to_list()) / sum_links * 100):.1f}")
    else:
        threshold = 100
        n_links = n_linksT
    
    if n_links < 11:
        size2 = 4
    elif 10 < n_links < 25:
        size2 = n_links//3
    else:
        size2 = 9
    
    my_dpi = 192
    fig, ax = plt.subplots(figsize=(8, size2), dpi=my_dpi)
    
    data_scoring4plot.sort_values(by="Link", ascending=True, inplace=True)
    confidence_colors = list(data_scoring4plot["colors"])
    
    data_scoring4plot.plot.barh(ax=ax, x="event-event link", y="Link", color=confidence_colors)
    ax.set_title(f"Distribution of event-event links according to the {len(data_scoring4plot)} most common links,\nrepresenting {threshold}% of the total data set ({n_linksT} different links and {sum_links} total links)", 
                  fontsize=12, pad=20)
    ax.grid(axis="x", linestyle=":", linewidth=0.2, alpha=0.8, color="grey")
    
    ax.bar_label(ax.containers[0], padding=2)
    ax.set_xlabel("count of occurences", fontsize=12, labelpad=12)
    ax.set_ylabel("event-event link", fontsize=12, labelpad=12)
    plt.xticks(fontsize=11, rotation=0)
    plt.yticks(fontsize=11, rotation=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_axisbelow(True)
    
    low_patch = mpatches.Patch(color=low_color, label="Low")
    quite_low_patch = mpatches.Patch(color=quite_low_color, label="Quite Low")
    moderate_patch = mpatches.Patch(color=moderate_color, label="Moderate")
    high_patch = mpatches.Patch(color=high_color, label="High")
    very_high_patch = mpatches.Patch(color=very_high_color, label="Very High")
    
    legend = ax.legend(handles=[very_high_patch, high_patch, moderate_patch, quite_low_patch, low_patch], 
                        loc="center right", prop={'size': 11}, bbox_to_anchor=(1.25, 0.5))
    legend.set_title("Confidence score", prop={'size':12})
    
    output_name = "distribution_scores.pdf"
    output_path = os.path.join(pp_directory, output_name)
    
    fig.savefig(output_path, bbox_inches="tight", format="pdf")
    return None