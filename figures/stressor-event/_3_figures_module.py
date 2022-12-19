#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 09:38:45 2022

@author: <Thibaut Coustillet>, <Florence Jornod>
__lab__: SysTox Team (U1124 | T3S) FRANCE
__version__: 1.0.0
"""

####################

import csv
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import upsetplot
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import ticker
from pathlib import Path
import sys

####################


# Get the parent directory of the directory containing this script.
pp_directory = sys.argv[4]


def plot_years(file_name):
    """Plot the number of articles published per year for a given data set.
    
    Args:
    -----
        file_name: A .csv file with AOPhF data: "all_stressors.csv".
    
    Returns:
    --------
        A barplot "distribution_years.pdf" showing the number of articles published per year.
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data = pd.read_csv(file_name, sep=";")
    data["pubdate"] = data["pubdate"].astype(str)
    data["pubdate"] = data["pubdate"].str.replace("-(.*)", "", regex=True)
    data = data.drop_duplicates(subset=["PMID"])
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
    plt.xticks(fontsize=8, rotation=90, ha='center')
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
    
    output_name = "figures/distribution_years.pdf"
    output_path = os.path.join(pp_directory, output_name)
    
    fig.savefig(output_path, bbox_inches="tight", format="pdf")
    plt.close(fig)
    
    return None


def roundup(x):
    """To round a number up to the next hundred.
    """
    y = x if x % 100 == 0 else x + 100 - x % 100
    return y


def plot_stressors(file_name):
    """Plot the number of stressor-event links found according to the stressor.
    
    Args:
    -----
        file_name: A .csv file with AOPhF data: "all_stressors.csv".
    
    Returns:
    --------
        A barplot "distribution_stressors.pdf" showing the number of 
        stressor-event links found according to the stressor.
    
    Note:
    -----
        If the number of stressors is > 30, then the figure displays only the first 30.
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data = pd.read_csv(file_name, sep=";")
    count_stressors = data["stressor"].value_counts()
    
    n_stressors = len(count_stressors)
    n_links = len(data)
    ind = np.arange(n_stressors)
    
    percentage = [0 for k in ind]
    tot_p = 0
    for s in ind:
        tot_p += count_stressors[s] / n_links * 100
        percentage[s] = tot_p
    
    diff_percentage = [percentage[0]]
    diff_percentage += [j-i for i,j in zip(percentage[:-1], percentage[1:])]
    
    # threshold = 90                                                                      # pour ne garder que les 90% 1ers
    if n_stressors > 30:
        # eX = bisect.bisect_left(percentage, threshold)                                  # pour ne garder que les 90% 1ers
        # count_stressors = count_stressors[:eX+1]                                        # pour ne garder que les 90% 1ers
        count_stressors = count_stressors[:30]                                            # pour ne garder que les 30 1ers
        threshold = float(f"{(np.sum(count_stressors.to_list()) / n_links * 100):.1f}")  # pour ne garder que les 30 1ers
    else:
        threshold = 100
    
    my_dpi = 192
    width = 0.55
    
    fig, ax1 = plt.subplots(figsize=(12, 4), dpi=my_dpi)
    ax1.set_title(f"Distribution of stressor-event links according to the {len(count_stressors)} most common stressors,\nrepresenting {threshold}% of the total data set ({n_stressors} stressors and {n_links} links)", 
                  fontsize=12, pad=20)
    
    #14203E
    count_stressors.plot.bar(ax=ax1, color="lightgrey", width=width, label="count of occurrences")
    
    ax1.grid(axis="x", linestyle=":", linewidth=0.2, alpha=0.8, color="grey")
    ax1.axis(ymin=0, ymax=roundup(max(count_stressors)))
    
    ax1.set_xlabel("stressor", fontsize=12, labelpad=12)
    ax1.set_ylabel("count of occurrences", fontsize=12, labelpad=12)
    plt.xticks(fontsize=9, rotation=90)
    plt.yticks(fontsize=10, rotation=0)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["bottom"].set_visible(True)
    ax1.spines["left"].set_visible(True)
    ax1.set_axisbelow(True)
    
    ax2 = ax1.twinx()
    
    ax2.plot(ind, [i for i in percentage], c="#F7A617", zorder=2)
    ax2.scatter(ind, [i for i in percentage], c="#F7A617", zorder=2, label="% of total links")
    
    ax2.axis(ymin=0, ymax=110)
    ax2.set_ylabel("% of total links", fontsize=12, labelpad=12)
    plt.yticks(fontsize=10, rotation=0, c="#F7A617")
    
    ax2.set_yticks(np.linspace(0, 100, 11))
    my_colors = ["#F7A617" for i in range(11)]
    for ticklabel, tickcolor in zip(ax2.get_yticklabels(), my_colors):
        ticklabel.set_color(tickcolor)
    
    y_text_dict = [0, 10, 20, 30, 40, "50%", 60, 70, 80, 90, 100]
    ax2.set_yticklabels(y_text_dict)
    for i in range(len(percentage)):
        p = str("{:.1f}%".format(diff_percentage[i]))
        ax2.annotate(p, (i, percentage[i]+4), c="darkorange", fontsize=7, ha="center", 
                     bbox={"facecolor": "#E5E5E5", "boxstyle": "round", 
                           "pad": 0.2, "edgecolor": "none", "alpha": 0.4})
    
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(True)
    ax2.spines["bottom"].set_visible(True)
    ax2.spines["left"].set_visible(False)
    
    ax3 = ax1.twinx()
    
    ax3.axis(ymin=0, ymax=110)
    ax3.axhline(y=50, color="#F7A617", linestyle="dotted", lw=0.75)
    
    ax3.set_yticks([])
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(True)
    ax3.spines["bottom"].set_visible(True)
    ax3.spines["left"].set_visible(False) 
    
    ax4 = ax1.twinx()
    
    ax4.bar_label(ax1.containers[0], fontsize=8)
    ax4.axis(ymin=0, ymax=roundup(max(count_stressors)))
    ax4.set_yticks([])
    
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)
    ax4.spines["bottom"].set_visible(False)
    ax4.spines["left"].set_visible(False)
    
    ax1.set_zorder(2)
    ax2.set_zorder(3)
    ax3.set_zorder(1)
    ax4.set_zorder(4)
    ax1.patch.set_visible(False)
    ax2.patch.set_visible(False)
    ax3.patch.set_visible(False)
    ax4.patch.set_visible(False)
    
    fig.legend(fontsize=9)
    
    output_name = "figures/distribution_stressors.pdf"
    output_path = os.path.join(pp_directory, output_name)
    
    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    
    return None


def plot_events(file_name):
    """Plot the number of stressor-event links found according to the event.
    
    Args:
    -----
        file_name: A .csv file with AOPhF data: "all_stressors.csv".
    
    Returns:
    --------
        A barplot "distribution_events.pdf" showing the number of stressor-event links 
        found according to the event.
    
    Note:
    -----
        If the number of events is > 30, then the figure displays only the first 30.
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data = pd.read_csv(file_name, sep=";")
    count_events = data["event"].value_counts()
    
    n_links = len(data)
    n_distinct_events = data["event"].nunique()
    
    ind = np.arange(len(count_events))
    percentage = [0 for k in ind]
    tot_p = 0
    for s in ind:
        tot_p += count_events[s] / n_links * 100
        percentage[s] = float("{:.1f}".format(tot_p))
    
    diff_percentage = [percentage[0]]
    diff_percentage += [j-i for i,j in zip(percentage[:-1], percentage[1:])]
    
    # threshold = 80                                                                   # pour ne garder que les 80% 1ers
    if n_distinct_events > 30:
        # eX = bisect.bisect_left(percentage, threshold)                               # pour ne garder que les 80% 1ers
        # count_events = count_events[:eX+1]                                           # pour ne garder que les 80% 1ers
        count_events = count_events[:30]                                               # pour ne garder que les 30 1ers
        threshold = float(f"{(np.sum(count_events.to_list()) / n_links * 100):.1f}")  # pour ne garder que les 30 1ers
    else:
        threshold = 100
    
    my_dpi = 192
    fig, ax1 = plt.subplots(figsize=(15, 4), dpi=my_dpi)
    ax1.set_title(f"Distribution of stressor-event links according to the {len(count_events)} most common events,\nrepresenting {threshold}% of the total data set ({n_distinct_events} distinct events and {n_links} links)", 
                  fontsize=12, pad=20)
    
    #14203E
    count_events.plot.bar(ax=ax1, color="lightgrey", label="count of occurrences")
    ax1.grid(axis="x", linestyle=":", linewidth=0.2, alpha=0.8, color="grey")
    ax1.axis(ymin=0, ymax=roundup(max(count_events)))
    
    ax1.set_xlabel("event", fontsize=12, labelpad=12)
    ax1.set_ylabel("count of occurrences", fontsize=12, labelpad=12)
    ID = [f"#{i+1}" for i in range(len(count_events))]
    plt.xticks(range(len(ID)), ID, fontsize=10, rotation=0)
    
    # ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
    # ax1.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    # ax1.xaxis.set_minor_formatter(ticker.IndexFormatter(ID))
    # ax1.tick_params(axis="x", which="minor", length=16, width=1.5)
    # ax1.tick_params(axis="x", which="both", color="lightgrey")
    # ID2 = ["#{}".format(2*i) for i in range(len(ID[1::2])+1)]
    # ax1.set_xticklabels(ID2, c="k", fontsize=10, minor=True)
    
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["bottom"].set_visible(True)
    ax1.spines["left"].set_visible(True)
    ax1.set_axisbelow(True)
    
    ax2 = ax1.twinx()
    
    ax2.plot(ind, [i for i in percentage], c="#F7A617", zorder=2)
    ax2.scatter(ind, [i for i in percentage], c="#F7A617", zorder=2, label="% of total links")
    
    ax2.axis(ymin=0, ymax=110)
    ax2.set_ylabel("% of total links", fontsize=12, labelpad=12)
    plt.yticks(fontsize=10, rotation=0, c="#F7A617")
    
    ax2.set_yticks(np.linspace(0, 100, 11))
    my_colors = ["#F7A617" for i in range(11)]
    for ticklabel, tickcolor in zip(ax2.get_yticklabels(), my_colors):
        ticklabel.set_color(tickcolor)
    
    y_text_dict = [0, 10, 20, 30, 40, "50%", 60, 70, 80, 90, 100]
    ax2.set_yticklabels(y_text_dict)
    for i in range(len(percentage)):
        p = str("{:.1f}%".format(diff_percentage[i]))
        ax2.annotate(p, (i, percentage[i]+4), c="darkorange", fontsize=7, ha="center", 
                     bbox={"facecolor": "#E5E5E5", "boxstyle": "round", 
                           "pad": 0.2, "edgecolor": "none", "alpha": 0.4})
    
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(True)
    ax2.spines["bottom"].set_visible(True)
    ax2.spines["left"].set_visible(False)
    
    ax3 = ax1.twinx()
    
    ax3.axis(ymin=0, ymax=110)
    ax3.axhline(y=50, color="#F7A617", linestyle="dotted", lw=0.75)
    
    ax3.set_yticks([])
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(True)
    ax3.spines["bottom"].set_visible(True)
    ax3.spines["left"].set_visible(False) 
    
    ax4 = ax1.twinx()
    
    ax4.bar_label(ax1.containers[0], fontsize=8)
    ax4.axis(ymin=0, ymax=roundup(max(count_events)))
    ax4.set_yticks([])
    
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)
    ax4.spines["bottom"].set_visible(False)
    ax4.spines["left"].set_visible(False)
    
    ax1.set_zorder(2)
    ax2.set_zorder(3)
    ax3.set_zorder(1)
    ax4.set_zorder(4)
    ax1.patch.set_visible(False)
    ax2.patch.set_visible(False)
    ax3.patch.set_visible(False)
    ax4.patch.set_visible(False)
    
    fig.legend(fontsize=9)
    patch = []
    for e, i in zip(count_events.index, ID):
        patch += [mpatches.Patch(color="w", label=f"{i} : {e}")]
    
    plt.legend(handles=patch, 
               loc="lower left",
                bbox_to_anchor=(0, [-len(count_events)/10 if len(count_events)<10 else -1][0]), 
               ncol=((len(count_events)-1) // 10)+1, 
               handlelength=0, handletextpad=0, fancybox=0)
    
    output_name = "figures/distribution_events.pdf"
    output_path = os.path.join(pp_directory, output_name)
    
    fig.savefig(output_path, bbox_inches="tight", format="pdf")
    plt.close(fig)
    
    return None


def plot_distribution4artciles(file_name):
    """Plot the number of articles with X stressors found and Y events found.
    
    Args:
    -----
        file_name: A .csv file with AOPhF data: "all_stressors.csv".
    
    Returns:
    --------
        A plot "distribution_articles.pdf" showing the number of articles found 
        for a couple (#stressors, #events).
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data = pd.read_csv(file_name, sep=";")
    col_names = ["PMID", "stressor", "event"]
    data = data[col_names]
    n_articles = data["PMID"].nunique()
    
    output_df = []
    for pmid in data["PMID"].drop_duplicates():
        data_1p = data.loc[data["PMID"] == pmid]
        stressors_1p = data_1p["stressor"].nunique()
        events_1p = data_1p["event"].nunique()
        
        d1 = pd.DataFrame([[pmid, stressors_1p, events_1p]], columns=col_names)
        output_df += [d1]
    
    output_df = pd.concat(output_df)
    output_df["Total"] = output_df["stressor"] + output_df["event"]
    output_df.sort_values(by="Total", ascending=False, inplace=True)
    output_df["PMID"] = output_df["PMID"].astype(str)
    output_df.reset_index(drop=True, inplace=True)
    
    data4plot = pd.crosstab(output_df["event"], output_df["stressor"])
    data4plot = data4plot.reindex(index=data4plot.index[::-1])
    
    add_zeros_st = sorted(list(set(range(output_df["stressor"].max()+1)) - set(list(data4plot.columns))))[1::]
    for c in add_zeros_st:
        data4plot.insert(c-1, c, [0 for i in range(len(list(data4plot.index)))])
    
    add_zeros_en = sorted(list(set(list(range(max(data4plot.index)+1))) - set(list(data4plot.index))))[1::]
    for c in add_zeros_en:
        data4plot.loc[c] = [0 for i in range(len(list(data4plot.columns)))]
    
    data4plot.sort_index(inplace=True, ascending=False)
    
    #####fig, ax = plt.subplots(figsize=(5, 5))
    fig, ax = plt.subplots(figsize=(len(data4plot.columns), len(data4plot.columns)))
    ax.set_title(f"Distribution of the {n_articles} articles according to\nthe number of stressors and events detected in each", 
                 fontsize=12, pad=20)
    
    # cbar_kws = {"label": "number of articles", "pad": 0.22, "orientation" : "horizontal"}
    # colors2 = sns.color_palette("YlOrRd", as_cmap=True)
    
    sns.heatmap(data4plot, cmap="icefire", 
                annot=True, fmt='g', square=True,
                cbar=False, 
                # cbar_kws=cbar_kws, 
                linewidths=0.1, 
                linecolor="gainsboro",
                mask=(data4plot == 0))
    
    ax.plot((1), (output_df["event"].max()), ls="", marker=">", ms=5, color="k",
            transform=ax.get_yaxis_transform(), clip_on=False)
    ax.plot((0), (1), ls="", marker="^", ms=5, color="k",
            transform=ax.get_xaxis_transform(), clip_on=False)
    
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.figure.axes[-1].yaxis.label.set_size(12)
    
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.xlabel('count of stressors', fontsize=11, labelpad=15)
    plt.ylabel('count of events', fontsize=11, labelpad=15)
    
    output_name = "figures/distribution_articles.pdf"
    output_path = os.path.join(pp_directory, output_name)
    
    fig.savefig(output_path, bbox_inches="tight", format="pdf")
    plt.close(fig)
    
    return None


####################

low_color = "#fe9929"
quite_low_color = "#fec44f"
moderate_color = "#e6f598"
high_color = "#91cf60"
very_high_color = "#1a9850"

####################


def plot_stressor_event_links(file_name):
    """Plot the stressor-event links found according to the number of papers.
    
    Args:
    -----
        file_name: A .csv file with AOPhF scoring data: "scoring_stressor-event.csv".
    
    Returns:
    --------
        A barplot "distribution_scores.pdf" showing the most common stressor-event links.
    
    Note:
    -----
        If the number of stressor-event links is > 30, 
        then the figure displays only the first 30.
        The color indicates the confidence score.
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data_scoring = pd.read_csv(file_name, sep="\t")
    data_scoring.insert(2, "stressor-event link", (data_scoring["Stressor"].astype(str) + " + " + data_scoring["Event"].astype(str)))
    data_scoring.sort_values(by="Link", ascending=False, inplace=True)
    
    def make_confidence_colors(row):
        """To match a color to each confidence level.
        """
        
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
    
    #####n_links = len(data_scoring)
    n_linksT = len(data_scoring)
    sum_links = np.sum(data_scoring["Link"])
    #####ind = np.arange(n_links)
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
    #####if n_links > cut_off:
    if n_linksT > cut_off:
        data_scoring4plot = data_scoring4plot.head(cut_off)
        n_links = len(data_scoring4plot)
        threshold = float(f"{(np.sum(data_scoring4plot['Link'].to_list()) / sum_links * 100):.1f}")
    else:
        threshold = 100
        n_links = n_linksT #####
    
    #####if n_links < 25:
    #####    size2 = n_links//2
   ##### else:
     #####   size2 = 9

    if n_links < 8:
        size2 = 3
    elif 8 <= n_links < 15:
        size2 = 5
    else:
        size2 = 9
    
    my_dpi = 192
    fig, ax = plt.subplots(figsize=(8, size2), dpi=my_dpi)
    
    data_scoring4plot.sort_values(by="Link", ascending=True, inplace=True)
    confidence_colors = list(data_scoring4plot["colors"])
    
    data_scoring4plot.plot.barh(ax=ax, x="stressor-event link", y="Link", color=confidence_colors)
    ax.set_title(f"Distribution of stressor-event links according to the {len(data_scoring4plot)} most common links,\nrepresenting {threshold}% of the total data set ({n_links} different links and {sum_links} total links)", 
                  fontsize=12, pad=20)
    ax.grid(axis="x", linestyle=":", linewidth=0.2, alpha=0.8, color="grey")
    
    ax.bar_label(ax.containers[0], padding=2, fontsize=10)
    ax.set_xlabel("count of occurences", fontsize=12, labelpad=12)
    ax.set_ylabel("stressor-event link", fontsize=12, labelpad=12)
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
    
    output_name = "figures/distribution_scores.pdf"
    output_path = os.path.join(pp_directory, output_name)
    
    fig.savefig(output_path, bbox_inches="tight", format="pdf")
    plt.close(fig)
    
    return None


####################
###### HEATMAP #####
####################


def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp


def fig_heatmap_cut(df_stressor_event, figsize, x, colors):
    """Plot the stressor-event links found according to the number of papers.
    
    Args:
    -----
        df_stressor_event: A .csv file with the number of links for each stressor-event pair. 
        = output of 'make_file4heatmap' function (_1_preprocessing_module.py).
        = "file4heatmap.csv".
        figsize: size of the figure.
        x: maximum number of events for each heatmap.
        colors: color scale for the heatmap.
    
    Returns:
    --------
        A heatmap every x events.
    
    Note:
    -----
        Colors are AOP-helpFinder official colors and the color scale is based
        to the number of papers (=number of links between 1 event and 1 stressor).
        The location (i,j) holds the number of papers linking the event i to the stressor j.
    """
    
    flat_list = [item for sdf_np in df_stressor_event.to_numpy() for item in sdf_np]
    vmax = np.quantile(flat_list, 0.99)
    df_stressor_event = df_stressor_event.replace(0, np.nan)
    
    i = 0
    end = df_stressor_event.shape[1]
    j = end
    if end > x:
        j=x
    
    fig = plt.figure(figsize=figsize)
    sns.set(font_scale=1.25)
    htmp = sns.heatmap(df_stressor_event.T[i:j], annot=True, fmt="g", linecolor="whitesmoke", 
                       annot_kws={'verticalalignment':'center_baseline'}, 
                       yticklabels=True, vmax=vmax, linewidths=.5, cmap=get_continuous_cmap(colors), 
                       cbar_kws={"shrink": 1, "label": "count of links (= count of abstracts)"})
    htmp.set_facecolor('xkcd:white')
    htmp.tick_params(left=True, bottom=True)
    htmp.set_xlabel("stressor", fontsize=22, labelpad=20)
    htmp.set_ylabel("event", fontsize=22, labelpad=20)
    
    output_name = f"figures/heatmap_{i}.pdf"
    output_path = os.path.join(pp_directory, output_name)
    
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    while(j<end):
        i=j+1
        j=j+x
        
        fig, ax  = plt.subplots(figsize=figsize)
        sns.set(font_scale=1.25)
        htmp = sns.heatmap(df_stressor_event.T[i-1:j], annot=True, fmt="g", linecolor="whitesmoke", 
                           annot_kws={'verticalalignment':'center_baseline'}, ax=ax,
                           yticklabels=True, vmax=vmax, linewidths=.5, cmap=get_continuous_cmap(colors), 
                           cbar_kws={"shrink": 1, "label": "count of links (= count of abstracts)"})
        
        htmp.set_facecolor('xkcd:white')
        htmp.tick_params(left=True, bottom=True)
        htmp.set_xlabel("stressor", fontsize=22, labelpad=25)
        htmp.set_ylabel("event", fontsize=22, labelpad=25)
        
        output_name2 = f"figures/heatmap_{i}.pdf"
        output_path2 = os.path.join(pp_directory, output_name2)
        
        plt.savefig(output_path2, bbox_inches='tight', dpi=300)
        plt.close(fig)
    
    return None


def fig_heatmap_tot(df_stressor_event, figsize, colors):
    """Plot the stressor-event links found according to the number of papers.
    
    Args:
    -----
        df_stressor_event: A .csv file with the number of links for each stressor-event pair. 
        = output of 'make_file4heatmap' function (_1_preprocessing_module.py).
        = "file4heatmap.csv".
        figsize: size of the figure.
        colors: color scale for the heatmap.
    
    Returns:
    --------
        A global heatmap (all stressor-event pairs).
    
    Note:
    -----
        Colors are AOP-helpFinder official colors and the color scale is based
        to the number of papers (=number of links between 1 event and 1 stressor).
        The location (i,j) holds the number of papers linking the event i to the stressor j.
    """
    
    flat_list = [item for sdf_np in df_stressor_event.to_numpy() for item in sdf_np]
    vmax = np.quantile(flat_list, 0.99)
    df_stressor_event = df_stressor_event.replace(0, np.nan)
    n_links = df_stressor_event.values.sum()


    fig, ax = plt.subplots(figsize=figsize)
    htmp = sns.heatmap(df_stressor_event.T, annot=True,fmt="g", linecolor="whitesmoke", 
                       annot_kws={'verticalalignment':'center_baseline'}, ax=ax,
                       yticklabels=True, vmax=vmax, linewidths=.5, cmap=get_continuous_cmap(colors))
    
    cbar = htmp.collections[0].colorbar
    cbar.set_label("count of links (= count of abstracts)", labelpad=18)
    htmp.tick_params(left=True, bottom=True)
    htmp.set_xlabel("stressor", labelpad=18)
    htmp.set_ylabel("event", labelpad=18)
    ax.set_title(f"Distribution of the {n_links} stressor-event links according to stressors and events", 
                  pad=10)

    output_name = "figures/heatmap_total.pdf"
    output_path = os.path.join(pp_directory, output_name)
    
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    return None


cs_colors = (low_color, quite_low_color, moderate_color, high_color, very_high_color)
cmap_confidence = LinearSegmentedColormap.from_list('Custom', cs_colors, len(cs_colors))


def fig_heatmap_cut_cs(file_name, figsize, x):
    """Plot the stressor-event links found according to the confidence score.
    
    Args:
    -----
        file_name: A .csv file with the number of links for each stressor-event pair
        and with AOPhF scoring data = "scoring_stressor-event.csv".
        figsize: size of the figure.
        x: maximum number of events for each heatmap.
    
    Returns:
    --------
        A heatmap every x events.
    
    Note:
    -----
        Color scale is based to the confidence score.
        The location (i,j) holds the number of papers linking the event i to the stressor j.
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data = pd.read_csv(file_name, sep="\t")
    data = data.replace(0, np.nan)
    
    def confidence_to_number(row):
        """
        """
        
        val = -1
        if row["Confidence"] == "Low":
            val = 0
        elif row["Confidence"] == "Quite low":
            val = 1
        elif row["Confidence"] == "Moderate":
            val = 2
        elif row["Confidence"] == "High":
            val = 3
        elif row["Confidence"] == "Very High":
            val = 4
        
        return val
    
    data["confidence_number"] = data.apply(confidence_to_number, axis=1)
    df4heatmap = data.pivot(index="Stressor", columns="Event", values="confidence_number")
    labels = data.pivot(index="Stressor", columns="Event", values="Link").to_numpy()
    df4heatmap.index.name = "stressor \ event"
    df4heatmap.columns.name = ""
    
    i = 0
    end = df4heatmap.shape[1]
    j = end
    if end > x:
        j=x
    
    fig = plt.figure(figsize=figsize)
    sns.set(font_scale=1.25)
    htmp = sns.heatmap(df4heatmap.T[i:j], annot=labels.T[i:j], fmt="g", linecolor="whitesmoke", 
                       annot_kws={'verticalalignment':'center_baseline'}, 
                       yticklabels=True, linewidths=.5, cmap=cmap_confidence, 
                       vmin=0, vmax=4)
    
    cbar = htmp.collections[0].colorbar
    cbar.set_ticks([0.4, 1.2, 2, 2.8, 3.6])
    cbar.set_ticklabels(["Low", "Quite Low", "Moderate", "High", "Very High"])
    cbar.set_label("Confidence score", labelpad=18)
    htmp.set_facecolor('xkcd:white')
    htmp.tick_params(left=True, bottom=True)
    htmp.set_xlabel("stressor", fontsize=22, labelpad=20)
    htmp.set_ylabel("event", fontsize=22, labelpad=20)
    
    output_name = f"figures/heatmap_{i}_cs.pdf"
    output_path = os.path.join(pp_directory, output_name)
    
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    while(j<end):
        i=j+1
        j=j+x
        
        fig, ax  = plt.subplots(figsize=figsize)
        sns.set(font_scale=1.25)
        htmp = sns.heatmap(df4heatmap.T[i-1:j], annot=labels.T[i-1:j], fmt="g", linecolor="whitesmoke", 
                           annot_kws={'verticalalignment':'center_baseline'}, ax=ax,
                           yticklabels=True, linewidths=.5, cmap=cmap_confidence, 
                           vmin=0, vmax=4)
        
        cbar = htmp.collections[0].colorbar
        cbar.set_ticks([0.4, 1.2, 2, 2.8, 3.6])
        cbar.set_ticklabels(["Low", "Quite Low", "Moderate", "High", "Very High"])
        cbar.set_label("Confidence score", labelpad=18)
        htmp.set_facecolor('xkcd:white')
        htmp.tick_params(left=True, bottom=True)
        htmp.set_xlabel("stressor", fontsize=22, labelpad=25)
        htmp.set_ylabel("event", fontsize=22, labelpad=25)
        
        output_name2 = f"figures/heatmap_{i}_cs.pdf"
        output_path2 = os.path.join(pp_directory, output_name2)
        
        plt.savefig(output_path2, bbox_inches='tight', dpi=300)
        plt.close(fig)
    
    return None


def fig_heatmap_tot_cs(file_name, figsize):
    """Plot the stressor-event links found according to the confidence score.
    
    Args:
    -----
        file_name: A .csv file with the number of links for each stressor-event pair
        and with AOPhF scoring data = "scoring_stressor-event.csv".
        figsize: size of the figure.
    
    Returns:
    --------
        A global heatmap (all stressor-event pairs).
    
    Note:
    -----
        Color scale is based to the confidence score.
        The location (i,j) holds the number of papers linking the event i to the stressor j.
    """
    
    file_name = os.path.join(pp_directory, file_name)
    data = pd.read_csv(file_name, sep="\t")
    data = data.replace(0, np.nan)
    
    def confidence_to_number(row):
        """
        """
        
        val = -1
        if row["Confidence"] == "Low":
            val = 0
        elif row["Confidence"] == "Quite low":
            val = 1
        elif row["Confidence"] == "Moderate":
            val = 2
        elif row["Confidence"] == "High":
            val = 3
        elif row["Confidence"] == "Very High":
            val = 4
        
        return val
    
    data["confidence_number"] = data.apply(confidence_to_number, axis=1)
    df4heatmap = data.pivot(index="Stressor", columns="Event", values="confidence_number")
    labels = data.pivot(index="Stressor", columns="Event", values="Link").to_numpy()
    df4heatmap.index.name = "stressor \ event"
    df4heatmap.columns.name = ""
    n_links = sum(data["Link"])

    fig, ax = plt.subplots(figsize=figsize)
    htmp = sns.heatmap(df4heatmap.T, annot=labels.T, fmt="g", linecolor="whitesmoke", 
                       annot_kws={'verticalalignment':'center_baseline'}, ax=ax,
                       yticklabels=True, linewidths=.5, cmap=cmap_confidence, 
                       vmin=0, vmax=4)
    
    cbar = htmp.collections[0].colorbar
    cbar.set_ticks([0.4, 1.2, 2, 2.8, 3.6])
    cbar.set_ticklabels(["Low", "Quite Low", "Moderate", "High", "Very High"])
    cbar.set_label("Confidence score", labelpad=18)
    
    htmp.tick_params(left=True, bottom=True)
    htmp.set_xlabel("stressor", labelpad=18)
    htmp.set_ylabel("event", labelpad=18)
    ax.set_title(f"Distribution of the {n_links} stressor-event links according to stressors and events", 
                  pad=10)
    output_name = "figures/heatmap_total_cs.pdf"
    output_path = os.path.join(pp_directory, output_name)
    
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    return None


def plot_heatmaps(file_name, scoring_file):
    """To plot all heatmaps.
    
    Args:
    -----
        file_name: A .csv file with the number of links for each stressor-event pair. 
        = output of 'make_file4heatmap' function (_1_preprocessing_module.py).
        = "file4heatmap.csv".
        scoring_file: A .csv file with the number of links for each stressor-event pair
        and with AOPhF scoring data = "scoring_stressor-event.csv".
    
    Returns:
    --------
        2 types of heatmaps (color scale based to number of links and confidence score).
        The location (i,j) holds the number of papers linking the event i to the stressor j.
    """
    
    file_name = os.path.join(pp_directory, file_name)
    df_stressor_event = pd.read_csv(file_name, sep=";", index_col="stressor \ event")
    n_stressors = len(df_stressor_event)
    n_events = len(df_stressor_event.columns)
    
    AOPhF_colors = ["#EAEAEA", "#FCA616", "#14203E"]
    cut_off = 40
    
    if n_events > 30 or n_stressors > 30:
        fig_heatmap_tot(df_stressor_event, (25,25), AOPhF_colors)
        fig_heatmap_tot_cs(scoring_file, (25,25))
    elif n_events < 6 and n_stressors < 6:
         fig_heatmap_tot(df_stressor_event, (5,5), AOPhF_colors)
         fig_heatmap_tot_cs(scoring_file, (5,5))
    else:
        # fig_heatmap_tot(df_stressor_event, (n_stressors, n_events//2), AOPhF_colors)
        # fig_heatmap_tot_cs(scoring_file, (n_stressors, n_events//2))
        fig_heatmap_tot(df_stressor_event, (n_stressors, (n_events//2)+1), AOPhF_colors)
        fig_heatmap_tot_cs(scoring_file, (n_stressors, (n_events//2)+1))
    
    if df_stressor_event.shape[1] > cut_off:
        fig_heatmap_cut(df_stressor_event, (20,20), cut_off, AOPhF_colors)
        fig_heatmap_cut_cs(scoring_file, (20,20), cut_off)
    
    return None
    