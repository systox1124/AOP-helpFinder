#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 09:39:19 2022

@author: <Thibaut Coustillet>
__lab__: SysTox Team (U1124 | T3S) FRANCE
__version__: 1.0.0
"""

####################

import _3_figures_module as fm
import sys
import time

####################


if __name__ == "__main__":
    all_stressors = sys.argv[1]
    file4heatmap = sys.argv[2]
    scoring_file = sys.argv[3]
    t0 = time.time()
    fm.plot_years(all_stressors)
    fm.plot_stressors(all_stressors)
    fm.plot_events(all_stressors)
    fm.plot_distribution4artciles(all_stressors)
    fm.plot_stressor_event_links(scoring_file)
    fm.plot_heatmaps(file4heatmap, scoring_file)
    print(f"\u2023 Figures generated in {(time.time() - t0):.3f} seconds.")
    