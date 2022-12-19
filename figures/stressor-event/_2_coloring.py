#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:36:34 2022

@author: <Thibaut Coustillet>
__lab__: SysTox Team (U1124 | T3S) FRANCE
__version__: 1.0.0
"""

####################

import _2_coloring_module as cm
import sys
import time

####################

if __name__ == "__main__":
    all_stressors = sys.argv[1]
    output_AOPhF = sys.argv[2]
    file4heatmap = sys.argv[3]
    data_scoring_stressor_event = sys.argv[4]
    t0 = time.time()
    cm.coloring_csv(all_stressors)
    cm.coloring_alls(output_AOPhF)
    cm.coloring_stressor_event_links(file4heatmap, data_scoring_stressor_event)
    print(f"\u2023 Files coloring: done in {(time.time() - t0):.3f} seconds.")
    