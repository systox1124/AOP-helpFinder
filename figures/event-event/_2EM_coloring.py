#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:36:34 2022

@author: <Thibaut Coustillet>
__lab__: SysTox Team (U1124 | T3S) FRANCE
__version__: 1.0.0
"""

####################

import _2EM_coloring_module as cm
import sys
import time

####################

if __name__ == "__main__":
    event_event_file = sys.argv[1]
    scoring_file = sys.argv[2]    
    t0 = time.time()
    cm.coloring_csv(event_event_file)
    cm.coloring_event_event_links(scoring_file)
    print(f"\u2023 Files coloring: done in {(time.time() - t0):.3f} seconds.")
    