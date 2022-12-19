#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 09:39:19 2022

@author: <Thibaut Coustillet>
__lab__: SysTox Team (U1124 | T3S) FRANCE
__version__: 1.0.0
"""

####################

import _3EM_figures_module as fm
import sys
import time

####################

if __name__ == "__main__":
    event_event_file = sys.argv[1]
    score_file = sys.argv[2]
    t0 = time.time()
    fm.plot_years(event_event_file)
    fm.plot_event_event_links(score_file)
    print(f"\u2023 Figures generated in {(time.time() - t0):.3f} seconds.")
    