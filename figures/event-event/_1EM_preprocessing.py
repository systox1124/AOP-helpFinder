#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 10:07:38 2022

@author: <Thibaut Coustillet>
__lab__: SysTox Team (U1124 | T3S) FRANCE
__version__: 1.0.0
"""

####################

import _1EM_preprocessing_module as ppm
import os
import sys
import time
from pathlib import Path

####################

# Get the parent directory of the directory containing this script.
pp_directory = sys.argv[2]

if __name__ == "__main__":
    
    data_scoring = sys.argv[1]
    t0 = time.time()
    ppm.preprocess(data_scoring)
    print(f"\u2023 Preprocess event-event file: done in {(time.time() - t0):.3f} seconds.")
    