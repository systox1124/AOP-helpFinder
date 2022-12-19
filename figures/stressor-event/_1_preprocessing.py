#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:14:10 2022

@author: <Thibaut Coustillet>
__lab__: SysTox Team (U1124 | T3S) FRANCE
__version__: 1.0.0
"""

####################

import _1_preprocessing_module as ppm
import os
import sys
import time
from pathlib import Path

####################

# Get the parent directory of the directory containing this script.
pp_directory = sys.argv[2]
#pp_directory = Path(__file__).resolve().parents[1].as_posix()


if __name__ == "__main__":
    
    # Creation of the directory to store the results.
    resv2 = os.path.join(pp_directory, "")
    figv2 = os.path.join(pp_directory, "figures/")
    #os.mkdir(resv2)
    os.mkdir(figv2)
    
    output_aophf = sys.argv[1]
    t0 = time.time()
    ppm.add_stressors(output_aophf)
    ppm.make_file4heatmap(output_aophf)
    print(f"\u2023 Create table with all stressors and file for heatmap: done in {(time.time() - t0):.3f} seconds.")
    