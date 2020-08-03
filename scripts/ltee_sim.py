#!interpreter [optional-arg]
# -*- coding: utf-8 -*-

"""
Description
"""
#
import os
import sys
from datetime import datetime
from dataclasses import dataclass
#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
#
from params import *
from ltee_model import *
#
#############################
bacteria = []
bacteria.append(x)

# Create directory
today = datetime.now()
data_dir = ("./data_" + today.strftime('%Y-%m-%d'))
os.mkdir(data_dir)

# LTEE

for t, ti in enumerate(exp_time[:-1]):
    
    if not ti%day and ti:
        # Save data to file
        store_data(t, cycle_len,bacteria, res, data_dir)
        # Transfer
        print ("Transfer at: ", ti)
        res, bacteria = transfer(res,bacteria, d)
    sol = simulate(t, exp_time, res, bacteria)
    res, bacteria = update_res_od(sol.t, sol.y, bacteria, res) 
    bacteria = life_and_death(bacteria, ti) # prev just calling the fun
    if ti == exp_time[-2]: # Last element of the iteration
        store_data(t, cycle_len,bacteria, res, data_dir)
        print('*********')
    extintion = alive_strains(bacteria)
    if not extintion:
        store_data(t, cycle_len, bacteria, res, data_dir)
        print('There are not more strains alive')
        break