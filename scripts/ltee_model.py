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

#
@dataclass
class Bacteria:
    '''A class that stores bacteria values'''
    id: int
    od: np.ndarray
    ps: list
    mother: int
    alive: bool
        
def get_ivp(t, exp_time, res, bac):
    ''' Create a 1D array of time, resource and densities to use them
        to solve a ODE System'''
    time = [exp_time[t], exp_time[t+1]]
    get_res = res[[-1], [-1]] 
    get_bac = np.array([])
    get_bac = [b.od[1,-1] for b in bac] 
    y0 = np.append(get_res,get_bac)
    
    return time, y0#, get_ps

def batch(t, x, bac):
    dS = x[0]
    N = len(bac)
    dx = np.empty(N+1)
    
    for i, p in enumerate(bac):
        m, a = p.ps
        #xi = x[i+1]
        #
        fSi = (m * x[0]) / (a + x[0])
        dS = dS + (fSi * x[i+1])
        #dx[0] -= dx[0] + (fSi * x[i+1]) # this should be x+1?
        dx[i+1] = x[i+1] * fSi    
    dx[0] = - dS
    return dx

def simulate(t, exp_time, res, bac):
    this_time, this_y0 = get_ivp(t, exp_time, res, bac)
    t_eval = np.linspace(this_time[0],this_time[1], 2)
    sol = solve_ivp(batch, this_time, this_y0, args=[bac], dense_output=True, t_eval=t_eval)
    return sol

def mutate(b,len_b):
    rand_p = np.random.default_rng().uniform(0,1)
    p_mut = 5e-3 # was 5e-2 # 5e-3
    this_od = b.od[1,-1]
    ismut = (rand_p < (p_mut * this_od))
    if ismut:
        sigma = 0.02
        # Parameters for new strain 
        m_mut = np.random.default_rng().normal(b.ps[0], sigma, 1)
        a_mut = np.random.default_rng().normal(b.ps[1], sigma, 1)
        ps_mut = np.append(m_mut,a_mut)
        # density
        od_len = len(b.od[0])
        od_fill = np.zeros(od_len)
        new_od = np.vstack((b.od[0], od_fill))
        new_od[1,-1] = 0.033 # Density of new strains
        # Create new strain here
        b_mut = Bacteria(len_b+1, new_od, ps_mut, b.id, True)
    else:
        b_mut = None
    return b_mut

def alive_strains(bac):
    num_types = 0
    for i, b in enumerate(bac):
        if b.alive:
            num_types+=1
    return num_types
        
def update_res_od(sol_t, sol_y, bac, res):
    '''Update resource and optic densities from the ODE solution'''
    sol_res = np.stack((sol_t[1:],sol_y[0][1:]))
    new_res = np.append(res, sol_res, axis=1) # This works
    sol_od = np.delete(sol_y[1:],0,1)
    # Array with time in axis 0 and without the first elements 
    #     for each array of densities
    new_od = np.vstack((sol_t[1:], sol_od))
    for i, x in enumerate(bac):
        x.od = np.append(x.od, new_od[[0,i+1],:], axis=1)        
    return new_res, bac

def transfer(res,bac, d):
    for b in bac:
        b.od[1,-1] *= d
    res[1,-1] = 1.0
    return res, bac

def life_and_death(bac, ti):
    death = 0.033
    num_types = alive_strains(bac)
    for b in bac:
        if b.alive:
            if b.od[1,-1] < death: #works
            #if b.od[1,-1] <= death:
                print('Dead id: ', b.id, 'at ', ti)
                b.alive = False
                b.od[1,-1] = 0.0
            if b.alive and num_types<20:
                bac_mut = mutate(b,len(bac))
                if bac_mut:
                    bac.append(bac_mut)
                    print('Born: ', bac_mut.id, 'from: ', b.id, 'at ', ti)
    return bac
                
def cycle_df(bac, res, t0, t1, cycle_counter, data_dir):
    '''Create dataframe and store it as cvs'''
    len_bac = len(bac)
    len_od = len(bac[0].od[0][t0:t1])
    arr_ps = np.empty([len_bac, 4])
    result_array = np.empty([len_bac+1, len_od])
    result_array[0] = bac[0].od[0][t0:t1] # Time array
    col_name_bac = ['id', 'm', 'a', 'mother']
    bac_names = [f'x_{i}' for i, b in enumerate(bac)]
    bac_names.insert(0,'t')
    res_cycle = res[:,t0:t1]
    #test_result_array = [b.od[1] for i, b in enumerate(bacteria)]
    for i, b in enumerate(bac): # Use comprehension list instead (?)
        result_array[i+1] = b.od[1][t0:t1]
        arr_ps[i] = [b.id, b.ps[0], b.ps[1],b.mother]
    df_ods = pd.DataFrame(result_array.T)
    df_res = pd.DataFrame(res_cycle.T)
    df_ps = pd.DataFrame(arr_ps)
    # Naming
    df_ods.columns = bac_names
    df_res.columns = ['t', 'S']
    df_ps.columns = col_name_bac
    # Store
    df_ods.to_csv(data_dir + '/od_' + cycle_counter +'.csv', index=False)
    df_res.to_csv(data_dir + '/res_' + cycle_counter +'.csv', index=False)
    df_ps.to_csv(data_dir + '/ps_' + cycle_counter +'.csv', index=False)

def store_data(t, cycle_len, bac, res, data_dir):
    if t%10:# This would not work with odd sized cycles
        t+=2
    cycle_counter = str(int(t/cycle_len))
    cycle_counter =  cycle_counter.zfill(4)
    cycle_df(bac, res, t-cycle_len, t, cycle_counter, data_dir)