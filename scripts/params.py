import numpy as np
from ltee_model import Bacteria

res = np.array([[0.0],[1.0]])
od = np.array([[0.0], [0.033]])
m = 1.9
a = 0.02
ps = [m,a]
ps2 = [(1.3), (0.4)]
s = np.array([0, 1]) # (?)
x = Bacteria(0, od, ps, 0, True) # Ancestral strain
death = 0.033

d = 0.1
store = True
day = 10
cycles = 500
exp_time = np.arange(0, ((cycles)*day), 1e-3) #1e-2
cycle_len = int(len(exp_time)/cycles)


# Add exp parameters (d, sigma, etc)
