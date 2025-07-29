import numpy as np
import math
import matplotlib.pyplot as plt
import datetime
import random
import json
plt.rcParams.update({'font.size': 16})
##########################################################################################

import matplotlib.pyplot as plt

number_vs = 2 #1 for temperature, 2 for reaction time, 3 for initial concentration

# Replace this with your actual file path
# filename_1 = "4a.txt"
# filename_2 = "4b.txt"
# filename_3 = "4c.txt"
# number_vs = 3

# filename_1 = "8ba.txt"
# filename_2 = "8bb.txt"
# filename_3 = "8bc.txt"
# number_vs = 3

filename_1 = "9a_1.txt"
filename_2 = "9a_2.txt"
filename_3 = "9a_3.txt"
number_vs = 1
# 
# filename_1 = "9b_1.txt"
# filename_2 = "9b_2.txt"
# filename_3 = "9b_3.txt"
# number_vs = 2

# Initialize empty lists for the four columns
col1a, col1b, col2a, col2b, col3a, col3b = [], [], [], [], [], []

# Read the file and skip the first line
with open(filename_1, 'r') as f:
    next(f)  # Skip the header line
    for line in f:
        values = line.strip().split(',')  # Split by colon
        if len(values) >= 4:
            col1a.append(float(values[number_vs]))
            col1b.append(float(values[0]))
with open(filename_2, 'r') as f:
    next(f)  # Skip the header line
    for line in f:
        values = line.strip().split(',')  # Split by colon
        if len(values) >= 4:
            col2a.append(float(values[number_vs]))
            col2b.append(float(values[0]))
with open(filename_3, 'r') as f:
    next(f)  # Skip the header line
    for line in f:
        values = line.strip().split(',')  # Split by colon
        if len(values) >= 4:
            col3a.append(float(values[number_vs]))
            col3b.append(float(values[0]))

v_min = np.min(col1b)
v_max = np.max(col1b)
if v_max != v_min:
    col1b = (col1b - v_min) / (v_max - v_min)
else:
    col1b = np.zeros_like(col1b)

v_min = np.min(col2b)
v_max = np.max(col2b)
if v_max != v_min:
    col2b = (col2b - v_min) / (v_max - v_min)
else:
    col2b = np.zeros_like(col2b)
    
v_min = np.min(col3b)
v_max = np.max(col3b)
if v_max != v_min:
    col3b = (col3b - v_min) / (v_max - v_min)
else:
    col3b = np.zeros_like(col3b)

# Plot the vectors (against their index)
plt.figure(figsize=(10, 6))
# plt.scatter(col1a,col1b, label=r'$\sigma _i = 0.1$',c='r',s=70)
# plt.scatter(col2a,col2b, label=r'$\sigma _i = 0.5$',c='g',s=70)
# plt.scatter(col3a,col3b, label=r'$\sigma _i = 0.8$',c='b',s=70)
plt.scatter(col1a,col1b, label=r'$E_a = 26$',c='r',s=70)
plt.scatter(col2a,col2b, label=r'$E_a = 25$',c='g',s=70)
plt.scatter(col3a,col3b, label=r'$E_a = 24$',c='b',s=70)
#plt.xscale('log')
plt.xlabel("Concentration Kinetic Compressor (M)")
plt.xlabel("Temperature (°C)")
#plt.xlabel("Time (s)")
plt.ylabel("Number of Dimmers + Monomers (Normalized)")
plt.legend(loc=2)
plt.tight_layout()
plt.savefig(str(datetime.datetime.now)+'.svg', dpi=300)
plt.show()