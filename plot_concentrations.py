import numpy as np
import math
import matplotlib.pyplot as plt
import datetime
import random
import json
plt.rcParams.update({'font.size': 16})
##########################################################################################

time_plot = []
M_plot = [[] for _ in range(12)]		  # 12 lists for M_plot[i]
P2_plot = [[[] for _ in range(12)] for _ in range(12)]	# 12x12 lists for P2_plot[i][j]

with open("jjjj.txt", "r") as f:
	lines = f.readlines()

# Skip lines until the header (find the line that starts with "time,")
data_start_index = 0
for idx, line in enumerate(lines):
	if line.startswith("time,"):
		data_start_index = idx + 1
		break

# Parse the table lines
for line in lines[data_start_index:]:
	values = line.strip().split(",")
	
	# Read time
	time_plot.append(float(values[0]))
	
	# Read M_plot
	for i in range(12):
		M_plot[i].append(float(values[1 + i]))
	
	# Read P2_plot
	idx = 1 + 12
	for i in range(12):
		for j in range(12):
			P2_plot[i][j].append(float(values[idx]))
			idx += 1

##########################################################################################

# Parameters
threshold_detection = 0.0001
target_time = 10000.0  # <-- Change this to the time you want to check

# Find index closest to target_time
time_array = np.array(time_plot)
time_index = (np.abs(time_array - target_time)).argmin()

# Count species above threshold at the given time
num = len(P2_plot)
non_zero = 0

for i in range(num):
    if M_plot[i][time_index] > threshold_detection:  # Check monomer at time_index
        non_zero += 1
    for j in range(num):
        if P2_plot[i][j][time_index] > threshold_detection:  # Check dimer at time_index
            non_zero += 1

print('Species at time', time_plot[time_index], ':', non_zero)

##########################################################################################

transparent=0.4
for i in range(num):
#	if feeding_0==1 and i==0:
#		transparent=1.
#	elif feeding_1==1 and i==1:
#		transparent=1.
#	else:
#		transparent=.2
	if i==1:
		plt.plot(time_plot,M_plot[i],c='orange')
	else:
		plt.plot(time_plot,M_plot[i],c='r',alpha=transparent+0.2)
	for j in range(num):
		#plt.plot(time_plot,P2_plot[i][j],c='b',alpha=1.-float(i+1)/num)
		#plt.plot(time_plot,P2_plot[i][j],c='b',alpha=transparent)
		if i==1 or j==1:
			plt.plot(time_plot,P2_plot[i][j],c='g',alpha=transparent+0.2)
		else:
			plt.plot(time_plot,P2_plot[i][j],c='b',alpha=transparent)

plt.xscale('log',base=10)
plt.yscale('log',base=10)
#plt.ylim(1*0.000001, 10.00)
#plt.xlim(1*100, 120000)#0.85
#plt.xlim(0.02, 20000)#0.85
#plt.title('T='+str(T)+' time='+str(time_end)+' Q='+str(IC))
plt.axhline(y=threshold_detection, color='gray', linestyle='--',linewidth=1.)
plt.axvline(x=10000, color='gray', linestyle='--',linewidth=1.)
#plt.axvline(x=10000, color='gray', linestyle='--',linewidth=1.)
plt.savefig(str(datetime.datetime.now)+'.svg', dpi=300)
plt.show()

