import numpy as np
import math
import matplotlib.pyplot as plt
import datetime
import random
import json
plt.rcParams.update({'font.size': 16})	# Set default font size for plots

##########################################################################################
################################### Load Parameters	######################################
##########################################################################################

# Constants
kappa = 1.380649*10**(-23)
planck = 6.62607015*10**(-34)
R = 8.31446261815324  # Gas constant

# Load parameters from input JSON file
#def load_parameters(filename='input_parameters.json'):
def load_parameters(filename='Figure_try.json'):
	with open(filename, 'r') as f:
		params = json.load(f)
	return params
params = load_parameters()

# Load parameters
num = params['num_species']
connectivity = params['connectivity']
dominance_thermo = params['dominance_thermo']
dominance_kinetic = params['dominance_kinetic']
IC = params['IC']
feeding_0,feeding_1 = params['feeding_0'],params['feeding_1']
T = params['T']
time_end = params['time_end']

plotting = params["plotting"]
threshold_detection = params['threshold_detection']
target_time = params["target_time"]

IC_all_monomers = params['IC_all_monomers']
enthalpy = params['enthalpy']
entropy = params['entropy']
activation_energy = params['activation_energy']

precision_dt = params['precision_dt']
IC_eps = params['IC_eps']
water = params['water']
feeding_quantity = params['feeding_quantity']
randomnes_H = params['randomnes_H']
randomnes_Ea = params['randomnes_Ea']

# Simulation loop range
loop_start = params["loop_start"]
loop_end = params["loop_end"]
loop_increment = params["loop_increment"]

# Initialize result containers
number_species,temperatures_list,times_list,quantity_list = [],[],[],[]

##########################################################################################
################################# Main simulation loop	##################################
##########################################################################################

for num_loop in range(loop_start,loop_end,loop_increment):
	np.random.seed(10)	# Set random seed for reproducibility
	
	# SET LOOP CONDITION
	#T=num_loop*1.+params['T']
	#time_end = 1.*np.exp(0.1*num_loop)
	#target_time = time_end
	#IC = params['IC']+0.01*num_loop
	#precision_dt = params['precision_dt']/np.sqrt(np.abs(num_loop))

	# Initial concentrations for monomers and polymers
	M,M0 = np.full(shape=num, fill_value=IC_all_monomers),np.full(shape=num, fill_value=IC_all_monomers)
	P2,P20 = np.full(shape=(num,num),fill_value=IC_eps),np.zeros((num,num))
	P3,P30 = np.full(shape=(num,num,num),fill_value=IC_eps),np.zeros((num,num,num))
	P4,P40 = np.full(shape=(num,num,num,num),fill_value=IC_eps),np.zeros((num,num,num,num))
	P5,P50 = np.full(shape=(num,num,num,num,num),fill_value=IC_eps),np.zeros((num,num,num,num,num))
	W,W0 = water,water
	if enthalpy>dominance_thermo:
		M0[0],M[0]=IC,IC
		print ('Thermodynamic')
	if activation_energy>dominance_kinetic:
		M0[1],M[1]=IC*1,IC
		print ('Kinetic')

	# Define thermodynamic parameters
	HP2,SP2,GP2 =  np.zeros((num,num)),np.zeros((num,num)),np.zeros((num,num))
	for i in range(num):
		for j in range(num):
			# Add thermodynamic dominance and noise
			if i==0 or j==0:
				HP2[i][j] = dominance_thermo*4.184*1000.+np.random.uniform(-randomnes_H,randomnes_H,1)*4.184*1000.
			else:
				HP2[i][j] = enthalpy*4.184*1000.+np.random.uniform(-randomnes_H,randomnes_H,1)*4.184*1000.
			SP2[i][j] = entropy*4.184*1000.
			GP2[i][j] = HP2[i][j]-T*SP2[i][j]

	# Calculate ΔG and activation energies for M <-> P2 reactions
	H2,S2,G2,G2_act,Ea2_mp,Ea2_pm,k2_mp,k2_pm = np.zeros((num,num)),np.zeros((num,num)),np.zeros((num,num)),np.zeros((num,num)),np.zeros((num,num)),np.zeros((num,num)),np.zeros((num,num)),np.zeros((num,num))
	for i in range(num):
		for j in range(num):
			H2[i][j] = HP2[i][j]
			S2[i][j] = SP2[i][j]
			G2[i][j] = GP2[i][j]
			# Add kinetic dominance and noise
			if i==1 or j==1:
				G2_act[i][j] = dominance_kinetic*4.184*1000.+ np.random.uniform(-randomnes_Ea,randomnes_Ea,1)*4.184*1000.
			else:
				G2_act[i][j] = activation_energy*4.184*1000.+ np.random.uniform(-randomnes_Ea,randomnes_Ea,1)*4.184*1000.
			# Compute forward and reverse activation energies
			if G2[i][j]<0.:
				Ea2_mp[i][j] = G2_act[i][j]
				Ea2_pm[i][j] = np.abs(G2[i][j]) + G2_act[i][j]
			else:
				Ea2_mp[i][j] = G2[i][j] + G2_act[i][j]
				Ea2_pm[i][j] = G2_act[i][j]

			# Calculate reaction rate constants based on connectivity
			k2_mp[i][j] = (kappa*T/planck)*np.exp(-(Ea2_mp[i][j])/(R*T))
			k2_pm[i][j] = (kappa*T/planck)*np.exp(-(Ea2_pm[i][j])/(R*T))
			# Disconnect specific reactions if isolated
			if connectivity==0:
				if activation_energy>dominance_kinetic and ((i==1 and j!=1) or (j==1 and i!=1)):
					k2_mp[i][j] = 0.
					k2_pm[i][j] = 0.
				elif enthalpy>dominance_thermo and ((i==0 and j!=0) or (j==0 and i!=0)):
					k2_mp[i][j] = 0.
					k2_pm[i][j] = 0.
			# Disconnect specific reactions if partial
			if connectivity==1:
				if activation_energy>dominance_kinetic and ((i==1 and (j==9 or j==10 or j==11)) or (j==1 and (i==9 or i==10 or i==11))):
					k2_mp[i][j] = 0.
					k2_pm[i][j] = 0.
				elif enthalpy>dominance_thermo and ((i==0 and (j==9 or j==10 or j==11)) or (j==0 and (i==9 or i==10 or i==11))):
					k2_mp[i][j] = 0.
					k2_pm[i][j] = 0.

	# Set time step dt based on fastest rate
	max_k = max(k2_mp.max(),k2_pm.max())
	dt = precision_dt/max_k
	time = int(time_end/dt)
	if time<100:
		dt = time_end/10000.
		time = int(time_end/dt)
	if time > 50000000:
		time = 1000

	# Print key simulation parameters
	print ('K*dt',dt*max_k)
	print ('dt(s)',dt)
	print ('Time steps',"{:e}".format(time))
	print ('Reaction Time',"{:e}".format(time_end))
	print ('Temperature:',T-272)
	print ('Initial Concentration:',IC)

	# Initialize plots
	M_plot,P2_plot,W_plot,time_plot = [],[],[],[]
	for i in range(num):
		M_plot.append([])
		M_plot[i].append(M[i])
		P2_plot.append([])
		for j in range(num):
			P2_plot[i].append([])
			P2_plot[i][j].append(P2[i][j])
	W_plot.append(W)
	time_plot.append(0.000001*dt)

	# Run simulation
	next_print = 0 
	for t in range(time):
		progress = 100 * t // time
		if progress >= next_print:
			print(f"Simulation progress: {progress}%")# Print progress
			next_print += 10 # 10 to print every 10%

		W = water  # Reset water
		for l in range(num):
			increase_M_2,decrease_M_2 = 0.,0.
			for i in range(num):
				for j in range(num):
					# Calculate M production and consumption for monomer l
					if i ==l or j == l:
						if i==l and j == l:
							increase_M_2 += +2.*k2_pm[i][j]*P20[i][j]*W0
							decrease_M_2 += -2.*k2_mp[i][j]*M0[i]*M0[j]
						else:
							increase_M_2 += +k2_pm[i][j]*P20[i][j]*W0
							decrease_M_2 += -k2_mp[i][j]*M0[i]*M0[j]
					# Update P2 concentration
					if l==0:
						P2[i][j] = dt*(k2_mp[i][j]*M0[i]*M0[j]-k2_pm[i][j]*P20[i][j]*W0) + P20[i][j]
			
			# Apply feeding if enabled
			if (feeding_0==1 or feeding_1==1):						
				if (l==0 and feeding_0==1):
					M[l] = feeding_quantity
				elif (l==1 and feeding_1==1):
					M[l] = feeding_quantity
				else:
					M[l] = dt*(increase_M_2+decrease_M_2)+M0[l]
			else:
				M[l] = dt*(increase_M_2+decrease_M_2)+M0[l]
		
		# Update previous values
		M0,P20,P30,P40,P50,W0 = M,P2,P3,P4,P5,W

##########################################################################################
################################# SAVE DATA FOR PLOT #####################################
##########################################################################################
	
	# Save Evolution Concentrations
		if t<100:	
			W_plot.append(W)  # Save water concentration early steps
			time_plot.append((t+1)*dt)	# Save current time
			for i in range(num):
				M_plot[i].append(M[i])	# Save monomer concentrations
				for j in range(num):
					P2_plot[i][j].append(P2[i][j])	# Save dimer concentrations
		else:
			if t%100==0:  # Save less frequently after 100 steps
				W_plot.append(W)
				time_plot.append((t+1)*dt)
				for i in range(num):
					M_plot[i].append(M[i])
					for j in range(num):
						P2_plot[i][j].append(P2[i][j])
				
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
	
	number_species.append(non_zero)
	temperatures_list.append(T)
	times_list.append(time_end)
	quantity_list.append(IC)

	print('Number Species:',number_species)


# Collect simulation parameters for output file
simulation_parameters = {
	"connectivity": connectivity,
	"dominance_thermo": dominance_thermo,
	"dominance_kinetic": dominance_kinetic,
	"IC": IC,
	"feeding_0": feeding_0,
	"feeding_1": feeding_1,
	"T": T,
	"time_end": time_end,
	"num_species": num,
	"threshold_detection": threshold_detection, 
	"randomnes_H": randomnes_H,
	"randomnes_Ea": randomnes_Ea,
	"loop_start": loop_start,
	"loop_end": loop_end,
	"loop_increment": loop_increment,
	"precision_dt": precision_dt,
	"enthalpy": enthalpy,
	"entropy": entropy,
	"activation_energy": activation_energy,
	"feeding_quantity": feeding_quantity,
	"IC_eps": IC_eps,
	"water": water	
}

# Write parameters and concentration data to file
with open("Concentrations_time"+str(datetime.datetime.now)+".txt", "w") as f:
	# Write simulation parameters first
	f.write("Simulation Parameters:\n")
	for key, value in simulation_parameters.items():
		f.write(f"{key}: {value}\n")
	f.write("\n")  # Blank line before the data table

	# Prepare header row
	headers = ["time"] + [f"M_{i}" for i in range(num)] + [f"P_{i}_{j}" for i in range(num) for j in range(num)]
	f.write(",".join(headers) + "\n")

	# Write data rows for each saved time point
	for k in range(len(time_plot)):
		time_val = f"{time_plot[k]:.3e}"  # Format time
		m_vals = [f"{M_plot[i][k]:.3e}" for i in range(num)]  # Format monomers
		p2_vals = [f"{P2_plot[i][j][k]:.3e}" for i in range(num) for j in range(num)]  # Format dimers
		f.write(",".join([time_val] + m_vals + p2_vals) + "\n")	

# Write Number Species above threshold 
with open("Species_above_thershold"+str(datetime.datetime.now)+".txt", "w") as f:
	f.write(",".join(['Number Specie'] + ['Temperature']+ ['Time']+ ['Initial Concentraiton']) + "\n")
	for k in range(len(number_species)):
		column1 = f"{number_species[k]:.3e}"  
		column2 = f"{temperatures_list[k]:.3e}" 
		column3 = f"{times_list[k]:.3e}"
		column4 = f"{quantity_list[k]:.3e}" 
		f.write(",".join([column1] + [column2] + [column3] + [column4]) + "\n")	

					   
##########################################################################################
################################# MAKING PLOTS ###########################################
##########################################################################################

# Plot results
if plotting == 1:
	transparent=0.4
	for i in range(num):
		# Highlight monomer 1 in orange, others in red with transparency
		if enthalpy>dominance_thermo and i==0:
			plt.plot(time_plot,M_plot[i],c='orange')
		elif activation_energy>dominance_kinetic and i==1:
			plt.plot(time_plot,M_plot[i],c='orange')
		else:
			plt.plot(time_plot,M_plot[i],c='r',alpha=transparent+0.2)
		for j in range(num):
			# Highlight dimers involving species 1 in green, others in blue
			if enthalpy>dominance_thermo and (i==0 or j==0):
				plt.plot(time_plot,P2_plot[i][j],c='g',alpha=transparent+0.2)
			elif activation_energy>dominance_kinetic and (i==1 or j==1):
				plt.plot(time_plot,P2_plot[i][j],c='g',alpha=transparent+0.2)
			else:
				plt.plot(time_plot,P2_plot[i][j],c='b',alpha=transparent)
	
	plt.xscale('log',base=10)
	plt.yscale('log',base=10)
	plt.ylim(1*0.000001, 1.00)
	#plt.xlim(1*100, 120000)
	plt.axhline(y=threshold_detection, color='gray', linestyle='--',linewidth=1.)
	plt.axvline(x=10000, color='gray', linestyle='--',linewidth=1.)
	plt.savefig(str(datetime.datetime.now)+'.svg', dpi=300)
	plt.show()
