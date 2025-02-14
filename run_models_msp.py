#######################
#  run_optimizers.py  #
#######################

from utilities_msp import *
from functions_msp import *
from time import perf_counter

# Initial time in ns:
time_i = int(perf_counter())

# Get command line args:
args = parse_command_line_generic()


amode = args['amode'] # EvolutionaryAlgorithm/SimulatedAnnealing/TabuSearch/StochasticHillClimbing/SmartRunner/ExtremalOptimization

# Process generic 'non-array coordinates':

# Population size:
Npop = args['Npop']
assert type(Npop) is int and Npop > 1, "Error: provide an integer >1 for Npop!"

# Total number of moves:
ltot = args['ltot']
assert type(ltot) is int and ltot > 0, "Error: provide a positive integer for ltot!"

# Total number of runs per hyperparameter set:
if 'nruns' in args:
	nruns = args['nruns']
	assert type(nruns) is int and nruns > 0, "Error: provide a positive integer for nruns!"
else:
	nruns = 20

# Model type and moveset modes:
mode = args['mode'] # full/basic/reduced

# Output file:
fout = args['out']


if amode == "WrightFisher" or amode == "WF":

	# Process 'array coordinates':
	mu_rate = postprocess_val('mu_rate',float,args)
	sigma = postprocess_val('sigma',float,args)
	A_PD = postprocess_val('A_PD',float,args)
	r_PD = postprocess_val('r_PD',float,args)

	prm1 = mu_rate
	prm2 = sigma
	prm3 = A_PD
	prm4 = r_PD
	prm_labels = ["mu_rate","sigma","A","r"]

	#debug
	#print("mu_rate =",mu_rate)
	#print("sigma =",sigma)
	#print("A_PD =",A_PD)
	#print("r_PD =",r_PD)

elif amode == "ContinuousTime" or amode == "CT":

	# Process 'array coordinates':
	mu_rate = postprocess_val('mu_rate',float,args)
	sigma = postprocess_val('sigma',float,args)
	r_PD = postprocess_val('r_PD',float,args)
	mu_delta = postprocess_val('mu_delta',int,args)

	prm1 = mu_rate
	prm2 = sigma
	prm3 = r_PD
	prm4 = mu_delta
	prm_labels = ["mu_rate","sigma","r","mu_delta"]

	# Process auxiliary 'non-array coordinates':
	# Time step:
	if 'delta_t' in args:
		delta_t = args['delta_t']
		assert type(delta_t) is float and delta_t > 0, "Error: provide a positive float for delta_t!"
	else:
		delta_t = 0.005

	# Initial number of occupied genotypes:
	if 'Nu' in args:
		Nu = args['Nu']
		assert type(Nu) is int and (Nu > 0 and Nu <= Npop), "Error: provide a positive integer <= Npop for Nu!"
	else:
		Nu = Npop

else:
	print("Unknown amode = {}!".format(amode))
	sys.exit("\nTerminating ..")


l1 = len(prm1)
l2 = len(prm2)
l3 = len(prm3)
l4 = len(prm4)


# Set up auxiliary parameters related to data input/output:

# Save data option:
if 'sdata' in args:
	sdata = args['sdata']
	assert type(sdata) is int and (sdata == 0 or sdata == 1), "Error: provide a {0,1} integer for sdata!"
else:
	sdata = 0 # off by default

# Probabilities to cooperate:
rv_range = [0.01, 0.99] # determines the range of the values in the pc matrix

# Set up mode-dependent parameters:
if mode == "full":
	if 'pc_in' in args:
		####
		pc_in = args['pc_in']
		pc_set = np.loadtxt(pc_in, dtype=np.float64)
		assert pc_set.shape[0] == Npop and pc_set.shape[1] == Npop, "Error: pc matrix dim mismatch!"
	else:
		pc_set = np.array(np.random.uniform(size=(Npop,Npop),low=rv_range[0],high=rv_range[1]))

	if 'pc_arr_in' in args:
		print(f"WARNING: -pc_arr_in option is ignored in the {mode} mode!")

	if 'a_in' in args:
		print(f"WARNING: -a_in option is ignored in the {mode} mode!")

	if 'b_in' in args:
		print(f"WARNING: -b_in option is ignored in the {mode} mode!")

	if 'pc_arr_out' in args or 'a_out' in args or 'b_out' in args:
		print(f"WARNING: some output options are ignored in the {mode} mode!")

	params = (pc_set,)

elif mode == "basic":

	if 'pc_arr_in' in args:
		####
		pc_arr_in = args['pc_arr_in']
		pc_arr_set = np.loadtxt(pc_arr_in, dtype=np.float64)
		assert pc_arr_set.shape[0] == Npop, "Error: pc_arr vector dim mismatch!"
	else:
		pc_arr_set = np.array(np.random.uniform(size=Npop,low=rv_range[0],high=rv_range[1])) # p_c = 1 - p_d

	pc_set = np.repeat(pc_arr_set,Npop).reshape(Npop,Npop)

	if 'pc_in' in args:
		print(f"WARNING: -pc_in option is ignored in the {mode} mode!")

	if 'a_in' in args:
		print(f"WARNING: -a_in option is ignored in the {mode} mode!")

	if 'b_in' in args:
		print(f"WARNING: -b_in option is ignored in the {mode} mode!")

	# Save *initial* settings if requested:
	if 'pc_arr_out' in args:
		pc_arr_out = args['pc_arr_out']
		np.savetxt(pc_arr_out, pc_arr_set, delimiter=' ', fmt='%12.7f')

	if 'a_out' in args or 'b_out' in args:
		print(f"WARNING: some output options are ignored in the {mode} mode!")

	params = (pc_set, pc_arr_set,)

elif mode == "reduced":

	if 'a_in' in args:
		####
		a_in = args['a_in']
		a_set = np.loadtxt(a_in, dtype=np.float64)
		assert a_set.shape[0] == Npop, "Error: a vector dim mismatch!"
	else:
		a_set = np.array(np.random.uniform(size=Npop,low=rv_range[0],high=rv_range[1]))
	
	if 'b_in' in args:
		####
		b_in = args['b_in']
		b_set = np.loadtxt(b_in, dtype=np.float64)
		assert b_set.shape[0] == Npop, "Error: b vector dim mismatch!"
	else:
		b_set = np.array(np.random.uniform(size=Npop,low=rv_range[0],high=rv_range[1]))

	pc_set = np.outer(a_set,b_set) # outer product

	if 'pc_in' in args:
		print(f"WARNING: -pc_in option is ignored in the {mode} mode!")

	if 'pc_arr_in' in args:
		print(f"WARNING: -pc_arr_in option is ignored in the {mode} mode!")

	# Save *initial* settings if requested:
	if 'a_out' in args:
		a_out = args['a_out']
		np.savetxt(a_out, a_set, delimiter=' ', fmt='%12.7f')

	if 'b_out' in args:
		b_out = args['b_out']
		np.savetxt(b_out, b_set, delimiter=' ', fmt='%12.7f')

	if 'pc_arr_out' in args:
		print(f"WARNING: some output options are ignored in the {mode} mode!")

	params = (pc_set, a_set, b_set,)

else:
	print(f"Unknown mode = {mode}!")
	####
	sys.exit("\nTerminating ..")

# Save *initial* settings if requested:
if 'pc_out' in args:
	pc_out = args['pc_out']
	np.savetxt(pc_out, pc_set, delimiter=' ', fmt='%12.7f')

# Initialize population counts:
if amode == "ContinuousTime" or amode == "CT":

	if 'counts_init' in args:

		if 'Nu' in args:
			print(f"WARNING: -Nu option is ignored since initial counts are provided!")

		####
		counts_init = args['counts_init']
		counts_0 = np.loadtxt(counts_init, dtype=np.float64)
		assert counts_0.shape[0] == Npop, "Error: counts vector dim mismatch!"

	else:

	    n_arr_i = np.array(np.random.uniform(size=Nu,low=0.0,high=1.0))
	    n_arr_i = Npop*(n_arr_i/np.sum(n_arr_i)) # 'normalization' to N
	    n_arr_i = np.array(iteround.saferound(n_arr_i,0)) # snap counts to integers; preserves Npop
	    ####
	    counts_0 = np.zeros(Npop)
	    counts_0[:Nu] = n_arr_i

	# Save *initial* counts if requested:
	if 'counts_out' in args:
		counts_out = args['counts_out']
		np.savetxt(counts_out, counts_0, delimiter=' ', fmt='%4d')


#debug
#print("numpy_version =",np.__version__)

############################################################

prm1_glob = []
prm2_glob = []
prm3_glob = []
prm4_glob = []

Ffirst_glob = []
Flast_glob = []

run_cnt = 1
for m in range(l1):

	print("{}. {} = {}".format(m+1,prm_labels[0],prm1[m]))

	for k in range(l2):
	
		#print("  {}. {} = {:.2f}".format(k+1,prm_labels[1],prm2[k]))
		print("  {}. {} = {}".format(k+1,prm_labels[1],prm2[k]))

		for i in range(l3):
				    
			#print("     {}. {} = {:.2f}".format(i+1,prm_labels[2],prm3[i]))
			print("     {}. {} = {}".format(i+1,prm_labels[2],prm3[i]))

			for p in range(l4):
				    
				#print("     {}. {} = {:.2f}".format(i+1,prm_labels[2],prm3[i]))
				print("         {}. {} = {}".format(p+1,prm_labels[3],prm4[p]))

				print("           Starting {} {} run(s) ..".format(nruns,amode))

				step = 10

				for j in range(nruns):

					if j > 0 and j%step==0:
						print("        run {} ..".format(j))

					### Main optimizer function ###
					if amode == "WrightFisher" or amode == "WF":

						### Main function of the Wright-Fisher simulation ###
						[Fave_traj,Fstd_traj,Flower_traj,Fupper_traj,params_out] = \
							WrightFisher(prm1[m],prm2[k],prm3[i],prm4[p],ltot,Npop,mode,params)

					elif amode == "ContinuousTime" or amode == "CT":
						
						### Main function of the hybrid CT simulation ###
						[Fave_traj,Flower_traj,Fupper_traj,dF_recon_traj,dF_num_traj,Fs_var_traj,Fa_var_traj,Fsa_cov_traj,Nspecies_traj,params_out,n_arr] = \
							Hybrid_CT(prm1[m],prm2[k],prm3[i],prm4[p],ltot,Npop,mode,delta_t,params,counts_0)
					
					else:
						print("Unknown amode = {}!".format(amode))
						sys.exit("\nTerminating ..")

					print("******************")

					prm1_glob.append(prm1[m])
					prm2_glob.append(prm2[k])
					prm3_glob.append(prm3[i])
					prm4_glob.append(prm4[p])

					Ffirst_glob.append(float(Fave_traj[0]))
					Flast_glob.append(float(Fave_traj[-1]))

					# Output fitness trajectories:
					if sdata == 1:

						# Write the final state of the system:
						pc_cur = fout + '.pc.' + str(run_cnt)
						np.savetxt(pc_cur, params_out[0], delimiter=' ', fmt='%12.7f')

						if amode == "ContinuousTime" or amode == "CT":
							n_arr_cur = fout + '.n.' + str(run_cnt)
							np.savetxt(n_arr_cur, n_arr, delimiter=' ', fmt='%4d')

						if mode == "basic":
							pc_arr_cur = fout + '.pc_arr.' + str(run_cnt)
							np.savetxt(pc_arr_cur, params_out[1], delimiter=' ', fmt='%12.7f')

						if mode == "reduced":
							a_cur = fout + '.a.' + str(run_cnt)
							np.savetxt(a_cur, params_out[1], delimiter=' ', fmt='%12.7f')
							####
							b_cur = fout + '.b.' + str(run_cnt)
							np.savetxt(b_cur, params_out[2], delimiter=' ', fmt='%12.7f')

						# Write the trajectory data:
						ftraj_cur = fout + '.f.' + str(run_cnt)
						####
						comment_str = "# run = {}\n".format(run_cnt)
						#ltot_arr = list(range(1,len(Ftraj)+1)) # algorithms may terminate early sometimes ..
						#ltot_arr = list(range(len(Fave_traj))) # algorithms may terminate early sometimes ..
						ltot_arr = np.arange(len(Fave_traj)) + 1

						if amode == "WrightFisher" or amode == "WF":
							header_arr = ['ltot','Fave','Fstd','Flow','Fupp']
							data_arr = [ltot_arr,Fave_traj,Fstd_traj,Flower_traj,Fupper_traj]

						if amode == "ContinuousTime" or amode == "CT":
							header_arr = ['ltot','Fave','Flow','Fupp','dF_recon','dF_num','Fs_var','Fa_var','Fsa_cov','Nspecies']
							data_arr = [ltot_arr,Fave_traj,Flower_traj,Fupper_traj,dF_recon_traj,dF_num_traj, \
										Fs_var_traj,Fa_var_traj,Fsa_cov_traj,Nspecies_traj]

						write_data(ftraj_cur,comment_str,header_arr,data_arr)

					run_cnt += 1


	print("\nElapsed time =",convert(int(perf_counter())-time_i))
	print("\n")


print("=== Finished {} {} run(s) .. ===".format(run_cnt-1,amode))


#############
if amode == "WrightFisher" or amode == "WF":
	comment_str = "# amode = {}, mode = {}\n# {}x{}x{}x{}x{} array\n# nruns = {}, ltot = {}, Npop = {}, sdata = {}\n". \
				format(amode,mode,l1,l2,l3,l4,nruns,nruns,ltot,Npop,sdata)

	header_arr = ['mu',' sigma','  A','  r',' Finit','  Ffin']
elif amode == "ContinuousTime" or amode == "CT":
	comment_str = "# amode = {}, mode = {}\n# {}x{}x{}x{}x{} array\n# nruns = {}, ltot = {}, Npop = {}, dt = {}, sdata = {}\n". \
				format(amode,mode,l1,l2,l3,l4,nruns,nruns,ltot,Npop,delta_t,sdata)

	header_arr = ['mu',' sigma','  r',' mu_delta',' Finit','  Ffin']
else:
	print("Unknown amode = {}!".format(amode))
	sys.exit("\nTerminating ..")

data_arr = [prm1_glob,prm2_glob,prm3_glob,prm4_glob,Ffirst_glob,Flast_glob]

write_data(fout,comment_str,header_arr,data_arr)


# Final time in ns:
time_f = int(perf_counter())
print("\nTotal elapsed time =",convert(time_f-time_i))
print("\n")

