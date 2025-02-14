from math import *
import numpy as np
import copy
import iteround
import scipy.stats as stats
from scipy.stats import bernoulli
from scipy.stats import norm
from utilities_msp import *

#import warnings
#warnings.filterwarnings("ignore", category=RuntimeWarning)

#################################
### Function definitions (CT) ###
#################################

# Computes change in agent counts due to mutations:
def mutate_counts(n_old, mu, sigma, prms, n_zero_flat, n_pos_flat, MODE):
    
    # Probabilities to cooperate:


    if MODE == "full":
        pc = prms[0]
    elif MODE == "basic":
        pc = prms[0]
        pc_arr = prms[1]
    elif MODE == "reduced":
        pc = prms[0]
        a = prms[1]
        b = prms[2]
    else:
        print(f"Unknown mode = {MODE}!")
        ####
        sys.exit("\nTerminating ..")

    Ncur = pc.shape[0]
    sigmas = np.repeat(sigma,Ncur)
    two_sigmas = np.repeat(sigma,2)
    ones = np.ones(Ncur)
    
    # Potential mutant holders:
    #zero_cnt_args = list(np.argwhere(n_old == 0).flatten())
    zero_cnt_args = list(n_zero_flat)
    
    # Currently occupied slots:
    #nonzero_cnt_args = list(np.argwhere(n_old > 0).flatten())
    nonzero_cnt_args = list(n_pos_flat)
    
    n_new = copy.deepcopy(n_old)
    
    for ind in nonzero_cnt_args:
        for i in range(int(n_old[ind])):
            u = np.random.uniform(low=0.0,high=1.0)
            if u < mu: # mutate this agent
                #debug
                #print("Mutating agent {} at position {} ..".format(i+1,ind+1))
                if MODE == "full": # mutate a row
                    pc_row_new = stats.truncnorm.rvs(-pc[ind,:]/sigma,(ones-pc[ind,:])/sigma,\
                                                    loc=pc[ind,:],scale=sigmas)               
                    #pc_row_new = np.zeros(Ncur)
                    #for j in range(Ncur):
                    #    pc_row_new[j] = stats.truncnorm.rvs(-pc[ind,j]/sigma,(1.0-pc[ind,j])/sigma,\
                    #                                        loc=pc[ind,j],scale=sigma,size=1)
                    ######
                    if int(n_new[ind]) == 1: # mutate in situ
                        pc[ind,:] = pc_row_new
                    else:
                        new_ind = zero_cnt_args[-1]
                        ####
                        # Copy the impressions of others onto the new mutant's position:
                        pc[:,new_ind] = pc[:,ind]
                        pc[new_ind,:] = pc_row_new # this will overwrite the diag. elem., needs to be fixed
                        pc[new_ind,new_ind] = \
                            stats.truncnorm.rvs(-pc[ind,ind]/sigma,(1.0-pc[ind,ind])/sigma,\
                                                loc=pc[ind,ind],scale=sigma,size=1)
                        ####
                        n_new[ind] -= 1
                        n_new[new_ind] += 1
                        zero_cnt_args.pop() # remove the newly occupied position from the n=0 list
                elif MODE == "basic":
                    pc_new = stats.truncnorm.rvs(-pc_arr[ind]/sigma,(1.0-pc_arr[ind])/sigma,\
                                                 loc=pc_arr[ind],scale=sigma,size=1)
                    if int(n_new[ind]) == 1: # mutate in situ
                        pc_arr[ind] = pc_new
                        pc[ind,:] = np.repeat(pc_new,Ncur)
                    else:
                        new_ind = zero_cnt_args[-1]
                        ####
                        pc_arr[new_ind] = pc_new
                        # Copy impressions of others onto the new mutant's position:
                        pc[:,new_ind] = pc[:,ind]
                        pc[new_ind,:] = np.repeat(pc_new,Ncur)
                        ####
                        n_new[ind] -= 1
                        n_new[new_ind] += 1
                        zero_cnt_args.pop() # remove the newly occupied position from the n=0 list
                else: # MODE == "reduced"
                    lower = np.array([-a[ind]/sigma, -b[ind]/sigma])
                    upper = np.array([(1.0-a[ind])/sigma, (1.0-b[ind])/sigma])
                    loc = np.array([a[ind], b[ind]])
                    [a_new,b_new] = stats.truncnorm.rvs(lower,upper,loc=loc,scale=two_sigmas)
                    #a_new = stats.truncnorm.rvs(-a[ind]/sigma,(1.0-a[ind])/sigma,loc=a[ind],scale=sigma,size=1)
                    #b_new = stats.truncnorm.rvs(-b[ind]/sigma,(1.0-b[ind])/sigma,loc=b[ind],scale=sigma,size=1)
                    #debug
                    #print("n_old[ind] =",n_old[ind])
                    #print("a_old = {}, a_new = {}".format(a[ind],a_new))
                    #print("b_old = {}, b_new = {}".format(b[ind],b_new))
                    if int(n_new[ind]) == 1: # mutate in situ
                        a[ind] = a_new
                        b[ind] = b_new
                        pc[ind,:] = a_new * b
                        pc[:,ind] = b_new * a
                        pc[ind,ind] = a_new * b_new
                    else:
                        #debug
                        #print("zero_cnt_args =",zero_cnt_args)
                        new_ind = zero_cnt_args[-1]
                        #debug
                        #print("new_ind =",new_ind)
                        ####
                        a[new_ind] = a_new
                        b[new_ind] = b_new
                        # Copy impressions of others onto the new mutant's position:
                        pc[new_ind,:] = a_new * b
                        pc[:,new_ind] = b_new * a
                        pc[new_ind,new_ind] = a_new * b_new # fix the diag. elem.
                        ####
                        n_new[ind] -= 1
                        n_new[new_ind] += 1
                        zero_cnt_args.pop() # remove the newly occupied position from the n=0 list
    
    return n_new


# Computes approx. mean of the truncated Gaussian distribution (testing ONLY):
def getG(pc_cur,sigma):
    
    offset = 3.0
    ######
    if (pc_cur < offset * sigma) or (pc_cur > 1 - offset * sigma):
        Gval = stats.truncnorm.mean(-pc_cur/sigma, (1-pc_cur)/sigma, loc=0, scale=sigma)
    else:
        Gval = 0.0
        
    return Gval


# This function precomputes the means of the truncated Gaussian distribution:
def precomputeG(M,sigma):
    
    pc = np.linspace(0.0, 1.0, num=M)
    
    # This generates a runtime warning for some unclear reason:
    #means = stats.truncnorm.mean(-pc/sigma, (np.ones(M) - pc)/sigma, \
    #                             loc=0.0, scale=sigma)
    
    means = np.zeros(M)
    
    cnt = 0
    for x in pc:
        means[cnt] = stats.truncnorm.mean(-x/sigma, (1-x)/sigma, loc=0, scale=sigma)
        cnt += 1
    
    #debug
    #print(-pc/sigma)
    #print((np.ones(M) - pc)/sigma)
    
    return pc,means


# This function uses precomputed values of G to do a fast mapping onto the current 1D array or 2D matrix:
def get_all_G(pc_cur,precomp_means):
    
    pc_shape = pc_cur.shape
    Delta = 1.0/(precomp_means.size - 1)
    i_cur = np.round(pc_cur.ravel() / Delta).astype(int)
    
    return precomp_means[i_cur].reshape(pc_shape)

###############################
##### Aux. functions (WF) #####
###############################

# This function is used to create mutational perturbations:
def create_mut_mask(Ns, mu, sigma):
    mask = bernoulli.rvs(mu, size=Ns).astype(np.float64)
    #debug
    #print("Bernoulli mask =",mask)
    mask_pos = np.argwhere(mask == 1.0).ravel()
    ####
    #debug
    #print("mask_pos_size =",mask_pos.size)
    if mask_pos.size > 0:
        mask_delta = norm.rvs(loc=0.0, scale=sigma, size=mask_pos.size)
        #debug
        #print("Mutations =",mask_delta)
        #print("mask_delta_shape =",mask_delta.shape)
        #mask[mask_pos] = mask_delta
        np.put(mask, mask_pos, mask_delta)
    
    #debug
    #print("Final mask =",mask)
    
    return mask

# This function counts the number of identical rows in a matrix.
def count_identical_rows(matrix, digits = 6):
    """
    Counts the number of identical rows in a matrix.

    Args:
        matrix: 2D NumPy matrix.

    Returns:
        Sorted counts of identical rows.
    """
    matrix_r = np.around(matrix, decimals = digits) # rounding array elements
    row_counts = {}
    
    for row in matrix_r:
        row_tuple = tuple(row)  # Convert row to tuple
        #debug
        #print(row_tuple)
        if row_tuple in row_counts:
            row_counts[row_tuple] += 1
        else:
            row_counts[row_tuple] = 1
    
    # Sorted counts:
    cnt_srt = np.sort(list(row_counts.values()))
    #freq_srt = np.flip(cnt_srt)/matrix.shape[0]
    cnt_srt = np.flip(cnt_srt)
    
    return cnt_srt

##########################

##############################
##### Main function (WF) #####
##############################


# mu = mutation rate (prob. to mutate per generation)
# sigma = std_dev of the mutational Gaussian
# A = Fitness scale
# r = PD strength
# Ltot = total number of steps
# N = total population size
# MODE = "full"/"basic"/"reduced"
def WrightFisher(mu,sigma,A,r,Ltot,N,MODE,prms):

	Mstep = int(Ltot // 10) # intervals at which progress is reported

	#print(f"*** MODE = {MODE} ***")
	#print(f"N = {N}, Ltot = {Ltot}")
	#print(f"A = {A}, r = {r}, mu = {mu}, sigma = {sigma}")

	# Probabilities to cooperate:
	if MODE == "full":
		pc = prms[0]
	elif MODE == "basic":
		pc = prms[0]
		pc_arr = prms[1]
	elif MODE == "reduced":
		pc = prms[0]
		a = prms[1]
		b = prms[2]
	else:
		print(f"Unknown mode = {MODE}!")
	    ####
		sys.exit("\nTerminating ..")


	# Fixed genotype frequencies:
	x_arr = np.ones(N)/N
	#debug
	#print(x_arr)

	#debug
	#print("Initial pc Matrix:")
	#with np.printoptions(precision=4, suppress=False, formatter={'float': '{:0.4f}'.format}, linewidth=100):
	#    print(pc)
	#print("Initial a and b vectors:")
	#with np.printoptions(precision=4, suppress=False, formatter={'float': '{:0.4f}'.format}, linewidth=100):
	#    print(a)
	#    print(b)

	# Percentile bounds:
	LowerP = 0.25 #10
	UpperP = 0.75 #90

	Ftraj = np.zeros(Ltot) # mean fitness
	Ftraj_std = np.zeros(Ltot) # fitness std_dev
	Ftraj_lower_p = np.zeros(Ltot) # fitness lower_percentile
	Ftraj_upper_p = np.zeros(Ltot) # fitness upper_percentile

	#Fcorr = (1/(N-1)) * np.diagonal(pc) # finite-size fitness correction, DISABLED

	# Evolve the population:
	for l in range(Ltot):
	    ####
	    if Mstep > 0 and (l+1)%Mstep == 0:
	        print("{} steps ..".format(l+1))
	    
	    # Compute fitness for each organism:
	    pc_to_i = np.dot(pc.T,x_arr)
	    pc_i_to = np.dot(pc,x_arr)
	    
	    # Array of fitnesses (reduced):
	    F_arr_s = 0.5 * (pc_to_i + pc_i_to)
	    F_arr_a = (r + 0.5) * (pc_to_i - pc_i_to)
	    F_arr = A * (F_arr_s + F_arr_a) # A = overall scale factor
	    
	    # Finite-size correction:
	    #F_arr = (N/(N-1)) * F_arr - Fcorr # DISABLED
	    
	    Ftraj[l] = np.mean(F_arr)
	    
	    # Sample with fitnesses as weights (organisms are rows in pd_arr):
	    #debug
	    #print("l =",l)
	    #if l == Ltot - 1:
	    #    with np.printoptions(precision=4, suppress=False, formatter={'float': '{:0.4f}'.format}, linewidth=100):
	    #        print("Fitness:\n",F_arr)
	    #        print("pc_0:\n",pc[0,:])
	    #        print("pc_1:\n",pc[1,:])
	    #        print("pc_2:\n",pc[2,:])
	    #    print("========================")
	    
	    Fstd = np.std(F_arr)
	    Ftraj_std[l] = Fstd

	    [Ftraj_lower_p[l],Ftraj_upper_p[l]] = weighted_quantile(F_arr, [LowerP,UpperP], sample_weight=None)
	    #Ftraj_lower_p[l]= np.percentile(F_arr, LowerP)
	    #Ftraj_upper_p[l] = np.percentile(F_arr, UpperP)
	    
	    if Fstd > 1e-08: # unequal fitness values
	        Farr_exp = np.exp(F_arr)
	        Fnorm = Farr_exp/np.sum(Farr_exp)
	        #Farr_shifted = Farr - np.min(Farr)
	        #Fnorm = Farr_shifted/np.sum(Farr_shifted)
	        new_ind = np.sort(np.random.choice(N, size=N, replace=True, p=Fnorm))
	        #debug
	        #with np.printoptions(precision=4, suppress=False, formatter={'float': '{:0.4f}'.format}, linewidth=100):
	        #    print("Fnorm =",Fnorm)
	        #print("========================")
	    else: # [nearly] equal fitness values
	        new_ind = np.sort(np.random.choice(N, size=N, replace=True, p=None))
	    
	    #debug
	    #print("new_ind =",new_ind)
	    #print("========================")
	    
	    ###########################
	    # Build a new population:
	    ###########################
	    
	    ### Selection: ###
	    if MODE == "full":
	        pc_new = pc[np.ix_(new_ind,new_ind)]
	        
	        #pc_new = np.zeros((N, N))
	        #mutated = np.zeros((N, N), dtype=int) # mutation counter, to prevent multiple mutations

	        #for i in range(N):
	        #    for j in range(N):
	        #        pc_new[i,j] = pc[new_ind[i],new_ind[j]]
	        #        # Remove kin selection:
	        #        if i != j and new_ind[i] == new_ind[j]:
	        #            #pc_new[i,j] = np.random.uniform(low=0.9,high=1.0)
	        #            pc_new[i,j] = np.random.uniform()
	        #    ######
	        #    pc_new[i,i] = 0.0 # remove self-interactions
	            
	    elif MODE == "basic":
	        pc_arr_new = pc_arr[new_ind]
	        pc_new = np.repeat(pc_arr_new,N).reshape(N,N)
	    elif MODE == "reduced":
	        a_new = a[new_ind]
	        b_new = b[new_ind]
	        pc_new = np.outer(a_new,b_new) # outer product
	    
	    #debug
	    #print("pc after selection:")
	    #with np.printoptions(precision=4, suppress=False, formatter={'float': '{:0.4f}'.format}, linewidth=100):
	    #    print(pc_new)
	    #print("a/b after selection:")
	    #with np.printoptions(precision=4, suppress=False, formatter={'float': '{:0.4f}'.format}, linewidth=100):
	    #    print(a_new)
	    #    print(b_new)
	    #print("*======================*")
	    
	    ### Mutation: ###
	    if mu > 0.0:
	        # Create a Bernoulli "mask" for which entries to mutate, multiply by mutational 'deltas',
	        # and add to the old values.
	        if MODE == "full":
	            pc_new = pc_new + create_mut_mask(N*N, mu, sigma).reshape(N,N)
	            # Clip out-of-range values:
	            pc_new[pc_new < 0.0] = 0.0
	            pc_new[pc_new > 1.0] = 1.0
	            # Remove self-interactions:
	            #np.fill_diagonal(pc_new, 0.0) # so it does not affect the mean
	        elif MODE == "basic":
	            pc_arr_new = pc_arr_new + create_mut_mask(N, mu, sigma)
	            # Clip out-of-range values:
	            pc_arr_new[pc_arr_new < 0.0] = 0.0
	            pc_arr_new[pc_arr_new > 1.0] = 1.0
	            ####
	            pc_new = np.repeat(pc_arr_new,N).reshape(N,N)
	        elif MODE == "reduced":
	            a_new = a_new + create_mut_mask(N, mu, sigma)
	            b_new = b_new + create_mut_mask(N, mu, sigma)
	            # Clip out-of-range values:
	            a_new[a_new < 0.0] = 0.0
	            a_new[a_new > 1.0] = 1.0
	            ####
	            b_new[b_new < 0.0] = 0.0
	            b_new[b_new > 1.0] = 1.0
	            ####
	            pc_new = np.outer(a_new,b_new) # outer product
	        else:
	            print("Unknown MODE = {}!".format(MODE))
	            sys.exit("\nTerminating ..")
	    
	    #debug
	    #print("pc after mutation:")
	    #with np.printoptions(precision=4, suppress=False, formatter={'float': '{:0.4f}'.format}, linewidth=100):
	    #    print(pc_new)
	    #print("a/b after mutation:")
	    #with np.printoptions(precision=4, suppress=False, formatter={'float': '{:0.4f}'.format}, linewidth=100):
	    #    print(a_new)
	    #    print(b_new)
	    #print("*======================*")
	    
	    # Copy the new population into the old one:
	    pc = copy.deepcopy(pc_new)
	    if MODE == "basic":
	        pc_arr = copy.deepcopy(pc_arr_new)
	    elif MODE == "reduced":
	        a = copy.deepcopy(a_new)
	        b = copy.deepcopy(b_new)


	# Final population state:
	if MODE == "full":
		prms_out = (pc,)
	elif MODE == "basic":
		prms_out = (pc,pc_arr,)
	elif MODE == "reduced":
		prms_out = (pc,a,b)
	else:
		print(f"Unknown mode = {MODE}!")
		####
		sys.exit("\nTerminating ..")

	return [Ftraj,Ftraj_std,Ftraj_lower_p,Ftraj_upper_p,prms_out]



##############################
##### Main function (CT) #####
##############################


# mu = mutation rate (prob. to mutate per generation)
# sigma = std_dev of the mutational Gaussian
# r = PD strength
# mu_delta = determines the rate of rounding and mutations
# Ltot = total number of steps
# N = total population size
# MODE = "full"/"basic"/"reduced"
# delta_t = timestep for selection
def Hybrid_CT(mu,sigma,r,mu_delta,Ltot,N,MODE,delta_t,prms,n_arr_0):

	#Mstep = 10000 # intervals at which progress is reported
	Mstep = int(Ltot // 10) # intervals at which progress is reported

	#print(f"*** MODE = {MODE} ***")
	#print(f"N = {N}, Ltot = {Ltot} (t = {Ltot*delta_t}), r = {r}, mu = {mu}")
	#print(f"Ltot = {Ltot}: t = {Ltot*delta_t}, Nmut = {Ltot//mu_delta}")

	# Precompute G:
	M = 2000
	if mu > 0.0:
	    pc_precomp,means_precomp = precomputeG(M,sigma)

	# Create a truncated normal distribution object:
	# (lower_bound - mean) / sigma
	# (upper_bound - mean) / sigma
	#trunc_norm = stats.truncnorm(0.0, 1.0/sigma, loc=0.0, scale=sigma)

	# Probabilities to cooperate:
	if MODE == "full":
		pc = copy.deepcopy(prms[0])
	elif MODE == "basic":
		pc = copy.deepcopy(prms[0])
		pc_arr = copy.deepcopy(prms[1])
	elif MODE == "reduced":
		pc = copy.deepcopy(prms[0])
		a = copy.deepcopy(prms[1])
		b = copy.deepcopy(prms[2])
	else:
		print(f"Unknown mode = {MODE}!")
		####
		sys.exit("\nTerminating ..")

	
	# Initial genotype counts:
	n_arr = copy.deepcopy(n_arr_0)

	# Initial genotype frequencies:
	x_arr = n_arr/N

	#debug
	#print("n_arr =",n_arr)
	#print("x_arr =",x_arr)
	#print("a =",a)
	#print("b =",b)
	#print("pc:")
	#print(pc)

	# Init. "reduced" arrays:
	n_arr_pos = np.argwhere(n_arr > 0)
	n_arr_pos_flat = n_arr_pos.ravel()
	n_arr_zero_flat = np.argwhere(n_arr == 0).ravel()
	####
	x_red = x_arr[n_arr_pos].ravel()
	pc_red = pc[np.ix_(n_arr_pos_flat,n_arr_pos_flat)]

	#debug
	#print("x_red =",x_red)
	#print("pc_red:")
	#print(pc_red)
	#print("***** end_setup *****")
	#assert 0

	# Percentile bounds:
	LowerP = 0.25 #10
	UpperP = 0.75 #90

	Ftraj = np.zeros(Ltot) # average fitness
	#Ftraj_std = np.zeros(Ltot) # fitness std_dev
	Ftraj_lower_p = np.zeros(Ltot) # fitness lower_percentile
	Ftraj_upper_p = np.zeros(Ltot) # fitness upper_percentile

	dF_recon = np.zeros(Ltot) # d<F>/dt computed through var/cov

	# Array with all genotype counts:
	#n_all = np.zeros((N,Ltot+1),dtype=np.int32)
	#n_all[:,0] = np.copy(n_arr.astype(int))

	# Array with all genotype fitnesses:
	#F_all = np.zeros((N,Ltot))

	####
	Nspecies = np.zeros(Ltot,dtype=np.int64) # number of species with freqs above a threshold
	Fs_sigma_sqr = np.zeros(Ltot) # var of F_s
	Fa_sigma_sqr = np.zeros(Ltot) # var of F_a
	Fsa_cov = np.zeros(Ltot) # cov(F_s,F_a)

	#t_cur = 0.0
	######################################
	# Evolve the population:
	for l in range(Ltot):
	    
	    ####
	    if Mstep > 0 and (l+1)%Mstep == 0:
	        print("{} steps ..".format(l+1))
	    
	    # Project onto the subspace of non-zero counts:
	    #if l % mu_delta == 0: # note that multiple rounds of selection use the same notion of zero/non-zero slots
	    #    n_arr_pos = np.argwhere(n_arr > 0)
	    #    n_arr_pos_flat = n_arr_pos.ravel()
	    #    n_arr_zero_flat = np.argwhere(n_arr == 0).ravel()
	    #    ####
	    #    x_red = x_arr[n_arr_pos].ravel()
	    #    pc_red = pc[np.ix_(n_arr_pos_flat,n_arr_pos_flat)]
	        #debug
	        #print("shapes = {} {}".format(pc.shape,pc_red.shape))
	        #print("n_arr_pos_flat =",n_arr_pos_flat)
	        #print("n_arr_zero_flat =",n_arr_zero_flat)
	        #print("pc_red:")
	        #print(pc_red)
	        #print("x_red =",x_red)
	    
	    pc_to_i = np.dot(pc_red.T,x_red)
	    pc_i_to = np.dot(pc_red,x_red)
	    
	    # Array of fitnesses (reduced):
	    F_arr_s = 0.5 * (pc_to_i + pc_i_to)
	    F_arr_a = (r + 0.5) * (pc_to_i - pc_i_to)
	    F_arr = F_arr_s + F_arr_a
	    
	    #debug
	    #print("F_arr_shape =",F_arr.shape)
	    #print("F_arr =",F_arr)
	    
	    # Compute variances and covariances:
	    F_s_ave = np.dot(F_arr_s,x_red)
	    #F_a_ave = np.dot(F_arr_a,x_red) # this is =0.0 identically, from theory
	    
	    Fs_var = np.dot(np.square(F_arr_s),x_red) - F_s_ave*F_s_ave
	    #Fa_var = np.dot(np.square(F_arr_a),x_red) - F_a_ave*F_a_ave
	    Fa_var = np.dot(np.square(F_arr_a),x_red)
	    #F_sa_cov = np.dot(np.multiply(F_arr_s,F_arr_a),x_red) - F_s_ave*F_a_ave
	    F_sa_cov = np.dot(np.multiply(F_arr_s,F_arr_a),x_red)
	    
	    #debug
	    #print("F_s_ave = {}, F_a_ave = {}".format(F_s_ave,F_a_ave))
	    #print("dF_recon =",2*(Fs_var + F_sa_cov))
	    
	    # Selection:
	    dF_recon[l] = 2*(Fs_var + F_sa_cov)
	    Fs_sigma_sqr[l] = Fs_var
	    Fa_sigma_sqr[l] = Fa_var
	    Fsa_cov[l] = F_sa_cov

	    # Add the mutation contribution every delta_t_mut:
	    if mu > 0.0 and (l > 0 and l % mu_delta == 0):
	        ######
	        if MODE == "full":
	            Gvec = get_all_G(pc_red, means_precomp) @ x_red # should be a vector of <G> values, in reduced (n>0) subspace
	            Gvec = Gvec.ravel()
	            #debug
	            #print(F_arr[n_arr_pos])
	            #print("Gvec =",Gvec)
	            #for i in range(N):
	            #    if n_arr[i] == 0:
	            #        continue
	            #    ####
	            #    Gave = 0.0
	            #    for j in range(N):
	            #        if n_arr[j] > 0:
	            #            Gave += x_arr[j] * getG(pc[i,j],sigma)
	            #    ####
	            #    dF_mut += x_arr[i] * F_arr[i] * Gave
	        elif MODE == "basic":
	            Gvec = get_all_G(pc_arr[n_arr_pos].ravel(),means_precomp)
	            #for i in range(N):
	            #    if n_arr[i] == 0:
	            #        continue
	            #    ####
	            #    dF_mut += x_arr[i] * F_arr[i] * getG(pc_arr[i],sigma)
	        else: # MODE == "reduced"
	            a_red = a[n_arr_pos].ravel()
	            b_red = b[n_arr_pos].ravel()
	            ####
	            a_ave_red = np.dot(a_red,x_red)
	            b_ave_red = np.dot(b_red,x_red)
	            Gvec_a = get_all_G(a_red,means_precomp)
	            Gvec_b = get_all_G(b_red,means_precomp)
	            Gvec = np.add(b_ave_red * Gvec_a, a_ave_red * Gvec_b)
	            #for i in range(N):
	            #    if n_arr[i] == 0:
	            #        continue
	            #    ####
	            #    Gave = np.dot(b,x_arr)*getG(a[i],sigma) + np.dot(a,x_arr)*getG(b[i],sigma)
	            #    dF_mut += x_arr[i] * F_arr[i] * Gave
	        ######
	        dF_mut = np.dot(x_red, np.multiply(F_arr,Gvec))
	        dF_recon[l] += mu * dF_mut
	    
	    # Average fitness:
	    F_pop = np.dot(F_arr,x_red)
	    Ftraj[l] = F_pop # save the average fitness
	    #F_all[:,l] = np.copy(F_arr) # save all fitnesses

	    [Ftraj_lower_p[l],Ftraj_upper_p[l]] = weighted_quantile(F_arr, [LowerP,UpperP], sample_weight=x_red)
	    #Ftraj_lower_p[l]= np.percentile(F_arr, LowerP, method='inverted_cdf', weights=x_red)
	    #Ftraj_upper_p[l] = np.percentile(F_arr, UpperP, method='inverted_cdf', weights=x_red)
	    
	    ####
	    #if l==0:
	    #    print("Initial <F> = {:0.6f}".format(F_pop))
	    
	    # Number of species with >=1 individuals:
	    #Fnz_cur = F_arr[n_arr_pos] # array of fitnesses with at least one individual present
	    #Fnz_var[l] = np.var(Fnz_cur)
	    #Nspecies[l] = Fnz_cur.size
	    Nspecies[l] = F_arr[np.argwhere(x_red > 1./N)].size
	    
	    ########
	    #debug
	    #print("Before selection: x_arr =",x_arr)
	    # Compute new frequencies - selection:
	    x_red = x_red + delta_t * x_red * (F_arr - F_pop) # update non-zero freqs only
	    #for j in range(N):
	    #    x_arr[j] = x_arr[j] + delta_t * x_arr[j] * (F_arr[j] - F_pop)
	    
	    
	    if (l > 0 and l % mu_delta == 0) or l==Ltot-1:
	        
	        #debug
	        #print("After selection (l = {}): old x_red = {}".format(l,x_red))
	        
	        x_arr[n_arr_pos] = x_red.reshape(x_red.shape[0],1)
	        
	        #debug
	        #print("Before rounding: x_arr =",x_arr)
	        
	        n_arr[n_arr_pos] = x_arr[n_arr_pos] * N

	        #debug
	        #print("n_arr_pos =",n_arr_pos)
	        #print("Before rounding: n_arr =",n_arr)

	        n_arr = np.array(iteround.saferound(n_arr,0)) # snap counts to integers; preserves N
	        
	        #debug
	        #print("After rounding: x_arr =",n_arr/N)
	        
	        # Update zero-nonzero indices after rounding:
	        n_arr_pos = np.argwhere(n_arr > 0)
	        n_arr_pos_flat = n_arr_pos.ravel()
	        n_arr_zero_flat = np.argwhere(n_arr == 0).ravel()
	        
	        # Compute new counts - mutation:
	        if mu > 0.0 and (l > 0 and l % mu_delta == 0):
	            #debug
	            #print("n_arr_zero_flat =",n_arr_zero_flat)
	            #print("n_arr_pos_flat =",n_arr_pos_flat)
	            #print("old pc:")
	            #print(pc)
	            #print("**==**==**==**")

				# Final population state:
	            if MODE == "full":
	                prms_cur = (pc,)
	            elif MODE == "basic":
	                prms_cur = (pc,pc_arr,)
	            elif MODE == "reduced":
	                prms_cur = (pc,a,b)
	            else:
	                print(f"Unknown mode = {MODE}!")
	                ####
	                sys.exit("\nTerminating ..")

	            n_arr = mutate_counts(n_arr, mu, sigma, prms_cur, n_arr_zero_flat, n_arr_pos_flat, MODE)
	            ######
	            # Update zero-nonzero indices after mutation:
	            n_arr_pos = np.argwhere(n_arr > 0)
	            n_arr_pos_flat = n_arr_pos.ravel()
	            n_arr_zero_flat = np.argwhere(n_arr == 0).ravel()

	        x_arr = n_arr/N # update frequencies
	        
	        #debug
	        #print("After mutation: x_arr =",x_arr)
	        #print("new pc:")
	        #print(pc)
	        
	        #n_arr_pos = np.argwhere(n_arr > 0)
	        #n_arr_pos_flat = n_arr_pos.ravel()
	        #n_arr_zero_flat = np.argwhere(n_arr == 0).ravel()
	        ####
	        x_red = x_arr[n_arr_pos].ravel()
	        pc_red = pc[np.ix_(n_arr_pos_flat,n_arr_pos_flat)]
	        
	        #debug
	        #print("New x_red = {}".format(x_red))
	        #print("**************************")
	    
	    #t_cur += delta_t # update the time counter

	# Final population state:
	if MODE == "full":
		prms_out = (pc,)
	elif MODE == "basic":
		prms_out = (pc,pc_arr,)
	elif MODE == "reduced":
		prms_out = (pc,a,b)
	else:
		print(f"Unknown mode = {MODE}!")
		####
		sys.exit("\nTerminating ..")

	# Compute derivative of <F> by finite central differences:
	dF_num = np.zeros(Ltot)
	if Ltot > 2:
	    for k in range(1,Ltot-1):
	        dF_num[k] = (Ftraj[k+1] - Ftraj[k-1])/(2*delta_t)
	    ####
	    dF_num[0] = np.nan
	    dF_num[-1] = np.nan


	return [Ftraj,Ftraj_lower_p,Ftraj_upper_p,dF_recon,dF_num,Fs_sigma_sqr,Fa_sigma_sqr,Fsa_cov,Nspecies,prms_out,n_arr]

