########################
### utilities_min.py ###
########################

from math import *
import numpy as np
import sys
import re

#import matplotlib.pyplot as plt
#import seaborn as sns
#import pandas as pd

#import glob
#import os


# This function parses run_models_msp command line.
def parse_command_line_generic():

    args = {}

    argumentList = sys.argv[1:]

    if len(argumentList) == 0:
        print("No input arguments provided, please use '-h -amode=algorithm_name' or '--help -amode=algorithm_name' for help!")
        print("Allowed algorithm_name values: WrightFisher OR WF/ContinuousTime OR CT")
        sys.exit("\nTerminating ..")
    else:

        for arg in argumentList:

            if arg == '-h' or arg == '--help':
                continue

            toks = arg.split("=")
            assert len(toks) == 2, "Error: malformed token(s) in parse_command_line_generic!"
            assert type(toks[0]) is str and toks[0][0] == '-', "Error: malformed token(s) in parse_command_line_generic!"
            toks[0] = toks[0][1:] # get rid of the leading '-'

            if re.search('[a-zA-Z]', toks[1]) and not re.search('^[([]', toks[1]): # string with non-numeric entries, NOT an array
                args[toks[0]] = toks[1]
            else:
                try:
                    tok_val = eval(toks[1])
                except:
                    raise Exception("Error: malformed list token: {} in parse_command_line_generic!".format(toks[1]))
                ####
                if type(tok_val) is list or type(tok_val) is tuple:
                    if len(tok_val) == 3 or len(tok_val) == 4: # start,stop,num_vals
                        # Check stop > start:
                        assert tok_val[1] >= tok_val[0], \
                            "Error: start < stop in parse_command_line_generic!"
                        # Check num_vals:
                        assert type(tok_val[2]) is int and tok_val[2] > 0, \
                            "Error: provide a non-negative integer for num_vals in parse_command_line_generic!"
                    ######
                    if len(tok_val) == 3: # start,stop,num_vals
                        ####
                        if tok_val[2] == 1:
                            if tok_val[1] == tok_val[0]:
                                if type(tok_val[0]) is int and type(tok_val[1]) is int:
                                    tok_arr = [int(tok_val[0])]
                                else:
                                    tok_arr = [float(tok_val[0])]
                            else:
                                raise Exception("Unable to create an array with start={}, stop={}, num_vals={} in parse_command_line_generic!".\
                                    format(tok_val[0],tok_val[1],tok_val[2]))
                        else:
                            tok_arr = list(np.linspace(tok_val[0], tok_val[1], num=tok_val[2], dtype=np.float64))
                            step = (tok_arr[-1] - tok_arr[0])/(len(tok_arr)-1)
                            ####
                            if type(tok_val[0]) is int and type(tok_val[1]) is int and abs(step - int(step)) < 1e-10: # cast array to int
                                for i in range(len(tok_arr)):
                                    tok_arr[i] = int(tok_arr[i])
                            else: # cast array to float
                                for i in range(len(tok_arr)):
                                    tok_arr[i] = float(tok_arr[i])
                            
                    elif len(tok_val) == 4: # start,stop,num_vals,log_token
                        if type(tok_val[3]) == str and tok_val[3] == 'log10':
                            tok_arr = list(np.linspace(np.log10(tok_val[0]), np.log10(tok_val[1]), \
                                    num=tok_val[2], dtype=np.float64))
                            for i in range(len(tok_arr)):
                                tok_arr[i] = float(10**tok_arr[i])
                        else:
                            raise Exception("Unknown log_token={} in parse_command_line_generic!".format(tok_val[3]))
                    else:
                        raise Exception("Expected 3 or 4 entries in {}!".format(toks[1]))
                    ######
                    args[toks[0]] = tok_arr
                else: # int,float
                    args[toks[0]] = tok_val


        if '-h' in argumentList or '--help' in argumentList:

            if not 'amode' in args:
                print("Please provide -amode=algorithm_name with -h or --help for algorithm-specific help!")
                print("Allowed algorithm_name values: WrightFisher OR WF/ContinuousTime OR CT")
                sys.exit("\nTerminating ..")
            else:
                amode = args['amode']

                if amode == "WrightFisher" or amode == "WF":
                    print("python3 {} -amode=WF [-mu_rate=mu_val OR -mu_rate=(mu_min,mu_max,mu_num[,\"log10\"]) OR -mu_rate=[mu_min,mu_max,mu_num[,\"log10\"]]] \\".format(sys.argv[0]))
                    print("[-sigma=s_val OR -sigma=(s_min,s_max,s_num[,\"log10\"]) OR -sigma=[s_min,s_max,s_num[,\"log10\"]]] \\")
                    print("[-A_PD=A_val OR -A_PD=(A_min,A_max,A_num[,\"log10\"]) OR -A_PD=[A_min,A_max,A_num[,\"log10\"]]] \\")
                    print("[-r_PD=r_val OR -r_PD=(r_min,r_max,r_num[,\"log10\"]) OR -r_PD=[r_min,r_max,r_num[,\"log10\"]]] \\")
                    print("-ltot=ltot_val -Npop=N_val -mode=full/basic/reduced -out=outfilename \\")
                    print("[-nruns=nr_val -sdata={0,1}] \\")
                    print("[-pc_in=pc_in.dat -pc_out=pc_out.dat] \\")
                    print("[-pc_arr_in=pc_arr_in.dat -pc_arr_out=pc_arr_out.dat] \\")
                    print("[-a_in=a_in.dat -a_out=a_out.dat] \\")
                    print("[-b_in=b_in.dat -b_out=b_out.dat]")
                elif amode == "ContinuousTime" or amode == "CT":
                    print("python3 {} -amode=CT [-mu_rate=mu_val OR -mu_rate=(mu_min,mu_max,mu_num[,\"log10\"]) OR -mu_rate=[mu_min,mu_max,mu_num[,\"log10\"]]] \\".format(sys.argv[0]))
                    print("[-sigma=s_val OR -sigma=(s_min,s_max,s_num[,\"log10\"]) OR -sigma=[s_min,s_max,s_num[,\"log10\"]]] \\")
                    print("[-r_PD=r_val OR -r_PD=(r_min,r_max,r_num[,\"log10\"]) OR -r_PD=[r_min,r_max,r_num[,\"log10\"]]] \\")
                    print("[-mu_delta=mu_dt_val OR -mu_delta=(mu_dt_min,mu_dt_max,mu_dt_num[,\"log10\"]) OR -mu_delta=[mu_dt_min,mu_dt_max,mu_dt_num[,\"log10\"]]] \\")
                    print("-ltot=ltot_val -Npop=N_val -mode=full/basic/reduced -out=outfilename \\")
                    print("[-nruns=nr_val -Nu=Nu_val -delta_t=dt_val -sdata={0,1}] \\")
                    print("[-counts_init=counts_init.dat -counts_out=counts_out.dat] \\")
                    print("[-pc_in=pc_in.dat -pc_out=pc_out.dat] \\")
                    print("[-pc_arr_in=pc_arr_in.dat -pc_arr_out=pc_arr_out.dat] \\")
                    print("[-a_in=a_in.dat -a_out=a_out.dat] \\")
                    print("[-b_in=b_in.dat -b_out=b_out.dat]")
                else:
                    print("Unknown amode = {} in {}!".format(amode,sys.argv[0]))
                ####
                sys.exit("\nTerminating ..")

    return args


# This function reads data from a file output by run_models_msp.py
# into multiple Numpy matrices.
def read_data_generic(filename):

    f1 = open(filename, "r")

    header_cnt = 0
    line_cnt = 0

    data_lines = []

    for x in f1:
        toks = x.split()
        if len(toks) == 0:
            continue # skip empty lines
        
        if toks[0] == '#' or toks[0][0] == '#':

            if header_cnt == 0:
                assert toks[1] == "amode", "Error: unexpected header line in read_data_generic!"
                amode = toks[3]
                amode = amode[:-1] # get rid of the trailing ','
                assert toks[4] == "mode", "Error: unexpected header line in read_data_generic!"
                mode = toks[6]
                #mode = mode[:-1] # get rid of the trailing ','
            elif header_cnt == 1:
                dim = []
                for x in toks[1].split('x'):
                    dim.append(int(x))
                assert len(dim) == 5, "Error: dimension mismatch in read_data_generic!"
                ####
                dtot = 1
                for x in dim: # l1*l2*l3*l4*nruns
                    dtot *= x
                ####
                #D = int(toks[-1]) # tuple dim
                ####
                if amode == "WrightFisher" or amode == "WF":
                    prm1_arr = np.zeros(dtot,dtype=np.float64).reshape((dim[0],dim[1],dim[2],dim[3],dim[4]))
                    prm2_arr = np.zeros(dtot,dtype=np.float64).reshape((dim[0],dim[1],dim[2],dim[3],dim[4]))
                    prm3_arr = np.zeros(dtot,dtype=np.float64).reshape((dim[0],dim[1],dim[2],dim[3],dim[4]))
                    prm4_arr = np.zeros(dtot,dtype=np.float64).reshape((dim[0],dim[1],dim[2],dim[3],dim[4]))
                elif amode == "ContinuousTime" or amode == "CT":
                    prm1_arr = np.zeros(dtot,dtype=np.float64).reshape((dim[0],dim[1],dim[2],dim[3],dim[4]))
                    prm2_arr = np.zeros(dtot,dtype=np.float64).reshape((dim[0],dim[1],dim[2],dim[3],dim[4]))
                    prm3_arr = np.zeros(dtot,dtype=np.float64).reshape((dim[0],dim[1],dim[2],dim[3],dim[4]))
                    prm4_arr = np.zeros(dtot,dtype=np.int64).reshape((dim[0],dim[1],dim[2],dim[3],dim[4]))
                else:
                    raise Exception("Unknown amode: {} in read_data_generic!".format(amode))
                ####
                Ffirst_arr = np.zeros(dtot).reshape((dim[0],dim[1],dim[2],dim[3],dim[4]))
                Flast_arr = np.zeros(dtot).reshape((dim[0],dim[1],dim[2],dim[3],dim[4]))

            header_cnt += 1

            continue # skip commented lines

        else:
            assert len(toks) == 6, "Error: num_columns mismatch in read_data_generic!"

            if re.search('[a-df-zA-Z]', toks[-1]) or re.search('[a-df-zA-Z]', toks[-2]): # a string with non-numeric entries, supposed to be column names
                #debug
                #print(toks[-1])
                continue # this is a header
            else:
                data_lines.append(x)
                line_cnt += 1

    #debug
    #print("dtot =",dtot)
    #print("line_cnt =",line_cnt)
    assert line_cnt == dtot, "Error: data dimension mismatch in read_data_generic!"

    f1.close()

    ######
    cnt = 0
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                for l in range(dim[3]):
                    for m in range(dim[4]):
                        toks = data_lines[cnt].split()
                        if amode == "WrightFisher" or amode == "WF":
                            prm1_arr[i,j,k,l,m] = float(toks[0])
                            prm2_arr[i,j,k,l,m] = float(toks[1])
                            prm3_arr[i,j,k,l,m] = float(toks[2])
                            prm4_arr[i,j,k,l,m] = float(toks[3])
                        elif amode == "ContinuousTime" or amode == "CT":
                            prm1_arr[i,j,k,l,m] = float(toks[0])
                            prm2_arr[i,j,k,l,m] = float(toks[1])
                            prm3_arr[i,j,k,l,m] = float(toks[2])
                            prm4_arr[i,j,k,l,m] = int(toks[3])
                        else:
                            raise Exception("Unknown amode: {} in read_data_generic!".format(amode))
                        ####
                        Ffirst_arr[i,j,k,l,m] = float(toks[4])
                        Flast_arr[i,j,k,l,m] = float(toks[5])
                        
                        ####
                        cnt += 1

    return [dim,prm1_arr,prm2_arr,prm3_arr,prm4_arr,Ffirst_arr,Flast_arr,amode,mode]


# Converts integer seconds into HH:MM:SS format:
def convert(seconds):

    seconds = seconds % (24 * 3600)
    hour = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
     
    return "%dH:%02dM:%02dS" % (hour, minutes, seconds)


# Reads fixed-format data from a trajectory file:
def read_traj_data(datafile,amode):
    
    if amode == "WrightFisher" or amode == "WF":
        ltot = []
        Fave = []
        Fstd = []
        Flower = []
        Fupper = []
    elif amode == "ContinuousTime" or amode == "CT":
        ltot = []
        Fave = []
        Flower = []
        Fupper = []
        dF_recon = []
        dF_num = []
        Fs_var = []
        Fa_var = []
        Fsa_cov = []
        Nspecies =[]
    else:
        raise Exception(f"Unknown amode: {amode} in read_traj_data!") 

    f1 = open(datafile, "r")
    
    flag_header = 1
    for x in f1:
        toks = x.split()
        toks[-1] = toks[-1].strip('\n')
        #debug
        #print(toks," ",len(toks))
        if (len(toks) == 0):
            continue # skip empty lines
        elif (toks[0] == '' or toks[0] == '#' or toks[0][0] == '#'):
            continue # skip commented lines
        else:
            if flag_header == 1:
                # Header check:
                if amode == "WrightFisher" or amode == "WF":
                    assert len(toks) == 5, "Error: dimension mismatch in read_traj_data!"
                    assert toks[0] == "ltot" and toks[1] == "Fave" and toks[2] == "Fstd" and \
                            toks[3] == "Flow" and toks[4] == "Fupp", "Error: unexpected header in read_traj_data!"
                elif amode == "ContinuousTime" or amode == "CT":
                    assert len(toks) == 10, "Error: dimension mismatch in read_traj_data!"
                    assert toks[0] == "ltot" and toks[1] == "Fave" and toks[2] == "Flow" and toks[3] == "Fupp" and \
                            toks[4] == "dF_recon" and toks[5] == "dF_num" and toks[6] == "Fs_var" and \
                            toks[7] == "Fa_var" and toks[8] == "Fsa_cov" and toks[9] == "Nspecies", \
                            "Error: unexpected header in read_traj_data!"
                ####
                flag_header = 0
            else:
                if amode == "WrightFisher" or amode == "WF":
                    assert len(toks) == 5, "Error: dimension mismatch in read_traj_data!"
                    ltot.append(int(toks[0]))
                    Fave.append(float(toks[1]))
                    Fstd.append(float(toks[2]))
                    Flower.append(float(toks[3]))
                    Fupper.append(float(toks[4]))
                elif amode == "ContinuousTime" or amode == "CT":
                    assert len(toks) == 10, "Error: dimension mismatch in read_traj_data!"
                    ltot.append(int(toks[0]))
                    Fave.append(float(toks[1]))
                    Flower.append(float(toks[2]))
                    Fupper.append(float(toks[3]))
                    dF_recon.append(float(toks[4]))
                    dF_num.append(float(toks[5]))
                    Fs_var.append(float(toks[6]))
                    Fa_var.append(float(toks[7]))
                    Fsa_cov.append(float(toks[8]))
                    Nspecies.append(int(toks[9]))

    f1.close()

    if amode == "WrightFisher" or amode == "WF":
        data_out = [ltot,Fave,Fstd,Flower,Fupper]
    elif amode == "ContinuousTime" or amode == "CT":
        data_out = [ltot,Fave,Flower,Fupper,dF_recon,dF_num,Fs_var,Fa_var,Fsa_cov,Nspecies]
    else:
        raise Exception(f"Unknown amode: {amode} in read_traj_data!")

    return data_out


# Processes multiple output files, concatenates data arrays.
def read_multiple_files(filenames):

    print("Reading data from {} ..".format(filenames[0]))

    dim,prm1_arr,prm2_arr,prm3_arr,prm4_arr,Ffirst_arr,Flast_arr,amode,mode = \
            read_data_generic(filenames[0])

    # Keep the start and the end index of each original dataset in the final concatenated arrays:
    sind = []
    find = []
    sind.append(0)
    find.append(dim[4]-1)

    # Read in additional runs if provided:
    for j in range(1,len(filenames)):
        dim_cur,prm1_arr_cur,prm2_arr_cur,prm3_arr_cur,prm4_arr_cur,Ffirst_arr_cur,Flast_arr_cur,amode_cur,mode_cur = \
            read_data_generic(filenames[j])
        
        assert amode_cur == amode and mode_cur == mode, \
                "Error: mode mismatch, unable to merge datafiles!"
        
        assert dim_cur[0] == dim[0] and dim_cur[1] == dim[1] and dim_cur[2] == dim[2] and dim_cur[3] == dim[3], \
                "Error: dim mismatch, unable to merge datafiles!"
        
        # Concatenate data matrices:
        prm1_arr = np.concatenate((prm1_arr, prm1_arr_cur), axis=4)
        prm2_arr = np.concatenate((prm2_arr, prm2_arr_cur), axis=4)
        prm3_arr = np.concatenate((prm3_arr, prm3_arr_cur), axis=4)
        prm4_arr = np.concatenate((prm4_arr, prm4_arr_cur), axis=4)
        ######
        Ffirst_arr = np.concatenate((Ffirst_arr, Ffirst_arr_cur), axis=4)
        Flast_arr = np.concatenate((Flast_arr, Flast_arr_cur), axis=4)
        
        # Keep the start and the end index of each original dataset in the final concatenated arrays:
        sind.append(find[j-1]+1)
        find.append(sind[j]+dim_cur[4]-1)
        
        # Augment the dim array:
        dim[4] += dim_cur[4]
        
        print("Concatenated data from {} ..".format(filenames[j]))

    ######
    return [dim,sind,find,prm1_arr,prm2_arr,prm3_arr,prm4_arr,Ffirst_arr,Flast_arr,amode,mode]


# "Smart" rounding which handles very small and very large numbers correctly.
def precision_round(number, digits=3):
    digits -= 1 # this is to make the rounding more intuitive: e.g. 0.123 -> 0.123 with digits=3, 0.123 -> 0.12 with digits=2
    power = "{:e}".format(number).split('e')[1]
    return round(number, -(int(power) - digits))


# This aux. function post-processes 'array' arguments:
def postprocess_val(tok_string,tok_type,args): # tok_type = int/float/str

    if tok_string in args:
        tok = args[tok_string]
        #debug
        #print(tok_type)
        #print("tok =",tok," type =",type(tok))
        if type(tok) is tok_type: # single value, convert to array
            tok_arr = []
            tok_arr.append(tok)
            tok = tok_arr
            #tok = eval('['+str(tok)+']')
        elif type(tok) is list:
            for x in tok:
                if not type(x) is tok_type:
                    raise Exception("Invalid type of {} value: {}".format(tok_string,type(x)))
        else:
            raise Exception("Invalid value of {}: {}".format(tok_string,tok))

    else:
        raise Exception("Please pass {} as an {} argument!".format(tok_string,sys.argv[0]))

    return tok


# This function outputs run data to a file in a flexible format.
def write_data(filename,comment_str,header_arr,data_arr):
    assert len(header_arr) == len(data_arr), "Error: input arrays must have equal lengths!"

    header_str = ''
    for j in range(len(header_arr)-1):
        header_str += "{:<4}         ".format(header_arr[j])
    header_str += "{:<4}\n".format(header_arr[-1])

    f2 = open(filename,"w")
    f2.write(comment_str)
    f2.write(header_str)
    
    for j in range(len(data_arr[0])):
        cur_str = ''
        for k in range(len(data_arr)):
            if type(data_arr[k][j]) is int or type(data_arr[k][j]) is np.int64:
                cur_str += "{:<6}    ".format(data_arr[k][j])
            elif type(data_arr[k][j]) is float or type(data_arr[k][j]) is np.float64:
                cur_str += "{:<.4e}    ".format(data_arr[k][j])
            else:
                raise Exception("Unknown data type: {} in write_data!".format(type(data_arr[k][j])))

        ####
        cur_str += "\n"
        #cur_str += "{}\n".format(data_arr[-1][j]) # this slot is for the best_state tuple
        f2.write(cur_str)
    
    f2.close()


def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)

    if sample_weight is None:
        sample_weight = np.ones(len(values))

    sample_weight = np.array(sample_weight)

    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'Quantiles should be in the [0, 1] range in weighted_quantile!'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight

    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)

    return np.interp(quantiles, weighted_quantiles, values)


