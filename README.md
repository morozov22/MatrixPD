# MatrixPD Manual

## Developer: Alexandre V. Morozov (morozov@physics.rutgers.edu)

## Paper: A.V. Morozov and A. Feigel, "Emergence of cooperation due to opponent-specific responses in Prisoner's Dilemma"

<br/>

## MatrixPD code:

* <span style="color:blue"><font size="4"> __run_models_msp.py__ (high-level Python3 script for running evolutionary simulations in batch mode) </font></span>

* <span style="color:blue"><font size="4"> __functions_msp.py__ (evolutionary simulations library) </font></span>

* <span style="color:blue"><font size="4"> __utilities_msp.py__ (auxiliary functions library) </font></span>

`run_models_msp.py` is a high-level script for running either Wright-Fisher evolutionary simulations (non-overlapping generations) with selection and mutations (the WF algorithm) or solving a system of continuous-time replicator equations augmented with mutations (the CT algorithm). In both cases, mutations act directly on the individual's probabilities to cooperate.

<br/>

Run this command:
```
python3 run_models_msp.py -h
```
or
```
python3 run_models_msp.py --help
```
to obtain the list of allowed algorithm names.

<br/>

Run this command:
```
python3 run_models_msp.py -h -amode=algorithm_name
```
or
```
python3 run_models_msp.py --help -amode=algorithm_name
```
where `algorithm_name = {WrightFisher || WF, ContinuousTime || CT}` to
obtain an algorithm-specific list of available options.
Most of these options are self-explanatory and the notation corresponds to that used in the manuscript.
All the options are described in more detail below.

<br/>

Please note that some options accept a single value or an array of
values either on a linear or a $\log_{10}$ scale. For example, `-r_PD=0.3` will set $r=0.3$ in all runs,
`-r_PD=\[0.1,1.0,10\]` will perform a scan over 10 uniformly spaced values of the Prisonner's Dilemma (PD)
strength parameter: $r = \{0.1, 0.2, \dots , 1.0\}$,
while `-r_PD=\[0.1,10.0,3,\"log10\"\]` will result in $r = \{0.1, 1.0, 10.0\}$.
In all array-type options, parentheses `()` can be used instead of brackets `[]`. The library that handles input options tries to guess intelligently whether the resulting array values
should be *int* or *float*, although all MatrixPD array-type input parameters should be of the *float* type.
__Note:__ The backslashes in the array-type options were tested on a Mac and should work on any Linux box. For Windows machines, try `-r_PD="[0.1,10.0,3,'log10']"` instead - some experimentation may be required.

<br/>

The `-mode=full/basic/reduced` option determines the type of the cooperativity matrix:
`-mode=full` specifies the most general model characterized by the opponent-specific probabilities to cooperate, $p_{c, i \to j}$, `-mode=basic` corresponds to the opponent-blind model, $p_{c, i \to j} = p_{c,i}$, and `-mode=reduced` invokes the low-rank approximation to the cooperativity matrix (not used in the paper): $p_{c, i \to j} = a_i b_j$.

The `-out=outfilename` option specifies an output filename with a summary of the run.

The `-nruns=nr_val` option specifies the number of runs per unique combination of input parameter settings (the default is $20$). Multiple runs with the same parameter settings are used to collect statistics in stochastic simulations.

The `-ltot=ltot_val` option is used to provide the number of generations in the WF algorithm and the number
of timesteps in the CT algorithm. Each timestep advances the CT simulation by $\delta t$ provided via the
`-delta_t=dt_val` option (the default is $\delta t = 0.005$).

The `-Npop=N_val` option sets the total number of individuals in the population. In the CT algorithm,
the `-Nu=Nu_val` option sets the initial number of genotypes (novel genotypes can be created
only through mutations). Note that $N_u \le N_\text{pop}$; the defauilt is $N_u = N_\text{pop}$.

In the WF algorithm, the array-type options are `-mu_rate`, `-sigma`, `-A_PD` and `-r_PD` (the mutation rate $\mu$, the truncated Gaussian scale parameter $\sigma$, the PD scale factor $A$ and the PD strength
parameter $r$, respectively). In the CT algorithm, the array-type options are `-mu_rate`, `-sigma`, `-r_PD`
and `-mu_delta` (the mutation rate $\mu$, the truncated Gaussian scale parameter $\sigma$,
the PD strength parameter $r$ and the number of timesteps $\Delta_\mu$ between subsequent rounding/mutation
events, respectively).

The `-pc_in=pc_in.dat` option is restricted to the `full` mode and allows the user to provide an $N_\text{pop} \times N_\text{pop}$ matrix of the probabilities to cooperate, $\mathbf{P}_c$: $(\mathbf{P}_c)_{ij} = p_{c, i \to j}$. If the input matrix is not provided, it is generated randomly by sampling from the uniform distribution.

The `-pc_arr_in=pc_arr_in.dat` option is restricted to the `basic` mode and allows the user to provide
a vector of the probabilities of cooperate of length $N_\text{pop}$, $\mathbf{p}_c$: $(\mathbf{p}_c)_{i} = p_{c,i}$. If the input vector is not provided, it is generated randomly by sampling from the uniform distribution.

The `-a_in=a_in.dat` and `-b_in=b_in.dat` options are restricted to the `reduced` mode and allow the user to provide two vectors of length $N_\text{pop}$, $\mathbf{a}$ and $\mathbf{b}$, that are used to construct the matrix of the probabilities to cooperate: $p_{c, i \to j} = a_i b_j$. If the input vectors is not provided, they are generated randomly by sampling from the uniform distribution.

The corresponding `-pc_out=pc_out.dat`, `-pc_arr_out=pc_arr_out.dat`, `-a_out=a_out.dat` and
`-b_out=b_out.dat` options allow the user to output the __initial__ $\mathbf{P}_c$, $\mathbf{p}_c$,
$\mathbf{a}$ and $\mathbf{b}$, respectively. These initial probabilites are used to initialize each run in
the current batch of simulations. The `-*_out` options are useful when the data structures created
inside the code need to be saved.

In the CT algorithm, the `-counts_init=counts_init.dat` option is used to
provide a vector of length $N_\text{pop}$ which contains the initial genotype counts in the population.
The `-Nu=Nu_val` option is ignored in this case. The `-counts_out=counts_out.dat` is used to save the __initial__ genotype counts.

If the `-sdata={0,1}` option is set to $1$ (the default is $0$), trajectory data will be output in a separate file for each run in the series. In addition, the __final__ matrix of the probabilities to cooperate $\mathbf{P}_c$ and the corresponding mode-dependent vectors used to construct $\mathbf{P}_c$:
$\mathbf{p}_c$, $\mathbf{a}$ and $\mathbf{b}$ are output to separate files. Finally,
a file with the __final__ genotype counts is output if the algorithm type is CT (there are no
genotype counts in WF because each individual in the population is stored separately).
All filenames are created automatically based on the user-provided output filename
(via the `-out=outfilename` option).

<br/><br/>


The main output file of a typical `run_models_msp.py` looks as follows:
```
# amode = CT, mode = full
# 1x1x3x1x1 array
# nruns = 1, ltot = 500000, Npop = 500, dt = 0.005, sdata = 1
mu            sigma           r           mu_delta          Finit           Ffin
2.0000e-01    2.0000e-02    3.0000e-01    10000        5.0051e-01    7.7596e-01    
2.0000e-01    2.0000e-02    5.0000e-01    10000        5.0051e-01    8.6551e-01    
2.0000e-01    2.0000e-02    7.0000e-01    10000        5.0051e-01    9.3531e-01 
```

This is a series of 3 CT runs in the `full` mode, each containing $5 \times 10^5$ timesteps with
a time increment $\delta t = 0.005$, for $N_\text{pop} = 500$. The file contains the initial
fitness `Finit` and the final fitness `Ffin` for each combination of array-type hyperparameters
`mu`, `sigma`, `r` and `mu_delta` (three distinct sets of hyperparameters in this case,
corresponding to $r = (0.3,0.5,0.7)$). 
These user-provided input parameter ranges have resulted in a series of $1 \times 1 \times 3 \times 1 \times 1 = 3$ runs, where the last dimension of the 5D array was provided via the `-nruns=1` option.


The output file was generated using the following command (available as `run_CT_1.cmd` in the `examples/` folder):
```
python3 ../run_models_msp.py -amode=CT -mu_rate=0.2 -sigma=0.02 -r_PD=\(0.3,0.7,3\) -mu_delta=10000 -ltot=500000 -Npop=500 -mode=full -out=run_CT_1.dat -nruns=1 -Nu=500 -delta_t=0.005 -sdata=1 -counts_out=run_CT_1.n.init.dat
```

Note that rounding of genotype frequencies and mutations occur every $10^4$ steps (since `-mu_delta=10000`), so that evolutionary dynamics consists of solving a system of replicator equations with selection only, interspersed by mutational events that change the counts of individuals to create novel phenotypes.
Setting `mu_delta > ltot` would result in disabling both genotype frequency rounding and mutations,
so that evolutionary dynamics would consist of solving a system of replicator equations with selection only
(genotype frequencies are rounded once in the end to produce final genotype counts).

Since the `-sdata` option was set to $1$, the code has produced several output files, all available in the `examples/` folder. In particular, the `run_CT_1.dat.f.[1-3]` trajectory files contain the following data for each step of the run:


| Variable    | Meaning |
| -------- | ------- |
| ltot    | current timestep |
| Fave    | average fitness $\langle F \rangle$ |
| Flow    | 25% percentile of population fitness |
| Fupp    | 75% percentile of population fitness |
| dF_recon    | $d \langle F \rangle /  dt$ reconstructed using variance-covariance theorem|
| dF_num    | $d \langle F \rangle /  dt$ computed numerically using finite differences|
| Fs_var    | $\text{Var} (F^s)$ |
| Fa_var    | $\text{Var} (F^a)$ |
| Fsa_cov    | $\text{Cov} (F^s, F^a)$ |
| Nspecies    | number of genotypes with frequencies $>\!1/N_\text{pop}$|


The `examples/` folder contains the total of four representative runs, two CT and two WF.
The output of each run was produced using the Python3 command found in the corresponding `run_*.cmd` file.
The `run_CT_1` run is discussed above.
The `run_CT_2` run shows how the descendants of a populist quickly take over the population under selective forces.

The `run_WF_1` and `run_WF_2` runs are Wright-Fisher simulations with selection and mutation
in the `basic` and `full` modes, respectively. In both cases, $5$ independent runs starting from the same initial condition were carried out, as requested via the `-nruns=5` option.
Since the `-sdata` option was set to $1$ in both runs, the code has produced several output files, all available in the `examples/` folder. In particular, the `run_WF_1.dat.f.[1-5]` and `run_WF_2.dat.f.[1-5]` trajectory files contain the following data for each generation:


| Variable    | Meaning |
| -------- | ------- |
| ltot    | current generation number |
| Fave    | average fitness $\langle F \rangle$ |
| Fstd    | standard deviation $\sigma_F$ |
| Flow    | 25% percentile of population fitness |
| Fupp    | 75% percentile of population fitness |


<br/><br/>                                                                             
                                                                               
The `utilities_msp.py` library contains several functions designed to work with the `run_models_msp.py` output:

* <span style="color:green"><font size="3"> __read_data_generic(*filename*)__ - reads a single output file into multi-D Numpy arrays </font></span>

* <span style="color:green"><font size="3"> __read_multiple_files(*filenames*)__ - reads multiple output files and concatenates the data into multi-D Numpy arrays; typically used to combine output of multiple parallel runs </font></span>

* <span style="color:green"><font size="3"> __read_traj_data(*datafile*, *amode*)__ - reads a single file with the trajectory data; __amode__ refers to CT or WF </font></span>

These functions are provided in order to enable subsequent processing and plotting of the output data.

