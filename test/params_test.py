import os
## Parameters used in script_test.py ##

mcmc=True# if False, plot with linear fit. If True, plots with mcmc

########## Definition of the output file ##########

output_file_mcmc='./results_fit_sed_mcmc'
output_file_linear='./result_fit_sed_mat'
#output_file_interpolation='./results_interpolation'

plots=False

if mcmc==True:
    output=output_file_mcmc
else:
    output=output_file_linear
if os.path.exists(output):
    print('the output directory, '+output + ' exists already')
else:
    os.mkdir(output)

########## Definition of the object parameters ##########

z =0
distance_modulus=33.51
explosion_date=0 #in principle, this has to be in jd. For the test data it is in jd-explosion dates
EBV=0.035
csm=False

data_file='./data_files/data_13dqy_formatted_for_package.txt' #must have at least the following fiels, and a header with the fields written as such: jd,flux,fluxerr,filter.
# filter can be one of the following: ('UVW1','UVW2','UVM2','u_swift','v_swift','b_swift','g_sdss','r_sdss','i_sdss','z_sdss'
# ,'r_cousin','i_cousin','h_2mass','j_2mass','u_johnson','b_johnson','v_johnson')
dates_file='./data_files/13dqy_int_dates.txt'

lower_limit_on_flux=1e-40
filters_directory='../PhotoFit/Filters' #put the path to the Filters directory here

# Interpolation step
#already_run_interp_errors=dict() #don't touch this line
#already_run_interp_errors['UVW1']=            False
#already_run_interp_errors['UVW2']=            False
#already_run_interp_errors['UVM2']=            False
#already_run_interp_errors['u_swift']=         False
#already_run_interp_errors['b_swift']=         False
#already_run_interp_errors['v_swift']=         False
#already_run_interp_errors['r_sdss']=          False
#already_run_interp_errors['g_sdss']=          False
#already_run_interp_errors['i_sdss']=          False
#already_run_interp_errors['r_p48']=           False
#already_run_interp_errors['g_p48']=           False
#already_run_interp_errors['z_sdss']=          False
#already_run_interp_errors['u_johnson']=       False
#already_run_interp_errors['v_johnson']=       False
#already_run_interp_errors['b_johnson']=       False
#already_run_interp_errors['i_cousin']=        False
#already_run_interp_errors['r_cousin']=        False
#already_run_interp_errors['j_2mass']=         False
#already_run_interp_errors['h_2mass']=         False
#

# In case you fit with a linear-fitting algorythm
already_run_matrix=False
num_iterations=100

# In case you fit with mcmc
already_run_mcmc=True
nwalkers=80
num_steps=350

# In both case: either None or a list of Boolean of size the number of epochs, where True is for epochs already ran and False is for epoch left to run
already_run_fit=None

excluded_bands=[]

#Setting the priors on T and R:
priors=True
lowrad=[0.5e14,0.5e14,0.5e14,0.5e14,0.5e14,0.5e14,1e14,1e14,1e14,1e14,1e14,1e14,1e14,1e14,1e14,8e14,8e14,8e14,8e14,8e14,8e14,8e14,8e14,8e14,8e14]
hirad=[5e14,5e14,5e14,5e14,5e14,5e14,1e15,1e15,1e15,1e15,1e15,1e15,1e15,1e15,1e15,2e15,2e15,2e15,2e15,2e15,2e15,2e15,2e15,2e15,2e15]
lowtemp=[1e4,1e4,1e4,1e4,1e4,1e4,1e4,1e4,5e3,5e3,5e3,5e3,5e3,5e3,5e3,5e3,5e3,5e3,5e3,5e3,5e3,4e3,4e3,4e3,4e3]
hitemp=[4e4,4e4,3e4,3e4,3e4,3e4,3e4,3e4,1.8e4,1.8e4,1.8e4,1.8e4,1.8e4,1.8e4,1.8e4,1.8e4,1.8e4,1.8e4,1e4,1e4,1e4,1e4,1e4,1e4,1e4]



# In case you want to compare your R and T results with existing results from a file
data_compare='./data_files/Yaron_2017_results.txt'#file with column#1: time from explosion, column#2: temperature (K), column#3:radius (cm)



