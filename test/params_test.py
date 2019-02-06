## Parameters used in script_test.py ##

mcmc=False# if False, plot with linear fit. If True, plots with mcmc

########## Definition of the output file ##########

output_file_mcmc='./results_fit_sed_mcmc'
output_file_linear='./result_fit_sed_mat'
output_file_interpolation='./results_interpolation'
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
data_file='./data_files/data_13dqy_formatted_for_package.txt' #must have at least the following fiels, and a header with the fields written as such: jd,flux,fluxerr,filter.
# filter can be one of the following: ('UVW1','UVW2','UVM2','u_swift','v_swift','b_swift','g_sdss','r_sdss','i_sdss','z_sdss'
# ,'r_cousin','i_cousin','h_2mass','j_2mass','u_johnson','b_johnson','v_johnson')
dates_file='./data_files/13dqy_int_dates.txt'

lower_limit_on_flux=1e-40
filters_directory='/Users/maayanesoumagnac/AstroCosmo/GitHub_repositories/PhotoFit/PhotoFit/Filters' #put the path to the Filters directory here

# Interpolation step
already_run_interp_errors=dict() #don't touch this line
already_run_interp_errors['UVW1']=True
already_run_interp_errors['UVW2']=True
already_run_interp_errors['UVM2']=True
already_run_interp_errors['u_swift']=False
already_run_interp_errors['b_swift']=False
already_run_interp_errors['v_swift']=False
already_run_interp_errors['r_sdss']=True
already_run_interp_errors['g_sdss']=True
already_run_interp_errors['i_sdss']=True
already_run_interp_errors['r_p48']=True
already_run_interp_errors['g_p48']=True
already_run_interp_errors['z_sdss']=True
already_run_interp_errors['u_johnson']=True
already_run_interp_errors['v_johnson']=True
already_run_interp_errors['b_johnson']=True
already_run_interp_errors['i_cousin']=True
already_run_interp_errors['r_cousin']=True
already_run_interp_errors['j_2mass']=True
already_run_interp_errors['h_2mass']=True


# In case you fit with a linear-fitting algorythm
already_run_matrix=True
num_iterations=100

# In case you fit with mcmc
already_run_mcmc=False
nwalkers=100
num_steps=350

# In case you want to compare your R and T results with existing results from a file
data_compare='./data_files/Yaron_2017_results.txt'#file with column#1: time from explosion, column#2: temperature (K), column#3:radius (cm)
