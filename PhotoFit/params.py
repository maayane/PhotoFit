## Parameters used in script.py ##
mcmc=False# if False, plot with linear fit. If True, plots with mcmc

output_file_mcmc='./test/results_fit_sed_mcmc'
output_file_linear='./test/result_fit_sed_mat'
output_file_interpolation='./test/results_interpolation'
z =0#0.011855 #2018fif: 0.017189 #redshift to correct for.
distance_modulus=33.51#2018fif: 34.31
explosion_date=0# 2018fif: 2458351.653729907237
EBV=0.035#2018fif: 0.13564
data_file='/Users/maayanesoumagnac/PostDoc/projects/2018fif/comparision_with_others/data_13dqy_formatted_for_package.txt'
#2018fif: '/Users/maayanesoumagnac/PostDoc/projects/2018fif/data/data_ZTF_Swift_3_subtracted.txt' #must have at least the following fiels, and a header with the fields written as such: jd,flux,fluxerr,filter.
# filter can be one of the following: ('UVW1','UVW2','UVM2','u_swift','v_swift','b_swift','g_sdss','r_sdss','i_sdss','z_sdss'
# ,'r_cousin','i_cousin','h_2mass','j_2mass','u_johnson','b_johnson','v_johnson')
dates_file='/Users/maayanesoumagnac/PostDoc/projects/2018fif/comparision_with_others/13dqy_int_dates.txt'
#2018fif:'/Users/maayanesoumagnac/PostDoc/projects/2018fif/dates.txt'
lower_limit_on_flux=1e-40
filters_directory='/Users/maayanesoumagnac/Maayane_Astro_python_library/data/filters'

# Interpolation step
already_run_interp_errors=dict() #don't touch this line
already_run_interp_errors['UVW1']=False # set to False if the interpolation for this band has NOT been done yet. Otherwise, set to False to save time.
already_run_interp_errors['UVW2']=True#
already_run_interp_errors['UVM2']=False#
already_run_interp_errors['u_swift']=False # 
already_run_interp_errors['b_swift']=False # 
already_run_interp_errors['v_swift']=False # 
already_run_interp_errors['r_sdss']=True #
already_run_interp_errors['g_sdss']=True #
already_run_interp_errors['i_sdss']=True #
already_run_interp_errors['r_p48']=False # 
already_run_interp_errors['g_p48']=False #
already_run_interp_errors['z_sdss']=True
already_run_interp_errors['u_johnson']=False
already_run_interp_errors['v_johnson']=True
already_run_interp_errors['b_johnson']=False
already_run_interp_errors['i_cousin']=False
already_run_interp_errors['r_cousin']=False
already_run_interp_errors['j_2mass']=False
already_run_interp_errors['h_2mass']=False


# In case you fit with a linear-fitting algorythm
already_run_matrix=False
num_iterations=100

# In case you fit with mcmc
already_run_mcmc=False
nwalkers=100
num_steps=350