import os 

## Parameters used in script_test.py ##
mcmc=False# if False, plot with linear fit. If True, plots with mcmc
dynesty=True
########## Definition of the output file ##########

output_file_mcmc='./results_fit_sed_mcmc'
output_file_dynesty='./results_fit_sed_dynesty'
output_file_linear='./result_fit_sed_mat'
output_file_interpolation='./results_interpolation'

########## Definition of the object parameters ##########
excluded_bands=['b_swift']
z = 0.005224
distance_pc=None
distance_modulus=30.39#31.8
n=13 #number of epochs
explosion_date=2458975.0 #in principle, this has to be in jd. For the test data it is in jd-explosion dates
#EBV=0.035
EBV=0.02
csm=False
plots=False


data_file='/home/user/path/to/object/Photofit/data_file.txt' #must have at least the following fiels, and a header with the fields written as such: jd,flux,fluxerr,filter.
# filter can be one of the following: ('UVW1','UVW2','UVM2','u_swift','v_swift','b_swift','g_sdss','r_sdss','i_sdss','z_sdss'
# ,'r_cousin','i_cousin','h_2mass','j_2mass','u_johnson','b_johnson','v_johnson')
dates_file='/home/user/path/to/object/interpolation_dates.txt'

lower_limit_on_flux=1e-40
filters_directory='/home/user/path/to/Filters/' #put the path to the Filters directory here


# in case we are running with dynesty:
dyn_kwargs={'maxiter':1e5, #maximum number of bounding interations used in dynesty 
	    'maxcall':5e5, #maxmimum number of calls for log-likelihood function 
	    'nlive':250 #maxmimum number of live particles used in prior sampling
		}   



# In case you fit with a linear-fitting algorythm
already_run_matrix=False
num_iterations=500   #density of grid in log space 
already_run_fit=[False]*n #Bollian list of length # of epochs, specifying which epochs have already run. This is for both linear and mcmc

# In case you fit with mcmc
already_run_mcmc=False 
already_run_dynesty=False
#
already_run_fit=[False]*n
#already_run_fit[0]=False
#already_run_fit[1]=False
#already_run_fit[2]=True
#already_run_fit[3]=False
#already_run_fit[4]=False
#already_run_fit[5]=True
#already_run_fit[6]=True
#already_run_fit[7]=True
#already_run_fit[8]=True
#already_run_fit[9]=True
#already_run_fit[10]=True
#already_run_fit[11]=True
#already_run_fit[12]=True
#already_run_fit[13]=True

    

nwalkers=200 #density of mcmc parameters. see emcee documentation
num_steps=2000 #

# In case you want to compare your R and T results with existing results from a file
data_compare='./test/data_files/Yaron_2017_results.txt'#file with column#1: time from explosion, column#2: temperature (K), column#3:radius (cm)

show_plot=False #Show underlying plots



priors=True
lowrad=[1e13]*13#the lower limit on the radius prior
hirad=[4e15]*13#the upper limit on the radius prior
lowtemp=[5e3]*13
hitemp=[4e4]*13


number_of_plot=9 #number of sed plots can be 9/16

####################### DON'T TOUCH BELOW THIS POINT #######################
if mcmc==True:
    output=output_file_mcmc
elif dynesty==True:
    output=output_file_dynesty
else:
    output=output_file_linear
if os.path.exists(output):
    print('the output directory, '+output + ' exists already')
else:
    os.mkdir(output)