
from PhotoFit import PhotoFit_fun
import params_photofit as params
import os
import pylab

########## Definition of the output file ##########
output_file_mcmc=params.output_file_mcmc
output_file_linear=params.output_file_linear
output_file_dynesty=params.output_file_dynesty

if params.mcmc==True:
    output=output_file_mcmc
elif params.dynesty==True:
    output=output_file_dynesty
else:
    output=output_file_linear
if os.path.exists(output):
    print('the output directory, '+output + ' exists already')
else:
    os.mkdir(output)

########################################### Read and process the photometric data ############################################


Best=PhotoFit_fun.calculate_T_and_R_in_time(data_file=params.data_file,dates_file=params.dates_file,
                                            already_run_mcmc=params.already_run_mcmc,already_run_matrix=params.already_run_matrix,
                                            num_iterations=params.num_iterations,show_underlying_plots=params.plots,verbose=False,redshift=params.z,
                                            distance_modulus=params.distance_modulus,explosion_date=params.explosion_date,EBV=params.EBV,
                                            output=output,filters_directory=params.filters_directory,mcmc=params.mcmc, distance_pc=params.distance_pc,
                                            
                                            lower_limit_on_flux=params.lower_limit_on_flux,csm=params.csm,num_steps=params.num_steps,
                                            nwalkers=params.nwalkers,excluded_bands=params.excluded_bands,already_run_fit=params.already_run_fit,priors=params.priors,lowrad=params.lowrad,hirad=params.hirad,
                                            lowtemp=params.lowtemp,hitemp=params.hitemp,dynesty=params.dynesty,already_run_dynesty=params.already_run_dynesty,**params.dyn_kwargs)

PhotoFit_fun.plot_T_and_R_in_time(Best,compare=False,output=output)

PhotoFit_fun.plot_L_in_time(Best,data_file=params.data_file,lower_limit_on_flux=params.lower_limit_on_flux,dates_file=params.dates_file,
                            error_lum_ran=False,explosion_date=params.explosion_date,output=output,mcmc=params.mcmc,excluded_bands=params.excluded_bands)

PhotoFit_fun.plot_SEDs(Best,data_file=params.data_file,lower_limit_on_flux=params.lower_limit_on_flux,
                       dates_file=params.dates_file,
                       number_of_plot=params.number_of_plot,redshift=params.z,distance_modulus=params.distance_modulus,explosion_date=params.explosion_date,
                       output=output,filters_directory=params.filters_directory,
                       EBV=params.EBV,excluded_bands=params.excluded_bands)


