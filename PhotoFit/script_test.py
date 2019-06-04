from PhotoFit import PhotoFit_fun
import params_test as params
import os
import pylab

########## Definition of the output file ##########
output_file_mcmc=params.output_file_mcmc
output_file_linear=params.output_file_linear
output_file_interpolation=params.output_file_interpolation

if params.mcmc==True:
    output=output_file_mcmc
else:
    output=output_file_linear
if os.path.exists(output):
    print('the output directory, '+output + ' exists already')
else:
    os.mkdir(output)

########################################### Read and process the photometric data ############################################

Best=PhotoFit_fun.calculate_T_and_R_in_time(data_file=params.data_file,dates_file=params.dates_file,already_run_interp_errors_from_param=params.already_run_interp_errors,
                                            already_run_mcmc=params.already_run_mcmc,already_run_matrix=params.already_run_matrix,
                                            num_iterations=params.num_iterations,show_underlying_plots=True,verbose=True,redshift=params.z,
                                            distance_modulus=params.distance_modulus,explosion_date=params.explosion_date,EBV=params.EBV,
                                            output=output,filters_directory=params.filters_directory,mcmc=params.mcmc,
                                            output_file_interpolation=params.output_file_interpolation,
                                            lower_limit_on_flux=params.lower_limit_on_flux,csm=params.csm,already_run_fit=params.already_run_fit,num_steps=params.num_steps,nwalkers=params.nwalkers)
PhotoFit_fun.plot_T_and_R_in_time(Best,data_compare=params.data_compare,compare=False,label_comparision='PTF 13dqy',output=output)
pylab.show()

#PhotoFit_fun.plot_L_in_time(Best,data_file=params.data_file,lower_limit_on_flux=params.lower_limit_on_flux,dates_file=params.dates_file,error_lum_ran=False,explosion_date=params.explosion_date,output=output,mcmc=params.mcmc,output_file_interpolation=params.output_file_interpolation)

PhotoFit_fun.plot_SEDs(Best,already_interpolated=True,data_file=params.data_file,lower_limit_on_flux=params.lower_limit_on_flux,dates_file=params.dates_file,already_run_interp_errors_from_params=params.already_run_interp_errors, number_of_plot=9,redshift=params.z,distance_modulus=params.distance_modulus,explosion_date=params.explosion_date,output=output,filters_directory=params.filters_directory,output_file_interpolation=params.output_file_interpolation,EBV=params.EBV)

