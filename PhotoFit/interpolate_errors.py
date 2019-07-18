## NOTE: if the sampling does not work it may be due to zero errors.
import pdb

import numpy as np
import logging
import os
from . import fitter_general
import pylab
from scipy.interpolate import interp1d
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def interpolate_errors(data,x_on_which_to_interpolate,output_path=None,already_run=False,show_plot=True,title=None,verbose=False,band_name='unknown'):
    """Description: given a value, select the line of array with the first colomn closest to the value
    Input  :- N-3 array with the known x position, y positions, errors
            - an M-1 array with the x positions on which to interpolate
            -already_run: if true, the errors have already been ran and stored and the code only reads them from a file.
        Output : a M-5 array with the x positions on which to interpolate, interp1d, best b, lower sigma, higher sigma
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example: if you want to use the error bar closest to a given date MJD:
        error_bar_photo_UV=nearest_values_in_arrays.line_of_closest_first_value(MJD_basis_interpolation_UV,LC_UV)[0,2]
        Reliable:  """
    #print('we will interpolate {0} on {1}'.format(data,x_on_which_to_interpolate)
    #print(data)
    #pdb.set_trace()

    pylab.figure()
    pylab.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], color='red',label='data')
    #print('data[:,0] is',data[:,0])
    #print('x is',x_on_which_to_interpolate)
    #pdb.set_trace()
    for i,j in enumerate(x_on_which_to_interpolate):
        pylab.axvline(j, linestyle='--')
    pylab.axvline(x_on_which_to_interpolate[0],label='interpolation dates')
    if title is not None:
        pylab.title(title)
    pylab.legend()
    pylab.savefig(output_path+'/data_and_interpolation_dates.png',
                  facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
                  bbox_inches=None, pad_inches=0.5)
    pylab.close()
    #pylab.show()
    #print('data is',data)
    Results_array=np.empty((np.shape(x_on_which_to_interpolate)[0],5))
    # calculate the interpolated values with interp1d
    interpol_1d= interp1d(data[:,0], data[:,1])
    #print('test')
    #print(interpol_1d(x_on_which_to_interpolate[0]))
    #print(np.shape(data[:,0]))
    #print(np.shape(data[:, 1]))
    #print(x_on_which_to_interpolate)
   # print('the shape of x is',np.shape(x_on_which_to_interpolate))
    values_with_interp1d=np.zeros(np.shape(x_on_which_to_interpolate))
    for i,j in enumerate(values_with_interp1d):#je sais pas pourquoi, en python3 c est necessaire
        values_with_interp1d[i]=interpol_1d(x_on_which_to_interpolate[i])
    already_run=already_run
    if already_run==False:
        for i,j in enumerate(x_on_which_to_interpolate):
            print('epoch #{0}'.format(i+1))
            #print('i is {0} and j is {1}'.format(i,j))
            #print('I am interpolating on ',x_on_which_to_interpolate)
            #print('the data x I am interpolating is',data[:,0])
            #pdb.set_trace()
            if os.path.exists(output_path+'/error_calc_'+str(i)+'th_point'):
                print(output_path+'/error_calc_'+str(i)+'th_point exists')
            else:
                #print('bla')
                os.mkdir(output_path + '/error_calc_' + str(i) + 'th_point')
            #		logger.info('output_path/txt_files did exist')
            #		#shutil.rmtree('./outputs_from_fit_black_body_flux_filters_function')
            x_lower=np.max(data[data[:,0]<j,0])#le plus grand jour plus petit que j
            x_higher=np.min(data[data[:,0]>j,0])#le plus petit jour plus grans que j
            y_lower=data[data[:,0]==x_lower,1]#resoudre le cas u plusieurs #y quand le jour est =
            y_higher=data[data[:,0]==x_higher,1]#resoudre le cas u plusieurs
            err_lower=data[data[:,0]==x_lower,2]#resoudre le cas u plusieurs
            err_higher=data[data[:,0]==x_higher,2]#resoudre le cas u plusieurs
            #print('the pont on which I am now interpolating is',j)
            #print('x_lower is,',x_lower)
            if verbose==True:
                print('y_lower is',y_lower)
            if np.shape(y_lower)[0] > 1:
                print(
                'warining, there is a repeated y point in the data file, at x point {0}. solve this and continue'.format(
                    x_lower))
                pdb.set_trace()

            #print('err_lower is',err_lower)
            #print('x_higher is,',x_higher)
            #print(np.shape(y_higher))
            if verbose == True:
                print('y_higher is',y_higher)
            if np.shape(y_higher)[0]>1:
                print('warining, there is a repeated y point in the data file, at x point {0}. solve this and continue'.format(x_higher))
                pdb.set_trace()
            #print('err_higher is',err_higher)
            #pdb.set_trace()
            my_data=np.empty((2,3))
            my_data[0,:]=x_lower,y_lower,err_lower
            my_data[1,:]=x_higher,y_higher,err_higher
    #        print(my_data
            '''
            pylab.figure()
            pylab.errorbar(my_data[:,0],my_data[:,1],yerr=my_data[:,2],color='red')
            pylab.axvline(j,linestyle='--')
            '''
            #print(y_higher)
            #print(y_lower)
            #if y_higher==y_lower:
            #    a_ini=(y_higher+0.001-y_lower)/(x_higher-x_lower)
            #else:
            a_ini=(y_higher-y_lower)/(x_higher-x_lower)
            b_ini=(y_higher+y_lower)/2

            #print('b_ini is', b_ini
            #pdb.set_trace()
            #print('a_ini is {0}'.format(a_ini))
            #print('b_ini is {0}'.format(b_ini))
            if a_ini==0.:#CE CAS MARCHE MAL, preferer une alternative
                prior_a=np.array([-0.001,0.001])
                prior_b = np.sort(np.array([y_lower + 0.* y_lower, y_higher - 0.1 * y_higher]))

            else:
                prior_a=np.array([0.1*a_ini,10*a_ini])
                prior_bx=np.sort(np.array([y_lower,y_higher]))
                prior_b = np.array([0.1*np.min(prior_bx),1.9*np.max(prior_bx)])
            #print('the prior on b is {0}'.format(prior_b))
            #fit a my_data a(x-j)+b avec meme:
            class model_linear_x_fixed(object):  # given tref,A(t-tref])**n
                def __init__(self, a, b, t):
                    self.a = a
                    self.b = b
                    self.t = t

                def model_array(self):
                    g = np.zeros((len(self.t), 2))
                    g[:, 0] = self.t[:]
                    g[:, 1] = self.a*(self.t[:]-j)+self.b
                    return g
            #print(my_data[:,0]-2458300)
            #print(my_data[:, 1])
            #pdb.set_trace()
            samples = fitter_general.emcee_n_param(ndim=2, model_nparam=model_linear_x_fixed,
                                                           prior_param=[prior_a, prior_b], data=my_data[:,0:2],
                                                           uncertainties=my_data[:,2], initial_conditions=[a_ini, b_ini],
                                                           flatchain_path=output_path+'/error_calc_'+str(i)+'th_point/flatchain.txt',
                                                           already_run=False,nwalkers=70,num_steps=500,verbose=verbose)

            best = fitter_general.calc_best_fit_n_param(ndim=2, model_nparam=model_linear_x_fixed,
                                                                flatchain_path=output_path+'/error_calc_'+str(i)+'th_point/flatchain.txt',
                                                                data=my_data, uncertainties=my_data[:,2], winners=50,
                                                                output_file_path=output_path+'/error_calc_'+str(i)+'th_point',
                                                                bounds=[prior_a, prior_b], already_run_calc_all_chis=False,
                                                                show_plots=False)

            bests = best[:-1]

            '''
            plots_opt_fit = fitter_general.plot_opt_fit_n_param(ndim=2, model_nparam=model_linear_x_fixed,
                                                                        bests=bests, data=my_data,
                                                                        flatchain_path=output_path+'/error_calc_'+str(i)+'th_point/flatchain.txt',
                                                                        uncertainties=my_data[:,2],
                                                                        output_file_path=output_path+'/error_calc_'+str(i)+'th_point',
                                                                        xlabel='x', ylabel='y')

            triangle = fitter_general.plot_2D_distributions(
                        flatchain_path=output_path+'/error_calc_'+str(i)+'th_point/flatchain.txt', bests=bests,
                        title='test', output_file_path=output_path+'/error_calc_'+str(i)+'th_point', parameters_labels=['a', 'b'])
            '''

            histos = fitter_general.plot_1D_marginalized_distribution(
                        flatchain_path=output_path+'/error_calc_'+str(i)+'th_point/flatchain.txt', bests=bests,
                output_pdf_file_path=output_path+'/error_calc_'+str(i)+'th_point',output_txt_file_path=output_path+'/error_calc_'+str(i)+'th_point', parameters_labels=['a', 'b'], number_bins=2000,
                title='generated for the interpolation, epoch #{0}, band {1}'.format(i,band_name))
            #pylab.show()


            #def plot_1D_marginalized_distribution(flatchain_path, bests=None, output_png_file_path='.',
            #                                      output_txt_file_path='.', parameters_labels=None, number_bins=None):


            # plot le best et l erreur.
            lower_sigma=np.genfromtxt(output_path+'/error_calc_'+str(i)+'th_point/1sigma.txt',skip_header=1)[1,0]
            upper_sigma=np.genfromtxt(output_path+'/error_calc_'+str(i)+'th_point/1sigma.txt',skip_header=1)[1,1]
            Results_array[i,:]=j,values_with_interp1d[i],bests[1],lower_sigma,upper_sigma
        #print('Results_array is',Results_array)
        np.savetxt(output_path+'/Results_array.txt',Results_array,header='days on which we interpolated, interp value with interp1d, best b [from a(x-xint)+b], lower sigma on b, upper sigma on b')


    else:
        #print('I am looking into',output_path+'/Results_array.txt')
        Results_array=np.genfromtxt(output_path+'/Results_array.txt',skip_header=1)
        #print('Result is',Results_array)
    #print('Results_array is',Results_array)

    pylab.figure()
    pylab.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], color='red',label='data')
    for i,j in enumerate(x_on_which_to_interpolate):
        #pylab.axvline(j, linestyle='--')
        #print(np.shape(Results_array))
        #print(Results_array)
        #print(i)
        if Results_array.ndim>1:
            #print(Results_array[i,3])
            #print(Results_array[i,4])
            pylab.plot(j,Results_array[i,1],'bo')#,label='interpolation with interp1d')
            pylab.plot(j,Results_array[i,2],'go')#,label='best b in linear fit')
            #pylab.errorbar(j,Results_array[i.T,color='green')
            pylab.vlines(j,Results_array[i,3],Results_array[i,4],color='green')
        else:
            #print(Results_array[i, 3])
            #print(Results_array[i, 4])
            pylab.plot(j, Results_array[1], 'bo')  # ,label='interpolation with interp1d')
            pylab.plot(j, Results_array[2], 'go')  # ,label='best b in linear fit')
            # pylab.errorbar(j,Results_array[i.T,color='green')
            pylab.vlines(j, Results_array[3], Results_array[4], color='green')

    if Results_array.ndim > 1:
        pylab.plot(x_on_which_to_interpolate[0], Results_array[0, 1], 'bo', label='interpolation with interp1d',alpha=0.5)
        pylab.plot(x_on_which_to_interpolate[0], Results_array[0, 2], 'go', label='best b in mcmc fit',alpha=0.5)
    else:
        pylab.plot(x_on_which_to_interpolate[0], Results_array[1], 'bo', label='interpolation with interp1d',alpha=0.5)
        pylab.plot(x_on_which_to_interpolate[0], Results_array[2], 'go', label='best b in mcmc fit',alpha=0.5)

    pylab.legend()
    pylab.savefig(output_path+'/Plot_w_interpolated_errors.png',
        facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False,
        bbox_inches=None, pad_inches=0.5)
    if title is not None:
        pylab.title(title)
    if show_plot==True:
        pylab.show()
        #pylab.plot()

    return Results_array


    #savetxt
    #return [x,errors]