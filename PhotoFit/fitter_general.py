#! /usr/bin/env python


import numpy as np
from . import class_chi2
import emcee
import pylab
import corner
import scipy.optimize as op
import os


# simple power lax ax^n or a(x-xref)^n
'''
def emcee_1_param(model_1param,prior_param1,data,uncertainties,initial_conditions,flatchain_path=None,already_run=False):
    """Description: given data, uncertainties, and priors on a and b, run emcee to fit the data with the model a(x-xref)^n (0 for x smaller than xref)
    Input  :-
    Output :-
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example:     Reliable:  """
    if initial_conditions[0] < np.min(prior_param1) or initial_conditions[0]>np.max(prior_param1):
        print('ERROR: you need to give initial condition for the one parameter which are within the prior range on a'
        exit()

    def lnprior(theta):
        param1 = theta
        if np.min(prior_param1) < param1 < np.max(prior_param1) :
            return 0.0
        return -np.inf

    def myloglike(theta):
        param1 = theta
        my_model = model_1param(param1, data[:,0]).model_array()
        my_objective = class_chi2.objective_with_uncertainties(my_model, data,
                                                                   sigmas=uncertainties)
        loglik = my_objective.chi_square_value()
        return -0.5 * loglik

    def lnprob(theta):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + myloglike(theta)

    if already_run!=True:
        print('*** EMCEE run ***'
        ndim=1
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
        pos = np.zeros((nwalkers, ndim))
        print(pos)
        pos[:, 0] = [initial_conditions[0] + min(np.average(prior_param1),100.) * 1e-4*np.random.randn() for j in range(nwalkers)]

        #print(pos)
        #print('the min of pos is', np.min(pos[:, 2])
        #print('the min of pos is', np.max(pos[:, 2])
        #print('prior is', prior_xref
        if np.min(pos[:, 0]) <= np.min(prior_param1) or np.max(pos[:, 0]) >= np.max(prior_param1):
            print('ERROR! the pos is beyond the prior for the first param'
            exit()
        #pdb.set_trace()
        #pos = [[2.1,2.2]+ 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        pos, prob, state = sampler.run_mcmc(pos, 100)
        sampler.reset()
        sampler.run_mcmc(pos,1000)
        if flatchain_path==None:
            np.savetxt('./flatchain_fitter_powerlaw.txt',sampler.flatchain)
        else:
            np.savetxt(flatchain_path, sampler.flatchain)
        print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))
#print("Autocorrelation time:", sampler.get_autocorr_time())
        samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    else:
        if flatchain_path==None:
            samples=np.genfromtxt('./flatchain_fitter_powerlaw.txt')
        else:
            samples = np.genfromtxt(flatchain_path)
    return samples
'''
def emcee_n_param(ndim,model_nparam, prior_param, data, uncertainties, initial_conditions, nwalkers=100,num_steps=200,flatchain_path=None,already_run=False,verbose=False):
    """Description: given data, uncertainties, and priors on a and b, run emcee to fit the data with the model a(x-xref)^n (0 for x smaller than xref)
    Input  :-model_nparam a class like in the models.py file, that has a model_array method
            -a tuple where every element is a 2-1 array with the prior
            -num_steps goes into sampler.run_mcmc(pos,num_steps), default is 200.
    Output :-
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example: samples=fitter_general.emcee_n_param(ndim=2,model_nparam=model_lo_to,prior_param=[prior_alpha,prior_tref],data=my_data,uncertainties=errors,initial_conditions=[alpha_true,tref_true],flatchain_path='output_test_fitter_general/flatchain_test.txt',already_run=False)
    Reliable:  """
    if np.shape(initial_conditions)[0]!=ndim:
        #print np.shape(initial_conditions)[0]
        print('you did not give initial conditions for every parameter')
        exit()
    for i, j in enumerate(initial_conditions):
        if j < np.min(prior_param[i]) or j > np.max(prior_param[i]):
            print('ERROR: you need to give initial condition for parameter number {0} which are within the prior range'.format(i+1))
            exit()

    def lnprior(theta):
        if ndim==1:
            param0 = theta
            if np.min(prior_param[0]) < param0 < np.max(prior_param[0]):
                return 0.0
        elif ndim==2:
            param0,param1=theta
            if np.min(prior_param[0]) < param0 < np.max(prior_param[0]) and np.min(prior_param[1]) < param1 < np.max(prior_param[1]):
                return 0.0
        elif ndim==3:
            param0,param1,param2=theta
            if np.min(prior_param[0]) < param0 < np.max(prior_param[0]) and np.min(prior_param[1]) < param1 < np.max(prior_param[1]) and np.min(prior_param[2]) < param2 < np.max(prior_param[2]):
                return 0.0
        elif ndim==4:
            param0,param1,param2,param3=theta
            if np.min(prior_param[0]) < param0 < np.max(prior_param[0]) and np.min(prior_param[1]) < param1 < np.max(prior_param[1]) and np.min(prior_param[2]) < param2 < np.max(prior_param[2]) and np.min(prior_param[3]) < param3 < np.max(prior_param[3]):
                return 0.0
        elif ndim==5:
            param0, param1, param2, param3,param4 = theta
            if np.min(prior_param[0]) < param0 < np.max(prior_param[0]) and np.min(prior_param[1]) < param1 < np.max(
                    prior_param[1]) and np.min(prior_param[2]) < param2 < np.max(prior_param[2]) and np.min(
                    prior_param[3]) < param3 < np.max(prior_param[3]) and np.min(
                    prior_param[4]) < param4 < np.max(prior_param[4]):
                return 0.0
        elif ndim==6:
            param0, param1, param2, param3,param4,param5 = theta
            if np.min(prior_param[0]) < param0 < np.max(prior_param[0]) and np.min(prior_param[1]) < param1 < np.max(
                    prior_param[1]) and np.min(prior_param[2]) < param2 < np.max(prior_param[2]) and np.min(
                    prior_param[3]) < param3 < np.max(prior_param[3]) and np.min(
                    prior_param[4]) < param4 < np.max(prior_param[4]) and np.min(
                    prior_param[5]) < param5 < np.max(prior_param[5]):
                return 0.0
        elif ndim==7:
            param0, param1, param2, param3,param4,param5,param6 = theta
            if np.min(prior_param[0]) < param0 < np.max(prior_param[0]) and np.min(prior_param[1]) < param1 < np.max(
                    prior_param[1]) and np.min(prior_param[2]) < param2 < np.max(prior_param[2]) and np.min(
                    prior_param[3]) < param3 < np.max(prior_param[3]) and np.min(
                    prior_param[4]) < param4 < np.max(prior_param[4]) and np.min(
                    prior_param[5]) < param5 < np.max(prior_param[5]) and np.min(
                    prior_param[6]) < param6 < np.max(prior_param[6]):
                return 0.0
        return -np.inf

    def myloglike(theta):
        if ndim==1:
            param1 = theta
            my_model = model_nparam(param1, data[:, 0]).model_array()
        elif ndim==2:
            param1,param2 = theta
            #print param1
            #print param2
            my_model = model_nparam(param1,param2, data[:, 0]).model_array()
        elif ndim==3:
            param1,param2,param3 = theta
            my_model = model_nparam(param1,param2,param3, data[:, 0]).model_array()
        elif ndim==4:
            param1,param2,param3,param4 = theta
            my_model = model_nparam(param1,param2,param3, param4,data[:, 0]).model_array()
        elif ndim==5:
            param1,param2,param3,param4,param5 = theta
            my_model = model_nparam(param1,param2,param3, param4,param5,data[:, 0]).model_array()
        elif ndim==6:
            param1,param2,param3,param4,param5,param6 = theta
            my_model = model_nparam(param1,param2,param3, param4,param5,param6,data[:, 0]).model_array()
        elif ndim == 7:
            param1, param2, param3, param4, param5, param6, param7 = theta
            my_model = model_nparam(param1, param2, param3, param4, param5, param6,param7,data[:, 0]).model_array()

        my_objective = class_chi2.objective_with_uncertainties(my_model, data,
                                                               sigmas=uncertainties)
        loglik = my_objective.chi_square_value()
        return -0.5 * loglik

    def lnprob(theta):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + myloglike(theta)

    if already_run != True:
        if verbose==True:
            print('*** EMCEE run ***')
        #ndim = 1
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
        pos = np.zeros((nwalkers, ndim))
        for i in range(ndim):
            #print min(np.average(prior_param[i]), 100.)
            #print initial_conditions[i]
            #print np.average(prior_param[i])
            #pos[:, i] = [initial_conditions[i] + (np.max(prior_param[i])-np.average(prior_param[i])) * 1e-1 * np.random.randn() for j in range(nwalkers)]
            pos[:, i] = [
                np.average(prior_param[i]) + (np.max(prior_param[i]) - np.min(prior_param[i])) * 1e-1 * np.random.randn()
                for j in range(nwalkers)]
            if np.min(pos[:, i]) <= np.min(prior_param[i]) or np.max(pos[:, i]) >= np.max(prior_param[i]):
                #print(pos)
                print('initial condition is',initial_conditions[i])
                print(' 0.5*(np.max(prior_param[i]) - np.min(prior_param[i])) * 1e-3  is', (np.max(prior_param[i]) - np.min(prior_param[i])) * 1e-3)
                print('the min of pos is',np.min(pos[:, i]))
                print('the max of pos is',np.max(pos[:, i]))
                print('the prior is',prior_param[i])
                print('ERROR! the pos is beyond the prior for the '+str(i+1)+'th param')
                exit()
        #print(pos)
        for i in range(ndim):
            if verbose ==True:
                print('the min of pos is', np.min(pos[:, i]))
                print('the max of pos is', np.max(pos[:, i]))
                print('the prior is', prior_param[i])
        #pos, prob, state = sampler.run_mcmc(pos, 100)
        #sampler.reset()
        sampler.run_mcmc(pos, num_steps)
        if flatchain_path == None:
            np.savetxt('./flatchain_fitter_general.txt', sampler.flatchain)
        else:
            np.savetxt(flatchain_path, sampler.flatchain)
        print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))
        # print("Autocorrelation time:", sampler.get_autocorr_time())
        samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    else:
        if flatchain_path == None:
            samples = np.genfromtxt('./flatchain_fitter_general.txt')
        else:
            samples = np.genfromtxt(flatchain_path)
    return samples

def calc_best_fit_n_param(ndim,model_nparam,flatchain_path,data,uncertainties,winners=50,output_file_path=None,bounds=None,already_run_calc_all_chis=False,show_plots=False,dilution_factor=0,verbose=False):
    """Description:
        -Input  :
        -Output : a numpy array with [best params, chi2]
        Plots and txt files:
        txt files:
       pdf files:
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example: best=fitter_general.calc_best_fit_n_param()
        Reliable:  """
    if verbose==True:
        print('*** Calculation of the maximum likelihood values ***')
    samples=np.genfromtxt(flatchain_path,delimiter=None,dtype=float)
    if output_file_path==None:
        output_file_path='.'
    file_path_to_all_chain_chis = output_file_path + '/all_param_and_chis_of_chain.txt'
    #### 1. calculate the chi2 of all the parameters combinations of the chain, store this in file_path_to_all_chain_chis/all_param_and_chis_of_chain.txt'
    if already_run_calc_all_chis != True:
        if dilution_factor==0:
            if verbose == True:
                print('I am calculating the chi value for each combination of parameters in the chain')
            chis = np.zeros((np.shape(samples)[0], ndim+1))
            for i, line in enumerate(samples):
                # print('i am a line {0}'.format(i)
                if ndim==1:
                    my_objective = class_chi2.objective_with_uncertainties(model_nparam(line[0],data[:, 0]).model_array(), data, sigmas=uncertainties)
                elif ndim==2:
                    my_objective = class_chi2.objective_with_uncertainties(model_nparam(line[0], line[1],data[:, 0]).model_array(),
                                                                           data, sigmas=uncertainties)
                elif ndim==3:
                    my_objective = class_chi2.objective_with_uncertainties(model_nparam(line[0],line[1],line[2], data[:, 0]).model_array(),
                                                                           data, sigmas=uncertainties)
                elif ndim==4:
                    my_objective=class_chi2.objective_with_uncertainties(model_nparam(line[0],line[1],line[2], line[3],data[:, 0]).model_array(),
                                                                           data, sigmas=uncertainties)
                elif ndim==5:
                    my_objective=class_chi2.objective_with_uncertainties(model_nparam(line[0],line[1],line[2], line[3],line[4],data[:, 0]).model_array(),
                                                                           data, sigmas=uncertainties)
                elif ndim==6:
                    my_objective=class_chi2.objective_with_uncertainties(model_nparam(line[0],line[1],line[2], line[3],line[4],line[5],data[:, 0]).model_array(),
                                                                           data, sigmas=uncertainties)
                elif ndim==7:
                    my_objective=class_chi2.objective_with_uncertainties(model_nparam(line[0],line[1],line[2], line[3],line[4],line[5],line[6],data[:, 0]).model_array(),
                                                                           data, sigmas=uncertainties)
                chis[i, 0:ndim] = samples[i, 0:ndim]
                chis[i, ndim] = my_objective.chi_square_value()
        else:
            if verbose == True:
                print('I am calculating the chi value for every {0} combination of parameters in the chain'.format(dilution_factor))
            samples_diluted=samples[::dilution_factor].copy()
            chis = np.zeros((np.shape(samples_diluted)[0], ndim + 1))
            for i, line in enumerate(samples_diluted):
                # print('i am a line {0}'.format(i)
                if ndim == 1:
                    my_objective = class_chi2.objective_with_uncertainties(
                        model_nparam(line[0], data[:, 0]).model_array(), data, sigmas=uncertainties)
                elif ndim == 2:
                    my_objective = class_chi2.objective_with_uncertainties(
                        model_nparam(line[0], line[1], data[:, 0]).model_array(),
                        data, sigmas=uncertainties)
                elif ndim == 3:
                    my_objective = class_chi2.objective_with_uncertainties(
                        model_nparam(line[0], line[1], line[2], data[:, 0]).model_array(),
                        data, sigmas=uncertainties)
                elif ndim == 4:
                    my_objective = class_chi2.objective_with_uncertainties(
                        model_nparam(line[0], line[1], line[2], line[3], data[:, 0]).model_array(),
                        data, sigmas=uncertainties)
                elif ndim == 5:
                    my_objective = class_chi2.objective_with_uncertainties(
                        model_nparam(line[0], line[1], line[2], line[3], line[4], data[:, 0]).model_array(),
                        data, sigmas=uncertainties)
                elif ndim == 6:
                    my_objective = class_chi2.objective_with_uncertainties(
                        model_nparam(line[0], line[1], line[2], line[3], line[4], line[5], data[:, 0]).model_array(),
                        data, sigmas=uncertainties)
                elif ndim == 7:
                    my_objective = class_chi2.objective_with_uncertainties(
                        model_nparam(line[0], line[1], line[2], line[3], line[4], line[5], line[6],
                                     data[:, 0]).model_array(),
                        data, sigmas=uncertainties)
                chis[i, 0:ndim] = samples[i, 0:ndim]
                chis[i, ndim] = my_objective.chi_square_value()
        np.savetxt(file_path_to_all_chain_chis, chis)
    else:
        chis = np.genfromtxt(file_path_to_all_chain_chis)
    #### 2. sort all the chis by increasing order
    param_and_chi_sorted = chis[np.argsort(chis[:, ndim])]
    #### 3. take the n=winners combinations with lowest chi2, and use them as initial point of the op.minimize algorythm
    result_op_min = np.zeros((winners, np.shape(param_and_chi_sorted)[1]))

    def lnlike(theta, data, uncert):
        if ndim==1:
            x= theta
            return -0.5 * class_chi2.objective_with_uncertainties(model_nparam(x,data[:, 0]).model_array(),
                                                              data, sigmas=uncert).chi_square_value()
        elif ndim==2:
            x,y= theta
            return -0.5 * class_chi2.objective_with_uncertainties(model_nparam(x,y,data[:, 0]).model_array(),
                                                              data, sigmas=uncert).chi_square_value()
        elif ndim==3:
            x,y,z = theta
            return -0.5 * class_chi2.objective_with_uncertainties(model_nparam(x,y,z, data[:, 0]).model_array(),
                                                              data, sigmas=uncert).chi_square_value()
        elif ndim==4:
            x,y,z,w = theta
            return -0.5 * class_chi2.objective_with_uncertainties(model_nparam(x,y,z,w, data[:, 0]).model_array(),
                                                              data, sigmas=uncert).chi_square_value()
        elif ndim==5:
            x,y,z,w,u = theta
            return -0.5 * class_chi2.objective_with_uncertainties(model_nparam(x,y,z,w,u, data[:, 0]).model_array(),
                                                              data, sigmas=uncert).chi_square_value()
        elif ndim==6:
            x,y,z,w,u,r = theta
            return -0.5 * class_chi2.objective_with_uncertainties(model_nparam(x,y,z,w,u,r, data[:, 0]).model_array(),
                                                              data, sigmas=uncert).chi_square_value()
        elif ndim==7:
            x,y,z,w,u,r,s = theta
            return -0.5 * class_chi2.objective_with_uncertainties(model_nparam(x,y,z,w,u,r,s, data[:, 0]).model_array(),
                                                              data, sigmas=uncert).chi_square_value()

    nll = lambda *args: -lnlike(*args) / 0.5
    bounds=[[np.min(bounds[i]),np.max(bounds[i])] for i in range(ndim)]
    #print bounds
    for i in range(winners):
        result_op_mini = op.minimize(nll, [param_and_chi_sorted[i, j] for j in range(ndim)],args=(data, uncertainties), bounds=bounds)
        besties = result_op_mini["x"]
        #print besties
        result_op_min[i, 0:np.shape(param_and_chi_sorted)[1] - 1] = besties
        if ndim==1:
            result_op_min[i, np.shape(param_and_chi_sorted)[1] - 1] = class_chi2.objective_with_uncertainties(
            model_nparam(besties[0],data[:, 0]).model_array(), data,
            sigmas=uncertainties).chi_square_value()
        elif ndim==2:
            result_op_min[i, np.shape(param_and_chi_sorted)[1] - 1] = class_chi2.objective_with_uncertainties(
            model_nparam(besties[0],besties[1],data[:, 0]).model_array(), data,
            sigmas=uncertainties).chi_square_value()
        elif ndim==3:
            result_op_min[i, np.shape(param_and_chi_sorted)[1] - 1] = class_chi2.objective_with_uncertainties(
            model_nparam(besties[0],besties[1],besties[2],data[:, 0]).model_array(), data,
            sigmas=uncertainties).chi_square_value()
        elif ndim==4:
            result_op_min[i, np.shape(param_and_chi_sorted)[1] - 1] = class_chi2.objective_with_uncertainties(
            model_nparam(besties[0],besties[1],besties[2],besties[3],data[:, 0]).model_array(), data,
            sigmas=uncertainties).chi_square_value()
        elif ndim==5:
            result_op_min[i, np.shape(param_and_chi_sorted)[1] - 1] = class_chi2.objective_with_uncertainties(
            model_nparam(besties[0],besties[1],besties[2],besties[3],besties[4],data[:, 0]).model_array(), data,
            sigmas=uncertainties).chi_square_value()
        elif ndim==6:
            result_op_min[i, np.shape(param_and_chi_sorted)[1] - 1] = class_chi2.objective_with_uncertainties(
            model_nparam(besties[0],besties[1],besties[2],besties[3],besties[4],besties[5],data[:, 0]).model_array(), data,
            sigmas=uncertainties).chi_square_value()
        elif ndim == 7:
            result_op_min[i, np.shape(param_and_chi_sorted)[1] - 1] = class_chi2.objective_with_uncertainties(
                model_nparam(besties[0], besties[1], besties[2], besties[3], besties[4], besties[5],besties[6],
                             data[:, 0]).model_array(), data,
                sigmas=uncertainties).chi_square_value()

    #### 4. sort and store the obtained best fit parameters and corresponding chi in file_path_to_all_chain_chis/winners_best_chis.txt
    best_param_and_chi_sorted = result_op_min[np.argsort(result_op_min[:, ndim])]
    # print best_param_and_chi_sorted
    np.savetxt(output_file_path + '/' + str(winners) + '_best_combination_and_chis.txt', best_param_and_chi_sorted)

    #### 5. find the three optimized best fit combinations with lowest chi2 and store them in best_com
    maxparam = np.zeros((np.shape(best_param_and_chi_sorted)))
    # print np.shape(maxparam)
    # pdb.set_trace()
    j = 0
    maxparam[0, :] = best_param_and_chi_sorted[0, 0:np.shape(best_param_and_chi_sorted)[1]]
    for i in range(winners - 1):
        # print('param number ', i
        # print best_param_and_chi_sorted[i, 0:np.shape(best_param_and_chi_sorted)[1] - 1]  # print the a and b of the line
        a = (
        best_param_and_chi_sorted[i + 1, 0:np.shape(best_param_and_chi_sorted)[1] - 1] != best_param_and_chi_sorted[
                                                                                          i,
                                                                                          0:np.shape(
                                                                                              best_param_and_chi_sorted)[
                                                                                                1] - 1])
        # print a
        # print a.all
        if a.all():  # if the param are not the same
            # print('the line {0} is not identical to the line {1}'.format(i, i + 1)
            j = j + 1
            maxparam[j, :] = best_param_and_chi_sorted[i + 1,
                             0:np.shape(best_param_and_chi_sorted)[1]]  # we add the line to the column
            # else:
            #    print('the line {0} is identical to the line {1}'.format(i, i + 1
    maxiparam_good_size = maxparam

    np.savetxt(output_file_path + '/best_optimized_combinations.txt', maxiparam_good_size,
               header='best parameters,chi^2')
    # np.savetxt(output_file_path+'/first_best_combintation.txt', maxiparam_good_size[0, 0:np.shape(param_and_chi_sorted)[1]],header='best a, best n, chi^2')
    # if np.shape(maxiparam_good_size)[1]>1:
    #    np.savetxt(output_file_path + '/second_best_combintation.txt', maxparam[1, 0:np.shape(param_and_chi_sorted)[1]],
    #           header='#best a, best n, chi^2')
    # if np.shape(maxiparam_good_size)[1]>2:
    #    np.savetxt(output_file_path + '/third_best_combintation.txt', maxparam[2, 0:np.shape(param_and_chi_sorted)[1]],
    #           header='#best a, best n, chi^2')
    if show_plots == True:
        pylab.figure()
        pylab.plot(param_and_chi_sorted[:, ndim], 'b', label=r'$\chi^2 of the combinations in the chain$')
        pylab.title(r'sorted $\chi^2$ values in the chain {0}'.format(flatchain_path))
        pylab.xlabel('index')
        pylab.ylabel('$\chi^2$')
        pylab.savefig(output_file_path + '/all_chis_plot.pdf', facecolor='w', edgecolor='w', orientation='portrait',
                      papertype=None, format='pdf',
                      transparent=False, bbox_inches=None, pad_inches=0.1)
        pylab.figure()
        pylab.plot(best_param_and_chi_sorted[:, ndim], 'r')
        pylab.title(r'sorted optimized $\chi^2$ values')
        pylab.xlabel('index')
        pylab.ylabel('$\chi^2$')
        pylab.savefig(output_file_path + '/opt_chis_plot.pdf', facecolor='w', edgecolor='w', orientation='portrait',
                      papertype=None, format='pdf',
                      transparent=False, bbox_inches=None, pad_inches=0.1)
    #print(maxiparam_good_size)
    if verbose==True:
        print('the best fit param are {0} and the corrisponding chi^2 is {1}'.format(maxiparam_good_size[0,:ndim],maxiparam_good_size[0, ndim]))
        print('the reduced chi^2 is {0}'.format(round(maxiparam_good_size[0, ndim] / (np.shape(data)[0] - ndim), 4)))
    #print('the denominator is',(np.shape(data)[0] - ndim)
    #print ndim
    #print np.shape(data)[0]
    #pdb.set_trace()
    return maxiparam_good_size[0, :]
    #pylab.show()

def calc_best_fit_n_param_fast(ndim,flatchain_path,output_file_path=None):
    """Description:
        -Input  :
        -Output : a numpy array with [best params, chi2]
        Plots and txt files:
        txt files:
       pdf files:
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example: best=fitter_general.calc_best_fit_n_param()
        Reliable:  """
    print('*** Calculation of the maximum likelihood values ***')
    samples=np.genfromtxt(flatchain_path,delimiter=None,dtype=float)
    if output_file_path==None:
        output_file_path='.'
    bests=np.zeros(np.shape(samples)[1])
    file_path_to_all_chain_chis = output_file_path + '/all_param_and_chis_of_chain.txt'
    #### 1. calculate the chi2 of all the parameters combinations of the chain, store this in file_path_to_all_chain_chis/all_param_and_chis_of_chain.txt'
    for i in range(ndim):
        bests[i]=np.median(samples[:,i])
    np.savetxt(output_file_path + '/' +'best_from_median.txt',bests)
    return bests

def plot_opt_fit_n_param(ndim, model_nparam,bests,data,flatchain_path=None,uncertainties=None, output_file_path='.',xlabel=None,ylabel=None):
    """Description:
        Input  :
        Output : 
        Plots and txt files:
            - 
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example:
        Reliable:  """

    print('*** Plots of the best fit solution ***')
    if ndim==1:
        best_fit=model_nparam(bests[0],data[:,0]).model_array()
        emcee_fit = model_nparam(bests[0],np.linspace(np.min(data[:, 0]), np.max(data[:, 0]),
                                                              1000)).model_array()
    elif ndim==2:
        best_fit=model_nparam(bests[0],bests[1],data[:,0]).model_array()
        #print np.linspace(np.min(data[:, 0]), np.max(data[:, 0]),1000)
        emcee_fit = model_nparam(bests[0], bests[1],np.linspace(np.min(data[:, 0]), np.max(data[:, 0]),
                                                              1000)).model_array()
    elif ndim==3:
        best_fit=model_nparam(bests[0],bests[1],bests[2],data[:,0]).model_array()
        emcee_fit = model_nparam(bests[0], bests[1], bests[2],np.linspace(np.min(data[:, 0]), np.max(data[:, 0]),
                                                              1000)).model_array()
    elif ndim==4:
        best_fit=model_nparam(bests[0],bests[1],bests[2],bests[3],data[:,0]).model_array()
        emcee_fit = model_nparam(bests[0], bests[1], bests[2],bests[3],np.linspace(np.min(data[:, 0]), np.max(data[:, 0]),
                                                              1000)).model_array()
    elif ndim==5:
        best_fit=model_nparam(bests[0],bests[1],bests[2],bests[3],bests[4],data[:,0]).model_array()
        emcee_fit = model_nparam(bests[0], bests[1], bests[2],bests[3],bests[4],np.linspace(np.min(data[:, 0]), np.max(data[:, 0]),
                                                              1000)).model_array()
    elif ndim==6:
        best_fit=model_nparam(bests[0],bests[1],bests[2],bests[3],bests[4],bests[5],data[:,0]).model_array()
        emcee_fit = model_nparam(bests[0], bests[1], bests[2],bests[3],bests[4],bests[5],np.linspace(np.min(data[:, 0]), np.max(data[:, 0]),
                                                              1000)).model_array()
    elif ndim==7:
        best_fit=model_nparam(bests[0],bests[1],bests[2],bests[3],bests[4],bests[5],bests[6],data[:,0]).model_array()
        emcee_fit = model_nparam(bests[0], bests[1], bests[2],bests[3],bests[4],bests[5],bests[6],np.linspace(np.min(data[:, 0]), np.max(data[:, 0]),
                                                              1000)).model_array()
    elif ndim==8:
        best_fit=model_nparam(bests[0],bests[1],bests[2],bests[3],bests[4],bests[5],bests[6],bests[7],data[:,0]).model_array()
        emcee_fit = model_nparam(bests[0], bests[1], bests[2],bests[3],bests[4],bests[5],bests[6],bests[7],np.linspace(np.min(data[:, 0]), np.max(data[:, 0]),
                                                              1000)).model_array()
    elif ndim==9:
        best_fit=model_nparam(bests[0],bests[1],bests[2],bests[3],bests[4],bests[5],bests[6],bests[7],bests[8],data[:,0]).model_array()
        emcee_fit = model_nparam(bests[0], bests[1], bests[2],bests[3],bests[4],bests[5],bests[6],bests[7],bests[8],np.linspace(np.min(data[:, 0]), np.max(data[:, 0]),
                                                              1000)).model_array()
    fig = pylab.figure()
    pylab.plot(best_fit[:,0],best_fit[:,1], 'bo', label=r'best-fit model')
    #pylab.plot(emcee_fit[:,0],emcee_fit[:,1],'b-')
    if uncertainties is None==False:
        pylab.errorbar(data[:,0],data[:,1],yerr=0.5 * uncertainties, fmt='r.', label=r'data')
    else:
        pylab.plot(data[:, 0], data[:, 1],'r.', label=r'data')
    #pylab.title(r'Data and best-fit model, $\chi^2/dof={0}$'.format(round(class_chi2.objective_with_uncertainties(best_fit,data,uncertainties).chi_square_value() / (np.shape(data)[0] - ndim), 3)))
    ax = pylab.gca()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    pylab.legend(loc=1)
    my_objective=class_chi2.objective_with_uncertainties(best_fit,data, sigmas=uncertainties)
    chi=my_objective.chi_square_value()
    dof=np.shape(data)[0]-ndim
    pylab.title(r'Best fit, $\chi=$'+str(chi)+', $\chi/dof=$'+str(chi/dof))
    pylab.savefig(output_file_path + '/best_fit.pdf', facecolor='w', edgecolor='w',
                  orientation='portrait', papertype=None, format='pdf', transparent=False, bbox_inches=None,
                  pad_inches=0.1)

    if flatchain_path!=None:
        samples = np.genfromtxt(flatchain_path, comments='#')
        fig = pylab.figure()
        if ndim==1:
            for a in samples[np.random.randint(len(samples), size=15)]:
                emcee_fit_indata=model_nparam(a,np.linspace(np.min(data[:,0]),np.max(data[:,0]),1000))
                pylab.plot(emcee_fit_indata.model_array()[:, 0], emcee_fit_indata.model_array()[:, 1], 'b-', alpha=0.5)
        elif ndim==2:
            for a, b, in samples[np.random.randint(len(samples), size=15)]:
                emcee_fit_indata = model_nparam(a, b, data[:,0])
                emcee_fit_indata = model_nparam(a, b,np.linspace(np.min(data[:, 0]), np.max(data[:, 0]),
                                                                                 1000))
                pylab.plot(emcee_fit_indata.model_array()[:, 0], emcee_fit_indata.model_array()[:, 1], 'b-', alpha=0.5)
        elif ndim==3:
            for a,b,c in samples[np.random.randint(len(samples), size=15)]:
                emcee_fit_indata=model_nparam(a,b,c,np.linspace(np.min(data[:,0]),np.max(data[:,0]),1000))
                pylab.plot(emcee_fit_indata.model_array()[:, 0], emcee_fit_indata.model_array()[:, 1], 'b-', alpha=0.5)
        elif ndim==4:
            for a,b,c,d in samples[np.random.randint(len(samples), size=15)]:
                emcee_fit_indata=model_nparam(a,b,c,d,np.linspace(np.min(data[:,0]),np.max(data[:,0]),1000))
                pylab.plot(emcee_fit_indata.model_array()[:, 0], emcee_fit_indata.model_array()[:, 1], 'b-', alpha=0.5)
        elif ndim==5:
            for a,b,c,d,e in samples[np.random.randint(len(samples), size=15)]:
                emcee_fit_indata=model_nparam(a,b,c,d,e,np.linspace(np.min(data[:,0]),np.max(data[:,0]),1000))
                pylab.plot(emcee_fit_indata.model_array()[:, 0], emcee_fit_indata.model_array()[:, 1], 'b-', alpha=0.5)
        elif ndim == 6:
            for a, b, c, d, e,f in samples[np.random.randint(len(samples), size=15)]:
                emcee_fit_indata = model_nparam(a, b, c, d, e,f,
                                                np.linspace(np.min(data[:, 0]), np.max(data[:, 0]), 1000))
                pylab.plot(emcee_fit_indata.model_array()[:, 0], emcee_fit_indata.model_array()[:, 1], 'b-', alpha=0.5)
        elif ndim == 7:
            for a, b, c, d, e,f,g in samples[np.random.randint(len(samples), size=15)]:
                emcee_fit_indata = model_nparam(a, b, c, d, e,f,g,
                                                np.linspace(np.min(data[:, 0]), np.max(data[:, 0]), 1000))
                pylab.plot(emcee_fit_indata.model_array()[:, 0], emcee_fit_indata.model_array()[:, 1], 'b-', alpha=0.5)
        elif ndim == 8:
            for a, b, c, d, e,f,g,h in samples[np.random.randint(len(samples), size=15)]:
                emcee_fit_indata = model_nparam(a, b, c, d, e,f,g,h,
                                                np.linspace(np.min(data[:, 0]), np.max(data[:, 0]), 1000))
                pylab.plot(emcee_fit_indata.model_array()[:, 0], emcee_fit_indata.model_array()[:, 1], 'b-', alpha=0.5)
        elif ndim == 9:
            for a, b, c, d, e,f,g,h,k in samples[np.random.randint(len(samples), size=15)]:
                emcee_fit_indata = model_nparam(a, b, c, d, e,f,g,h,k,
                                                np.linspace(np.min(data[:, 0]), np.max(data[:, 0]), 1000))
                pylab.plot(emcee_fit_indata.model_array()[:, 0], emcee_fit_indata.model_array()[:, 1], 'b-', alpha=0.5)

        pylab.plot(best_fit[:, 0], best_fit[:, 1], 'ro', label=r'Best-fit model')
        #pylab.plot(emcee_fit[:, 0], emcee_fit[:, 1], 'r-')

        #uncomment
        '''
        if uncertainties != 'None':
            pylab.errorbar(data[:, 0], data[:, 1], yerr=0.5 * uncertainties, fmt='r.', label=r'data')
        else:
            pylab.plot(data[:, 0], data[:, 1], fmt='r.', label=r'data')
        '''
        #pylab.title(r'Data and best-fit model, $\chi^2/dof={0}$'.format(round(
        #    class_chi2.objective_with_uncertainties(best_fit, data, uncertainties).chi_square_value() / (
        #    np.shape(data)[0] - ndim), 3)))
        ax = pylab.gca()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        pylab.legend(loc=1)
        pylab.savefig(output_file_path + '/best_fit_with_chain.pdf', facecolor='w', edgecolor='w',
                  orientation='portrait', papertype=None, format='pdf', transparent=False, bbox_inches=None,
                  pad_inches=0.1)
    #pylab.show()

def plot_2D_distributions(flatchain_path,bests=None,output_file_path='.',parameters_labels=None,title=None,verbose=False):
    """Description: Given a chain, plots the 2-D marginalized distributions.
    Input  :-flatchain_Path: the path to the chain.
            -bests (optional):a numpy array with values of the parameters to be shown on the plots (ex: calculated best-fit values, true known values, etc)
            -output_file_path (optional): the path to the outputs folders. Default is '.'
            -parameters_labels (optional): a tuple with the names of the parameters to be added to the plot. Default is ['param 1','param 2']
    Output :- no output. plots.
    pdf files:
             - 2D-distributions.png : a corner plot with the 2D marginalized distributions.
    Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
    Example: samples=fitter_powerlaw.fitter_powerlaw(np.array([0.,5.]),np.array([0.,5.]),my_data,uncertainties=errors*np.ones(50),flatchain_path='./flatchain_2_test.txt',already_run=True)
    Reliable:  """
    if verbose == True:
        print('*** Plots of the 2D distributions ***')
    samples=np.genfromtxt(flatchain_path)
    #if output_file_path==None:
    #    output_file_path='.'
    #if parameters_labels==None:
    #    labels=['param 1','param 2']
    #else:
    #    labels=parameters_labels
        #print labels
    if title!=None:
        big_title=title
    if bests is None:
        #fig = triangle.corner(samples, labels=parameters_labels,
         #                quantiles=[0.159, 0.5, 0.841],show_titles=False,big_title=big_title)
        fig = corner.corner(samples, labels=parameters_labels,
                         quantiles=[0.159, 0.5, 0.841],show_titles=False,big_title=big_title)
        for ax in fig.get_axes():
        #    ax.tick_params(axis='both', labelsize=14)
            #ax.set_xlabel(fontsize=14)
            ax.xaxis.label.set_size(20)
            ax.yaxis.label.set_size(20)
    else:
        #fig = triangle.corner(samples,labels=parameters_labels,
        #                      show_titles=False,
        #                      truths=bests[:],quantiles=[0.159, 0.5, 0.841],big_title=big_title
        fig = corner.corner(samples,labels=parameters_labels,
                              show_titles=False,
                              truths=bests[:],quantiles=[0.159, 0.5, 0.841],big_title=big_title)
        for ax in fig.get_axes():
            #ax.tick_params(axis='both', labelsize=14)
            #ax.set_xlabel(fontsize=14)
            ax.yaxis.label.set_size(20)
            ax.xaxis.label.set_size(20)
            #ax.get_ylabel().set_fontsize(60)
    if os.path.isdir(output_file_path)==True:
        pylab.tight_layout()
        fig.savefig(output_file_path+'/2D-distributions.png')
    else:
        print('output_file_path is', output_file_path)
        pylab.tight_layout()
        fig.savefig(output_file_path)

def plot_1D_marginalized_distribution(flatchain_path, bests=None,output_pdf_file_path='.',output_txt_file_path='.', parameters_labels=None,number_bins=None,verbose=False,title='marginalized posterior'):
        """Description: Given a chain, plots the 1-D marginalized distributions.
        Input  :-flatchain_Path: the path to the chain.
                -bests (optional):a numpy array with values of the parameters to be shown on the plots (ex: calculated best-fit values, true known values, etc)
                -output_file_path (optional): the path to the outputs folders. Default is '.'
                -parameters_labels (optional): a tuple with the names of the parameters to be added as xlabels. Default is Nonex
                -number of bins in the histograms. Default is 100
        Output :- no output. plots.
        pdf files:
                 - histo_param_' + parameters_labels(j) + 'pdf' histograms.
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example: histos=fitter_powerlaw.plot_1D_marginalized_distribution(flatchain_path='flatchain_2_test.txt',bests=bests,output_file_path='./output_test_fitter_powerlaw',parameters_labels=['a','b'],number_bins=50)
        Reliable:  """
        #if output_file_path == None:
        #    output_file_path = '.'
        #if number_bins == None:
        #    number_bins = 100
        if verbose == True:
            print('*** Plots of the 1D distributions ***')
        if isinstance(flatchain_path,np.ndarray)==True:
            if verbose==True:
                print('you gave flatchain_path as a numpy array')
            if len(np.shape(flatchain_path)) < 2:
                if verbose==True:
                    print('1D array')
                    print(np.shape(flatchain_path))
                sample=np.zeros((np.shape(flatchain_path)[0],1))
                sample[:,0]=flatchain_path
            else:
                sample=flatchain_path
            if verbose == True:
                print (np.shape(sample))
                print(sample)
        elif isinstance(flatchain_path,str):
            if verbose == True:
                print('you gave flatchain_path as an actual path')
            sample = np.loadtxt(flatchain_path)
        if bests is None:
            sigma = np.zeros((np.shape(sample)[1], 4))
        else:
            sigma = np.zeros((np.shape(sample)[1], 5))
        for j in range(np.shape(sample)[1]):
            if parameters_labels is None:
                if verbose == True:
                    print('I am plotting the histogram for the parameter number {0}'.format(j))
            else:
                if verbose == True:
                    print('I am plotting the histogram for the parameter {0}'.format(parameters_labels[j]))

            pylab.figure()
            n, bin, patches = pylab.hist(sample[:, j], bins=number_bins, alpha=0.5, color='orange')
            #               pylab.axvline(bn_mcmc2, color='r', linestyle='dashed', linewidth=3,label=r'$\chi^2={0}$'.format(maxparam[1,np.shape(maxparam)[1]-1]))
            #               pylab.axvline(bn_mcmc3, color='b', linestyle='dashed', linewidth=3,label=r'$\chi^2={0}$'.format(maxparam[2,np.shape(maxparam)[1]-1]))
            #               pylab.legend()
            bin2 = np.zeros((np.shape(n)))
            #print np.shape(bin)
            #pdb.set_trace()
            for i in range(np.shape(n)[0]):
                bin2[i] = bin[i] + (bin[i + 1] - bin[i]) / 2
            for i, b in enumerate(bin):
                if float(sum(n[:i])) / sum(n) >= 0.16:
                    #print np.shape(bin)
                    #print i
                    #print np.shape(n)
                    #print sum(n[:i])
                    #print bin2[i]
                    if verbose == True:
                        print('i is ',i)
                        print('bin2 is of shape',np.shape(bin2))
                        print('sun(n[:i]) is',sum(n[:i]))
                    nlower = bin[i]
                    if verbose == True:
                        print('at the {0}th bin, i.e. at b={1}, the ration float(sum(n[:i]))/sum(n) is {2} and the 16percent limit has been reached'.format(
                            i, bin2[i], float(sum(n[:i])) / sum(n)))
                    break
            for i, b in enumerate(bin):
                if float(sum(n[:i])) / sum(n) >= 0.84:
                    nupper = bin[i]
                    if verbose == True:
                        print('the rached ratio is {0}'.format(float(sum(n[:i])) / sum(n)))
                    #print('i is ',i
                    #print('bin2 is of shape',np.shape(bin2)
                    #print bin2[i]
                    #print('sun(n[:i]) is',sum(n[:i])
                    if verbose == True:
                        print('at the {0}th bin, i.e. at b={1}, the ration float(sum(n[:i]))/sum(n) is {2} and the 16percent upper limit has been reached'.format(
                        i, bin2[i], float(sum(n[:i])) / sum(n)))
                    break
            for i, b in enumerate(bin):
                if float(sum(n[:i])) / sum(n) >= 0.5:
                    nmedian = bin[i]
                    if verbose == True:
                        print('the rached ratio is {0}'.format(float(sum(n[:i])) / sum(n)))
                    #print('i is ',i
                    #print('bin2 is of shape',np.shape(bin2)
                    #print bin2[i]
                    #print('sun(n[:i]) is',sum(n[:i])
                        print('at the {0}th bin, i.e. at b={1}, the ration float(sum(n[:i]))/sum(n) is {2} and the 50percent upper limit has been reached'.format(
                        i, bin2[i], float(sum(n[:i])) / sum(n)))
                    break
            pylab.plot(bin2, n, color='k')
            if bests is not None:
                pylab.axvline(bests[j], color='k', linestyle='dashed', linewidth=3,label=r'Best fit')#. $\chi^2={0}$'.format(maxparam[0, np.shape(maxparam)[1] - 1]))
            pylab.axvline(nlower, color='b', linestyle='dashed', linewidth=3, label=r'$1\sigma$-range lower limit')
            pylab.axvline(nmedian, color='c', linestyle='dashed', linewidth=3, label=r'median')
            pylab.axvline(nupper, color='r', linestyle='dashed', linewidth=3, label=r'$1\sigma$-range upper limit')
            pylab.legend(fontsize = 'x-small')
            pylab.ylabel('Marginalized Posterior Distribution')
            pylab.title(title)
            if parameters_labels!=None:
                pylab.xlabel(parameters_labels[j])
                pylab.savefig(output_pdf_file_path + '/histo_param_' + parameters_labels[j] + '.pdf', facecolor='w', edgecolor='w',
                              orientation='portrait', papertype=None, format='pdf', transparent=False, bbox_inches=None,
                              pad_inches=0.1)
            else:
                pylab.savefig(output_pdf_file_path + '/histo_param_'+str(j)+'.pdf', facecolor='w', edgecolor='w',
                          orientation='portrait', papertype=None, format='pdf', transparent=False, bbox_inches=None,
                          pad_inches=0.1)
            pylab.close()
            if bests is None:
                sigma[j, :] = nlower, nupper, float(nupper - nlower) / 2,nmedian
                np.savetxt(output_txt_file_path + '/1sigma.txt', sigma, header='lower 1sigma, upper 1sigma, diff/2,median')
            else:
                sigma[j, :] = nlower, nupper, float(nupper - nlower) / 2, nmedian,bests[j]
                np.savetxt(output_txt_file_path + '/1sigma.txt', sigma, header='lower 1sigma, upper 1sigma, diff/2,median, best-fit')


        #pylab.show()