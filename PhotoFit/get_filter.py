#! //anaconda/bin/python



import pdb
import numpy as np
import pyphot

# It would make sens to reorganize this as a class Filter, with attributes: family_name, name,
# and methods effective wavelength, tranmission curve, filter half width, see for example Eran's get_filter

def filter_effective_wavelength(filter_family,filtername):
    """Description: inspired by Eran's finction of the same name. Given a filter name, finds the effective wavelength in microns
    Input  :-Filter family among:
                Johnson, sdss, poss
            -Filter name among:
                Johnson: 'U','B','V','R','I','J','H','K' or in minuscules
                sdss: 'u','g','r','i','z'
                poss: 'O','E'
    Output  :-effective wavelength in microns (kind of a hybrid between wavelength-weighted average and frequency-weighted average. See reference)
    Reference: Fukugita et al. 1996, Table 2a.
    Tested : ?
    By : Maayane T. Soumagnac Nov 2016
    URL :
    Example:
    Reliable:  """
    if filter_family.lower()=='sdss':
        if filtername.lower()=='u':
            W=0.3557
        elif filtername.lower()=='g':
            W=0.4825
        elif filtername.lower()=='r':
            W=0.6261
        elif filtername.lower()=='i':
            W=0.7672
        elif filtername.lower()=='z':
            W=0.9097
        else: print('unknown SDSS filter name')
    elif filter_family.lower()=='johnson':
        if filtername.lower()=='u':
            W=0.367#?
        elif filtername.lower()=='b':
            W=0.436#?
        elif filtername.lower()=='v':
            W=0.545
        elif filtername.lower()=='r':
            W=0.638
        elif filtername.lower()=='i':
            W=0.797
        elif filtername.lower()=='j':
            W=1.220
        elif filtername.lower()=='h':
            W=1.630
        elif filtername.lower()=='k':
            W=2.190
        else: print('unknown Johnson filter name')
    elif filter_family.lower()=='poss':
        if filtername.lower()=='o':
            W=0.420
        elif filtername.lower()=='e':
            W=0.640
    else:
        print('unknown family filter name. You can choose among [sdss, johnson, poss] and it is case insensitive')
        pdb.set_trace()
    return W


def make_filter_object(Filter_vector,filters_directory,central=True,verbose=False):
    """Description: from a filter vector where each element is a couple [filter family,filter name], create a filter object P as in pyphoy
            Input  :- a filter vector: a numpy array where each element is a tuple of strings [filter family,filter name]
                    - central. If true, gives pyphot .cl wavelength, which corresponds to Eran's AstFilter.get('family','band').eff_wl
                                else, gives phyphot .eff wavelength, which I am not sure what it is..
            Output :- filter vector, with additional column : P object
                    - a np.array with the effective wavelength of each filter
            with the corresponding data
            Tested : ?
                 By : Maayane T. Soumagnac Nov 2016
                URL :
            Example:Filter_vector = np.empty([2, 2], dtype=object)
                    Filter_vector[0] = [str('GALEX'), str('GALEX_NUV')]
                    Filter_vector[1]=[str('ptf_p48'),str('r')]
                    [P, wav]=make_filter_object(Filter_vector)
            Reliable:  """
    wavelength_filter_effective = dict()#np.empty(np.shape(Filter_vector)[0])
    wavelength_filter_central=dict()
    wavelength_filter_pivot = dict()
    P_vector=dict()
    if  isinstance(Filter_vector, (list,)):
        P_vector['filter_family']=Filter_vector[0][0]
        P_vector['filter_name']=Filter_vector[0][1]
        P_vector['filter_object'] = []
    else:
        P_vector['filter_family']=Filter_vector[:, 0]
        P_vector['filter_name']=Filter_vector[:,1]
        P_vector['filter_object'] = []
    #print("P_vector['filter_object'] is",P_vector['filter_object'])
    #pdb.set_trace()
    for i, j in enumerate(Filter_vector):
        '''
        if j[0].lower() == 'ptf_p48':
            #print(j[1])
            if j[1].lower() == 'g_p48_ptf':
                if verbose == True:
                    print('You gave the G filter of the PTF_P48 family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(
                    '/Users/maayanesoumagnac/Maayane_Astro_python_library/data/filters/PTF_G.rtf', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='PTF_P48_G', dtype='photon', unit='Angstrom')
            elif j[1].lower() == 'r_p48_ptf':
                if verbose==True:
                    print('You gave the R filter of the PTF_P48 family')
                Transmission = np.genfromtxt(
                    '/Users/maayanesoumagnac/Maayane_Astro_python_library/data/filters/P48_R_T.rtf', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='PTF_P48_R', dtype='photon', unit='Angstrom')
        '''
        if j[0].lower() == 'ztf_p48':
            #print(j[1])
            if j[1].lower() == 'g_p48':
                if verbose == True:
                    print('You gave the G filter of the ztf_P48 family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(filters_directory+'/ZTF_g_fromgithub_AA.txt', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='ZTF_P48_G', dtype='photon', unit='Angstrom')
            elif j[1].lower() == 'r_p48':
                if verbose==True:
                    print('You gave the R filter of the ZTF_P48 family')
                Transmission = np.genfromtxt(filters_directory+'/ZTF_r_fromgithub_AA.txt', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='ZTF_P48_R', dtype='photon', unit='Angstrom')
            elif j[1].lower() == 'i_p48':
                if verbose==True:
                    print('You gave the I filter of the ZTF_P48 family')
                Transmission = np.genfromtxt(filters_directory+'/ZTF_i_fromgithub_AA_reorder.txt', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='ZTF_P48_I', dtype='photon', unit='Angstrom')
        elif j[0].lower() == 'swift':
            if j[1].lower() == 'uvw1':
                if verbose == True:
                    print('You gave the uvw1 filter of the swift family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(filters_directory+'/Swift_UVW1.rtf', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_UVW1', dtype='photon',
                                  unit='Angstrom')
        #elif j[0].lower() == 'swift':
            elif j[1].lower() == 'uvw2':
                if verbose == True:
                    print('You gave the uvw2 filter of the swift family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/Swift_UVW2.rtf', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_UVW2', dtype='photon',
                                  unit='Angstrom')
        #elif j[0].lower() == 'swift':
            elif j[1].lower() == 'uvm2':
                if verbose == True:
                    print('You gave the uvm2 filter of the swift family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/Swift_UVM2.rtf', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_UVM2', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'u_swift':
                if verbose == True:
                    print('You gave the u filter of the swift family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/Swift_u.rtf', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_u', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'v_swift':
                if verbose == True:
                    print('You gave the v filter of the swift family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/Swift_V.rtf', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_V', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'b_swift':
                if verbose == True:
                    print('You gave the b filter of the swift family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/Swift_B.rtf', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='Swift_B', dtype='photon',
                                  unit='Angstrom')
        elif j[0].lower() == 'galex':
            if j[1].lower() == 'galex_nuv':
                if verbose == True:
                    print('You gave the nuv filter of the galex family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/GALEX_NUV.dat', delimiter=None)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='galex_nuv', dtype='photon',
                                  unit='Angstrom')
        elif j[0].lower() == 'sdss':
            if verbose == True:
                print('You gave the sdss family')
            #print(j[1].lower())
            if j[1].lower() == 'r_sdss':
                if verbose == True:
                    print('You gave the r filter of the sdss family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/SDSS_r.txt', delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='r_sdss', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'g_sdss':
                if verbose == True:
                    print('You gave the g filter of the sdss family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/SDSS_g.txt', delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='g_sdss', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'i_sdss':
                if verbose == True:
                    print('You gave the i filter of the sdss family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/SDSS_i.txt', delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='i_sdss', dtype='photon',
                                  unit='Angstrom')
            elif j[1].lower() == 'u_sdss':
                if verbose == True:
                    print('You gave the u filter of the sdss family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/SDSS_u.txt', delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='u_sdss', dtype='photon',
                                  unit='Angstrom')

            elif j[1].lower() == 'z_sdss':
                if verbose == True:
                    print('You gave the z filter of the sdss family')
                #pdb.set_trace()
                Transmission = np.genfromtxt(filters_directory+'/SDSS_z.txt', delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='z_sdss', dtype='photon',
                                  unit='Angstrom')

        elif j[0].lower() == '2mass':
            if verbose == True:
                print('You gave the 2MASS family')
                print(j[1].lower())
            if j[1].lower() == 'j_2mass':
                if verbose == True:
                    print('You gave the J filter of the 2MASS family')
            # pdb.set_trace()
                Transmission = np.genfromtxt(filters_directory+'/2MASS_J.txt',
                                         delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='2MASS_J', dtype='photon',
                              unit='Angstrom')
            elif j[1].lower() == 'h_2mass':
                if verbose == True:
                    print('You gave the h filter of the 2MASS family')
            # pdb.set_trace()
                Transmission = np.genfromtxt(filters_directory+'/2MASS_H.txt',
                                         delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='h_2mass', dtype='photon',
                              unit='Angstrom')
            elif j[1].lower() == 'k_2mass':
                if verbose == True:
                    print('You gave the k filter of the 2MASS family')
                    # pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/2MASS_K.txt',
                    delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='k_2mass', dtype='photon',
                                  unit='Angstrom')
        elif j[0].lower() == 'cousin':
            if verbose == True:
                print('You gave the cousin family')
                #print(j[1].lower())
            if j[1].lower() == 'i_cousin':
                if verbose == True:
                    print('You gave the i filter of the cousin family')
            # pdb.set_trace()
                Transmission = np.genfromtxt(filters_directory+'/cousin_i.txt',
                                         delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='i_counsin', dtype='photon',
                              unit='Angstrom')
            elif j[1].lower() == 'r_cousin':
                if verbose == True:
                    print('You gave the r filter of the cousin family')
                Transmission = np.genfromtxt(filters_directory+'/cousin_r.txt',
                                         delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='r_cousin', dtype='photon',
                              unit='Angstrom')
        elif j[0].lower() == 'johnson':
            if verbose == True:
                print('You gave the johnson family')
                print(j[1].lower())
            if j[1].lower() == 'u_johnson':
                if verbose == True:
                    print('You gave the u filter of the johnson family')
            # pdb.set_trace()
                Transmission = np.genfromtxt(filters_directory+'/johnson_u.txt',
                                         delimiter=',')
                #print('the shape of transission is',Transmission)
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='u_johnson', dtype='photon',
                              unit='Angstrom')
            elif j[1].lower() == 'b_johnson':
                if verbose == True:
                    print('You gave the b filter of the johnson family')
            # pdb.set_trace()
                Transmission = np.genfromtxt(filters_directory+'/johnson_b.txt',
                                         delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='b_johnson', dtype='photon',
                              unit='Angstrom')
            elif j[1].lower() == 'v_johnson':
                if verbose == True:
                    print('You gave the v filter of the johnson family')
                    # pdb.set_trace()
                Transmission = np.genfromtxt(
                    filters_directory+'/johnson_v.txt',
                    delimiter=',')
                P = pyphot.Filter(Transmission[:, 0], Transmission[:, 1], name='v_johnson', dtype='photon',
                                  unit='Angstrom')

        else:
            print('I HAVE NOT RECOGNIZE THE FILTER')
            lib = pyphot.get_library()
            f = lib.find(j[0].lower())  # filter family
            if verbose == True:
                for name in f:
                    lib[name].info(show_zeropoints=True)
            P = lib[j[1]]  # filter name
        #print(P_vector)
        #print(P)
        P_vector['filter_object'].append(P)
        wavelength_filter_effective[P_vector['filter_name'][i]]=P.leff.item()
        wavelength_filter_central[P_vector['filter_name'][i]]=P.cl.item()
        if isinstance(Filter_vector, (list,)):
            P_vector['filter_object']=P_vector['filter_object'][0]
    if central==True:
        return P_vector, wavelength_filter_central#same as Eran eff_wl
    else:
        return P_vector,wavelength_filter_effective

        #return 'hi'