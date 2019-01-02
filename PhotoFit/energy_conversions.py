"""*******************************************************
This module has functions converting energies
*****************************************************
"""
#print __doc__

def convert_energy(Energy, Input_unit, Output_unit):#converts distance modulus to distance in parsec
    """Description: ocnverts energy from Input_unit to Output_unit.
    Input  :- Energy to be converted
    		- Input unit. Can be:
    		'erg'   - ergs
    		'J'     - Jouls
    		'Hz'    - Frequency [1/s]
    		'A','Ang'- Wavelength [Ang]
    		'cm'    - Wavelength [cm]
    		'nm'    - Wavelength [nm]
    		'm'     - Wavelength [m]
    		'eV'    - Electron volts [h nu/q]
    		'keV'   - kilo Electron volts [h nu/q]
    		'MeV'   - Mega Electron volts [h nu/q]
    		'GeV'   - Giga Electron volts [h nu/q]
    		'T'     - Temperature [K]
    		'me'    - Electron mass [E/m_e]
    		'mp'    - Proton mass [E/m_p]
    		'cal'   - calorie (4.184 J)
    		'Btu'   - (1.055x10^3 J)
    		'kWh'   - kilowatt-hour (3.6x10^6 J)
    		'TNT'   - one ton of TNT (4.2x10^9 J)
    		'gr'    - Energy equivalent of 1 gram of matter (9x10^13 J)
    		- requested output unit. One of the above.
    Output :- Energy in output system
     Tested : ?
         By : Maayane T. Soumagnac Nov 2016
        URL :
     Example:
     Reliable:  """
    Erg2Erg = 1.
    Erg2J = 1e-7# 1erg=10**-7J
    #Erg2Hz = 1.5092e26#?
    #Erg2A = 1.9864e-8 #?
    Erg2eV = 6.2415e11
    #Erg2T = 7.2430e15#?
    Erg2me = 1.2214e6
    Erg2mp = 665.214577
    Erg2cal = Erg2J/4.184
    Erg2Btu = Erg2J/1.055e3
    Erg2kWh = Erg2J/3.6e6
    Erg2TNT = Erg2J/4.2e9
    #Erg2gr = get_constant('c', 'cgs'). ^ -2

    relation='inv'
    if Input_unit.lower()=='erg':
        ConvFactor = Erg2Erg
    elif Input_unit.lower()=='j':
        ConvFactor=Erg2J
    #elif Input_unit.lower()=='hz':
    #    ConvFactor=Erg2Hz
    #elif Input_unit.lower() in ('a','ang'):
    #    relation='inv'
    #    ConvFactor=Erg2A
    '''
    elif Input_unit.lower()=='cm':
        relation='inv'
        ConvFactor=Erg2A*1e-8
    elif Input_unit.lower()=='nm':
        relation='inv'
        ConvFactor=Erg2A*1e-4
    #elif Input_unit.lower()=='m':
    #    relation='inv'
    #    ConvFactor=Erg2A*1e-10
    elif Input_unit.lower()=='ev':
        ConvFactor = Erg2eV
    elif Input_unit.lower()=='kev':
        ConvFactor = Erg2eV*1e-3
    elif Input_unit.lower()=='mev':
        ConvFactor = Erg2eV*1e-6
    elif Input_unit.lower()=='gev':
        ConvFactor = Erg2eV* 1e-9
    #elif Input_unit.lower()=='t':
    #    ConvFactor = Erg2T
    elif Input_unit.lower()=='me':
        ConvFactor = Erg2me
    elif Input_unit.lower()=='mp':
        ConvFactor = Erg2mp
    elif Input_unit.lower()=='cal':
        ConvFactor = Erg2cal
    elif Input_unit.lower()=='btu':
        ConvFactor = Erg2Btu
    elif Input_unit.lower()=='kwh':
        ConvFactor = Erg2kWh
    elif Input_unit.lower()=='tnt':
        ConvFactor = Erg2TNT
    #elif Input_unit =='gr':
    #    ConvFactor = Erg2gr
    else: print 'error: unknown InUnit option'
    '''
    if relation=='lin':
        ErgE=Energy*ConvFactor
    elif relation=='inv':
        ErgE=Energy/ConvFactor
    else: print('unknown relation')
    
    relation = 'lin'
    if Output_unit.lower()=='erg':
        ConvFactor = Erg2Erg
    elif Output_unit.lower()== 'j':
        ConvFactor = Erg2J
    #elif Output_unit.lower()== 'hz':
    #    ConvFactor = Erg2Hz
    #elif Output_unit.lower()== {'a','ang'}:
    #    relation   = 'inv'
    #    ConvFactor = Erg2A
    #elif Output_unit.lower()== 'nm':
    #    relation   = 'inv'
    #    ConvFactor = Erg2A*1e-4
    #elif Output_unit.lower()== 'cm':
    #    relation   = 'inv'
    #     ConvFactor = Erg2A*1e-8
    #elif Output_unit.lower()== 'm':
    #    relation = 'inv'
    #    ConvFactor = Erg2A*1e-10
    '''
    elif Output_unit.lower()== 'ev':
        ConvFactor = Erg2eV
    elif Output_unit.lower()== 'kev':
        ConvFactor = Erg2eV*1e-3
    elif Output_unit.lower()== 'mev':
        ConvFactor = Erg2eV*1e-6
    elif Output_unit.lower()== 'gev':
        ConvFactor = Erg2eV*1e-9
    #elif Output_unit.lower()== 't':
    #    ConvFactor = Erg2T
    elif Output_unit.lower()== 'me':
        ConvFactor = Erg2me
    elif Output_unit.lower()== 'mp':
        ConvFactor = Erg2mp
    elif Output_unit.lower()== 'cal':
        ConvFactor = Erg2cal
    elif Output_unit.lower()== 'btu':
        ConvFactor = Erg2Btu
    elif Output_unit.lower()== 'kwh':
        ConvFactor = Erg2kWh
    elif Output_unit.lower()== 'tnt':
        ConvFactor = Erg2TNT
    #elif Output_unit.lower()== 'gr':
    #    ConvFactor = Erg2gr
    else: print 'Unknown InUnit Option'
'''
    if relation == 'lin':
        new_energy = ErgE * ConvFactor
    elif relation == 'inv':
        new_energy= ErgE/ConvFactor
    else:
        print('unknown relation')
    return new_energy