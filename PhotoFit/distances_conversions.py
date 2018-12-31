
"""*******************************************************
This module has functions converting distances
*****************************************************
"""
#print __doc__


def DM_to_pc(x):#converts distance modulus to distance in parsec
    psec=10**(x/5+1)
    return psec

def pc_to_cm(x):#convert pc to cm
    cm=x*3.0857e18
    return cm
def pc_to_m(x):#convert pc to cm
    m=x*3.0857e16
    return m
def cm_to_pc(x):
    ps=x/(3.0857e18)
    return pc
def solar_radius_to_m(x):#convert pc to cm
    m=x*6.957e8
    return m
def solar_radius_to_pc(x):#convert pc to cm
    m=x*2.25461e-8
    return m
def cm_to_solar_radius(x):
    sr=x/6.957e10
    return sr
def solar_radius_to_cm(x):
    cm=x*6.957e10
    return cm