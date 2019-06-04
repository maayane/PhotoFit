import numpy as np
import pylab

ZTF_i_github=np.genfromtxt('ZTF_i_fromgithub_AA.txt')
ZTF_g_github=np.genfromtxt('ZTF_g_fromgithub_AA.txt')
ZTF_r_github=np.genfromtxt('ZTF_r_fromgithub_AA.txt')
ZTF_g_Eran=np.genfromtxt('PTF_G.rtf')
ZTF_r_Eran=np.genfromtxt('P48_R_T.rtf') 

pylab.figure()
pylab.plot(ZTF_i_github[:,0],ZTF_i_github[:,1],'b-',label='ZTF i')
pylab.plot(ZTF_g_github[:,0],ZTF_g_github[:,1],'g-',label='ZTF g')
pylab.plot(ZTF_g_Eran[:,0],ZTF_g_Eran[:,1],'g--',label='PTF g')
pylab.plot(ZTF_r_Eran[:,0],ZTF_r_Eran[:,1],'r--',label='PTF r')
pylab.plot(ZTF_r_github[:,0],ZTF_r_github[:,1],'r-',label='ZTF r')
pylab.legend()
pylab.show()

