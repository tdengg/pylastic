from pylastic.thermo import plot_thermal

import os
import json
import collections
import matplotlib.pyplot as plt
import numpy as np
path = '/home/t.dengg/results/W/Cij_T_exp/T_expansion_PBEsol'

###### alpha(T), l_min(T) ########
inst = plot_thermal.PLOT_THERMAL(path)
import numpy as np
inst.free_ene(np.linspace(0,1000,101))
alpha = inst.alpha
#print inst.minl


CV = {}
###### C_v(T,V) ########
for d in os.listdir(path):
	if 'scale_' in d:
		l = float(d.lstrip('scale_'))
				
		f = open(path + '/' + d+ '/thermal_properties')
		lines = f.readlines()
		CV[l] = [float(line.split()[-1]) for line in lines]
		f.close
	
	else:
		continue
		
pCV = []
for i in range(len(np.linspace(0,1000,101))):
	C_VatT = []
	latt = []
	for key in CV.keys():
		latt.append(key)
		C_VatT.append(CV[key][i])
	coeff = np.polyfit(latt,C_VatT,5)
	pCV.append(np.poly1d(coeff))

	#plt.plot(np.linspace(min(latt),max(latt),101),pCV(np.linspace(min(latt),max(latt),101)))

####### C_ij ###########
with open('/home/t.dengg/results/W/Cij_T_exp/Cij_0/FT_out.json') as data_file:    
    Cij_dict = json.load(data_file)

with open('/home/t.dengg/results/W/Cij_T_exp/T_expansion/FT_out.json') as data_file:    
    lMin_T_dict = json.load(data_file)
    
Cij_dict_sorted = collections.OrderedDict(sorted(Cij_dict.items()))

lMin_T_dict_sorted = collections.OrderedDict(sorted(lMin_T_dict.items()))
C11 = []
C12 = []
C44 = []
l = []
for key in Cij_dict_sorted.keys(): 
    C11.append( Cij_dict_sorted[key]['SM'][0] )
    C12.append( Cij_dict_sorted[key]['SM'][1] )
    C44.append( Cij_dict_sorted[key]['SM'][-1] )
for key in Cij_dict_sorted.keys(): l.append(float(key))

#### C11(T)
coeff11 = np.polyfit(l,C11,5)
pC11 = np.poly1d(coeff11)
polyCx = np.linspace(min(l),max(l),1000)

#### C12(T)
coeff12 = np.polyfit(l,C12,5)
pC12 = np.poly1d(coeff12)
polyCx = np.linspace(min(l),max(l),1000)

#### C44(T)
coeff44 = np.polyfit(l,C44,5)
pC44 = np.poly1d(coeff44)
polyCx = np.linspace(min(l),max(l),1000)

F_min = lMin_T_dict_sorted['F_min']
l_min = lMin_T_dict_sorted['l_min']
T = lMin_T_dict_sorted['T']


ax1 = plt.subplot(231)
ax1.plot(l, C11,'+')
ax1.plot(polyCx, pC11(polyCx))
ax1.plot(l_min,pC11(l_min),'o')

ax2 = plt.subplot(232)
ax2.plot(T,pC11(l_min))
ax2.plot(T,pC12(l_min))
ax2.plot(T,pC44(l_min))

ax3 = plt.subplot(233)

pCV_T = []
for i in range(len(np.linspace(0,1000,101))):
	pCV_T.append(pCV[i](l_min[i]))
	
## Calculate Adiabatic correction term
correction = []
contribution_alpha = []
contribution_C = []
contribution_CV = []
contribution_V = []
for i in range(len(T)-1):
	
	c_alpha = (alpha[i]*10**(-6.))**2.
	c_C = ((pC11(l_min[i])+2.*pC12(l_min[i]))*10**9)**2.
	c_CV = pCV_T[i]
	c_V = (l_min[i]*10**(-10))**3./2.
	contribution_alpha.append(c_alpha)
	contribution_C.append(c_C)
	contribution_CV.append(c_CV)
	contribution_V.append(c_V)
	correction.append((  c_V  * c_alpha * c_C * T[i] / c_CV * 6.022*10**(23.) )*10**(-9))

	#correction.append((   ( (l_min[i]*10**(-10))**3./2.) * (alpha[i]*10**(-6.))**2. * ((pC11(l_min[i])+2.*pC12(l_min[i]))*10**9)**2. * T[i] / pCV_T[i] * 6.022*10**(23.) )*10**(-9))
T.pop()


C11_T_adiabatic = []
C12_T_adiabatic = []
C44_T_adiabatic = []
for i in range(len(T)):
	C11_T_adiabatic.append(pC11(l_min[i]) + correction[i])
	C12_T_adiabatic.append(pC12(l_min[i]) + correction[i])
	C44_T_adiabatic.append(pC44(l_min[i]))
	

ax3.plot(T,C11_T_adiabatic)
ax3.plot(T,C12_T_adiabatic)
ax3.plot(T,C44_T_adiabatic)

ax4 = plt.subplot(234)
ax4.plot(T,contribution_alpha)
ax5 = plt.subplot(235)
ax5.plot(T,contribution_C)
ax6 = plt.subplot(236)
ax6.plot(T,contribution_CV)

plt.show()
