import numpy as np
from pylastic.tools.CVS_2nd_deriv import ANALYTICS
from pylastic.elatoms import Structures, ElAtoms
from pylastic.postprocess import ECs
from pylastic.tools.eos import Birch
import matplotlib.pyplot as plt

from pylastic.thermo.fullqh import Postprocess

import pickle
import json
import os


Cij_dic = {}
GV=[]
GR=[]
GH=[]
BV=[]
BR=[]
BH=[]



ec = ECs('vasp', thermo=True, thmod='structures_phonons')
ec.set_structures()
ec.set_gsenergy()
gs_struct = ec.get_structures()


f=open('structures_phonons.pkl')
dic = pickle.load(f)
f.close()
# Fit vib free energgy:

ph_struct = dic.get_structures()


#ec.fitorder = [4,4,6]
#ec.set_analytics()

#ec.set_ec(['0.04','0.04','0.05'])




f=open('structures_phonons.pkl')
dic = pickle.load(f)
f.close()
# Fit vib free energgy:

ph_struct = dic.get_structures()
Vols = []
for a in sorted(ph_struct.keys()):
	if not round(ph_struct[a].V0,3) in Vols: Vols.append(round(ph_struct[a].V0,3))

etas = [-0.05,-0.025,0,0.025,0.05]

for v in Vols:
	for d in ['01','08','23']:
		Fenergy = []
		for i,eta in enumerate(etas):
			Fenergy.append(ph_struct[(d,eta,v)].fenergy)
			Temps = ph_struct[(d,eta,v)].T
		Fenergy = np.transpose(np.array(Fenergy)/125.)
		
		p=[]
		
		for i,T in enumerate(Temps):
			
			coeff = np.polyfit(etas,Fenergy[:][i],3)
			
			p.append(np.poly1d(coeff))
		
		for ETA in [-0.05,-0.045,-0.04,-0.035,-0.03,-0.025,-0.02,-0.015,-0.01,-0.005,0.,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05]:
			phenergy = []	
			for i,T in enumerate(Temps):
				phenergy.append(p[i](ETA))
			
		
			gs_struct[d,ETA,v].phenergy = phenergy
			gs_struct[d,ETA,v].T = Temps

print gs_struct[('01',-0.05,Vols[0])].phenergy, gs_struct[('01',-0.05,Vols[0])].gsenergy

post = Postprocess()
post.set_structures(gs_struct)
post.rewrite_structures()
post.ECs()
