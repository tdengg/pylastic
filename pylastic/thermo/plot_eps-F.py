import matplotlib.pyplot as plt
import os
import lxml.etree as et
import numpy as np
import copy
import threading

import pylastic.tools.CVS_2nd_deriv as CVS2

class AnalyzeF(object):

	def __init__(self,distortion='Dst02',scale = 'scale_3.17226',rootdir='/home/t.dengg/results/W/Cij_T_F/Cij_T',path = '/',etamin = -0.05,etamax = 0.05,n = 21,T_index=100,fitorder = 6,delta_x_deriv = 0.001):

               	_e     = 1.602176565e-19 	     # elementary charge
      	       	Bohr   = 5.291772086e-11 	     # a.u. to meter
      	       	Ryd2eV = 13.605698066		     # Ryd to eV
      	       	Angstroem = 1.e-10		     # Angstroem to meter
      	       	ToGPa  = (_e*Ryd2eV)/(1e9*Bohr**3)    # Ryd/[a.u.^3] to GPa
      	       	vToGPa  = (_e)/(1e9*Angstroem**3)   # eV/[Angstroem^3] to GPa




		plt.figure()
      	       	ax1=plt.subplot(221)
      	       	ax2=plt.subplot(222)
      	       	ax3=plt.subplot(223)
      	       	ax4=plt.subplot(224)


      	       	##########################################
      	       	############## MANUAL INPUT ##############
      	       	#distortion='Dst02'
      	       	#scale = 'scale_3.17226'
      	       	#rootdir='/home/t.dengg/results/W/Cij_T_F/Cij_T'
      	       	#path = rootdir+'/'+scale+'/'+distortion+'/'
      	       	#etamin = -0.05
      	       	#etamax = 0.05
      	       	#n = 21
		#
      	       	#T_index=100
		#
      	       	#fitorder = 6
		#
      	       	#delta_x_deriv = 0.001
      	        ##########################################
      	       	##########################################


		colors = ['b','g','r','c','m','k','y']


               	F = []
               	F1 = []
               	F2 = []
               	eps = []
               	volume=[]
               	i=0
               	for d in sorted(os.listdir(path)):
                       	if os.path.isdir(path+d):
                       
                       
                	       	vasprun = et.parse(path+d+'/'+'vasprun.xml')
                	       	volume.append(vasprun.xpath("//i[@name='volume']")[0].text)
                       
                       
                	       	f=open(path+d+'/'+'F_TV_back')
                	       	data = f.readlines()
                	       	F.append(float(data[T_index].split()[1])/96.47244)
                	       	F1.append(float(data[10].split()[1])/96.47244)
                	       	F2.append(float(data[100].split()[1])/96.47244)
                	       	f.close()
                	       	eps.append(i)
                	       	i+=1
                       	else:
                	      	 continue


               	eps = map(float,eps)
               	eta = np.linspace(etamin,etamax,n)

               	coeff = np.polyfit(eta,F,fitorder)
               	p = np.poly1d(coeff)
               	coeffl = np.polyfit(eta,F,1)
               	pl = np.poly1d(coeffl)
               	polyx = np.linspace(min(eta),max(eta),1000)


               	r=0
               	for i in range(len(F)):
               		r = r + (F[i]-p(eta[i]))**2.
               	rms1 = np.sqrt(r/float(len(F)))
               	r=0
               	for i in range(len(F)):
               	        r = r + (F[i]-pl(eta[i]))**2.
               	rms2 = np.sqrt(r/float(len(F)))
               	print "##############################################"
               	print "# --> Calculated root mean square error (rms):"
               	print "#---------------------------------------------"
               	print '# rms (fitorder: %s)'%fitorder,rms1
               	print '# rms (linear fit)',rms2
               	print "##############################################"

		m=0
               	for j in [2,3,4,6,8]:
                       	CVS=[]
                       	ETA=[]
                       	DERIV=[]
                       	DERIVx=[]
                       	strain = list(copy.copy(eta))
                       	energy = list(copy.copy(F))
                       	while (len(strain) > j+1): 
                	       	emax = max(strain)
                	       	emin = min(strain)
                	       	emax = max(abs(emin),abs(emax))
                	       	coeff = np.polyfit(strain,energy,j)
                	       	DERIVx.append(np.polyval(self.deriv(coeff,j),delta_x_deriv)*vToGPa*2)
                	       	DERIV.append(np.polyfit(strain,energy,j)[j-2]*2.*vToGPa*2)
                	       	S = 0
                	       	for k in range(len(strain)):
                		       	Y      = energy[k]
                		       	etatmp = []
                		       	enetmp = []

                		       	for l in range(len(strain)):
                			       	if (l==k): pass
                			       	else:		
                				       	etatmp.append(strain[l])
                				       	enetmp.append(energy[l])

                		       	Yfit = np.polyval(np.polyfit(etatmp,enetmp, j), strain[k])
                		       	S    = S + (Yfit-Y)**2
                       
                	       	CV = np.sqrt(S/len(strain))
                       
                       
                	       	CVS.append(CV)
                	       	ETA.append(emax)
                	       	if (abs(strain[0]+emax) < 1.e-7):
                		       	strain.pop(0)
                		       	energy.pop(0)
                	       	if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                		       	strain.pop()
                		       	energy.pop()
                       	ax3.plot(ETA,CVS, '%s'%colors[m] ,label='%s'%j)
                       	if j%2==0:
                	       	ax4.plot(ETA,DERIV,'%s'%colors[m],label='%s'%j)
			if j==8:
                       		ax4.plot(ETA,DERIVx, '%s--'%colors[m],label='%s, dx=%s'%(j,delta_x_deriv))
			m+=1
		
		m=0	
		for j in [2,3,4,6,8]:
			cvs2 = CVS2.CVS(eta,F)
			CVS = cvs2.calc_cvs(j)
			A1 = [i[1] for i in CVS]
			A2 = [i[0] for i in CVS]
			ax2.plot(A1,A2, '%s'%colors[m])
			m+=1
			
               	ax3.legend(title='Fitorder:')
               	ax3.set_ylabel('CVS')
               	ax4.legend(title='Fitorder:',loc=4)
		ax4.set_ylabel(r'$d^2E/d\epsilon^2$   $Pa$')

               	ax1.plot(eta,F,'o')

               	ax1.plot(polyx,p(polyx), label='6')
               	ax1.plot(polyx,pl(polyx), label='Linear fit')
		ax1.legend(title='Fitorder:')
               	#ax2.plot(eta,volume)

               	#plt.plot(eps,F1)
               	#plt.plot(eps,F2)
		ax1.set_title('Vibrational Free Energy')
		ax3.set_xlabel(r'$\eta_{max}$')
		ax2.set_ylabel(r'CVS($d^2F_{vib}/d\eta^2$)')
		ax1.set_ylabel(r'$F_{vib}$    $eV$')
		ax4.set_xlabel(r'$\eta_{max}$')
              	plt.draw()
	       
	def deriv(self, coeff,order):
		dc=[]
		for i in range(len(coeff)-2):
			dc.append(((order-i)*(order-1-i))*coeff[i])
		return dc

class AnalyzeFE(object):

	def __init__(self,distortion='Dst02',scale = 'scale_3.17226',rootdir='/home/t.dengg/results/W/Cij_T_F/Cij_T',path = '/',etamin = -0.05,etamax = 0.05,n = 21,T_index=100,fitorder = 6,delta_x_deriv = 0.001):

               	_e     = 1.602176565e-19 	     # elementary charge
      	       	Bohr   = 5.291772086e-11 	     # a.u. to meter
      	       	Ryd2eV = 13.605698066		     # Ryd to eV
      	       	Angstroem = 1.e-10		     # Angstroem to meter
      	       	ToGPa  = (_e*Ryd2eV)/(1e9*Bohr**3)    # Ryd/[a.u.^3] to GPa
      	       	vToGPa  = (_e)/(1e9*Angstroem**3)   # eV/[Angstroem^3] to GPa




		plt.figure()
      	       	ax1=plt.subplot(221)
      	       	ax2=plt.subplot(222)
      	       	ax3=plt.subplot(223)
      	       	ax4=plt.subplot(224)


      	       	##########################################
      	       	############## MANUAL INPUT ##############
      	       	#distortion='Dst02'
      	       	#scale = 'scale_3.17226'
      	       	#rootdir='/home/t.dengg/results/W/Cij_T_F/Cij_T'
      	       	#path = rootdir+'/'+scale+'/'+distortion+'/'
      	       	#etamin = -0.05
      	       	#etamax = 0.05
      	       	#n = 21

      	       	#T_index=100

      	       	#fitorder = 6

      	       	#delta_x_deriv = 0.001
      	        ##########################################
      	       	##########################################


		colors = ['b','g','r','c','m','k','y']


               	F = []
               	F1 = []
               	F2 = []
               	eps = []
               	volume=[]
               	i=0
               	for d in sorted(os.listdir(path)):
                       	if os.path.isdir(path+d):
                       
                       		
				
				cell = 125.#supercell size
                    
                    		vaspout = et.parse(path+d+'/'+'vasprun.xml')
                    		elem = vaspout.xpath("//scstep/energy/i[@name = 'e_fr_energy']")
                    		allengys = []
                    		for k in elem:
                        		try:
                            			allengys.append(float(k.text))
                        		except:
                            			allengys.append(0.)
                    
                    		trueengys = []
                    		for engy in allengys:
                        		if engy < -1000. and engy > -3000.: trueengys.append(engy)
                    		gsenergy = trueengys[-1]/cell
				
				
				
                	       	vasprun = et.parse(path+d+'/'+'vasprun.xml')
                	       	volume.append(vasprun.xpath("//i[@name='volume']")[0].text)
                       
                       
                	       	f=open(path+d+'/'+'F_TV_back')
                	       	data = f.readlines()
                	       	F.append(float(data[T_index].split()[1])/96.47244+gsenergy)
                	       	F1.append(float(data[10].split()[1])/96.47244)
                	       	F2.append(float(data[100].split()[1])/96.47244)
                	       	f.close()
                	       	eps.append(i)
                	       	i+=1
                       	else:
                	      	 continue


               	eps = map(float,eps)
               	eta = np.linspace(etamin,etamax,n)

               	coeff = np.polyfit(eta,F,fitorder)
               	p = np.poly1d(coeff)
               	coeffl = np.polyfit(eta,F,4)
               	pl = np.poly1d(coeffl)
               	polyx = np.linspace(min(eta),max(eta),1000)


               	r=0
               	for i in range(len(F)):
               		r = r + (F[i]-p(eta[i]))**2.
               	rms1 = np.sqrt(r/float(len(F)))
               	r=0
               	for i in range(len(F)):
               	        r = r + (F[i]-pl(eta[i]))**2.
               	rms2 = np.sqrt(r/float(len(F)))
               	print "##############################################"
               	print "# --> Calculated root mean square error (rms):"
               	print "#---------------------------------------------"
               	print '# rms (fitorder: %s)'%fitorder,rms1
               	print '# rms (linear fit)',rms2
               	print "##############################################"

		m=0
               	for j in [2,3,4,6,8]:
                       	CVS=[]
                       	ETA=[]
                       	DERIV=[]
                       	DERIVx=[]
                       	strain = list(copy.copy(eta))
                       	energy = list(copy.copy(F))
                       	while (len(strain) > j+1): 
                	       	emax = max(strain)
                	       	emin = min(strain)
                	       	emax = max(abs(emin),abs(emax))
                	       	coeff = np.polyfit(strain,energy,j)
                	       	DERIVx.append(np.polyval(self.deriv(coeff,j),delta_x_deriv)*vToGPa*2)
                	       	DERIV.append(np.polyfit(strain,energy,j)[j-2]*2.*vToGPa*2)
                	       	S = 0
                	       	for k in range(len(strain)):
                		       	Y      = energy[k]
                		       	etatmp = []
                		       	enetmp = []

                		       	for l in range(len(strain)):
                			       	if (l==k): pass
                			       	else:		
                				       	etatmp.append(strain[l])
                				       	enetmp.append(energy[l])

                		       	Yfit = np.polyval(np.polyfit(etatmp,enetmp, j), strain[k])
                		       	S    = S + (Yfit-Y)**2
                       
                	       	CV = np.sqrt(S/len(strain))
                       
                       
                	       	CVS.append(CV)
                	       	ETA.append(emax)
                	       	if (abs(strain[0]+emax) < 1.e-7):
                		       	strain.pop(0)
                		       	energy.pop(0)
                	       	if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                		       	strain.pop()
                		       	energy.pop()
                       	ax3.plot(ETA,CVS,'%s'%colors[m], label='%s'%j)
                       	if j%2==0:
                	       	ax4.plot(ETA,DERIV,'%s'%colors[m] ,label='%s'%j)
			if j==8:
                       		ax4.plot(ETA,DERIVx, '%s--'%colors[m],label='%s, dx=%s'%(j,delta_x_deriv))
			m+=1
			
               	ax3.legend(title='Fitorder:')
               	ax3.set_ylabel('CVS')
               	ax4.legend(title='Fitorder:', loc=2)
		ax4.set_ylabel(r'$d^2F/d\eta^2$    $Pa$')

               	ax1.plot(eta,F,'o')

               	ax1.plot(polyx,p(polyx), label='6')
               	ax1.plot(polyx,pl(polyx), label='4')
		ax1.legend(title='Fitorder:')
               	#ax2.plot(eta,volume)
	
               	#plt.plot(eps,F1)
               	#plt.plot(eps,F2)
		ax1.set_title('Free Energy')
		ax3.set_xlabel('$\eta_{max}$')
		ax2.set_ylabel(r'CVS($d^2F/d\eta^2$)')
		ax1.set_ylabel('$F$    $eV$')
		ax4.set_xlabel('$\eta_{max}$')
              	plt.draw()
		
		m=0	
		for j in [2,3,4,6,8]:
			cvs2 = CVS2.CVS(eta,F)
			CVS = cvs2.calc_cvs(j)
			A1 = [i[1] for i in CVS]
			A2 = [i[0] for i in CVS]
			ax2.plot(A1,A2, '%s'%colors[m])
			m+=1
	       
	def deriv(self, coeff,order):
		dc=[]
		for i in range(len(coeff)-2):
			dc.append(((order-i)*(order-1-i))*coeff[i])
		return dc

class AnalyzeFE_seperated_fit(object):

	def __init__(self,distortion='Dst02',scale = 'scale_3.17226',rootdir='/home/t.dengg/results/W/Cij_T_F/Cij_T',path = '/',etamin = -0.05,etamax = 0.05,n = 21,T_index=100,fitorder = 6,delta_x_deriv = 0.001):

               	_e     = 1.602176565e-19 	     # elementary charge
      	       	Bohr   = 5.291772086e-11 	     # a.u. to meter
      	       	Ryd2eV = 13.605698066		     # Ryd to eV
      	       	Angstroem = 1.e-10		     # Angstroem to meter
      	       	ToGPa  = (_e*Ryd2eV)/(1e9*Bohr**3)    # Ryd/[a.u.^3] to GPa
      	       	vToGPa  = (_e)/(1e9*Angstroem**3)   # eV/[Angstroem^3] to GPa




		plt.figure()
      	       	ax1=plt.subplot(221)
      	       	ax2=plt.subplot(222)
      	       	ax3=plt.subplot(223)
      	       	ax4=plt.subplot(224)


      	       	##########################################
      	       	############## MANUAL INPUT ##############
      	       	#distortion='Dst02'
      	       	#scale = 'scale_3.17226'
      	       	#rootdir='/home/t.dengg/results/W/Cij_T_F/Cij_T'
      	       	#path = rootdir+'/'+scale+'/'+distortion+'/'
      	       	#etamin = -0.05
      	       	#etamax = 0.05
      	       	#n = 21

      	       	#T_index=100

      	       	#fitorder = 6

      	       	#delta_x_deriv = 0.001
      	        ##########################################
      	       	##########################################


		colors = ['b','g','r','c','m','k','y']


               	F = []
               	F1 = []
               	F2 = []
               	eps = []
               	volume=[]
               	i=0
               	for d in sorted(os.listdir(path)):
                       	if os.path.isdir(path+d):
                       
                       		
				
				cell = 125.#supercell size
                    
                    		vaspout = et.parse(path+d+'/'+'vasprun.xml')
                    		elem = vaspout.xpath("//scstep/energy/i[@name = 'e_fr_energy']")
                    		allengys = []
                    		for k in elem:
                        		try:
                            			allengys.append(float(k.text))
                        		except:
                            			allengys.append(0.)
                    
                    		trueengys = []
                    		for engy in allengys:
                        		if engy < -1000. and engy > -3000.: trueengys.append(engy)
                    		gsenergy = trueengys[-1]/cell
				
				
				
                	       	vasprun = et.parse(path+d+'/'+'vasprun.xml')
                	       	volume.append(vasprun.xpath("//i[@name='volume']")[0].text)
                       
                       
                	       	f=open(path+d+'/'+'F_TV')
                	       	data = f.readlines()
                	       	F.append(float(data[T_index].split()[1])/96.47244+gsenergy)
                	       	F1.append(float(data[10].split()[1])/96.47244)
                	       	F2.append(float(data[100].split()[1])/96.47244)
                	       	f.close()
                	       	eps.append(i)
                	       	i+=1
                       	else:
                	      	 continue


               	eps = map(float,eps)
               	eta = np.linspace(etamin,etamax,n)

               	coeff = np.polyfit(eta,F,fitorder)
               	p = np.poly1d(coeff)
               	coeffl = np.polyfit(eta,F,4)
               	pl = np.poly1d(coeffl)
               	polyx = np.linspace(min(eta),max(eta),1000)


               	r=0
               	for i in range(len(F)):
               		r = r + (F[i]-p(eta[i]))**2.
               	rms1 = np.sqrt(r/float(len(F)))
               	r=0
               	for i in range(len(F)):
               	        r = r + (F[i]-pl(eta[i]))**2.
               	rms2 = np.sqrt(r/float(len(F)))
               	print "##############################################"
               	print "# --> Calculated root mean square error (rms):"
               	print "#---------------------------------------------"
               	print '# rms (fitorder: %s)'%fitorder,rms1
               	print '# rms (linear fit)',rms2
               	print "##############################################"

		m=0
               	for j in [2,3,4,6,8]:
                       	CVS=[]
                       	ETA=[]
                       	DERIV=[]
                       	DERIVx=[]
                       	strain = list(copy.copy(eta))
                       	energy = list(copy.copy(F))
                       	while (len(strain) > j+1): 
                	       	emax = max(strain)
                	       	emin = min(strain)
                	       	emax = max(abs(emin),abs(emax))
                	       	coeff = np.polyfit(strain,energy,j)
                	       	DERIVx.append(np.polyval(self.deriv(coeff,j),delta_x_deriv)*vToGPa*2)
                	       	DERIV.append(np.polyfit(strain,energy,j)[j-2]*2.*vToGPa*2)
                	       	S = 0
                	       	for k in range(len(strain)):
                		       	Y      = energy[k]
                		       	etatmp = []
                		       	enetmp = []

                		       	for l in range(len(strain)):
                			       	if (l==k): pass
                			       	else:		
                				       	etatmp.append(strain[l])
                				       	enetmp.append(energy[l])

                		       	Yfit = np.polyval(np.polyfit(etatmp,enetmp, j), strain[k])
                		       	S    = S + (Yfit-Y)**2
                       
                	       	CV = np.sqrt(S/len(strain))
                       
                       
                	       	CVS.append(CV)
                	       	ETA.append(emax)
                	       	if (abs(strain[0]+emax) < 1.e-7):
                		       	strain.pop(0)
                		       	energy.pop(0)
                	       	if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                		       	strain.pop()
                		       	energy.pop()
                       	ax3.plot(ETA,CVS,'%s'%colors[m], label='%s'%j)
                       	if j%2==0:
                	       	ax4.plot(ETA,DERIV,'%s'%colors[m],label='%s'%j)
			if j==8:
                       		ax4.plot(ETA,DERIVx, '%s--'%colors[m],label='%s, dx=%s'%(j,delta_x_deriv))
			m+=1
		
		m=0	
		for j in [2,3,4,6,8]:
			cvs2 = CVS2.CVS(eta,F)
			CVS = cvs2.calc_cvs(j)
			A1 = [i[1] for i in CVS]
			A2 = [i[0] for i in CVS]
			ax2.plot(A1,A2, '%s'%colors[m])
			m+=1

               	ax3.legend(title='Fitorder:')
               	ax3.set_ylabel('CVS')
               	ax4.legend(title='Fitorder:', loc=2)
		ax4.set_ylabel(r'$d^2F/d\epsilon^2$    $Pa$')

               	ax1.plot(eta,F,'o')

               	ax1.plot(polyx,p(polyx), label='6')
               	ax1.plot(polyx,pl(polyx), label='4')
		ax1.legend(title='Fitorder:')
               	#ax2.plot(eta,volume)
	
               	#plt.plot(eps,F1)
               	#plt.plot(eps,F2)
		ax1.set_title('Free Energy - 2stage fitting')
		ax3.set_xlabel(r'$\eta_{max}$')
		ax2.set_ylabel(r'CVS($d^2F/d\eta^2$)')
		
		ax1.set_ylabel('$F$    $eV$')
		ax4.set_xlabel(r'$\eta_{max}$')
              	plt.draw()
	       
	def deriv(self, coeff,order):
		dc=[]
		for i in range(len(coeff)-2):
			dc.append(((order-i)*(order-1-i))*coeff[i])
		return dc


class AnalyzeE(object):

	def __init__(self,distortion='Dst02',scale = 'scale_3.17226',rootdir='/home/t.dengg/results/W/Cij_T_F/Cij_T',path = '/',etamin = -0.05,etamax = 0.05,n = 21,T_index=100,fitorder = 6,delta_x_deriv = 0.001):

               	_e     = 1.602176565e-19 	     # elementary charge
      	       	Bohr   = 5.291772086e-11 	     # a.u. to meter
      	       	Ryd2eV = 13.605698066		     # Ryd to eV
      	       	Angstroem = 1.e-10		     # Angstroem to meter
      	       	ToGPa  = (_e*Ryd2eV)/(1e9*Bohr**3)    # Ryd/[a.u.^3] to GPa
      	       	vToGPa  = (_e)/(1e9*Angstroem**3)   # eV/[Angstroem^3] to GPa



		plt.figure()

      	       	ax1=plt.subplot(221)
      	       	ax2=plt.subplot(222)
      	       	ax3=plt.subplot(223)
      	       	ax4=plt.subplot(224)


      	       	##########################################
      	       	############## MANUAL INPUT ##############
      	       	#distortion='Dst02'
      	       	#scale = 'scale_3.17226'
      	       	#rootdir='/home/t.dengg/results/W/Cij_T_F/Cij_T'
      	       	#path = rootdir+'/'+scale+'/'+distortion+'/'
      	       	#etamin = -0.05
      	       	#etamax = 0.05
      	       	#n = 21
#
      	       	#T_index=100

      	       	#fitorder = 6

      	       	#delta_x_deriv = 0.001
      	        ##########################################
      	       	##########################################

		colors = ['b','g','r','c','m','k','y']



               	F = []
               	F1 = []
               	F2 = []
               	eps = []
               	volume=[]
               	i=0
               	for d in sorted(os.listdir(path)):
                       	if os.path.isdir(path+d):
                       
                       		
				
				cell = 125.#supercell size
                    
                    		vaspout = et.parse(path+d+'/'+'vasprun.xml')
                    		elem = vaspout.xpath("//scstep/energy/i[@name = 'e_fr_energy']")
                    		allengys = []
                    		for k in elem:
                        		try:
                            			allengys.append(float(k.text))
                        		except:
                            			allengys.append(0.)
                    
                    		trueengys = []
                    		for engy in allengys:
                        		if engy < -1000. and engy > -3000.: trueengys.append(engy)
                    		gsenergy = trueengys[-1]/cell
				
				
				
                	       	vasprun = et.parse(path+d+'/'+'vasprun.xml')
                	       	volume.append(vasprun.xpath("//i[@name='volume']")[0].text)
                       
                       
                	       	f=open(path+d+'/'+'F_TV_back')
                	       	data = f.readlines()
                	       	F.append(gsenergy)
                	       	F1.append(float(data[10].split()[1])/96.47244)
                	       	F2.append(float(data[100].split()[1])/96.47244)
                	       	f.close()
                	       	eps.append(i)
                	       	i+=1
                       	else:
                	      	 continue


               	eps = map(float,eps)
               	eta = np.linspace(etamin,etamax,n)

               	coeff = np.polyfit(eta,F,fitorder)
               	p = np.poly1d(coeff)
               	coeffl = np.polyfit(eta,F,4)
               	pl = np.poly1d(coeffl)
               	polyx = np.linspace(min(eta),max(eta),1000)


               	r=0
               	for i in range(len(F)):
               		r = r + (F[i]-p(eta[i]))**2.
               	rms1 = np.sqrt(r/float(len(F)))
               	r=0
               	for i in range(len(F)):
               	        r = r + (F[i]-pl(eta[i]))**2.
               	rms2 = np.sqrt(r/float(len(F)))
               	print "##############################################"
               	print "# --> Calculated root mean square error (rms):"
               	print "#---------------------------------------------"
               	print '# rms (fitorder: %s)'%fitorder,rms1
               	print '# rms (linear fit)',rms2
               	print "##############################################"

		m=0
               	for j in [2,3,4,6,8]:
                       	CVS=[]
                       	ETA=[]
                       	DERIV=[]
                       	DERIVx=[]
                       	strain = list(copy.copy(eta))
                       	energy = list(copy.copy(F))
                       	while (len(strain) > j+1): 
                	       	emax = max(strain)
                	       	emin = min(strain)
                	       	emax = max(abs(emin),abs(emax))
                	       	coeff = np.polyfit(strain,energy,j)
                	       	DERIVx.append(np.polyval(self.deriv(coeff,j),delta_x_deriv)*vToGPa*2)
                	       	DERIV.append(np.polyfit(strain,energy,j)[j-2]*2.*vToGPa*2)
                	       	S = 0
                	       	for k in range(len(strain)):
                		       	Y      = energy[k]
                		       	etatmp = []
                		       	enetmp = []

                		       	for l in range(len(strain)):
                			       	if (l==k): pass
                			       	else:		
                				       	etatmp.append(strain[l])
                				       	enetmp.append(energy[l])

                		       	Yfit = np.polyval(np.polyfit(etatmp,enetmp, j), strain[k])
                		       	S    = S + (Yfit-Y)**2
                       
                	       	CV = np.sqrt(S/len(strain))
                       
                       
                	       	CVS.append(CV)
                	       	ETA.append(emax)
                	       	if (abs(strain[0]+emax) < 1.e-7):
                		       	strain.pop(0)
                		       	energy.pop(0)
                	       	if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                		       	strain.pop()
                		       	energy.pop()
                       	ax3.plot(ETA,CVS,'%s'%colors[m], label='%s'%j)
                       	if j%2==0:
                	       	ax4.plot(ETA,DERIV,'%s'%colors[m],label='%s'%j)
			if j==8:
                       		ax4.plot(ETA,DERIVx, '%s--'%colors[m],label='%s, dx=%s'%(j,delta_x_deriv))
			m+=1
			
		
		m=0	
		for j in [2,3,4,6,8]:
			cvs2 = CVS2.CVS(eta,F)
			CVS = cvs2.calc_cvs(j)
			A1 = [i[1] for i in CVS]
			A2 = [i[0] for i in CVS]
			ax2.plot(A1,A2, '%s'%colors[m])
			m+=1
		
               	ax3.legend(title='Fitorder:')
		ax3.set_ylabel('CVS')
               	ax4.legend(title='Fitorder:', loc=2)
		ax4.set_ylabel(r'$d^2U/d\eta^2$   $Pa$')

               	ax1.plot(eta,F,'o')

               	ax1.plot(polyx,p(polyx), label='fitorder 6')
               	ax1.plot(polyx,pl(polyx), label='fitorder 4')
		ax1.legend(title='Fitorder:')
               	#ax2.plot(eta,volume)

               	#plt.plot(eps,F1)
               	#plt.plot(eps,F2)
		ax1.set_title('Total Energy')
		ax3.set_xlabel(r'$\eta_{max}$')
		ax2.set_ylabel(r'CVS($d^2U/d\eta^2$)')
		ax1.set_ylabel('$U$   $eV$')
		ax4.set_xlabel(r'$\eta_{max}$')
              	plt.draw()
		
	       	
	def deriv(self, coeff,order):
		dc=[]
		for i in range(len(coeff)-2):
			dc.append(((order-i)*(order-1-i))*coeff[i])
		return dc

if __name__=='__main__':

	distortion='Dst02'
	scale = 'scale_3.15226'
	rootdir='/home/t.dengg/results/W/Cij_T_F/Cij_T'
	path = rootdir+'/'+scale+'/'+distortion+'/'
	etamin = -0.06
	etamax = 0.06
	n = 21
	T_index=100
	fitorder = 6
	delta_x_deriv = 0.001
	
	AnalyzeF(distortion=distortion,scale=scale,rootdir=rootdir,path=path,etamin=etamin,etamax=etamax,n=n,T_index=T_index,fitorder=fitorder,delta_x_deriv=delta_x_deriv)
	AnalyzeE(distortion=distortion,scale=scale,rootdir=rootdir,path=path,etamin=etamin,etamax=etamax,n=n,T_index=T_index,fitorder=fitorder,delta_x_deriv=delta_x_deriv)
	AnalyzeFE(distortion=distortion,scale=scale,rootdir=rootdir,path=path,etamin=etamin,etamax=etamax,n=n,T_index=T_index,fitorder=fitorder,delta_x_deriv=delta_x_deriv)
	AnalyzeFE_seperated_fit(distortion=distortion,scale=scale,rootdir=rootdir,path=path,etamin=etamin,etamax=etamax,n=n,T_index=T_index,fitorder=fitorder,delta_x_deriv=delta_x_deriv)
	plt.show()
