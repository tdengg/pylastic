import matplotlib.pyplot as plt
import os
import collections
import numpy as np
import lxml.etree as et

class PLOT_THERMAL(object):
	def __init__(self, dir, style='-'):
		self.dir = dir
		self._lname = None
		self.style = style
		self.conv = 96.47244
		self.__E0 = []
		self.__cellsize = 125.
	def get_data(self):
		f = []
		l = []
		dic = {}
		
		for d in sorted(os.listdir(self.dir)):
			
			if 'scale' in d and os.path.isdir(d):
				
				l = float(d.lstrip('scale_'))
				
				try:
					vrun = et.parse(self.dir+'/'+d+'/vasprun.xml')
					dic[l] = {}
					dic[l]['file'] = open(self.dir+'/'+d+'/F_TV').readlines()
					dic[l]['vasprun'] = vrun
					
					elem = vrun.xpath("//scstep/energy[i/@name='hartreedc']/i[@name = 'e_0_energy']")
					
					gsenergy = float(elem[1].text)
					dic[l]['E0']=gsenergy/self.__cellsize
					self.__E0.append(gsenergy/self.__cellsize)
					
					
				except:
					continue
				
				
			else:
				continue
		
		
		j = 0
		for out in sorted(dic):
			
			T = []
			F = []
			for i in dic[out]['file']:
				T.append(float(i.split()[0]))
				F.append(float(i.split()[1]))
				
			dic[out]['F'] = F
			dic[out]['T'] = T
			dic[out]['E0'] = float(dic[out]['vasprun'].xpath("//scstep/energy[i/@name='hartreedc']/i[@name = 'e_0_energy']")[1].text)
			#plt.plot(T,F)
			#plt.show()
			j+=1
		self.__dic = dic
		return dic
	
	def free_ene(self, trange):
		self.get_data()
		#self.get_E0()
		ax = plt.subplot(111)
		ndic = collections.OrderedDict(sorted(self.__dic.items()))	#sort dictionary
		minF = []
		minl = []
		for temp in trange:
			xdata = []
			ydata = []
			i=0
			
			for out in sorted(ndic):
				xdata.append(out)
				ind = self.__dic[out]['T'].index(temp)
				
				ydata.append(self.__dic[out]['F'][ind]/self.conv + self.__E0[i] + 13.5)
				i+=1
			#polyfit:
			coeff = np.polyfit(xdata,ydata,5)
			p = np.poly1d(coeff)
			polyx = np.linspace(min(xdata),max(xdata),1000)
			
			i=0
			for root in np.roots(p.deriv()):
				if np.imag(root) == 0. and min(xdata)*0.8 < np.real(root) < max(xdata)*1.03:
					minl.append(float(np.real(root)))
					minF.append(p(np.real(root)))
					break
				i+=1
			#minl.append(np.real(np.roots(p.deriv())[1]))
			#minF.append(p(np.roots(p.deriv())[1]))
			
			#polyfit F-T
			coeff = np.polyfit(minl,minF,21)
			p = np.poly1d(coeff)
			polyx = np.linspace(min(minl),max(minl),1000)
			
			#ax.plot(polyx,p(polyx))
			#ax.plot(minl,minF,'o')
			
		
		
		
		#polyfit thermal expansion:
		
		coeff = np.polyfit(trange,minl,3)
		p = np.poly1d(coeff)
		polyx = np.linspace(min(trange),max(trange),1000)
		
		#ax1 = plt.plot(polyx,p(polyx))
		#ax1 = plt.plot(trange,minl,'o')
		#plt.show()
		#thermal expansion/T
		
		alpha = []
		for i in range(len(trange)-1):
			alpha.append((minl[i+1]-minl[i])/(trange[i+1]-trange[i])/minl[0]*10**6.)
		plt.plot(trange[:-1], alpha, label = self._lname, lw=3., ls=self.style)
		self.__alpha = alpha
		self.__minl = minl
		self.__minF = minF
		
	def set_cellsize(self, cellsize):
		self.__cellsize = cellsize
		
	def get_cellsize(self):
		return self.__cellsize

	def set_minl(self, minl):
		self.__minl = minl
	def get_minl(self):
		return self.__minl
	
	def set_minF(self, minF):
		self.__minF = minF
	def get_minF(self):
		return self.__minF
	
	def set_alpha(self, alpha):
		self.__alpha = alpha
	def get_alpha(self):
		return self.__alpha
	
	@property
	def lname(self):
		return self._lname
	@lname.setter
	def lname(self, value):
		self._lname = value
		
	minl = property(fget=get_minl,fset=set_minl)
	minF = property(fget=get_minF,fset=set_minF)
	alpha = property(fget=get_alpha,fset=set_alpha)
	cellsize = property(fget=get_cellsize,fset=set_cellsize)
		
			
