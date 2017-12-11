import numpy as np
from copy import copy
try:
    import matplotlib.pyplot as plt
    mpl=True
except:
    mpl=False
import scipy.stats as stat

class CVS(object):
    def __init__(self, xdata, ydata):
        self.__strain = list(xdata)
        self.__energy = list(ydata)
    
    def calc_cvs(self, fitorder):
        
        
        
        strain = copy(self.__strain)
        energy = copy(self.__energy)
        
        self.__cvs = []
        
        while (len(strain) > fitorder+1):
            emax = max(strain)
            emin = min(strain)
            emax = max(abs(emin),abs(emax))
    
            S = 0
            for k in range(len(strain)):
                #Y      = energy[k]
                Y = np.polyfit(strain,energy, fitorder)[fitorder-2]*2.
                etatmp = []
                enetmp = []

                for l in range(len(strain)):
                    if (l==k): pass
                    else:
                        etatmp.append(strain[l])
                        enetmp.append(energy[l])
                
                Yfit = np.polyfit(etatmp,enetmp, fitorder)[fitorder-2]*2.
                #Yfit = np.polyval(np.polyfit(etatmp,enetmp, fitorder), strain[k])
                S    = S + (Yfit-Y)**2
            
            self.__cvs.append((np.sqrt(S/len(strain)),emax,fitorder))
            
    
            if (abs(strain[0]+emax) < 1.e-7):
                strain.pop(0)
                energy.pop(0)
            if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                strain.pop()
                energy.pop()
                
        
        
        return self.__cvs
    
class ANALYTICS(object):
    def __init__(self, strain, E):
        self.__strain = strain
        self.__E = E
        
        
    def pwcvs(self):
        if mpl:
            fig = plt.figure()
            ax1=fig.add_subplot(311)
            ax2=fig.add_subplot(312)
            ax3=fig.add_subplot(313)
        
        strain = self.__strain
        gsenergy = self.__E
        
        pwcvs = []
            
        j=0
        delta=0.
        m=0
        for fitorder in range(2,9,2):
            j=0
            marker=['+','o','d','s','<']
            ###########################################
            strain = copy(strain)
            energy = copy(gsenergy)
            CV_0=[]
            E2nd_0=[]
            while (len(strain) > fitorder+1):
                emax = max(strain)
                emin = min(strain)
                emax = max(abs(emin),abs(emax))
        
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
                    coeff = np.polyfit(etatmp,enetmp, fitorder)
                    Yfit = np.polyval(coeff, strain[k])
                    S    = S + (Yfit-Y)**2
                        
                E2nd_0.append((coeff[fitorder-2],emax,fitorder))   
                CV_0.append((np.sqrt(S/len(strain)),emax,fitorder))
                cv=(np.sqrt(S/len(strain)),emax,fitorder)
        
                if (abs(strain[0]+emax) < 1.e-7):
                    strain.pop(0)
                    energy.pop(0)
                if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                    strain.pop()
                    energy.pop()
            if mpl:
                ax1.plot([c[1] for c in E2nd_0],[c[0] for c in E2nd_0])
                ax2.plot([c[1] for c in CV_0],[c[0] for c in CV_0])
            ###################################################
            CVN=np.zeros((len(strain),100))
            CVmin= np.zeros(len(strain))
            while j < 100:
                
                CV = []
                E2nd=[]
                E_fluct = [e * (np.random.uniform(-1,1) * 0.000001*(1.+j/10) + 1.) for e in gsenergy]
                
                coeff = np.polyfit(strain, gsenergy, fitorder)
                poly = np.poly1d(coeff)
                dpoly = np.poly1d.deriv(poly)
                ddpoly = np.poly1d.deriv(dpoly)
        
            
                coeff_fluct = np.polyfit(strain, E_fluct, fitorder)
                poly_fluct = np.poly1d(coeff_fluct)
                dpoly_fluct = np.poly1d.deriv(poly_fluct)
                ddpoly_fluct = np.poly1d.deriv(dpoly_fluct)
            
                pvol = np.linspace(min(strain),max(strain),100)
                #print np.linalg.norm(poly_fluct(pvol)-poly(pvol))
                delta += np.linalg.norm(dpoly_fluct(pvol)-dpoly(pvol))
                
                
                strain = copy(strain)
                energy = copy(E_fluct)
                a=0
                while (len(strain) > fitorder+1):
                    
                    emax = max(strain)
                    emin = min(strain)
                    emax = max(abs(emin),abs(emax))
        
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
                        coeff = np.polyfit(etatmp,enetmp, fitorder)
                        Yfit = np.polyval(coeff, strain[k])
                        S    = S + (Yfit-Y)**2
                    E2nd.append((coeff[fitorder-2],emax,fitorder))
                    CV.append((np.sqrt(S/len(strain)),emax,fitorder))
                    cv=(np.sqrt(S/len(strain)),emax,fitorder)
        
                    if (abs(strain[0]+emax) < 1.e-7):
                        strain.pop(0)
                        energy.pop(0)
                    if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                        strain.pop()
                        energy.pop()
                              
                    CVN[a][j] = np.sqrt(S/len(strain))
                    a+=1
        
                j+=1
                if mpl:
                    ax1.plot([c[1] for c in E2nd],[c[0] for c in E2nd],marker[m])
                    ax2.plot([c[1] for c in CV],[c[0] for c in CV],marker[m])
            
            
                
            sumC=0.
            index=0
            for point in CVN:
                sumC+=point
                CVmin[index] = max(point)-min(point)
                
                index+=1
            
            m+=1
            
            
            
            CVmin = np.trim_zeros(CVmin)    
            array_var = np.array(CVmin)#list(reversed(CVmin)))
            array_CV = np.array([c[0] for c in CV_0])
            #print array_CV,array_var
            stability = array_CV**(1.)*array_var#/min(array_var)
            #ax3.plot(range(len(CVmin)),array_CV,':')
            #ax3.plot(range(len(CVmin)),array_var,'--')
            ax3.plot(range(len(CVmin)),stability)
        #plt.plot(pvol,)
        #plt.plot(pvol,poly(pvol),'--')
        
        
        
        #plt.show()

    def phist(self, Ndiv=100, Nsample=100, rand_dev=0.00001):
        if mpl:
            fig = plt.figure()
            ax1=fig.add_subplot(221)
            ax2=fig.add_subplot(222)
            ax3=fig.add_subplot(223)
            ax4=fig.add_subplot(224)
        
            fig2 = plt.figure()
            ax5=fig2.add_subplot(211)
            ax6=fig2.add_subplot(212)
        
        eta = self.__strain
        gsenergy = self.__E
        
        pwCVS = []
        pwCVS_val = []
        fordr = []
        fetamax = []
        
        stddev = []
        hist_E2nd=[]
        j=0
        delta=0.
        m=0
        N=0
        E2nd_all = []
        color=['b','g','r','c','m','k']
        for fitorder in range(2,9,2):
            j=0
            marker=['+','o','d','s','<']
            ###########################################
            strain = copy(eta)
            energy = copy(gsenergy)
            CV_0=[]
            E2nd_0=[]
            while (len(strain) > fitorder+1):
                emax = max(strain)
                emin = min(strain)
                emax = max(abs(emin),abs(emax))
        
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
                    coeff = np.polyfit(etatmp,enetmp, fitorder)
                    Yfit = np.polyval(coeff, strain[k])
                    S    = S + (Yfit-Y)**2
                        
                E2nd_0.append((coeff[fitorder-2],emax,fitorder))   
                E2nd_all.append((coeff[fitorder-2],emax,fitorder))
                
                CV_0.append((np.sqrt(S/len(strain)),emax,fitorder))
                cv=(np.sqrt(S/len(strain)),emax,fitorder)
        
                if (abs(strain[0]+emax) < 1.e-7):
                    strain.pop(0)
                    energy.pop(0)
                if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                    strain.pop()
                    energy.pop()
            if mpl:
                ax1.plot([c[1] for c in E2nd_0],[c[0] for c in E2nd_0],color=color[m],label='n=%s'%fitorder,lw=3)
                ax2.plot([c[1] for c in CV_0],[c[0] for c in CV_0],color=color[m],lw=3)
            
                ax5.plot(strain,energy, color=color[m],lw=3)
            ###################################################
            CVN=np.zeros((len(eta),Nsample))
            CVmin= np.zeros(len(eta))
            CVstd = np.zeros(len(eta))
            while j < Nsample:
                
                CV = []
                E2nd=[]
                E_fluct = [e * (np.random.uniform(-1,1) * rand_dev*(1.+j/10) + 1.) for e in gsenergy]
                
                coeff = np.polyfit(eta, gsenergy, fitorder)
                poly = np.poly1d(coeff)
                dpoly = np.poly1d.deriv(poly)
                ddpoly = np.poly1d.deriv(dpoly)
        
            
                coeff_fluct = np.polyfit(eta, E_fluct, fitorder)
                poly_fluct = np.poly1d(coeff_fluct)
                dpoly_fluct = np.poly1d.deriv(poly_fluct)
                ddpoly_fluct = np.poly1d.deriv(dpoly_fluct)
            
                pvol = np.linspace(min(eta),max(eta),100)
                #print np.linalg.norm(poly_fluct(pvol)-poly(pvol))
                delta += np.linalg.norm(dpoly_fluct(pvol)-dpoly(pvol))
                
                
                strain = copy(eta)
                energy = copy(E_fluct)
                a=0
                while (len(strain) > fitorder+1):
                    
                    emax = max(strain)
                    emin = min(strain)
                    emax = max(abs(emin),abs(emax))
        
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
                        coeff = np.polyfit(etatmp,enetmp, fitorder)
                        Yfit = np.polyval(coeff, strain[k])
                        S    = S + (Yfit-Y)**2
                    E2nd.append((coeff[fitorder-2],emax,fitorder))
                    hist_E2nd.append(coeff[fitorder-2])
                    N+=1
                    CV.append((np.sqrt(S/len(strain)),emax,fitorder))
                    cv=(np.sqrt(S/len(strain)),emax,fitorder)
        
                    if (abs(strain[0]+emax) < 1.e-7):
                        strain.pop(0)
                        energy.pop(0)
                    if (abs(strain[len(strain)-1]-emax) < 1.e-7):
                        strain.pop()
                        energy.pop()
                              
                    CVN[a][j] = np.sqrt(S/len(strain))
                    a+=1
                    
                j+=1
                if mpl:
                    ax1.plot([c[1] for c in E2nd],[c[0] for c in E2nd],marker[m],color='y')
                    ax2.plot([c[1] for c in CV],[c[0] for c in CV],marker[m],color='y')
                
                    ax5.plot(strain,energy,marker[m],color='y')
                    
                    ax2.set_xticks([0.01,0.02,0.03,0.04,0.05])
                    
                    ax1.set_xlabel(r'$\eta_{max}$')
                    ax1.set_ylabel(r'$d^2E/d\eta^2$     (GPa)')
                    
                    ax2.set_xlabel(r'$\eta_{max}$')
                    ax2.set_ylabel(r'CVS')
                    ax3.set_xlabel(r'$\eta_{max}$')
                    ax3.set_ylabel(r'p-CVS')
                    ax4.set_xlabel(r'$d^2E/d\eta^2$    (GPa)')
                    ax4.set_ylabel(r'density')
                
            sumC=0.
            index=0
            for point in CVN:
                sumC+=point
                CVmin[index] = (max(point)-min(point))**2.
                CVstd[index] = stat.tstd(point)
                index+=1
            
            
            print "Standard deviation of CVS of fitorder %s: %s"%(fitorder,stat.tstd(CVN))
            print "Standard deviation of d2E/d(eta)2 of fitorder %s: %s"%(fitorder,stat.tstd(hist_E2nd))
            self.standardev = stat.tstd(CVN)
            CVmin = np.trim_zeros(CVmin)   
            CVstd =  np.trim_zeros(CVstd) 
            array_var = np.array(CVmin)#list(reversed(CVmin)))
            array_CV = np.array([c[0] for c in CV_0])
            #print array_CV,array_var
            #stability = array_CV**(1.)*array_var#/min(array_var)
            stability = array_CV**(1.)*CVstd**(1.)
            #ax3.plot(range(len(CVmin)),array_CV,':')
            #ax3.plot(range(len(CVmin)),array_var,'--')
            if mpl:
                ax3.plot([c[1] for c in CV],stability,color=color[m],lw=3)
                ax3.set_xticks([0.01,0.02,0.03,0.04,0.05])
                #ax3.plot([c[1] for c in CV],array_CV*array_var,'--',color=color[m])
            
            stddev.append((fitorder,self.standardev))
            pwCVS.extend(stability)
            
            fordr.extend([fitorder for c in CV])
            fetamax.extend([c[1] for c in CV])
            
            m+=1
        
        
        hist_vect = self.make_hist(Ndiv,hist_E2nd)
        
        deltasum = 0
        for i in range(len(hist_vect)):
            
            if max(hist_vect[i-1],hist_vect[i])!=0: deltasum += abs(hist_vect[i-1]-hist_vect[i])/max(hist_vect[i-1],hist_vect[i])
        
        print deltasum/Ndiv
        
        """
        temp_vec=copy(hist_vect)
        for index in range(4):
            max1 = max(temp_vec)
            i1 = temp_vec.index(max1)
            temp_vec.pop(i1)
            max2 = max(temp_vec)
            i2 = temp_vec.index(max2)
            if abs(i1-i2)>len(hist_vect)/50 and 0.9<(max1/max2)<1.1:
                Ndiv = Ndiv*2
                hist_vect = self.make_hist(Ndiv,hist_E2nd)
                print 'Changeing divisions in histogram method to %s'%(Ndiv)
        """
        
        """
        maxi=[]
        ndiv=[]
        temp_vec=copy(hist_vect)
        max1 = max(temp_vec)
        i1 = temp_vec.index(max1)
        
        k=0
        while temp_vec[i1-1] < temp_vec[i1]*0.95 and temp_vec[i1+1] < temp_vec[i1]*0.95:
            Ndiv=int(Ndiv*1.2)
            temp_vec = self.make_hist(Ndiv,hist_E2nd)
            max1 = max(temp_vec)
            i1 = temp_vec.index(max1)
            print 'Changeing divisions in histogram method to %s'%(Ndiv)
            maxi.append(min(hist_E2nd)+i1*(max(hist_E2nd)-min(hist_E2nd))/(Ndiv+1))
            ndiv.append(Ndiv)
            if k>2 and (maxi[k]-maxi[k-1])<1. and (maxi[k-1]-maxi[k-2])<1.: break 
            k+=1
        ax6.plot(ndiv,maxi)
        if not 19*Ndiv/20<i1<Ndiv/20: 
            Ndiv=Ndiv*2
            temp_vec = self.make_hist(Ndiv,hist_E2nd)
        hist_vect=temp_vec
        
        """
        
        x= np.linspace(min(hist_E2nd),max(hist_E2nd),Ndiv+1)
        if mpl:
            ax4.plot(x, hist_vect,lw=3)
        
        
        #plt.plot(pvol,)
        #plt.plot(pvol,poly(pvol),'--')
        pw_prediction=[]
        for i in range(3):
            pw_ind = pwCVS.index(min(pwCVS))
            
            for el in E2nd_all:
                if fordr[pw_ind]==el[2] and fetamax[pw_ind]==el[1] :
                    pw_prediction.append(el)
                    
            pwCVS.pop(pw_ind)
                    
        hist_prediction  = x[hist_vect.index(max(hist_vect))]
        
        
        dist_x = np.linspace(stat.norm.ppf(0.001, loc=hist_prediction, scale = 0.1),stat.norm.ppf(0.999, loc=hist_prediction, scale = 0.1),100) 
        dist_y = stat.norm.pdf(dist_x, loc=hist_prediction, scale = 0.1)
        if mpl:
            ax4.plot(dist_x,[y*100. for y in dist_y],lw=3)
        
        
        
            ax1.hlines(hist_prediction, 0., max(eta), lw=4., alpha=0.5, color='g')
            ax1.hlines(pw_prediction[0][0], 0., max(eta), lw=4., alpha=0.6, color='r')
            #ax1.hlines(pw_prediction[1][0], 0., max(eta), lw=4., alpha=0.4, color='r')
            #ax1.hlines(pw_prediction[2][0], 0., max(eta), lw=4., alpha=0.2, color='r')
            
            ax4.vlines(hist_prediction, 0., max(hist_vect), lw=4., alpha=0.5, color='g')
            ax4.vlines(pw_prediction[0][0], 0., max(hist_vect), lw=4., alpha=0.6, color='r')
            #ax4.vlines(pw_prediction[1][0], 0., max(hist_vect), lw=4., alpha=0.4, color='r')
            #ax4.vlines(pw_prediction[2][0], 0., max(hist_vect), lw=4., alpha=0.2, color='r')
            
        best_prediction=0.
        for val in pw_prediction:
            if best_prediction==0. or abs(val[0]-hist_prediction)<(best_prediction-hist_prediction):
                best_prediction=val[0]
                best_eta = val[1]
                best_forder = val[2]
                
        print best_prediction,best_eta,best_forder
        print pw_prediction
        if mpl:
            ax1.legend(title='Fitorder')
            #plt.show()
        
        #return (best_eta,best_forder),pw_prediction, hist_prediction, best_prediction, fig, fig2
        return (best_eta,best_forder),pw_prediction, stddev, best_prediction, fig, fig2
                
    def make_hist(self, Ndiv, hist_E2nd):
        hist_vect=np.zeros(Ndiv+1)
        
        min_E2nd = min(hist_E2nd) 
        max_E2nd = max(hist_E2nd)
        d_E2nd = (max_E2nd - min_E2nd)/Ndiv
        hist_indices = [int((p-min_E2nd)/d_E2nd) for p in hist_E2nd]
        hist_vect = list(hist_vect)
        
        for ind in hist_indices:
            hist_vect[ind]+=1
        
        return hist_vect
    