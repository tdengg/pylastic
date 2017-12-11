
import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    import yaml
except ImportError:
    print("You need to install python-yaml.")
    sys.exit(1)
    
try:
    from yaml import CLoader as Loader
    from yaml import CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from phonopy.units import VaspToTHz

class PlotBand(object):

    def __init__(self, fnames, co=['r-','r--', 'b-', 'g-', 'c-', 'c--', 'm-', 'm--','k-', 'k--'], factor=[1.,1.,1.,1.,1.,1.,1.,1.], mode_gamma=False, a1=3.1356, a2=3.1556):
        self.filenames = fnames
        
        colors = co


        count = 0
        dist_i=0
        if mode_gamma:
            
            FREQ=[]
            DIST=[]
            V2 = (a1)**3./2.*10**(-10.)
            V1 = (a2)**3./2.*10**(-10.)
            delta = V1-V2
            plt.ylabel(r'$\gamma$')
        m=0
        for i, filename in enumerate(self.filenames):
            if mode_gamma:
                FREQ.append([])
                DIST.append([])
            (distances,
             frequencies,
             segment_nqpoint,
             labels) = self.read_band_yaml(filename)
            if i>0: distances=dist_i
            end_points = [0,]
            for nq in segment_nqpoint:
                end_points.append(nq + end_points[-1])
            end_points[-1] -= 1
            segment_positions = distances[end_points]
        
            if all(x is None for x in labels):
                labels_at_ends = None
            else:
                labels_at_ends = [r"$%s$"%labels[n] for n in end_points]
        
            
            
            q = 0
            for j, nq in enumerate(segment_nqpoint):
                if mode_gamma:
                    FREQ[i].append(frequencies[q:(q + nq)])
                    DIST[i].append(distances[q:(q + nq)])
                else:
                    if j == 0:
                        plt.plot(distances[q:(q + nq)],
                                 frequencies[q:(q + nq)]*factor[i],
                                 colors[i])
                    else:
                        plt.plot(distances[q:(q + nq)],
                                 frequencies[q:(q + nq)]*factor[i],
                                 colors[i])
                q += nq
        
            plt.xlim(distances[0], distances[-1])
            if not mode_gamma:
                plt.plot([],[],colors[i], label=filename)
                plt.ylim(ymin=0.)
                plt.ylabel('Frequency  (THz)')
            
            xticks = segment_positions
            plt.xticks(xticks, labels_at_ends)
            for v in segment_positions[1:-1]:
                plt.axvline(x=v, linewidth=1, color='k', ls='--')
                
            self.seg_pos = segment_positions
            
                
            if i==0: dist_i = distances
        
            if mode_gamma and (i+1)%2==0:
                
                vec_d = DIST
                
                vec_f1 = [V1/f1[:,0]*(f1[:,0]-f2[:,0])/(delta) for f1,f2 in zip(FREQ[1+2*m], FREQ[0+2*m])]
                for k,vec in enumerate(vec_f1[0]): 
                    if vec<=-10. or vec>=10.: 
                        try:
                            vec_f1[0][k]=vec_f1[0][k+1]
                        except:
                            vec_f1[0][k]=vec_f1[0][k-1]
                for k,vec in enumerate(vec_f1[2]): 
                    if vec<=-10. or vec>=10.: vec_f1[2][k]=vec_f1[2][k-1]
                for k,vec in enumerate(vec_f1[3]): 
                    if vec<=-10. or vec>=10.: vec_f1[3][k]=vec_f1[3][k+1]
                    
                vec_f2 = [V1/f1[:,1]*(f1[:,1]-f2[:,1])/(delta) for f1,f2 in zip(FREQ[1+2*m], FREQ[0+2*m])]
                for k,vec in enumerate(vec_f2[0]): 
                    if vec<=-10. or vec>=10.: 
                        try:
                            vec_f2[0][k]=vec_f2[0][k+1]
                        except:
                            vec_f2[0][k]=vec_f2[0][k-1]
                for k,vec in enumerate(vec_f2[2]): 
                    if vec<=-10. or vec>=10.: vec_f2[2][k]=vec_f2[2][k-1]
                for k,vec in enumerate(vec_f2[3]): 
                    if vec<=-10. or vec>=10.: vec_f2[3][k]=vec_f2[3][k+1]
                    
                vec_f3 = [V1/f1[:,2]*(f1[:,2]-f2[:,2])/(delta) for f1,f2 in zip(FREQ[1+2*m], FREQ[0+2*m])]
                for k,vec in enumerate(vec_f3[0]): 
                    if vec<=-10. or vec>=10.: 
                        try:
                            vec_f3[0][k]=vec_f3[0][k+1]
                        except:
                            vec_f3[0][k]=vec_f3[0][k-1]
                for k,vec in enumerate(vec_f3[2]): 
                    if vec<=-10. or vec>=10.: vec_f3[2][k]=vec_f3[2][k-1]
                for k,vec in enumerate(vec_f3[3]): 
                    if vec<=-10. or vec>=10.: vec_f3[3][k]=vec_f3[3][k+1]
                    
                plt.plot(DIST[0][0],vec_f1[0], colors[m])
                plt.plot(DIST[0][1],vec_f1[1], colors[m])
                plt.plot(DIST[0][2],vec_f1[2], colors[m])
                plt.plot(DIST[0][3],vec_f1[3], colors[m])
                
                plt.plot(DIST[0][0],vec_f2[0], colors[m])
                plt.plot(DIST[0][1],vec_f2[1], colors[m])
                plt.plot(DIST[0][2],vec_f2[2], colors[m])
                plt.plot(DIST[0][3],vec_f2[3], colors[m])
                
                plt.plot(DIST[0][0],vec_f3[0], colors[m])
                plt.plot(DIST[0][1],vec_f3[1], colors[m])
                plt.plot(DIST[0][2],vec_f3[2], colors[m])
                plt.plot(DIST[0][3],vec_f3[3], colors[m])
                
                m+=1
    
    def read_band_yaml(self, filename):
        data = yaml.load(open(filename), Loader=Loader)
        frequencies = []
        distances = []
        labels = []
        for j, v in enumerate(data['phonon']):
            if 'label' in v:
                labels.append(v['label'])
            else:
                labels.append(None)
            frequencies.append([f['frequency'] for f in v['band']])
            distances.append(v['distance'])
    
        return (np.array(distances),
                np.array(frequencies),
                data['segment_nqpoint'],
                labels)

    def read_dos_dat(self,filename):
        dos = []
        frequencies = []
        for line in open(filename):
            if line.strip()[0] == '#':
                continue
            ary = [float(x) for x in line.split()]
            frequencies.append(ary.pop(0))
            dos.append(ary)
        return np.array(frequencies), np.array(dos)











    




  
