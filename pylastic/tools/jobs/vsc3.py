import os
import time
import threading
import pickle
import Queue as queue
import subprocess

class threads(object):
    def __init__(self, structfile):
        self.__structfile = structfile
        self.__kstrpath = 'kstr'
        self.__kstrname = 'NiTi.prn'
        self.__shapepath = 'shape'
        self.__kgrnpath = 'kgrn'
        self.__kfcdpath = 'kfcd'
        self.__currpath = None
        return
    
    def submit_kstr(self):
        subprocess.Popen(['sbatch {0}/run_kstr'.format(self.__currpath)])
        return
    
    def checkstatus(self,path,subpath,fname,q):

    
        Finished = False
        
        
        
        while not Finished:
            time.sleep(5)
            if os.path.exists('{0}/prn/{1}'.format(self.__currpath,fname)) and os.stat('{0}/prn/{1}'.format(self.__currpath,fname)).st_size != 0:
                f=open('{0}/prn/{1}'.format(self.__currpath,fname))
                if 'Finished' in f.readlines()[-1].split(): 
                    Finished=True
                    print '{0} FINISHED'.format(self.__currpath)
                    
                f.close()
        q.put((path,'done'))
        return True
    
    def kstr(self):
        self.submit_kstr()
    
    def start_jobs(self):
        f = open(self.__structfile)
        structobj = pickle.load(f)
        f.close()
        
        paths=[]
        for struct in structobj.get_structures().values():
            paths.append(struct.path.split('/')[-2]+'/'+struct.path.split('/')[-1]+'/')
        
        q=queue.Queue() 
        t=[]   
        for path in paths:
            self.__currpath = '{0}/{1}'.format(path,self.__kstrpath)
            self.submit_kstr()
            
            t.append(threading.Thread(target=self.checkstatus, args=(path,self.__kstrpath, self.__kstrname, q)))
        for task in t: 
            task.deamon=True
            task.start()
        q.join()
        print 'FINISHED --> PROCEED'
            
        
#t = threading.Thread(target=subproc, args=(a,))
#t.start()


