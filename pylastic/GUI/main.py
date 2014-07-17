from ttk import Frame, Label, Style, Notebook, Treeview
from Tkinter import Toplevel, Label, Tk, TOP, Y, BOTH, Listbox, Entry, StringVar, END, Button, Radiobutton, IntVar, Scrollbar, NSEW, HORIZONTAL, VERTICAL, W, NS, EW 
from tkFileDialog import askopenfilename



import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import os
import sys
#sys.path.append('/home/t.dengg/git/pylastic/pylastic')
#sys.path.append('../pylastic')

import pickle
import numpy as np
import lxml.etree as et

from pylastic.vaspIO import POS
from pylastic.elatoms import Structures, ElAtoms
from pylastic.postprocess import ECs


class pylasticGUI(Frame):
  
    def __init__(self):
        
        self.__fileName = 'POSCAR'
        self.__NoP = 11
        self.__strain = 0.05
        self.__executable = '/home/t.dengg/bin/vasp/vasp.5.3/vasp'
        self.__status = None
        self.__finished = None
        self.__paths = None
        self.__structuresinst = None
        self.__p1 = None
        self.__p2 = None
        self.__workdir = os.getcwd()
        Frame.__init__(self, name = 'pylastic')   
        self.pack(expand=Y, fill=BOTH)
        self._create_panel()
            
        
    def _create_panel(self):
        
        panel = Frame(self, name='elastic')
        panel.pack(side=TOP, fill=BOTH, expand=Y)
        
        nb = Notebook(panel, name='notebook')
        
        nb.enable_traversal()
        nb.pack(fill=BOTH, expand=Y, padx=2, pady=3)
        self._create_setup_tab(nb)
        self._create_analyze_tab(nb)
        
    def _create_setup_tab(self, nb):
        frame1 = Frame(nb, name='entry1')
        
        frame = Frame(nb, name='setup')
        
        self.e1 = Entry(frame)
        self.e1.bind()
        self.e1.pack()
        self.e1.delete(0, END)
        self.e1.insert(0, '21')
        
        self.label1 = Label(frame1, text='Lagrangian strain.')
        self.label1.pack()
        
        self.e2 = Entry(frame)
        self.e2.bind()
        self.e2.pack()
        self.e2.delete(0, END)
        self.e2.insert(0, '0.05')
        
        
        #Choose inputfile (Button)
        b2 = Button(frame, text="Inputfile", width=10, command=self.select_file)
        b2.pack()
        
        b = Button(frame, text="set values", width=10, command=self.read_val)
        b.pack()
        
        #Read structure (Button)
        b1 = Button(frame, text="Read structure", width=15, command=self.read_POS)
        b1.pack()
        
        #Write structures (Button)
        b3 = Button(frame, text="Setup calculation", width=15, command=self.setup_calc)
        b3.pack()
        
        #Strart calculations
        b4 = Button(frame, text="Calculation status", width=15, command=self.check_calc)
        b4.pack()
        
        
        v = IntVar()
        Radiobutton(frame, text="Energy", variable=v, value=1).pack()
        Radiobutton(frame, text="Strain", variable=v, value=2).pack()
        
        w = IntVar()
        Radiobutton(frame, text="2nd", variable=w, value=1).pack()
        Radiobutton(frame, text="3rd", variable=w, value=2).pack()
        
        nb.add(frame, text='Setup', underline=0, padding=2)
        
        return
    
    
    def _create_analyze_tab(self, nb): 
        frame = Frame(nb, name='analyze')
        
        self._create_treeview(frame, )
        root = self.tree.insert('', END, text = 'calc', values=['text1','text2'])
        
        
        
        #self._populate_tree(root)
        
        #
        b1 = Button(frame, text="Refresh", width=15, command=lambda: self._populate_tree(root))
        b1.pack()    
        
        b2 = Button(frame, text="Analyze", width=15, command=lambda: self.analyze(frame))
        b2.pack()
        
        b2 = Button(frame, text="Result", width=15, command=lambda: self.create_window())
        b2.pack()
        
        
        if self.__status and '--------' in self.__status: self.__finished = False
        else: self.__finished = True
        
        if self.__status and self.__finished:
            ECs.set_analytics()
        
        
        

        
        
        nb.add(frame, text='Analyze', underline=0, padding=2)
        return
    
    def _populate_tree(self, root):
        
        print self.tree.get_children()
        if len(self.tree.get_children(self.tree.get_children())) > 1: 
            for i in self.tree.get_children(self.tree.get_children()): self.tree.delete(i)
        #self.tree.delete(self.tree.get_children())
        status, paths = self.check_calc()
        
        for key in status.keys():
            #print paths[key].lstrip(self.__workdir).split('/')
            self.tree.insert(root, END, text = status[key], values=['text3','text4', paths[key]])
            
        
    
    def _create_treeview(self, parent):
        f = Frame(parent)
        f.pack(side=TOP, fill=BOTH, expand=Y)
        
        # create the tree and scrollbars
        self.dataCols = ('fullpath', 'type', 'status')       
        self.tree = Treeview(columns=self.dataCols,
                                 displaycolumns='status')
        
        ysb = Scrollbar(orient=VERTICAL, command= self.tree.yview)
        xsb = Scrollbar(orient=HORIZONTAL, command= self.tree.xview)
        self.tree['yscroll'] = ysb.set
        self.tree['xscroll'] = xsb.set
        
        # setup column headings
        self.tree.heading('#0', text='Directory Structure', anchor=W)
        self.tree.heading('status', text='Status', anchor=W)
        self.tree.column('status', stretch=0, width=70)
        
        # add tree and scrollbars to frame
        self.tree.grid(in_=f, row=0, column=0, sticky=NSEW)
        ysb.grid(in_=f, row=0, column=1, sticky=NS)
        xsb.grid(in_=f, row=1, column=0, sticky=EW)
        
        # set frame resizing priorities
        f.rowconfigure(0, weight=1)
        f.columnconfigure(0, weight=1)
        
        # action to perform when a node is expanded
        self.tree.bind('<<TreeviewOpen>>', self._update_tree)
        
        self.tree.bind("<Double-1>", self.OnDoubleClick)
        
    def OnDoubleClick(self, event):
        item = self.tree.identify('item',event.x,event.y)
        if self.tree.item(item,"text") == 'calc':
            self.create_window()
        
        
    def _update_tree(self, event):
        # user expanded a node - build the related directory
        nodeId = self.tree.focus()      # the id of the expanded node
        
        if self.tree.parent(nodeId):    # not at root
            topChild = self.tree.get_children(nodeId)[0]
            
            # if the node only has a 'dummy' child, remove it and
            # build new directory; skip if the node is already
            # populated
            if self.tree.item(topChild, option='text') == 'dummy':
                self.tree.delete(topChild)
                path = self.tree.set(nodeId, 'fullpath')
                self._populate_tree(nodeId, path, os.listdir(path))
    
    def create_window(self):
        
        t = Toplevel(self)
        t.wm_title("Energy - Strain")
        #l = Label(t, text="This is window #%s" % self.counter)
        #l.pack(side="top", fill="both", expand=True, padx=100, pady=100)

    
    def read_val(self):
        self.__NoP = float(self.e1.get())
        self.__strain = float(self.e2.get())
        
    def read_POS(self):
        self.__poscar = POS(self.__fileName).read_pos()
        #print self.__poscar
        
    def read_strain(self):
        self.__strain = self.e2.get()
        
    def select_file(self):
        self.__fileName = askopenfilename()
        print self.__fileName
        self.__workdir = self.__fileName.rstrip('POSCAR')
        print self.__workdir
        
    def setup_calc(self):
        
        ###################### Create Structures instance: ###################
        self.__structuresinst = Structures()
        
        ## Generate distorted structures and add them to structures object: ##
        atom = ElAtoms()
        print self.__poscar
        atom.poscarToAtoms(self.__poscar)
        for etan in np.linspace(-self.__strain, self.__strain, self.__NoP):
            for strains in range(len(atom.strainList)):
                atom = ElAtoms()
                atom.poscarToAtoms(self.__poscar)
                atom.distort(eta=etan, strainType_index = strains)
                self.__structuresinst.append_structure(atom)
        
        ####################### Write vasp input files: #######################
        self.__structuresinst.workdir = self.__workdir
        self.__structuresinst.write_structures()
        
        #################### Start local vasp calculation: ####################
        if self.__executable:
            self.__structuresinst.executable = self.__executable
            self.__structuresinst.calc_vasp()
    
    def check_calc(self):
        self.__status = {}
        self.__paths = {}
        if not self.__structuresinst:
            try:
                with open(self.__workdir+'/structures.pkl', 'rb') as input:
                    atoms = pickle.load(input)
            except:
                raise Exception("Please do setup first!")
                
        else:
            atoms = self.__structuresinst.get_structures()
        
        for stype, eta in atoms:
            #print stype, eta, atoms[(stype,eta)].path
            
            try:
                et.parse(atoms[(stype,eta)].path+'/vasprun.xml')
                self.__status[(stype,eta)] = 'finished'
            except:
                self.__status[(stype,eta)] = '--------'
            self.__paths[(stype,eta)] = atoms[(stype,eta)].path
        
        return self.__status, self.__paths
        
    
    def analyze(self, frame):
        if self.__status and '--------' in self.__status: self.__finished = False
        else: self.__finished = True
        
        if self.__status and self.__finished:
        
            ec = ECs()
            ec.set_structures()
            ec.set_gsenergy()
            ec.set_analytics()
            
            self.__p1 = ec.plot_cvs()
            self.__p2 = ec.plot_2nd()
            canvas = FigureCanvasTkAgg(self.__p1, master=frame)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            canvas = FigureCanvasTkAgg(self.__p2, master=frame)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        


if __name__ == '__main__':
    pylasticGUI().mainloop()


