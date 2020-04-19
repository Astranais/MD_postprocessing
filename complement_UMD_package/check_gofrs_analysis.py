#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 16:06:22 2018

@author: akobsch

*** program to plot gofr and int of each pair of atoms of the selected file along with the corresponding coordinence and bond length ***
"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import glob
import os
import matplotlib.pyplot as plt
import numpy as np

def creation_plot(file):
    """     ********** Creation of the plot with two vertical axis  **********    """
    plt.close()
    fig, ax1 = plt.subplots(figsize=(7,7))
    #plt.title("g(r) and Int(g(r)) for a = %1.1f $\mathregular{\AA}$"% float(cell_size),fontsize=30)
    ax1.set_xlabel('Distance ($\AA$)',fontsize=11,fontweight='bold')
    for tl in ax1.get_xticklabels():
        tl.set_fontsize(11)
    for tl in ax1.get_yticklabels():
        tl.set_fontsize(11)
    ax1.set_ylabel("g(r)",fontsize=11,fontweight='bold')
    ax1.grid(True, axis = 'x')

    ax2 = ax1.twinx()
    ax2.set_ylabel("Int(g(r))",fontsize=11,fontweight='bold')
    for tl in ax2.get_yticklabels():
        tl.set_fontsize(11)
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.9, left = 0.1, hspace = 0, wspace = 0)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.title('Pair distribution function for '+file.split('.outcar.gofr.dat')[0] ,fontsize=13,fontweight='bold' )
    return fig,ax1,ax2            
                
            
def main(argv):
    """     ********* Main program *********     """
    colors = ['DarkGreen','MediumOrchid','DodgerBlue','Tomato','Olive','DarkKhaki','Silver' ]
    atoms = []
    allpairs = {}
    allpairs_ordered = []
    data = {}
    directory = os.curdir
    try:
        options,arg = getopt.getopt(argv,"hg:a:d:",["gofrsfile","atoms","directory"])
    except getopt.GetoptError:
        print("check_gofr_analysis.py -g <gofrs.txt filename> -a <pair of atoms>(ex: 'Ca-O,Ca-Ca') -d <directory where the gofr.dat files are located (option)>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('check_gofr_analysis.py program to plot g(r) and int of each pair of atoms of the selected file along with the corresponding coordinence and bond length')
            print("check_gofr_analysis.py -g <gofrs.txt filename> -a <pair of atoms>(ex: 'Ca-O,Ca-Ca') -d <directory where the gofr.dat files are located (option but required if you have different gofr.dat file with same name in different subfolders)>")
            print('')
            sys.exit()
        elif opt in ("-g", "--gofrsfile"):
            gofrfile = str(arg)
        elif opt in ("-a", "--atoms"):
            atoms = arg.split(',')                      #list of atom pairs we want to analyze here)
        elif opt in ('-d','--directory'):
            directory = str(arg)
    #********** 1st step: extraction of pairs and data from the gofrs.txt file
    #read the header and extract the pairs
    skip_head = 0
    with open(gofrfile, 'r') as f:
        while True:
            line = f.readline()
            if not line: break
            else:
                entry = line.split('\n')[0].split('\t')
                skip_head +=1
                if entry[0] == 'pair':
                    allpairs_ordered = [entry[1]]        #we initialize the pairs list
                    for ii in range(2,len(entry)):       #we fill the pairs list with atom pairs occuring only once
                        if entry[ii] != entry[ii-1]:
                            allpairs_ordered.append(entry[ii])
                    for ii in range(len(allpairs_ordered)):
                        allpairs[allpairs_ordered[ii]] = ii #and we add the index for each pair in the dictionnay allpairs with atom pairs as key
                elif entry[0] == 'file':
                    break
    #fill the data dictionnary with all the data from gofrs
    with open(gofrfile, 'r') as f:
        [f.readline() for ii in range(skip_head)]
        while True:
            line = f.readline()
            if not line: break
            else:
                entry = line.split('\n')[0].split('\t')
                data[entry[0].split('/')[-1]] = entry[1:]
    print('all pairs available are ', allpairs)
    print("we analyze", atoms)
    #********** 2st step: for each gofr.dat file in the directories below the directory indicated we test if they are in the gofr.txt file
    # if yes, then we plot the corresponding graph
    # if no we move to the next one
    for dirpath, dirnames, filenames in os.walk(directory):
        files = sorted(glob.glob(dirpath+'/*.gofr.dat')) #I list every gofr files in alphabetic order
        print('in folder ', dirpath, 'the matching files are', files)
        if files != []:
            for file in files:
                if file.split('/')[-1] in data:
                    print(file)
                    #creation of the plot
                    fig,ax1,ax2 = creation_plot(file)
                    #extraction of data
                    index = 0        #used for colors
                    for pair in allpairs:
                        for atom_pair in atoms:
                            if (pair == atom_pair):
                                print('plot', pair)
                                #we load the data of the gofr.dat files
                                if pair == 'O2':
                                    distance, gofr, int_gofr = np.loadtxt(file, usecols=(0, allpairs['O-O']*2 + 1, allpairs['O-O'] *2 + 2), skiprows=1, unpack=True) 
                                else:
                                    distance, gofr, int_gofr = np.loadtxt(file, usecols=(0, allpairs[pair]*2 + 1, allpairs[pair] *2 + 2), skiprows=1, unpack=True) 
                                distance = np.delete(distance,[len(distance)-1]) #we remove the last 0 value
                                gofr = np.delete(gofr,[len(gofr)-1])   #we remove the last 0 value
                                int_gofr = np.delete(int_gofr,[len(int_gofr)-1])   #we remove the last 0 value
                                #we load the computed data from the gofrs.txt file
                                try:
                                    xmax = float(data[file.split('/')[-1]][allpairs[pair]*5])
                                    ymax = float(data[file.split('/')[-1]][allpairs[pair]*5 + 1])
                                    xmin = float(data[file.split('/')[-1]][allpairs[pair]*5 + 2])
                                except ValueError:
                                    xmax = float('nan')
                                    ymax = float('nan')
                                    xmin = float('nan')
                                #we plot the data along with lines indicating the bond distances and coordination (obtained from the gofrs_XXXX.txt file)
                                ax1.plot(distance, gofr, marker = '.',linestyle = '-', linewidth=3, color = colors[index], label = pair)
                                ax1.vlines(x=xmax,ymin  = 0, ymax = ymax+ymax*0.1, linestyle = '-.',linewidth=2,  color = colors[index])
                                ax2.plot(distance,int_gofr, marker = 'None',linestyle = ':', linewidth=2, color = colors[index])
                                ax1.vlines(x=xmin,ymin  = 0, ymax = ymax+ymax*0.1, linestyle = '--', linewidth=2, color = colors[index])
                                if pair == 'O2':
                                    ax1.set_ylim([0,ymax+ymax*0.1])
                                    ax1.set_xlim([0,2.5])
                                index+=1      #used to change of colors
                    ax1.legend(loc = 'upper left',fontsize=15)      
                    plt.savefig(file.split('.outcar.gofr.dat')[0]+'.png', bbox_inches='tight', dpi=90)
    #                plt.show() # à mettre après savefig sinon on ne sauvegarde rien !!  

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])
 
