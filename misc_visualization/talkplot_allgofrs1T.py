#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 20 Sep 2018

@author: anais
Langage : Python3

                    ******* Plot every gofr of one compound *********
                                at one  T
"""
#*********** Importation of the packages and modules used here ************
import sys
import getopt
import re
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import crystallography as cr
import natsort
from matplotlib.lines import Line2D #useful to create a custom legend
from matplotlib.colors import LinearSegmentedColormap





def split_name(filename):
    """ Function to extract T and acell from filename """
    # *********** My filename are in the format CaAl2Si2O8_T3_nvt_a12.0.outcar.gofr.dat
    # ******* so I can split their name with _ and tale the acell and T from their name
    # **** Then change these two lines in harmony with your filename
    temperature = filename.split('_')[1]                                     #I extract the corresponding temperature
    acell = filename.split('.outcar.gofr.dat')[0].split('_')[3].strip('a')   #I extract the corresponding acell
    return temperature, acell



def creation_plot(temperature):
    """     ********** Creation of the plot  **********    """
    plt.close()
    fig, ax1 = plt.subplots(figsize=(14,9))

    ax1.set_xlabel('Distance r ($\AA$)',fontsize=25,fontweight='bold', labelpad = 10)
    ax1.set_ylabel('g(r)', fontsize=25, fontweight = 'bold', labelpad = 15)
    
#    ax1.grid(True, axis = 'x')
    major_xticks = np.arange(0, 10, 1)
    minor_xticks = np.arange(0, 10, 0.5)
    ax1.set_xticks(major_xticks)
    ax1.set_xticks(minor_xticks, minor=True)

    ax1.xaxis.set_tick_params(labelsize = 20, width = 2)
    ax1.yaxis.set_tick_params(labelsize = 20, width = 2)
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  #set one decimal to ticks label``        

    #addition of a zoom
    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    left, bottom, width, height = [0.5, 0.63, 0.25, 0.25]
    ax2 = fig.add_axes([left, bottom, width, height])    

    major_xticks_zoom = np.arange(0, 3, 0.5) 
    minor_xticks_zoom = np.arange(0, 3, 0.1)
    ax2.set_xticks(major_xticks_zoom)
    ax2.set_xticks(minor_xticks_zoom, minor=True)
    ax2.set_xticklabels(major_xticks_zoom)      #to choose to display only major tick labels
    major_yticks_zoom = np.arange(0, 0.6, 0.1) 
    minor_yticks_zoom = np.arange(0, 0.6, 0.05)
    ax2.set_yticks(major_yticks_zoom)
    ax2.set_yticks(minor_yticks_zoom, minor=True)                                           

    ax2.xaxis.set_tick_params(labelsize = 20, width = 2)
    ax2.yaxis.set_tick_params(labelsize = 20, width = 2)
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  #set one decimal to ticks label``        


    ax2.set_xlim([0.5,2.5])
    ax2.set_ylim([0,0.2])
#    ax2.grid(True, which='both',axis = 'x', )
    
    return fig,ax1,ax2        
                



def main(argv):
    atoms = []
    lines = {}   #dictionnary for the lines for legend
    Na=6.022*10**23  #avogadro constant
    #other dictionnaries and parameters for the figure
    colors_densities = {}
    """     ********* Main program *********     """
    try:
        options,arg = getopt.getopt(argv,"hm:a:t:",["mineralfile","atoms","temperature"])
    except getopt.GetoptError:
        print("plot_allgofrs1T.py  -m <mineralfile>  -a <1 couple of atoms>(ex: 'O-O') -t <temperature, ex: T4>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print("plot_allgofrs1T.py script to plot all the gofrs of the selected couple of atoms at every acell for 1 T with an inset zommed in the 0-2$\AA$ region")
            print("plot_allgofrs1T.py -m <mineralfile>  -a <1 couple of atoms>(ex: 'O-O') -t <temperature, ex: T4>")
            print(" adapted to show O2 molecules !")
            print("")
            sys.exit()
        elif opt in ("-a", "--atoms"):
            atoms = arg.split(',')                      #list of atom couples we want to analyze here
        elif opt in ("-m","--mineralfile"):
            mineralfile = str(arg)
        elif opt in ("-t","--temperature"):
            temperature = str(arg)
    #***** Calculation of the molecular mass
    with open(mineralfile,'r') as mf:
        entry = mf.readline()
        elements = entry.split()[1:]
        entry = mf.readline()
        number = entry.split()[1:]
    MN = 0
    for i in range(len(elements)):
        MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]  #I compute the molecular mass
    #***** Count of the number of densities (needed for the automatic color change)
    #loop throught all files of matching temperature
    files = sorted(glob.glob('T*/*_'+temperature+'*outcar.gofr.dat')) #I list every gofr files in alphabetic order
    print(files)
    files.reverse()
    #first store every density into the dictionary
    for file in files:
        #calculation of densities in g/cm3
        temperature, acell = split_name(file)  
        Volume = float(acell)**3
        Density = MN/(Na*Volume*10**(-24)) 
        #we store the densities into the dictionnary
        colors_densities[round(Density,1)] = []
    #now the dictionary has every possible densities of our mineral, we attribute the colors to each density
    name = 'lajolla'
    cm_data = np.loadtxt("/Users/akobsch/Dropbox/Recherche/PhD-MD-silicates/simulations/scripts/Analysis/ScientificColourMaps6/"+name+"/"+name+".txt")
    #cm_data = cm_data[::-1] #for reverse colors (hawaii,lajolla)
    new_map = LinearSegmentedColormap.from_list('new', cm_data[40:]) #cm_data[:-20] for devon, davos, oslo
    color = iter(new_map(np.linspace(0,1,len(colors_densities)))) #Creation of the color list  
    list_densities = []
    for key in natsort.natsorted(colors_densities):
        c = next(color)
        colors_densities[key] = c
        list_densities.append(key)   
    #now the dictionary has every possible densities of our mineral, we attribute the colors to each density
    #list_densities = []
    #color = iter(plt.cm.rainbow(np.linspace(0,1,len(colors_densities)))) #Creation of the color list
    #for key in natsort.natsorted(colors_densities):
    #    c = next(color)
    #    colors_densities[key] = c
    #    list_densities.append(key) #list usefull if we want to plot the value of densities near to the lines, see plot_speciation_r0
    print(list_densities, len(list_densities))
    #*****************************
    #*******************
    #*******
    #Read & plot 
    fig, ax1, ax2 = creation_plot(temperature)    
    firstfile = files[0]                                                      #I take the first file  
    temperature0, acell0 = split_name(firstfile)
    with open(firstfile, 'r') as f:
        couples = f.readline()
        couples = couples.strip('dist')
        couples = re.sub('(Int\([A-Za-z-]*\))', ' ',couples).split()        #I substitute all Int(...)  by a blankspace and split using spaces
    col = 0
    for couple in couples:
        for atom_couple in atoms:
            if (couple == atom_couple):
                for file in files:
                    temperature, acell = split_name(file) 
                    Volume = float(acell)**3
                    Density = MN/(Na*Volume*10**(-24)) #calculation of densities in g/cm3
                    distance, gofr = np.loadtxt(file, usecols=(0, couples.index(couple) + col + 1),
                                                skiprows=1, unpack=True) 
                    distance = distance[:-1]
                    gofr = gofr[:-1]
                    #ax.plot(distance,data_acell,intgofr, color = colors_acell[acell])
                    lines[str(round(Density,1))], = ax1.plot(distance,gofr, 
                          color = colors_densities[round(Density,1)], linewidth=2)
                    ax2.plot(distance,gofr, color = colors_densities[round(Density,1)],linewidth=2)

                    
                #
                ax1.set_title('Pair distribution function of '+couple ,fontsize=25,fontweight='bold' )
                #ax1.set_title('At '+temperature.strip('T')+'000 K',fontsize=25,fontweight='bold' )
                # Fine-tune figure
                plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.77, left = 0.1 , hspace = 0, wspace = 0)
                ax0 = fig.add_subplot(111, frameon=False)
                s = [(k, lines[k]) for k in sorted(lines.keys())] #sort the dictionnary
                legend = ax0.legend([v for k,v in s],[k for k,v in s], title = ' $\\bf{Density}$ \n (g.cm$^{-3}$)',
                                    bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize = 20)
                plt.setp(legend.get_title(),fontsize = 25)
                
                plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False,
                                labelleft=False, left=False, labelright = False, right=False)
                filename = 'gofrs_'+files[0].split('/')[-1].split('_')[0]+'_'+temperature+'_'+couple+'.pdf'
                print(filename,' is created')
                plt.savefig(filename, bbox_inches="tight", dpi=150)
        col +=1
    #plt.show()

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])
















