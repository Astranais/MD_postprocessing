#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  12 2017

@author: anais
Langage : Python3

                ****     Plot the mean square displacement of each atoms at all T and rho       ****
                            one T per row / one atom per column
                            every rho on each graph
                             *************  ARTICLE VERSION  *************

"""



#     ********* Importation of the packages and modules used here *********     """
import sys
import getopt
import glob
import numpy as np
import matplotlib.pyplot as plt
import crystallography as cr
import natsort
import re
from matplotlib.lines import Line2D #useful to create a custom legend



def split_name(filename):
    """ Function to extract T and acell from filename """
    # *********** My filename are in the format CaAl2Si2O8_T3_nvt_a12.0.outcar.gofr.dat
    # ******* so I can split their name with _ and tale the acell and T from their name
    # **** Then change these two lines in harmony with your filename
    temperature = filename.split('_')[1]                                     #I extract the corresponding temperature
    acell = filename.split('.outcar.msd.dat')[0].split('_')[3].strip('a')   #I extract the corresponding acell
    return temperature, acell


def main(argv):
    """     ********* Main program *********     """
    #parameters for the figures depending on the output format (presentation or article)
    size_fonts = {'oral':30,'article':12} #size in px
    size_fonts_ticks = {'oral':19,'article':10}
    size_figure = {'oral':(18.5,16),'article':(14,14)}
    size_markers = {'oral':11,'article':4}
    size_lines = {'oral':4,'article':2} #linewidth in pt
    shift_label = {'oral':45,'article':20}
    #other dictionnaries and parameters for the figure
    colors_densities = {}
    lines = {}   #dictionnary for the lines for legend (the keys are acell or densities, and values are the color of lines)
    limits = {} #dictionnary for y limits 
    AllTrho = [] #list for storing the temperature in order to sort them
    files = []
    #other parameters
    Na=6.022*10**23  #avogadro constant
    try:
        options,arg = getopt.getopt(argv,"hm:t:",["mineralfile","typefigure"])
    except getopt.GetoptError:
        print("plot_msd_allrhoT.py  -m <mineralfile with elements> -t <type figure>(oral or article)")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print("plot_msd_allrhoT.py script to plot all the msd at every acell and T (matrix with atom in columns and T in lines)")
            print("plot_msd_allrhoT.py  -m <mineralfile with elements>-t <type figure>(oral or article)")
            print('requires the file containing elements and number (in order to compute the densities)')
            sys.exit()
        elif opt in ("-m","--mineralfile"):
            mineralfile = str(arg)
        elif opt in ("-t","--typefigure"):
            size_fonts = size_fonts[str(arg)]  #in this step, the dictionnary became the float used for font size
            size_fonts_ticks = size_fonts_ticks[str(arg)]
            size_figure = size_figure[str(arg)]
            size_markers = size_markers[str(arg)]
            size_lines = size_lines[str(arg)]
            shift_label = shift_label[str(arg)]
    #***** Calculation of the molecular mass
    with open(mineralfile,'r') as mf:
        entry = mf.readline()
        elements = entry.split()[1:]
        entry = mf.readline()
        number = entry.split()[1:]
    MN = 0
    for i in range(len(elements)):
        MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]  
    #***** List the files in alphabetical order regarding to temperature (need to extract T and compute it because of the T4.5 listed before T4)
    initfiles = sorted(glob.glob('*outcar.msd.dat'))                             #I list every msd files in alphabetic order
    for file in initfiles:
        temperature, acell = split_name(file)
        temperature = int(float(temperature.strip('T'))*1000)
        AllTrho.append((temperature,acell))
    AllTrho = sorted(AllTrho)     #first I sort everything in alphabetical order           
    AllTrho = sorted(AllTrho,key=lambda tup: tup[0],reverse=True)        #then I reverse the order of the T only
    for Trho in AllTrho:
        T = Trho[0]/1000
        rem = Trho[0]%1000
        if rem == 0:
            T = int(T)
        pattern = '_T'+str(T)+'_nvt_a'+Trho[1]
        for file in initfiles:
            if re.search(pattern,file):
                files.append(file)
    #print(files)
    #***** Calculation of the number of figures (T * ntypat)
    firstfile = files[0]                                                      #I take the first file  
    temperature0, acell0 = split_name(firstfile)
    with open(firstfile, 'r') as f:
        atoms = f.readline()
        atoms = atoms.strip('time_(fs)').split()                            #I extract the atoms
    T_numsubplot = [temperature0] #list for the correspondance T / subplot number
    for file in files:
        temperature, acell = split_name(file)
        if temperature != temperature0:
            T_numsubplot.append(temperature)
            temperature0 = temperature
    temperature0, acell0 = split_name(firstfile)
    #***** Count of the number of densities (needed for the automatic color change)
    #first store every density into the dictionary
    for file in files:
        #calculation of densities in g/cm3
        temperature, acell = split_name(file)  
        Volume = float(acell)**3
        Density = MN/(Na*Volume*10**(-24)) 
        #we store the densities into the dictionnary
        colors_densities[round(Density,1)] = []
    #now the dictionary has every possible densities of our mineral, we attribute the colors to each density
    list_densities = []
    color = iter(plt.cm.rainbow(np.linspace(0,1,len(colors_densities)))) #Creation of the color list
    for key in natsort.natsorted(colors_densities):
        c = next(color)
        colors_densities[key] = c
        list_densities.append(key) #list usefull if we want to plot the value of densities near to the lines, see plot_speciation_r0
    print(list_densities, len(list_densities))
    #*****************************
    #*******************
    #*******
    #Read & plot 
    plt.close(1)
    fig = plt.figure(1,figsize=size_figure)
    plt.subplots_adjust(top = 0.97, bottom = 0.07, right = 0.89, left = 0.07, hspace = 0, wspace = 0)
    plt.subplot(len(T_numsubplot),len(atoms),1)
    ax = plt.gca()
    for atom in atoms:
        for file in files:
            #calculation of densities in g/cm3
            temperature, acell = split_name(file)  
            Volume = float(acell)**3
            Density = MN/(Na*Volume*10**(-24)) 
            #change of subplot
            if temperature != temperature0:
                plt.subplot(len(T_numsubplot),len(atoms),T_numsubplot.index(temperature)*4+1+atoms.index(atom))
                ax = plt.gca()
                temperature0 = temperature
            #importation of data
            Temps,DataAtom = np.loadtxt(file,skiprows=1,usecols=(0,atoms.index(atom)+1),unpack=True)
            #limitation of data along x 
            Xlimit = 10000
            ax.set_xlim(0,Xlimit)
            for i in range(len(Temps)): 
                if Temps[i] >= Xlimit:
                    Temps[i] = np.NaN
                    DataAtom[i] = np.NaN
            #plot
            lines[str(round(Density,1))], = ax.plot(Temps,DataAtom,'-', color = colors_densities[round(Density,1)], linewidth = size_lines)
            #limitation of data along y
            if atoms.index(atom) == 0:
                ax.autoscale(axis='y', tight=True)
                limits[temperature] = ax.get_ylim()
            ax.set_ylim(limits[temperature])
            #we add x and y labels outside the plot on the right columns and lines
            if T_numsubplot.index(temperature) == 0:
                ax.set_xlabel(atoms[atoms.index(atom)], fontweight = 'bold', fontsize = size_fonts)
                ax.xaxis.set_label_position('top')
            if atoms.index(atom) == len(atoms)-1:
                ax.set_ylabel(str(int(float(temperature.strip('T'))*1000))+' K', fontweight = 'bold', fontsize = size_fonts)
                ax.yaxis.set_label_position('right')
            #we make the graph prettier
            if  atoms.index(atom) != 0:
                plt.setp(ax.get_yticklabels(), visible=False)
            if T_numsubplot.index(temperature) != len(T_numsubplot)-1:
                plt.setp(ax.get_xticklabels(), visible=False)
            major_ticks = np.arange(0, Xlimit, 2000) 
            minor_ticks = np.arange(0, Xlimit, 500)
            ax.set_xticks(major_ticks)
            ax.set_xticks(minor_ticks, minor=True)
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.yaxis.set_tick_params(which = 'both', direction='inout')
            ax.xaxis.set_tick_params(which = 'both', direction='inout')
            ax.grid(axis = 'y', which = 'major', linestyle = '--', linewidth = size_lines/2, alpha = 0.5)
            #ax.get_xaxis().set_tick_params(direction='inout')
            #ax.get_yaxis().set_tick_params(direction='inout')
            #plt.setp(ax.get_xticklabels()[-2], visible=False)
            #plt.setp(ax.get_yticklabels()[-1], visible=False)
            ax.tick_params(which = 'both', labelsize = size_fonts_ticks, width = size_lines/2)
    # Fine-tune figure
    ax0 = fig.add_subplot(111, frameon=False)
    custom_lines = [Line2D([0],[0],color = colors_densities[key], ls = '-', marker = '', linewidth = size_lines/2) for key in natsort.natsorted(colors_densities)]
    #legend = plt.legend([line for line in custom_lines],[label for label in natsort.natsorted(colors_densities)],title = '$\\bf{Density}$ \n (g.cm$^{-3}$)', bbox_to_anchor=(1.05, 1), loc=2, fontsize = size_fonts, borderaxespad=0.)
    legend = plt.legend([line for line in custom_lines],[label for label in natsort.natsorted(colors_densities)],title = '$\\bf{Density}$ (g.cm$^{-3}$)', bbox_to_anchor=(0.5, 1.04), loc='lower center', fontsize = size_fonts, borderaxespad=0.,ncol=10)
    plt.setp(legend.get_title(),fontsize = size_fonts)

#    s = [(k, lines[k]) for k in sorted(lines.keys())]
#    ax0.legend([v for k,v in s],[k for k,v in s], title = ' $\\bf{density}$ \n (g.cm$^{-3}$)',bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    ax0.set_xlabel(r'Time (fs)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_label)
    ax0.set_ylabel(r"Mean square displacement ($\mathregular{\AA^2}$)",fontsize = size_fonts,fontweight='bold', labelpad = shift_label + shift_label/2)
    
    figurename = 'msd_'+file.split('_')[0]+'.pdf'
    plt.savefig(figurename, bbox_inches='tight', dpi = 300)
    print(figurename, 'is created')
    #plt.show()




#ax.zaxis.set_major_formatter(FormatStrFormatter('%.1f'))  #set one decimal to ticks label``
# sharex=True, sharey='row'

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])























