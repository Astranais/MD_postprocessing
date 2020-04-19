#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the coordination number       ****
                    One subplot per pairs of atoms
            
                    
"""



#     ********* Importation of the packages and modules used here *********
import sys
import getopt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend import Legend
from matplotlib.ticker import AutoLocator,AutoMinorLocator
import numpy as np
import crystallography as cr








def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./XX/CaAl2Si2O8_T3_nvt_a12.5.outcar.gofr.dat
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.split('/')[-1]
    temperature = filename.split('_')[1]
    acell = filename.split('.outcar.gofr.dat')[0].split('_')[3].strip('a')
    return temperature, acell


def main(argv):
    """     ********* Main program *********     """
    #parameters for the figures depending on the output format (presentation or article)
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 20}
    #other dictionnaries and parameters for the figure
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200','T5':'#ffcd01','T5.5':'#ff6e00','T6':'#ff0101','T6.5':'#ff00a2','T7':'#ff01de','T10':'#ffa6f4','T15':'#ffe86e','T20':'#ffbf90'}
    
    #initialization of parameters
    filename2 = ''
    #other parameters
    Na=6.022*10**23
    try:
        options,arg = getopt.getopt(argv,"hf:g:",["filename1","gfilename2"])
    except getopt.GetoptError:
        print("plot_coordinence_separated.py -f <gofr_compound.txt> -g <gofr_compound.txt 2> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_coordinence_separated.py program to plot coordinence as a function of density for each T and all cations with O')
            print("plot_coordinence_separated.py -f <gofr_compound.txt> -g <gofr_compound.txt 2>)")
            print("plot_coordinence_separated.py requires gofr txt file containing every bond distance and coordination of one compound for every acell and T with the number of elements in header (use analyze_gofrs script etc.)")
            print('')
            sys.exit()
        if opt in ('-f','--filename'):
            filename = str(arg)
        elif opt in ('-g','--gfilename'):
            filename2 = str(arg) 
    #************ creation of the figure to plot (with correct number of lines)
    if filename2 != '':
        files = [filename,filename2]
        figurename = 'coordinence_'+filename.split('.txt')[0].split('_')[1]+'-'+filename2.split('.txt')[0].split('_')[1]
        typeline = {filename : '*-', filename2 : 'x--'}
    else:
        files = [filename]
        figurename = 'coordinence_'+filename.split('.txt')[0]
        typeline = {filename : '*-'}
    atoms = ['Ca-O','K-O','Al-O','Si-O']
    nlines = 3
    size_figure = (6,2.3*nlines) #height of the figure depends on the number of pair we display
    offset = 0 #offset for printing text on Na-O K-O subplot
    plt.close(1)
    fig, (ax1,ax2,ax3) = plt.subplots(nrows = nlines, ncols = 1, sharex = True, sharey = False, figsize=size_figure)
    major_xticks = np.arange(0, 8, 0.5) 
    minor_xticks = np.arange(0, 8, 0.1)
    for ax in [ax1,ax2,ax3]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)                                           
        ax.xaxis.set_ticks_position('both')
        majorLocator = AutoLocator()
        minorLocator = AutoMinorLocator()                                                        
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_minor_locator(minorLocator)
        ax.set_xlim(0.9,5.58)
        plt.autoscale(axis='y')
        ax.tick_params(which = 'both', labelsize = plot_parameters["size_font_ticks"], width = plot_parameters["size_lines"]/2)   
    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    fig.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax0 = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    ax0.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], labelpad = plot_parameters["shift_labelpad"])
    ax0.set_ylabel(r'Coordination number', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], labelpad = plot_parameters["shift_labelpad"])
    #************** for each file we extract and plot data for each atom on the correct subplot
    for file in files:
        #***** Initialization of data dictionnaries
        pair_columns = {}
        data = {} #big dictionnary with inside coord  all coresponding to different T and atom pairs
        rho = {}  #idem for rho
        #***** Calculation of the molecular mass
        with open(file, 'r') as f:
            line = f.readline()
            entry=line.split()
            elements = entry[1:]
            line = f.readline()
            entry=line.split()
            number = entry[1:]
            line = f.readline()
            entry=line.split()
            pair = entry[1:]
        MN = 0
        for i in range(len(elements)):
            MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
        #******** extraction of data        
        pair_columns[pair[0]] =4
        for i in range(1,len(pair)):
            if pair[i] != pair[i-1]:
                pair_columns[pair[i]] = i+4
        print('All the pair availables in the file with their index are:', pair_columns)
        for atom_pair in pair_columns:
            for atom in atoms:
                if atom_pair == atom:
                    data[atom] = {}
                    rho[atom] = {}
        with open(file,'r') as f:
            [f.readline() for i in range(4)]
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    temperature, acell = split_name(entry[0]) 
                    for atom_pair in pair_columns:
                        for atom in atoms:
                            if atom_pair == atom:
                                try:
                                    data[atom][temperature].append( float(entry[pair_columns[atom]]) )
                                    rho[atom][temperature].append( MN/(Na*float(acell)**3*10**(-24)) )    #calculation density
                                except KeyError:
                                    data[atom][temperature] = [ float(entry[pair_columns[atom]]) ]
                                    rho[atom][temperature] = [ MN/(Na*float(acell)**3*10**(-24)) ]    #calculation density
        #***** Plot each atom pair on the correct subplot --> ax1 for Na/K-O, ax2 for Al-O, ax3 for Si-O
        for temperature in data['Si-O']:
            for atom in data:
                if atom == 'Ca-O' or atom == 'K-O':
                    ax1.plot(rho[atom][temperature],data[atom][temperature], marker = typeline[file][0], linestyle = typeline[file][1:], color=colors_T[temperature], markersize = plot_parameters["size_markers"], linewidth=plot_parameters["size_lines"] )
                    ax1.text( 0.02,0.91-offset, atom ,transform=ax1.transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
#                    plt.setp(ax.get_xticklabels(), visible=False)
                elif atom == 'Si-O':
                    ax2.plot(rho[atom][temperature],data[atom][temperature], marker = typeline[file][0], linestyle = typeline[file][1:], color=colors_T[temperature], markersize = plot_parameters["size_markers"], linewidth=plot_parameters["size_lines"] )
                    ax2.text( 0.02,0.91, atom ,transform=ax2.transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
                elif atom == 'Al-O':
                    ax3.plot(rho[atom][temperature],data[atom][temperature], marker = typeline[file][0], linestyle = typeline[file][1:], color=colors_T[temperature], markersize = plot_parameters["size_markers"], linewidth=plot_parameters["size_lines"] )
                    ax3.text(0.02,0.91, atom ,transform=ax3.transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
        offset = offset + 0.10        #we update the offset for printing text on Na-O K-O subplot after finishing the first file              
    #save the figure
    figurename = figurename + '.pdf'
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    
    print(figurename, '   created')
    #plt.show()
   
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])
