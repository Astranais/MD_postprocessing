#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the coordination number       ****
                    One subplot per pairs of atoms
                    3 lines and 2 columns
                    Na feldspar as a common reference
            
                    
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
        options,arg = getopt.getopt(argv,"hf:g:j:",["filename1","gfilename2","jfilename3"])
    except getopt.GetoptError:
        print("plot_coordinence_3x2.py -f <gofr_compound.txt> -g <gofr_compound.txt 2>  -j <gofr_compound.txt 3> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_coordinence_3x2.py program to plot coordinence as a function of density for each T and all cations with O')
            print("plot_coordinence_3x2.py -f <gofr_compound.txt> -g <gofr_compound.txt 2>  -j <gofr_compound.txt 3>")
            print("plot_coordinence_3x2.py requires gofr txt file containing every bond distance and coordination of one compound for every acell and T with the number of elements in header (use analyze_gofrs script etc.)")
            print('')
            sys.exit()
        if opt in ('-f','--filename1'):
            filename = str(arg)
        elif opt in ('-g','--gfilename2'):
            filename2 = str(arg) 
        elif opt in ('-j','--jfilename3'):
            filename3 = str(arg)
    #************ creation of the figure to plot (with correct number of lines)
    files = [filename,filename2, filename3]
    filename_index = {filename:0,filename2:0,filename3:1}
    figurename = 'coordinence_all'
    typeline = {filename : '--', filename2 : '-', filename3: ':'}
    typemarker = {filename : 'x', filename2 : None, filename3: '+'}
    atoms = ['Ca-O','K-O','Na-O','Al-O','Si-O']
    atoms_index = {'Ca-O':0,'K-O':0,'Na-O':0,'Al-O':2,'Si-O':1}
    nlines = 3
    size_figure = (10,2.3*nlines) #height of the figure depends on the number of pair we display
    plt.close(1)
    fig, axes = plt.subplots(nrows = nlines, ncols = 2, sharex = True, sharey = 'row', figsize=size_figure)
    major_xticks = np.arange(0, 8, 1) 
    minor_xticks = np.arange(0, 8, 0.5)
    for i in range(nlines):
        for j in range(2):
            ax = axes[i,j]
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
                axes[atoms_index[atom],filename_index[file]].plot(rho[atom][temperature],data[atom][temperature], marker = typemarker[file], linestyle = typeline[file], color=colors_T[temperature], markersize = plot_parameters["size_markers"], linewidth=plot_parameters["size_lines"] )
                if atom == 'Na-O':
                    axes[atoms_index[atom],filename_index[file]].text( 0.02,0.81, atom ,transform=axes[atoms_index[atom],filename_index[file]].transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
                else:
                    axes[atoms_index[atom],filename_index[file]].text( 0.02,0.91, atom ,transform=axes[atoms_index[atom],filename_index[file]].transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
                if file == filename2:
                    axes[atoms_index[atom],1].plot(rho[atom][temperature],data[atom][temperature], marker = typemarker[file], linestyle = typeline[file], color=colors_T[temperature], markersize = plot_parameters["size_markers"], linewidth=plot_parameters["size_lines"] )
                    if atom == 'Na-O':
                        axes[atoms_index[atom],1].text( 0.02,0.81, atom ,transform=axes[atoms_index[atom],1].transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
    #legend
    legend_labels = {}
    compo = {filename:r'KAlSi$_{3}$O$_{8}$',filename2:r'NaAlSi$_{3}$O$_{8}$',filename3:r'CaAl$_{2}$Si$_{2}$O$_{8}$'}
    for file in files: 
        legend_labels[compo[file]] =  plt.Line2D((0,1),(0,0), color='k', linestyle = typeline[file], markersize = plot_parameters["size_markers"], marker=typemarker[file])
    s = [(k, legend_labels[k]) for k in sorted(legend_labels.keys(),reverse = False)]      
    legend1 = ax0.legend([v for k,v in s],[k for k,v in s], loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize = plot_parameters["size_font_ticks"], borderaxespad=0.)
    #save the figure
    figurename = figurename + '.pdf'
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    
    print(figurename, '   created')
    #plt.show()
   
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])
