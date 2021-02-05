#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the coordination number       ****
                    One subplot per pairs of atoms
            
                    
"""



#     ********* Importation of the packages and modules used here *********
import sys, glob
import getopt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend import Legend
from matplotlib.ticker import AutoLocator,AutoMinorLocator
from matplotlib.gridspec import GridSpec
import numpy as np
import re
import natsort
import crystallography as cr
from matplotlib.colors import LinearSegmentedColormap




def creation_plot_P(plot_parameters):
    nlines = 3
    size_figure = (6,2.3*nlines) #height of the figure depends on the number of pair we display    
    plt.close(1)

    fig = plt.figure(figsize = (size_figure[0],size_figure[1]))
    
    gs = GridSpec(nlines,2,width_ratios=[1,2],top = 0.88, bottom = 0.12, 
                  right = 0.95, left = 0.15, wspace = 0, hspace = 0.03)
    ax1 = fig.add_subplot(gs[0,1])
    ax2 = fig.add_subplot(gs[1,1])    
    ax3 = fig.add_subplot(gs[2,1])
    
    ax11 = fig.add_subplot(gs[0,0],sharey=ax1)
    ax22 = fig.add_subplot(gs[1,0],sharey=ax2)
    ax33 = fig.add_subplot(gs[2,0],sharey=ax3)
     
    major_xticks = np.arange(-5,1.5,1)
    minor_xticks = np.arange(-5,1.5,0.2)
    #apply setup    
    for ax in [ax11,ax22,ax33]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.set_xlim(-1,1)
        ax.set_facecolor((1,1,1,0))

    for ax in [ax1,ax2,ax3]:
        ax.set_xlim(1,275)
        ax.set_xscale('log')
        plt.setp(ax.get_xticklabels()[1], visible=False) 
        ax.tick_params(labelleft=False)
        
    for ax in [ax11,ax1,ax22,ax2]:
        ax.tick_params(labelbottom=False)

    for ax in [ax1,ax2,ax3,ax11,ax22,ax33]:
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = plot_parameters["size_font_ticks"],
                       width = plot_parameters["size_lines"]/2)
        majorLocator = AutoLocator()
        minorLocator = AutoMinorLocator()                                                        
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_minor_locator(minorLocator)
        plt.autoscale(axis='y')
        ax.tick_params(which = 'both', labelsize = plot_parameters["size_font_ticks"], 
                       width = plot_parameters["size_lines"]/2)   
                        
    ax0 = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax0.set_xlabel(r'Pressure (GPa)', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                   labelpad = plot_parameters["shift_labelpad"])
    ax0.set_ylabel(r'Coordination number', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                   labelpad = plot_parameters["shift_labelpad"])
    return fig, ax1, ax2, ax3, ax11, ax22, ax33, ax0


def creation_plot(plot_parameters):
    nlines = 3
    size_figure = (6,2.3*nlines) #height of the figure depends on the number of pair we display    
    plt.close(1)
    fig, (ax1,ax2,ax3) = plt.subplots(nrows = nlines, ncols = 1, sharex = True, 
         sharey = False, figsize=size_figure)
    major_xticks = np.arange(0, 8, 0.5) 
    minor_xticks = np.arange(0, 8, 0.1)
    for ax in [ax1,ax2,ax3]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.set_xlim(0.9,5.58)
        ax.xaxis.set_ticks_position('both')
        majorLocator = AutoLocator()
        minorLocator = AutoMinorLocator()                                                        
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_minor_locator(minorLocator)
        
        plt.autoscale(axis='y')
        ax.tick_params(which = 'both', labelsize = plot_parameters["size_font_ticks"], 
                       width = plot_parameters["size_lines"]/2)   
    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    fig.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax0 = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax0.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', fontsize = plot_parameters["size_fonts"],
                   labelpad = plot_parameters["shift_labelpad"])
    ax0.set_ylabel(r'Coordination number', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                   labelpad = plot_parameters["shift_labelpad"])
    return fig, ax1, ax2, ax3, ax0


def create_colors():
    """ function to create the colors_T dictionnary from color map """
    colors_T = {}
    cm_data = np.loadtxt("/Users/akobsch/Dropbox/Recherche/PhD-MD-silicates/simulations/scripts/Analysis/ScientificColourMaps6/romaO/romaO.txt")
    cm_data = cm_data[::-1] #for reverse colors
    new_map = LinearSegmentedColormap.from_list('new', cm_data)
    temperatures = ['T2','T2.5','T3','T3.5','T4','T4.5','T5','T5.5','T6','T6.5','T7','T7.5','T7.7']
    color = iter(new_map(np.linspace(0,1,len(temperatures)))) #Creation of the color list    
    for T in temperatures:
        c = next(color)
        colors_T[T] = c    
    return colors_T


def format1label(label):
    """formatage of label with removing _1 """
    i=0
    while True:
        if i <= len(label)-1:
            if re.match('_',label[i]):
                num=0
                try:
                    while re.match('[0-9]',label[i+1+num]):
                        num +=1
                except IndexError: #index error when we arrive at the end of the cluster name
                    pass
                    #print('end of the cluster')
                if label[i+1:i+1+num] == '1':
                    label = label[:i] + label[i+1+num:]
                else:
                    label = label[:i]+'$_{'+label[i+1:i+1+num]+'}$' + label[i+1+num:]
                    i = i+5
            i = i+1
        else:break
    return label


def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./XX/CaAl2Si2O8_T3_nvt_a12.5.outcar.gofr.dat
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.split('/')[-1]
    temperature = filename.split('_')[1]
    acell = filename.split('.outcar.gofr.dat')[0].split('_')[3].strip('a')
    return temperature, acell

def extract_TP(thermofile, column_number, TP, addtxt):
    """extract temperature and pressure from thermofile"""
    with open(thermofile, 'r') as f:
        [f.readline() for i in range(3)]
        #extract data
        while True:
            line = f.readline()
            if not line: break
            else:    
                entry=line.split('\n')[0].split('\t')
                TP[entry[0].split('outcar.umd.dat')[0].split('/')[-1]+addtxt] = (int(entry[column_number['T']]),float(entry[column_number['P']]))
    return TP

def main(argv):
    """     ********* Main program *********     """
    #parameters for the figures depending on the output format (presentation or article)
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),
                       "size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 20}
    #other dictionnaries and parameters for the figure
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200',
                'T5':'#ffcd01','T5.5':'#ff6e00','T6':'#ff0101','T6.5':'#ff00a2',
                'T7':'#ff01de','T10':'#ffa6f4','T15':'#ffe86e','T20':'#ffbf90'}
    #colors_T = create_colors()
    allT = {}
    #initialization of parameters
    filename2 = ''
    xvariable = 'rho'
    #dictionnary for fullaverages file in full version
    column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,
                     'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,
                     'Cvm':25,'stdev_Cvm':26,'testCv':31,'stdev_testCv':32}
    #other parameters
    Na=6.022*10**23
    try:
        options,arg = getopt.getopt(argv,"hf:g:a:x:",["filename1","gfilename2","atoms",'xvariable'])
    except getopt.GetoptError:
        print("plot_coordinence_separated.py -f <gofr_compound.txt> -g <gofr_compound.txt 2> -a <atomic pair list:'Ca-O,K-O,Al-O,Si-O'> -x <xvariable (P or rho)>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_coordinence_separated.py program to plot coordinence as a function of density for each T and all cations with O')
            print("plot_coordinence_separated.py -f <gofr_compound.txt> -g <gofr_compound.txt 2> -a <atomic pair list:'Ca-O,K-O,Al-O,Si-O'> -x <xvariable (P or rho)>")
            print("plot_coordinence_separated.py requires gofr txt file containing every bond distance and coordination of one compound for every acell and T with the number of elements in header (use analyze_gofrs script etc.)")
            print('')
            print('For plots as function of P, make sure the files from fullaverages.py are in the current folder')
            sys.exit()
        if opt in ('-f','--filename'):
            filename = str(arg)
        elif opt in ('-g','--gfilename'):
            filename2 = str(arg) 
        elif opt in ('-a','--atoms'):
            atoms = arg.split(',')
        elif opt in ('-x','--xvariable'):
            xvariable = str(arg)
    #************ creation of the figure to plot (with correct number of lines)
    if filename2 != '':
        files = [filename,filename2]
        figurename = 'coordinence_'+filename.split('.txt')[0].split('_')[1]+'-'+filename2.split('.txt')[0].split('_')[1]+'_'+xvariable
        typeline = {filename : '--', filename2 : '-'}
        typemarker = {filename : 'x', filename2 : ''}
        compo = {filename2:format1label(filename2.split('_')[1]),filename:format1label(filename.split('_')[1])}
    else:
        files = [filename]
        figurename = 'coordinence_'+filename.split('.txt')[0]+'_'+xvariable
        typeline = {filename : '-'}
        typemarker = {filename : '*'}
        compo = {filename:format1label(filename.split('_')[1])}
    
    if xvariable == 'P':
        fig, ax1,ax2,ax3,ax11,ax22,ax33, ax0 = creation_plot_P(plot_parameters)
        axis1 = [ax1,ax11]
        axis2 = [ax2,ax22]
        axis3 = [ax3,ax33]
    else:
        fig, ax1,ax2,ax3,ax0 = creation_plot(plot_parameters)
        axis1 = [ax1]
        axis2 = [ax2]        
        axis3 = [ax3]
    offset = 0 #offset for printing text on Na-O K-O subplot
    #************** for each file we extract and plot data for each atom on the correct subplot
    for file in files:
        #***** Initialization of data dictionnaries
        pair_columns = {}
        data = {} #big dictionnary with inside coord  all coresponding to different T and atom pairs
        X = {}  #idem for rho or P
        #******* Extract P and T for each thermo file
        if xvariable == 'P':
            TP = {}
            mineralname = file.split('_')[1]
            thermofiles = sorted(glob.glob('thermo_'+mineralname+'*.txt'))
            for thermofile in thermofiles:
                TP = extract_TP(thermofile, column_number, TP,'')
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
                    X[atom] = {}
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
                                    if xvariable == 'P':
                                        X[atom][temperature].append(TP[entry[0].split('outcar.gofr.dat')[0].split('/')[-1]][1])
                                    else:
                                        X[atom][temperature].append( MN/(Na*float(acell)**3*10**(-24)) )    #calculation density
                                except KeyError:
                                    data[atom][temperature] = [ float(entry[pair_columns[atom]]) ]
                                    if xvariable == 'P':
                                        X[atom][temperature] = [TP[entry[0].split('outcar.gofr.dat')[0].split('/')[-1]][1]]
                                    else:
                                        X[atom][temperature] = [ MN/(Na*float(acell)**3*10**(-24)) ]    #calculation density
        #***** Plot each atom pair on the correct subplot --> ax1 for Na/K-O, ax2 for Al-O, ax3 for Si-O
        for temperature in data['Si-O']:
            allT[temperature] = ''
            for atom in data:
                if atom == 'Na-O' or atom == 'K-O' or atom == 'Ca-O' :
                    for ax in axis1:
                        ax.plot(X[atom][temperature],data[atom][temperature], 
                                marker = typemarker[file], linestyle = typeline[file], 
                                color=colors_T[temperature], markersize = plot_parameters["size_markers"],
                                linewidth=plot_parameters["size_lines"] )
                    ax1.text( 0.02,0.91-offset, atom ,transform=ax1.transAxes, 
                             fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
#                    plt.setp(ax.get_xticklabels(), visible=False)
                elif atom == 'Si-O':
                    for ax in axis3:
                        ax.plot(X[atom][temperature],data[atom][temperature], 
                                marker = typemarker[file], linestyle = typeline[file], 
                                color=colors_T[temperature], markersize = plot_parameters["size_markers"], 
                                linewidth=plot_parameters["size_lines"] )
                    ax3.text( 0.02,0.91, atom ,transform=ax3.transAxes, 
                             fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
                elif atom == 'Al-O':
                    for ax in axis2:
                        ax.plot(X[atom][temperature],data[atom][temperature], 
                                marker = typemarker[file], linestyle = typeline[file],
                                color=colors_T[temperature], markersize = plot_parameters["size_markers"], 
                                linewidth=plot_parameters["size_lines"] )
                    ax2.text(0.02,0.91, atom ,transform=ax2.transAxes, 
                             fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
        offset = offset + 0.10        #we update the offset for printing text on Na-O K-O subplot after finishing the first file              
    
    #************ legend
    Tlist = [label for label in natsort.natsorted(allT)]
    custom_patch = [mpatches.Patch(color=colors_T[key]) for key in Tlist]
    legend = ax0.legend([col for col in custom_patch],[str(int(float(label.strip('T'))*1000)) for label in Tlist],
                        title = '$\\bf{Temperature~(K)}$', bbox_to_anchor=(0.5, 1.07), 
                        loc="lower center", fontsize = plot_parameters["size_font_ticks"], 
                        borderaxespad=0., ncol=len(Tlist))
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])
    
    legend_labels = {}
    for file in files: 
        legend_labels[compo[file]] =  plt.Line2D((0,1),(0,0), color='k', 
                     linestyle = typeline[file], markersize = plot_parameters["size_markers"]+2, 
                     marker=typemarker[file])
    s = [(k, legend_labels[k]) for k in sorted(legend_labels.keys(),reverse = False)]      
    legend1 = ax0.legend([v for k,v in s],[k for k,v in s], loc='lower center', 
                         bbox_to_anchor=(0.5, 1.02), ncol=len(files), 
                         fontsize = plot_parameters["size_font_ticks"], borderaxespad=0.)
    
    ax0.add_artist(legend)
    
    #save the figure
    figurename = figurename + '.pdf'
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    
    print(figurename, '   created')
    #plt.show()
   
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])
