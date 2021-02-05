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
import sys, glob
import getopt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend import Legend
from matplotlib.ticker import AutoLocator,AutoMinorLocator
import matplotlib.gridspec as gridspec
import numpy as np
import natsort
import crystallography as cr
from matplotlib.colors import LinearSegmentedColormap


def creation_plot_P(plot_parameters):
    nlines = 3
    size_figure = (10,2.3*nlines) #height of the figure depends on the number of pair we display    
    plt.close(1)

    fig = plt.figure(figsize = (size_figure[0],size_figure[1]))
    
    #first gridspec to make the two columns
    gs0 = gridspec.GridSpec(1,2,top = 0.88, bottom = 0.12, 
                  right = 0.95, left = 0.12, wspace = 0.02, hspace = 0.03)
    
    #same gridspec as in plot_coordinence_separated (one column) for the 1st column
    gs00 = gridspec.GridSpecFromSubplotSpec(nlines,2,width_ratios=[1,2],
                  wspace = 0, hspace = 0.03,
                  subplot_spec=gs0[0])
    
    ax1 = fig.add_subplot(gs00[0,1])
    ax2 = fig.add_subplot(gs00[1,1])    
    ax3 = fig.add_subplot(gs00[2,1])
    
    ax11 = fig.add_subplot(gs00[0,0],sharey=ax1)
    ax22 = fig.add_subplot(gs00[1,0],sharey=ax2)
    ax33 = fig.add_subplot(gs00[2,0],sharey=ax3)

    #same gridspec as in plot_coordinence_separated (one column) for the 2nd column
    gs01 = gridspec.GridSpecFromSubplotSpec(nlines,2,width_ratios=[1,2],
                  wspace = 0, hspace = 0.03,
                  subplot_spec=gs0[1])
    
    ax4 = fig.add_subplot(gs01[0,1],sharey=ax1)
    ax5 = fig.add_subplot(gs01[1,1],sharey=ax2)    
    ax6 = fig.add_subplot(gs01[2,1],sharey=ax3)
    
    ax44 = fig.add_subplot(gs01[0,0],sharey=ax4)
    ax55 = fig.add_subplot(gs01[1,0],sharey=ax5)
    ax66 = fig.add_subplot(gs01[2,0],sharey=ax6)


     
    major_xticks = np.arange(-5,1.5,1)
    minor_xticks = np.arange(-5,1.5,0.2)
    #apply setup    
    for ax in [ax11,ax22,ax33,ax44,ax55,ax66]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.set_xlim(-1,1)
        ax.set_facecolor((1,1,1,0))

    for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
        ax.set_xlim(1,275)
        ax.set_xscale('log')
        plt.setp(ax.get_xticklabels()[1], visible=False) 
        ax.tick_params(labelleft=False)

    for ax in [ax44,ax55,ax66]:
        ax.tick_params(labelleft=False)
    for ax in [ax1,ax11,ax2,ax22,ax4,ax44,ax5,ax55]:
        ax.tick_params(labelbottom=False)
    
    for ax in [ax1,ax2,ax3,ax4,ax5,ax6,ax11,ax22,ax33,ax44,ax55,ax66]:
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
    return fig, ax1, ax2, ax3, ax4,ax5,ax6, ax11, ax22, ax33, ax44,ax55,ax66, ax0


def creation_plot(plot_parameters):
    nlines = 3
    size_figure = (10,2.3*nlines) #height of the figure depends on the number of pair we display    
    plt.close(1)
    fig, axes = plt.subplots(nrows = nlines, ncols = 2, sharex = True, sharey = 'row', figsize=size_figure)
    major_xticks = np.arange(0, 8, 1) 
    minor_xticks = np.arange(0, 8, 0.2)
    for i in range(nlines):
        for j in range(2):
            ax=axes[i,j]
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
    ax1,ax2,ax3,ax4,ax5,ax6 = axes[0,0],axes[1,0],axes[2,0],axes[0,1],axes[1,1],axes[2,1]
    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    fig.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax0 = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax0.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', fontsize = plot_parameters["size_fonts"],
                   labelpad = plot_parameters["shift_labelpad"])
    ax0.set_ylabel(r'Coordination number', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                   labelpad = plot_parameters["shift_labelpad"])
    return fig, ax1,ax2,ax3,ax4,ax5,ax6, ax0



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
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 20}
    #other dictionnaries and parameters for the figure
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200','T5':'#ffcd01','T5.5':'#ff6e00','T6':'#ff0101','T6.5':'#ff00a2','T7':'#ff01de','T10':'#ffa6f4','T15':'#ffe86e','T20':'#ffbf90'}
    #colors_T = create_colors()
    allT = {}
    #initialization of parameters
    filename2 = ''
    xvariable = 'rho'
    #dictionnary for fullaverages file in full version
    column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,'Cvm':25,'stdev_Cvm':26,'testCv':31,'stdev_testCv':32}
    #other parameters
    Na=6.022*10**23
    try:
        options,arg = getopt.getopt(argv,"hf:g:j:x:",["filename1","gfilename2","jfilename3",'xvariable'])
    except getopt.GetoptError:
        print("plot_coordinence_3x2.py -f <gofr_compound.txt> -g <gofr_compound.txt 2>  -j <gofr_compound.txt 3>  -x <xvariable (P or rho)>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_coordinence_3x2.py program to plot coordinence as a function of density for each T and all cations with O')
            print("plot_coordinence_3x2.py -f <gofr_compound.txt> -g <gofr_compound.txt 2>  -j <gofr_compound.txt 3> -x <xvariable (P or rho)>")
            print("plot_coordinence_3x2.py requires gofr txt file containing every bond distance and coordination of one compound for every acell and T with the number of elements in header (use analyze_gofrs script etc.)")
            print('')
            print('For plots as function of P, make sure the files from fullaverages.py are in the current folder')
            sys.exit()
        if opt in ('-f','--filename1'):
            filename = str(arg)
        elif opt in ('-g','--gfilename2'):
            filename2 = str(arg) 
        elif opt in ('-j','--jfilename3'):
            filename3 = str(arg)
        elif opt in ('-x','--xvariable'):
            xvariable = str(arg)
    #************ creation of the figure to plot (with correct number of lines)
    files = [filename,filename2, filename3]
    filename_index = {filename:0,filename2:0,filename3:1}
    figurename = 'coordinence_all_'+xvariable
    typeline = {filename : '--', filename2 : '-', filename3: ':'}
    typemarker = {filename : 'x', filename2 : None, filename3: '+'}
    atoms = ['Ca-O','K-O','Na-O','Al-O','Si-O']
    atoms_index = {'Ca-O':0,'K-O':0,'Na-O':0,'Al-O':2,'Si-O':1}
    if xvariable == 'P':
        fig, ax1,ax2,ax3,ax4,ax5,ax6, ax11,ax22,ax33,ax44,ax55,ax66, ax0 = creation_plot_P(plot_parameters)
        axes = [ [[ax1,ax11],[ax4,ax44]],
                 [[ax2,ax22],[ax5,ax55]],
                 [[ax3,ax33],[ax6,ax66]] ]
    else:
        fig, ax1,ax2,ax3, ax4,ax5,ax6, ax0 = creation_plot(plot_parameters)
        axes = [ [[ax1],[ax4]],
                 [[ax2],[ax5]],
                 [[ax3],[ax6]] ]
    
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
                for ax in axes[atoms_index[atom]][filename_index[file]]:
                    ax.plot(X[atom][temperature],data[atom][temperature], marker = typemarker[file], linestyle = typeline[file], color=colors_T[temperature], markersize = plot_parameters["size_markers"], linewidth=plot_parameters["size_lines"] )
                if atom == 'Na-O':
                    axes[atoms_index[atom]][filename_index[file]][0].text( 0.02,0.81, atom ,transform=axes[atoms_index[atom]][filename_index[file]][0].transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
                else:
                    axes[atoms_index[atom]][filename_index[file]][0].text( 0.02,0.91, atom ,transform=axes[atoms_index[atom]][filename_index[file]][0].transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
                if file == filename2:
                    for ax in axes[atoms_index[atom]][1]:
                        ax.plot(X[atom][temperature],data[atom][temperature], marker = typemarker[file], linestyle = typeline[file], color=colors_T[temperature], markersize = plot_parameters["size_markers"], linewidth=plot_parameters["size_lines"] )
                    if atom == 'Na-O':
                        axes[atoms_index[atom]][1][0].text( 0.02,0.81, atom ,transform=axes[atoms_index[atom]][1][0].transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  )
    #************ legend
    Tlist = [label for label in natsort.natsorted(allT)]
    custom_patch = [mpatches.Patch(color=colors_T[key]) for key in Tlist]
    legend = ax0.legend([col for col in custom_patch],[str(int(float(label.strip('T'))*1000)) for label in Tlist],title = '$\\bf{Temperature~(K)}$', bbox_to_anchor=(0.5, 1.07), loc="lower center", fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol=len(Tlist))
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])
    
    legend_labels = {}
    compo = {filename:r'KAlSi$_{3}$O$_{8}$',filename2:r'NaAlSi$_{3}$O$_{8}$',filename3:r'CaAl$_{2}$Si$_{2}$O$_{8}$'}
    for file in files: 
        legend_labels[compo[file]] =  plt.Line2D((0,1),(0,0), color='k', linestyle = typeline[file], 
                     markersize = plot_parameters["size_markers"], marker=typemarker[file])
    s = [(k, legend_labels[k]) for k in sorted(legend_labels.keys(),reverse = False)]      
    legend1 = ax0.legend([v for k,v in s],[k for k,v in s], loc='lower center', 
                         bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize = plot_parameters["size_font_ticks"], borderaxespad=0.)
    
    ax0.add_artist(legend)
    
    #save the figure
    figurename = figurename + '.pdf'
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    
    print(figurename, '   created')
    #plt.show()
   
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])
