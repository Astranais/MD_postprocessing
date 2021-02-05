#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: akobsch

                    Plot of gofrs at 2 selected T for all densities or 2 selected rho for all T
                    + addition of a zoom for O-O
                    *************  ARTICLE VERSION  *************
"""

#*********** Importation of the packages and modules used here ************
import sys
import getopt
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter,AutoLocator,AutoMinorLocator
import crystallography as cr
import natsort
from matplotlib.colors import LinearSegmentedColormap




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
    """ Function to extract T and acell from filename """
    # *********** My filename are in the format CaAl2Si2O8_T3_nvt_a12.0.outcar.gofr.dat
    # ******* so I can split their name with _ and tale the acell and T from their name
    # **** Then change these two lines in harmony with your filename
    temperature = filename.split('_')[1]      #I extract the corresponding temperature
    acell = filename.split('.outcar.gofr.dat')[0].split('_')[3].strip('a')   #I extract the corresponding acell
    return temperature, acell


def creation_plot(plot_parameters,atom_couple):
    """     ********** Creation of the plot  **********    """
    print("I use creation plot without insert")
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    #plot
    plt.close()
    fig, (ax1, ax2) = plt.subplots(1,2, sharex = True, sharey = False, figsize=size_figure)
    #Adjustment of ticks
    major_xticks = np.arange(0, 9, 1) 
    minor_xticks = np.arange(0, 8.5, 0.5)
    for ax in [ax1,ax2]:
        majorLocator = AutoLocator()
        minorLocator = AutoMinorLocator()
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_minor_locator(minorLocator) 
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  #set one decimal to ticks label``  
        plt.autoscale(enable=True,axis='y',tight=False)
        ax.yaxis.set_ticks_position('both')
        ax.set_ylabel(r'g$_{\mathrm{\bf'+atom_couple+'}}$(r)', fontsize=size_fonts, fontweight = 'bold', labelpad = shift_labelpad)
        
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True) 
        ax.set_xlim([0,8])
        ax.set_xlabel('Distance r ($\AA$)',fontsize=size_fonts,fontweight='bold', labelpad = shift_labelpad)
         
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)         
    #Fine-tune figure
    plt.subplots_adjust(top = 0.97, bottom = 0.12, right = 0.89, left = 0.07, hspace = 0, wspace = 0.27)
    return fig, ax1, ax2


def creation_plot_insert(plot_parameters,atom_couple):
    """     ********** Creation of the plot  **********    """
    print("I use creation plot with insert")
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    #plot
    plt.close()
    fig, (ax1, ax2) = plt.subplots(1,2, sharex = True, sharey = False, figsize=size_figure)
    #Adjustment of ticks
    major_xticks = np.arange(0, 9, 1) 
    minor_xticks = np.arange(0, 8.5, 0.5)
    for ax in [ax1,ax2]:
        majorLocator = AutoLocator()
        minorLocator = AutoMinorLocator()
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_minor_locator(minorLocator) 
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  #set one decimal to ticks label``  
        plt.autoscale(enable=True,axis='y',tight=False)
        ax.yaxis.set_ticks_position('both')
        ax.set_ylabel(r'g$_{\mathrm{\bf'+atom_couple+'}}$(r)', fontsize=size_fonts, fontweight = 'bold', labelpad = shift_labelpad)
        
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True) 
        ax.set_xlim([0,8])
        ax.set_xlabel('Distance r ($\AA$)',fontsize=size_fonts,fontweight='bold', labelpad = shift_labelpad)
         
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)         
    
    
    #addition of a zoom
    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    left1, bottom1, width1, height1 = [0.265, 0.63, 0.15, 0.32]
    ax11 = fig.add_axes([left1, bottom1, width1, height1])    
    left2, bottom2, width2, height2 = [0.725, 0.63, 0.15, 0.32]
    ax22 = fig.add_axes([left2, bottom2, width2, height2])  
    
    major_xticks_zoom = np.arange(0, 3, 0.5) 
    minor_xticks_zoom = np.arange(0, 3, 0.1)
    major_yticks_zoom = np.arange(0, 0.6, 0.1) 
    minor_yticks_zoom = np.arange(0, 0.6, 0.05)
    
    for ax in [ax11,ax22]:
        ax.set_xticks(major_xticks_zoom)
        ax.set_xticks(minor_xticks_zoom, minor=True)
        ax.set_xticklabels(major_xticks_zoom)      #to choose to display only major tick labels
        
        ax.set_yticks(major_yticks_zoom)
        ax.set_yticks(minor_yticks_zoom, minor=True)                                           
    
        ax.xaxis.set_tick_params(labelsize = size_font_ticks, width = size_lines/2)
        ax.yaxis.set_tick_params(labelsize = size_font_ticks, width = size_lines/2)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  #set one decimal to ticks label``        
        
        ax.set_xlim([0.5,2.5])
        ax.set_ylim([0,0.2])    
    
    #Fine-tune figure
    plt.subplots_adjust(top = 0.97, bottom = 0.12, right = 0.89, left = 0.07, hspace = 0, wspace = 0.27)
    return fig, ax1, ax2, ax11, ax22


def format_T(temperature):
    T = str(round(float(temperature.strip('T'))*1000))              
    return T

def extraction(file, indexcol):
    """ extraction of data """
    distance, gofr = np.loadtxt(file, usecols=(0, indexcol), skiprows=1, unpack=True) 
    distance = distance[:-1]
    gofr = gofr[:-1]
    return distance, gofr

def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (10,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 5}
    #other dictionnaries and parameters for the figure
    lines = {}   #dictionnary for the lines for legend
    Na=6.022*10**23  #avogadro constant
    #other dictionnaries and parameters for the figure
    colors_densities = {}
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200','T5':'#ffcd01','T6':'#ff0101','T6.5':'#ff00a2','T7':'#ff01de'}
    #colors_T = create_colors()
    allT = {} #dictionnary for the legend
    letters = ['a','b']
    try:
        options,arg = getopt.getopt(argv,"hm:a:v:f:g:l:",["mineralfile","atoms","variable","folder1","gfolder2","letters"])
    except getopt.GetoptError:
        print("plot_allgofrs_2subplots.py  -m <mineralfile>  -a <1 couple of atoms>(ex: 'O-O') -v <variable for colors (rho,T)> -f <folder 1> -g <folder 2> -l <two letters>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print("plot_allgofrs_2subplots.py script to plot all the gofrs of the selected couple of atoms at every acell for 2 T or at every T for 2 acell")
            print("plot_allgofrs_2subplots.py for O-O, an insert is added  to zoom in the 0-2$\AA$ region")
            print("WARNING!!!!!")
            print("plot_allgofrs_2subplots.py requires to have separated folders of T and acell in order to plot all the files from the 2 T folders or from the 2 acell folders")
            print("")
            print("plot_allgofrs_2subplots.py -m <mineralfile>  -a <1 couple of atoms>(ex: 'O-O') -v <variable for colors (rho,T)>  -f <folder 1> -g <folder 2> -l <two letters>")
            print("")
            sys.exit()
        elif opt in ("-a", "--atoms"):
            atom_couple = str(arg)                     #list of atom couples we want to analyze here
        elif opt in ("-m","--mineralfile"):
            mineralfile = str(arg)
        elif opt in ("-v","--variable"):
            variable = str(arg)
        elif opt in ("-f","--folder1"):
            folder1 = str(arg).split('/')[0]
        elif opt in ("-g","--gfolder2"):
            folder2 = str(arg).split('/')[0]
        elif opt in ("-l","--letters"):
            letters = arg.split(',')
    #***** Calculation of the molecular mass
    with open(mineralfile,'r') as mf:
        entry = mf.readline()
        elements = entry.split()[1:]
        entry = mf.readline()
        number = entry.split()[1:]
    MN = 0
    for i in range(len(elements)):
        MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]  #I compute the molecular mass
    #***** Count of the number of densities or temperatures (needed for the automatic color change)
    #loop throught all files of the selected folder
    files1 = sorted(glob.glob(folder1+'/*outcar.gofr.dat')) #I list every gofr files from folder1 in alphabetic order
    files2 = sorted(glob.glob(folder2+'/*outcar.gofr.dat')) #I list every gofr files from folder1 in alphabetic order
    files = files1[:]
    files.extend(files2)
    print(files)
    files.reverse()
    if variable == 'rho':
        #we store every density into the dictionary colors_densities
        for file in files:
            #calculation of densities in g/cm3
            temperature, acell = split_name(file)  
            Volume = float(acell)**3
            Density = MN/(Na*Volume*10**(-24)) 
            #we store the densities into the dictionnary
            colors_densities[round(Density,1)] = []
    #now the dictionary has every possible densities or temperatures of our mineral, we attribute the colors to each density/temperature
    list_densities = []
    if variable == 'rho':
        color = iter(plt.cm.rainbow(np.linspace(0,1,len(colors_densities)))) #Creation of the color list
        for key in natsort.natsorted(colors_densities):
            c = next(color)
            colors_densities[key] = c
            list_densities.append(key) #list usefull if we want to plot the value of densities near to the lines, see plot_speciation_r0
        print(list_densities, "densities list for color bar. If you don't want densities but temperatures, please change the variable option to T")
    else:
        print([key for key in colors_T], "temperature list for color bar. If you don't want temperatures but densities, please change the variable option to rho")
    #*****************************
    #*******************
    #*******
    #Read & plot
    firstfile = files[0]                                                      #I take the first file  
    temperature0, acell0 = split_name(firstfile)
    with open(firstfile, 'r') as f:
        couples = f.readline()
        couples = couples.strip('dist')
        couples = re.sub('(Int\([A-Za-z-]*\))', ' ',couples).split()        #I substitute all Int(...)  by a blankspace and split using spaces
    col = 0
    for couple in couples:
        if (couple == atom_couple):
            if atom_couple == 'O-O':
                fig, ax1, ax2, ax11, ax22 = creation_plot_insert(plot_parameters,atom_couple)
                indexcol = couples.index(couple) + col + 1
                if variable == 'rho':
                    for file in files1:
                        #extraction of data
                        distance, gofr = extraction(file, indexcol)
                        #extraction of key for color legend 
                        temperature, acell = split_name(file)
                        Volume = float(acell)**3
                        Density = MN/(Na*Volume*10**(-24)) #calculation of densities in g/cm3
                        #plot
                        lines[str(round(Density,1))], = ax1.plot(distance,gofr, color = colors_densities[round(Density,1)], linewidth=plot_parameters['size_lines'])
                        ax11.plot(distance,gofr, color = colors_densities[round(Density,1)],linewidth=plot_parameters['size_lines'])
                    for file in files2:
                        #extraction of data
                        distance, gofr = extraction(file, indexcol)
                        #extraction of key for color legend 
                        temperature, acell = split_name(file)
                        Volume = float(acell)**3
                        Density = MN/(Na*Volume*10**(-24)) #calculation of densities in g/cm3
                        #plot
                        lines[str(round(Density,1))], = ax2.plot(distance,gofr, color = colors_densities[round(Density,1)], linewidth=plot_parameters['size_lines'])
                        ax22.plot(distance,gofr, color = colors_densities[round(Density,1)],linewidth=plot_parameters['size_lines'])
                else: #variable == T
                    for file in files1:
                        #extraction of data
                        distance, gofr = extraction(file, indexcol)
                        #extraction of key for color legend 
                        temperature, acell = split_name(file)
                        #plot
                        lines[format_T(temperature)], = ax1.plot(distance,gofr, color = colors_T[temperature], linewidth=plot_parameters['size_lines'])
                        ax11.plot(distance,gofr, color = colors_T[temperature],linewidth=plot_parameters['size_lines'])
                    for file in files2:
                        #extraction of data
                        distance, gofr = extraction(file, indexcol)
                        #extraction of key for color legend 
                        temperature, acell = split_name(file)
                        #plot
                        lines[format_T(temperature)], = ax2.plot(distance,gofr, color = colors_T[temperature], linewidth=plot_parameters['size_lines'])
                        ax22.plot(distance,gofr, color = colors_T[temperature],linewidth=plot_parameters['size_lines'])
            else: #atom_couple != O-O
                fig, ax1, ax2 = creation_plot(plot_parameters,atom_couple)
                indexcol = couples.index(couple) + col + 1
                if variable == 'rho':
                    for file in files1:
                        #extraction of data
                        distance, gofr = extraction(file, indexcol)
                        #extraction of key for color legend 
                        temperature, acell = split_name(file)
                        Volume = float(acell)**3
                        Density = MN/(Na*Volume*10**(-24)) #calculation of densities in g/cm3
                        #plot
                        lines[str(round(Density,1))], = ax1.plot(distance,gofr, color = colors_densities[round(Density,1)], linewidth=plot_parameters['size_lines'])
                    for file in files2:
                        #extraction of data
                        distance, gofr = extraction(file, indexcol)
                        #extraction of key for color legend 
                        temperature, acell = split_name(file)
                        Volume = float(acell)**3
                        Density = MN/(Na*Volume*10**(-24)) #calculation of densities in g/cm3
                        #plot
                        lines[str(round(Density,1))], = ax2.plot(distance,gofr, color = colors_densities[round(Density,1)], linewidth=plot_parameters['size_lines'])
                else: #variable == T
                    for file in files1:
                        #extraction of data
                        distance, gofr = extraction(file, indexcol)
                        #extraction of key for color legend 
                        temperature, acell = split_name(file)
                        #plot
                        lines[format_T(temperature)], = ax1.plot(distance,gofr, color = colors_T[temperature], linewidth=plot_parameters['size_lines'])
                    for file in files2:
                        #extraction of data
                        distance, gofr = extraction(file, indexcol)
                        #extraction of key for color legend 
                        temperature, acell = split_name(file)
                        #plot
                        lines[format_T(temperature)], = ax2.plot(distance,gofr, color = colors_T[temperature], linewidth=plot_parameters['size_lines'])
        col +=1        
    #article letter
    ax1.text(0.012,0.985, letters[0] , transform=ax1.transAxes, verticalalignment = 'top', horizontalalignment = 'left', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='none', edgecolor='k', pad=3.0))                
    ax2.text(0.012,0.985, letters[1] , transform=ax2.transAxes, verticalalignment = 'top', horizontalalignment = 'left', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='none', edgecolor='k', pad=3.0))    
    # Legend
    ax0 = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    s = [(k, lines[k]) for k in sorted(lines.keys())] #sort the dictionnary
    if variable == 'rho':
        title_legend=' $\\bf{Density}$ \n (g.cm$^{-3}$)'
    else:
        title_legend=' $\\bf{Temperature}$ (K)'    
        #title_legend=' $\\bf{Temperature}$ \n           (K)'    
    ncol = 8
    legend = ax0.legend([v for k,v in s[:ncol]],[k for k,v in s[:ncol]], title = title_legend, bbox_to_anchor=(0.5, 1.05),  loc="lower center", borderaxespad=0., fontsize = plot_parameters['size_font_ticks'], ncol = ncol)
    #legendbis = ax0.legend([v for k,v in s[ncol:]],[k for k,v in s[ncol:]], bbox_to_anchor=(0.5, 1.05),  loc="upper center", borderaxespad=0., fontsize = plot_parameters['size_font_ticks'], ncol = ncol)
    #legend = ax0.legend([v for k,v in s],[k for k,v in s], title = title_legend, bbox_to_anchor=(1.05, 1), loc='upper left", borderaxespad=0., fontsize = plot_parameters['size_font_ticks'])
    plt.setp(legend.get_title(),fontsize = plot_parameters['size_fonts'])
    ax0.add_artist(legend)
    #Save fig    
    filename = 'gofrs_'+files[0].split('/')[-1].split('_')[0]+'_'+folder1+'_'+folder2+'_'+atom_couple+'.pdf'
    plt.savefig(filename, bbox_inches="tight", dpi=300)
    print(filename,' is created')
    #plt.show()

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])





















