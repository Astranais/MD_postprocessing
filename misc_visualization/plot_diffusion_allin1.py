#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the self diffusion coefficient of 4 different atoms       ****
                      *************  ARTICLE VERSION  *************

"""


#     ********* Importation of the packages and modules used here *********
import glob
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend import Legend
import crystallography as cr
import re
import natsort


def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.outcar.msd.dat
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.strip('./')
    temperature = filename.split('_')[1]
    acell = filename.split('.outcar.msd.dat')[0].split('_')[3].strip('a')
    return temperature, acell



def creation_plot_3(plot_parameters,xvariable):
    """     ********** Creation of the plot  **********    """
    print("I use creation plot 3")
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    #plot
    plt.close()
    f, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, sharey=True, 
       figsize = (size_figure[0]*2.5,size_figure[1]))
    
    if xvariable == 'rho':
        major_xticks = np.arange(0, 4.5, 0.5) 
        minor_xticks = np.arange(0, 4.1, 0.1) 
        for ax in [ax1,ax2,ax3]:
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True)
        ax1.set_xlim(1.0,4.1)
        label = r'Density (g.cm$^{-3}$)'
    else:
        ax1.set_xlim(1,275)
        ax1.set_xscale('log')
        label = r'Pressure (GPa)'
        
    for ax in [ax1,ax2,ax3]:
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_yscale('log')                                          
        ax.yaxis.set_ticks_position('both')
    for ax in [ax1,ax2]:
        plt.setp(ax.get_xticklabels()[-1], visible=False) 
        #ax.grid(True, axis = 'x')
    
    ax1.set_ylim(4e-10,3e-7)

    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, 
                      hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax.set_xlabel(label, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    ax.set_ylabel(r'Diffusion coefficient (m$^2$.s$^{-1}$)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad*3+shift_labelpad/2)
    return f, ax1, ax2, ax3, ax


def creation_plot_2(plot_parameters,xvariable):
    """     ********** Creation of the plot  **********    """
    print("I use creation plot 2")
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    #plot
    plt.close()
    f, (ax1, ax2) = plt.subplots(2,1, sharex=True, sharey=True, figsize = (size_figure[0],size_figure[1]*2))
    
    if xvariable == 'rho':
        major_xticks = np.arange(0, 4.5, 0.5) 
        minor_xticks = np.arange(0, 4.1, 0.1) 
        for ax in [ax1,ax2]:
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True)
        ax1.set_xlim(1.0,4.1)
        label = r'Density (g.cm$^{-3}$)'
    else:
        ax1.set_xlim(1,275)
        ax1.set_xscale('log')
        label = r'Pressure (GPa)'

    for ax in [ax1,ax2]:
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_yscale('log')                                          
        ax.yaxis.set_ticks_position('both')
    for ax in [ax1]:
        plt.setp(ax.get_yticklabels()[-1], visible=False) 
        #ax.grid(True, axis = 'x')
    ax1.set_ylim(4e-10,3e-7)

    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, 
                      hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax.set_xlabel(label, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    ax.set_ylabel(r'Diffusion coefficient (m$^2$.s$^{-1}$)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad*3+shift_labelpad/2)
    return f, ax1, ax2, ax

def creation_plot(plot_parameters,xvariable):
    """     ********** Creation of the plot  **********    """
    print("I use creation plot")
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    #plot
    plt.close()
    fig, ax = plt.subplots(1,1, figsize = size_figure, linewidth = size_lines /2)
    
    if xvariable == 'rho':
        major_xticks = np.arange(0, 7, 0.5) 
        minor_xticks = np.arange(0, 7, 0.1)    
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.set_xlim(1.0,4.1)
        label = r'Density (g.cm$^{-3}$)'
    else:
        ax.set_xlim(1,275)
        ax.set_xscale('log')
        label = r'Pressure (GPa)'
        
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
    ax.set_yscale('log')                                          
    ax.yaxis.set_ticks_position('both') 
    #ax.grid(True, axis = 'x')
    ax.set_ylim(4e-10,3e-7)
    #plt.autoscale()
    
    ax.set_xlabel(label, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    ax.set_ylabel(r'Diffusion coefficient (m$^2$.s$^{-1}$)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad)
    return fig, ax


def format1label(label):
    """formatage of compound label """
    i=0
    while True:
        if i <= len(label)-1:
            if re.match('[0-9]',label[i]):
                num=0
                try:
                    while re.match('[0-9]',label[i+1+num]):
                        num +=1
                except IndexError: #index error when we arrive at the end of the cluster name
                    pass
                    #print('end of the cluster')
                label = label[:i]+'$_{'+label[i:i+1+num]+'}$' + label[i+1+num:]
                i = i+5
            i = i+1
        else:break
    return label


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
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),
                       "size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    #other dictionnaries and parameters for the figure
    style_markers = {'Ca':'','K':'','Na':'','Al':'x','Si':'*','O':''}
    style_lines = {'Ca':'--','K':'--','Na':'--','Al':'--','Si':'-','O':'-'}
    legend_labels = {}
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200',
                'T5':'#ffcd01','T6':'#ff0101','T6.5':'#ff00a2','T7':'#ff01de'}
    filename = 'all'
    filename2 = ''
    atoms = 'all'
    xvariable = 'rho'
    #variables nedded fot the plot
    data = {} #dictionnary containing the data for each element, initialized for each T
    stdev = {} #same for stdev
    #other dictionnaries and parameters for the figure
    ions = []
    compounds = []
    Temperatures = []
    #other parameters
    Na=6.022*10**23
    #dictionnary for fullaverages file in full version
    column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,
                     'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,
                     'Cvm':25,'stdev_Cvm':26,'testCv':31,'stdev_testCv':32}
    try:
        options,arg = getopt.getopt(argv,"hf:g:a:v:",["filename","gfilename",'atoms','variable'])
    except getopt.GetoptError:
        print("plot_diffusion_allin1.py -f <filename>(default = all) -g <filename 2 (if we want to plot 2 files)> -a <atoms list, default = 'all'> -v <variable x axis (rho or P)> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_diffusion_allin1.py program to plot diffusion coefficient as a function of density for each T and the selected minerals containing 4 elements')
            print("plot_diffusion_allin1.py -f <filename>(default = all files of every compound) -g <filename 2 (if we want to plot 2 files)>  -a <atoms list, default = 'all'> -v <variable x axis (rho or P)> ")
            print("plot_diffusion_allin1.py requires to be lauched from the folder containing every diffusivities file created by the script analyze_msd")
            print('')
            print('For plots as function of P, make sure the files from fullaverages.py are in the current folder')
            sys.exit()
        if opt in ('-f','--filename'):
            filename = str(arg)
        elif opt in ('-g','--gfilename'):
            filename2 = str(arg)
        elif opt in ('-a','--atoms'):
            atoms = str(arg)
        elif opt in ('-v','--variable'):
            xvariable = str(arg)
    if filename == 'all':
        fig,ax1,ax2,ax3, ax0 = creation_plot_3(plot_parameters,xvariable)
        print("axis are ax1",ax1,"ax2",ax2,"ax3",ax3,"ax0",ax0)
        files = sorted(glob.glob('diffusivities_*.txt'),reverse=True) #I list every diffusivities files
        figurename = 'diffusivities_all_'+atoms
    elif filename2 != '':
        fig,ax1,ax2, ax0 = creation_plot_2(plot_parameters,xvariable)
        print("axis are ax1",ax1,"ax2",ax2,"ax0",ax0)
        files = [filename, filename2] #we will plot 2 files
        figurename = filename.split('.txt')[0] +'_'+ filename2.split('.txt')[0].split('_')[1]+'_'+atoms
    else:
        fig, ax = creation_plot(plot_parameters,xvariable)
        files = [filename] #I take only the file we want
        figurename = filename.split('.txt')[0]+'_'+atoms
    #************ creation of arrays and plot at the same time
    for file in files:
        #********initialisations
        #**change of subplot
        if filename == 'all' or filename2 !='':
            if files.index(file) == 0:
                ax=ax1
                letter = 'a'
            if files.index(file) == 1:
                ax=ax2
                letter = 'b'
            if files.index(file) == 2:
                ax=ax3 
                letter = 'c'
            print("I plot on axis",ax)
        #**extraction compound
        compound = file.split('_')[1].split('.txt')[0]  
        compounds.append(format1label(compound))
        #******* Extract P and T for each thermo file
        if xvariable == 'P':
            TP = {}
            mineralname = file.split('_')[1].split('.txt')[0]
            thermofiles = sorted(glob.glob('thermo_'+mineralname+'*.txt'))
            for thermofile in thermofiles:
                TP = extract_TP(thermofile, column_number, TP,'')
            if TP == {}:
                print('ERROR!!!! TP dictionnary empty')
        #**creation of elements and number lists and initialization of T
        with open(file,'r') as f:
            skiplines = 0
            while True:
                line = f.readline()
                skiplines += 1
                entry=line.split()
                if entry[0] == 'elements':
                    elements = entry[1:]
                    ions.append(elements[0])
                elif entry[0] == 'number':
                    number = entry[1:]
                elif entry[0] == 'file':
                    line = f.readline()
                    entry=line.split()
                    temperature0, acell0 = split_name(entry[0])
                    break
        #**calculation of M*N nedded for the calculation of densities
        MN = 0
        for i in range(len(elements)):
            MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
        #**initialisation of data
        data = {} #big dictionnary with inside self diffusion coefficient or time of regime change, all coresponding to different T and atom pairs
        stdev = {} #idem for stdev
        X = {}  #idem for rho or P
        for i in range(len(elements)):
            elem = elements[i]
            stdev[elem] = {}
            data[elem] = {}
            X[elem] = {}
            #we store all the ions only once in ions list
            indicator = 0
            for j in range(len(ions)):
                if elements[i] == ions[j]:
                    indicator = 1
            if indicator == 0:
                ions.append(elements[i])
        #****** Extraction of data
        with open(file,'r') as f:
            [f.readline() for i in range(skiplines)]
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    temp, acell = split_name(entry[0]) 
                    for i in range(len(elements)): 
                        elem = elements[i]
                        try:
                            data[elem][temp].append(float(entry[i*6+1]))
                            stdev[elem][temp].append(float(entry[i*6+2]))    
                            if xvariable == 'P':
                                X[elem][temp].append(TP[entry[0].split('outcar.msd.dat')[0].split('/')[-1]][1])
                            else:
                                X[elem][temp].append( MN/(Na*float(acell)**3*10**(-24)) )    #calculation density
                        except KeyError:
                            data[elem][temp] = [ float(entry[i*6+1])]
                            stdev[elem][temp] = [ float(entry[i*6+2]) ]
                            if xvariable == 'P':
                                X[elem][temp] = [TP[entry[0].split('outcar.msd.dat')[0].split('/')[-1]][1]]
                            else:
                                X[elem][temp] = [ MN/(Na*float(acell)**3*10**(-24)) ]    #calculation density 
        #****** Selection of elements we plot
        if atoms != 'all':
            elements = atoms.split(',')
            print('now elements list is',elements)
        #****** Plot
        for elem in elements:
            for temperature in data[elem]:
                print("******** ", temperature)
                print(data[elem][temperature])
                Temperatures.append(temperature)
                ax.plot(X[elem][temperature],data[elem][temperature], 
                        style_markers[elem]+style_lines[elem], markersize = plot_parameters["size_markers"], 
                        color = colors_T[temperature], linewidth = plot_parameters["size_lines"])
        #ax.text(0.95,0.85, compound , transform=ax.transAxes, horizontalalignment = 'right', fontweight = 'bold', fontsize = plot_parameters["size_fonts"])
        if filename == 'all':
            ax.text(0.985,0.95, letter , transform=ax.transAxes, horizontalalignment = 'right',
                    fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                    bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
        elif filename2 != '':
            ax.text(0.985,0.95, letter , transform=ax.transAxes, horizontalalignment = 'right', 
                    fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                    bbox=dict(facecolor='none', edgecolor='k', pad=3.0))      
    #****** Create legend from custom artist/label lists
    #********* Legend               
    #Create legend from custom artist/label lists
    Temperatures = list(set(Temperatures)) #get elements only once in the list
    custom_patch = [mpatches.Patch(color=colors_T[key]) for key in natsort.natsorted(Temperatures)]
    
    if atoms != 'all':
        ions = elements
    print(ions) 
    for ion in ions:
        legend_labels[ion] = plt.Line2D((0,1),(0,0), color='k', markersize = plot_parameters["size_markers"], 
                     marker=style_markers[ion], linestyle=style_lines[ion], linewidth = plot_parameters["size_lines"])
    s = [(k, legend_labels[k]) for k in sorted(legend_labels.keys(),reverse = False)]
    
    if filename == 'all':
        ax=ax0
        legend = ax0.legend([col for col in custom_patch],[str(int(float(label.strip('T'))*1000)) for label in natsort.natsorted(Temperatures)],title = '$\\bf{Temperature~(K)}$', bbox_to_anchor=(0.5, 1.12), loc="lower center", fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol=len(Temperatures))    
        plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])
        plt.legend([v for k,v in s],[k for k,v in s], bbox_to_anchor=(0.5, 1), 
                   loc='lower center',fancybox=True, fontsize = plot_parameters["size_fonts"], ncol=len(style_lines))
        ax0.add_artist(legend)
    elif filename2 != '':
        ax=ax0
        legend = ax0.legend([col for col in custom_patch],[str(int(float(label.strip('T'))*1000)) for label in natsort.natsorted(Temperatures)],title = '$\\bf{Temperature~(K)}$', bbox_to_anchor=(0.5, 1.12), loc="lower center", fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol=len(Temperatures))    
        plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])
        plt.legend([v for k,v in s],[k for k,v in s], bbox_to_anchor=(0.5, 1), 
                   loc='lower center',fancybox=True, fontsize = plot_parameters["size_fonts"], ncol=len(style_lines))
        ax0.add_artist(legend)
    figurename = figurename + '.pdf'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












