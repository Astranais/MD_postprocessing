#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the self diffusion coefficient of 4 different atoms       ****
                                        for 3 minerals
            and compare with values obtained from vibr_spectrum_umd.py and vibr2diffusion.py scripts
                      *************  ARTICLE VERSION  *************

"""


#     ********* Importation of the packages and modules used here *********
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import crystallography as cr
import re
import natsort


def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.outcar.msd.dat
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.strip('./')
    temperature = filename.split('_')[1]
    acell = filename.split('.outcar')[0].split('_')[3].strip('a')
    return temperature, acell


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


def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (10,10),
                       "size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    #other dictionnaries and parameters for the figure
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200',
                'T5':'#ffcd01','T5.5':'#ff6e00','T6':'#ff0101','T6.5':'#ff00a2','T7':'#ff01de'}
    xvariable = 'rho'
    #variables nedded fot the plot
    #other dictionnaries and parameters for the figure
    Temperatures = []
    #other parameters
    Na=6.022*10**23
    #dictionnary for fullaverages file in full version
    column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,
                     'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,
                     'Cvm':25,'stdev_Cvm':26,'testCv':31,'stdev_testCv':32}
    try:
        options,arg = getopt.getopt(argv,"hv:",['variable'])
    except getopt.GetoptError:
        print("plot_diffusion_all_compvibr.py -v <variable x axis (rho or P)> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_diffusion_all_compvibr.py program to plot diffusion coefficient as a function of density for each T and the selected minerals containing 4 elements')
            print("plot_diffusion_all_compvibr.py -v <variable x axis (rho or P)> ")
            print("plot_diffusion_all_compvibr.py requires to be lauched from the folder containing every diffusivities file created by the script analyze_msd and vibr2diffusion")
            print('')
            print('For plots as function of P, make sure the files from fullaverages.py are in the current folder')
            sys.exit()
        if opt in ('-v','--variable'):
            xvariable = str(arg)
    figurename = 'diffusivities_all_compvibr_'+xvariable
    if xvariable == 'rho':
        major_xticks = np.arange(0, 4.5, 0.5) 
        minor_xticks = np.arange(0, 4.1, 0.1) 
        label = r'Density (g.cm$^{-3}$)'
    else:
        label = r'Pressure (GPa)' 
                
    minerals = ['NaAlSi3O8','KAlSi3O8','CaAl2Si2O8']
    

    plt.close(1)
    fig = plt.figure(1,figsize=plot_parameters['size_figure'])
    plt.subplots_adjust(top = 0.97, bottom = 0.07, right = 0.89, left = 0.07, hspace = 0, wspace = 0)
    plt.subplot(4,3,1) #4 inles 3 columns
    ax = plt.gca()
    for mineral in minerals:
        #************ Extract Data
        diffusion_msd_file = 'diffusivities_'+mineral+'.txt'
        diffusion_vibr_file = 'diffusivities-vibr_'+mineral+'.txt'
        #** creation of elements and number lists
        with open(diffusion_vibr_file,'r') as f:
            skiplines = 0
            while True:
                line = f.readline()
                skiplines += 1
                entry=line.split()
                if entry[0] == 'elements':
                    elements = entry[1:]
                elif entry[0] == 'number':
                    number = entry[1:]
                elif entry[0] == 'file':
                    break
        #** calculation of M*N nedded for the calculation of densities
        MN = 0
        for i in range(len(elements)):
            MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
        #** Extraction of data in diffusivities-vibr_mineral.txt file
        datavibr = {}
        Xvibr = {}
        #initialisation of data
        datavibr = {} #big dictionnary with inside self diffusion coefficient or time of regime change, all coresponding to different T and atom pairs
        Xvibr = {}  #idem for rho or P
        for i in range(len(elements)):
            elem = elements[i]
            datavibr[elem] = {}
            Xvibr[elem] = {}
        #fill dictionnaries
        with open(diffusion_vibr_file,'r') as f:
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
                            datavibr[elem][temp].append(float(entry[i+2]))
                            if xvariable == 'P':
                                Xvibr[elem][temp].append(float(entry[1]))
                            else:
                                Xvibr[elem][temp].append( MN/(Na*float(acell)**3*10**(-24)) )    #calculation density
                        except KeyError:
                            datavibr[elem][temp] = [ float(entry[i+2])]
                            if xvariable == 'P':
                                Xvibr[elem][temp] = [float(entry[1])]
                            else:
                                Xvibr[elem][temp] = [ MN/(Na*float(acell)**3*10**(-24)) ]    #calculation density         
        #** Extract P and T from thermofile
        if xvariable == 'P':
            TP = {}
            thermofile = 'thermo_'+mineral+'_all.txt'
            TP = extract_TP(thermofile, column_number, TP,'')
            if TP == {}:
                print('ERROR!!!! TP dictionnary empty')
        #** Extraction of data in diffusivities_mineral.txt file
        dataMSD = {}
        XMSD = {}
        #count number lines in header
        with open(diffusion_msd_file,'r') as f:
            skiplines = 0
            while True:
                line = f.readline()
                skiplines += 1
                entry=line.split()
                if entry[0] == 'file':
                    break
        #initialisation of data
        dataMSD = {} #big dictionnary with inside self diffusion coefficient or time of regime change, all coresponding to different T and atom pairs
        XMSD = {}  #idem for rho or P
        for i in range(len(elements)):
            elem = elements[i]
            dataMSD[elem] = {}
            XMSD[elem] = {}
        #fill dictionnaries
        with open(diffusion_msd_file,'r') as f:
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
                            dataMSD[elem][temp].append(float(entry[i*6+1]))
                            if xvariable == 'P':
                                XMSD[elem][temp].append(TP[entry[0].split('outcar.msd.dat')[0].split('/')[-1]][1])
                            else:
                                XMSD[elem][temp].append( MN/(Na*float(acell)**3*10**(-24)) )    #calculation density
                        except KeyError:
                            dataMSD[elem][temp] = [ float(entry[i*6+1])]
                            if xvariable == 'P':
                                XMSD[elem][temp] = [TP[entry[0].split('outcar.msd.dat')[0].split('/')[-1]][1]]
                            else:
                                XMSD[elem][temp] = [ MN/(Na*float(acell)**3*10**(-24)) ]    #calculation density 
        #************ Plot Data
        for elem in elements:
            #print("******** ", elem)
            #change of subplot
            plt.subplot(4,3,elements.index(elem)*3+1+minerals.index(mineral))
            ax = plt.gca()
            #plot
            for temp in datavibr[elem]:
                #print("******** ", temp)
                #put only the T you want in colors_T (no T10, 15) 
                try:
                    ax.plot(Xvibr[elem][temp],datavibr[elem][temp], ls = ':', 
                        marker = 'x', markersize = plot_parameters["size_markers"],
                        color = colors_T[temp], linewidth = plot_parameters["size_lines"])
                    ax.plot(XMSD[elem][temp],dataMSD[elem][temp], ls = '-', 
                        marker = '*', markersize = plot_parameters["size_markers"],
                        color = colors_T[temp], linewidth = plot_parameters["size_lines"])
                    Temperatures.append(temp)
                except KeyError:
                    pass
            #we add x and y labels outside the plot on the right columns and lines
            if elements.index(elem) == 0:
                ax.set_xlabel(format1label(mineral), fontweight = 'bold', fontsize = plot_parameters["size_fonts"])
                ax.xaxis.set_label_position('top')
            if minerals.index(mineral) == len(minerals)-1:
                if elements.index(elem) == 0:
                    ax.set_ylabel('Na, K, Ca', fontweight = 'bold', fontsize = plot_parameters["size_fonts"])
                else:
                    ax.set_ylabel(elem, fontweight = 'bold', fontsize = plot_parameters["size_fonts"])
                ax.yaxis.set_label_position('right')
            #limitation of data along x
            if xvariable == 'rho':
                ax.set_xticks(major_xticks)
                ax.set_xticks(minor_xticks, minor=True)
                ax.set_xlim(1.0,4.1)
                if minerals.index(mineral) != len(minerals)-1:
                    plt.setp(ax.get_xticklabels()[-1], visible=False)                
            else:
                ax.set_xscale('log')
                ax.set_xlim(1,275)               
            #limitation of data along y            
            ax.set_yscale('log')
            ax.set_ylim(4e-10,3e-7)
            #we remove unwanted axis tick labels
            if  minerals.index(mineral) != 0:
                plt.setp(ax.get_yticklabels(), visible=False)
            if elements.index(elem) != len(elements)-1:
                plt.setp(ax.get_xticklabels(), visible=False)
            #we make the graph prettier
            ax.xaxis.set_ticks_position('both')
            ax.yaxis.set_ticks_position('both')
            ax.yaxis.set_tick_params(which = 'both', direction='inout')
            ax.xaxis.set_tick_params(which = 'both', direction='inout')            
            ax.grid(axis = 'y', which = 'major', linestyle = '--', 
                    linewidth = plot_parameters["size_lines"]/2, alpha = 0.5)
            ax.tick_params(which = 'both', labelsize = plot_parameters["size_font_ticks"],
                           width = plot_parameters["size_lines"]/2)
            ax.set_facecolor((1,1,1,0))

            #plt.setp(ax.get_yticklabels()[-1], visible=False)
    
    #************ Fine tune figure
    ax0 = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax0.set_xlabel(label, fontweight = 'bold', fontsize = plot_parameters["size_fonts"],
                   labelpad = plot_parameters["shift_labelpad"]*2)
    ax0.set_ylabel(r'Diffusion coefficient (m$^2$.s$^{-1}$)', fontweight = 'bold', 
                   fontsize = plot_parameters["size_fonts"], 
                   labelpad = plot_parameters["shift_labelpad"]*3+plot_parameters["shift_labelpad"]/2)
    
    #************ Create legend from custom artist/label lists
    Temperatures = list(set(Temperatures)) #get elements only once in the list
    custom_patch = [mpatches.Patch(color=colors_T[key]) for key in natsort.natsorted(Temperatures)]
    
    legend_labels = {'D from MSD': plt.Line2D((0,1),(0,0), color='k', markersize = plot_parameters["size_markers"],
                                    marker='*', linestyle='-', linewidth = plot_parameters["size_lines"]),
                     'D from vibrational spectrum': plt.Line2D((0,1),(0,0), color='k', markersize = plot_parameters["size_markers"],
                                    marker='x', linestyle='--', linewidth = plot_parameters["size_lines"])}
    s = [(k, legend_labels[k]) for k in sorted(legend_labels.keys(),reverse = False)]
    
    legend = ax0.legend([col for col in custom_patch],[str(int(float(label.strip('T'))*1000)) for label in natsort.natsorted(Temperatures)],title = '$\\bf{Temperature~(K)}$', bbox_to_anchor=(0.5, 1.08), loc="lower center", fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol=len(Temperatures))    
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])
    plt.legend([v for k,v in s],[k for k,v in s], bbox_to_anchor=(0.5, 1.03), loc='lower center',fancybox=True, fontsize = plot_parameters["size_fonts"], ncol=len(legend_labels))
    ax0.add_artist(legend)
    
    
    figurename = figurename + '.pdf'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












