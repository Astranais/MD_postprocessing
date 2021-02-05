#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3


       program to  Plot abundance of coordinations for selected cation as a function of rho or T
       
      
Scientific color maps from:
    The software : Crameri, F. (2018a), Scientific colour-maps. Zenodo. http://doi.org/10.5281/zenodo.1243862
    The research : Crameri, F. (2018b), Geodynamic diagnostics, scientific visualisation and StagLab 3.0, Geosci. Model Dev., 11, 2541-2562, doi:10.5194/gmd-11-2541-2018
"""
#    ********* Importation of the packages and modules used here *********    
import sys, glob
import getopt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D #useful to create a custom legend
import crystallography as cr
import natsort
import re
from matplotlib.colors import LinearSegmentedColormap





def label_line(ax, line, label, halign, color='0.35', fs=12):
    """Add an annotation to the given line with appropriate placement and rotation.
    Based on code from:
        [How to rotate matplotlib annotation to match a line?]
        (http://stackoverflow.com/a/18800233/230468)
        User: [Adam](http://stackoverflow.com/users/321772/adam)
    Arguments
    ---------
    ax : `matplotlib.axes.Axes` object
        Axes on which the label should be added.
    line : `matplotlib.lines.Line2D` object
        Line which is being labeled.
    label : str
        Text which should be drawn as the label.
    ...
    Returns
    -------
    text : `matplotlib.text.Text` object
    """
    xdata, ydata = line.get_data()
#    x1 = xdata[0]
#    x2 = xdata[-1]
#    y1 = ydata[0]
#    y2 = ydata[-1]
    
    #formatage of label
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
                label = label[:i]+'$_{'+label[i+1:i+1+num]+'}$' + label[i+1+num:]
                i = i+5
            i = i+1
        else:break


    if halign.startswith('l'):
        x1 = xdata[0]
        x2 = xdata[1]
        y1 = ydata[0]
        y2 = ydata[1]
        xx = x1
        halign = 'left'
    elif halign.startswith('r'):
        x1 = xdata[-2]
        x2 = xdata[-1]
        y1 = ydata[-2]
        y2 = ydata[-1]
        xx = x2
        halign = 'right'
    elif halign.startswith('c'):        
        x1 = xdata[int(len(xdata)/2)-1]
        x2 = xdata[int(len(xdata)/2)]
        y1 = ydata[int(len(xdata)/2)-1]
        y2 = ydata[int(len(xdata)/2)]
        if ax.get_xscale() == 'log':
            xx = 10**(0.5*(np.log10(x1) + np.log10(x2)))
        else:
            xx = 0.5*(x1 + x2)
        halign = 'center'
    else:
        raise ValueError("Unrecogznied `halign` = '{}'.".format(halign))

    if ax.get_xscale() == 'log' and ax.get_yscale() == 'log':
        yy = 10**(np.interp(np.log10(xx), np.log10(xdata), np.log10(ydata)))
    elif ax.get_xscale() == 'log' and ax.get_yscale() != 'log':
        yy = np.interp(np.log10(xx), np.log10(xdata), ydata)
    else:
        yy = np.interp(xx, xdata, ydata)

    ylim = ax.get_ylim()
    # xytext = (10, 10)
    xytext = (0, 0)
    text = ax.annotate(label, xy=(xx, yy), xytext=xytext, textcoords='offset points',
                       size=fs, color=color, zorder=1,
                       horizontalalignment=halign, verticalalignment='center_baseline')

    sp1 = ax.transData.transform_point((x1, y1))
    sp2 = ax.transData.transform_point((x2, y2))

    rise = (sp2[1] - sp1[1])
    run = (sp2[0] - sp1[0])

    slope_degrees = np.degrees(np.arctan2(rise, run))
    text.set_rotation_mode('anchor')
    text.set_rotation(slope_degrees)
    ax.set_ylim(ylim)
    return text



def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.split('/')[-1]
    temperature = str(int(float(filename.split('_')[1].strip('T'))*1000))
    acell = filename.split('_')[3].strip('a')[:-1]
    #print(acell)
    return temperature, acell

def creation_plot2(variable, file1, file2, MN, plot_parameters):
    """     ********** Creation of the plot  **********    """
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    Na=6.022*10**23
    
    #plot
    plt.close()
    fig, (ax, ax2) = plt.subplots(1,2, sharey = True, sharex = True, figsize=size_figure)
    plt.subplots_adjust(top = 0.97, bottom = 0.12, right = 0.89, left = 0.07, hspace = 0, wspace = 0)
    ax.set_ylabel(r'Coordination polyhedra proportion', fontweight = 'bold', fontsize = size_fonts )

    for axis in [ax,ax2]:
        #Adjustment of ticks
        major_yticks = np.arange(0, 1.1, 0.1) 
        minor_yticks = np.arange(0, 1.05, 0.05)    
        axis.set_yticks(major_yticks)
        axis.set_yticks(minor_yticks, minor=True) 
        axis.yaxis.set_ticks_position('both')
        axis.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        axis.grid(True, which='major',axis = 'y', linestyle=':', linewidth=size_lines/2 )
        plt.autoscale(enable=True,axis='x',tight=False)
        axis.set_facecolor((1,1,1,0))
                
    ax.set_ylim(0,1)
    #ax.set_ylim(0,0.5)
    #Fine-tune figure and addition of a title per graph
    temp1, acell1 = split_name(file1)
    temp2, acell2 = split_name(file2)
    ax0 = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    if variable == 'rho':
        major_xticks = np.arange(0, 6.5, 0.5) 
        minor_xticks = np.arange(0, 6.5, 0.1)    
        for axis in [ax,ax2]:
            axis.set_xticks(major_xticks)
            axis.set_xticks(minor_xticks, minor=True)
            axis.xaxis.set_ticks_position('both')
            ax.set_xlim(0.9,6)
            ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        plt.setp(ax.get_xticklabels()[-1], visible=False)  
        ax0.set_xlabel(r"Density (g.cm$^{-3}$)", fontweight = 'bold', 
                       fontsize = size_fonts , labelpad = shift_labelpad*2)
        ax.set_xlabel(temp1+' K', fontweight = 'bold', fontsize = size_fonts )
        ax.xaxis.set_label_position('top')
        ax2.set_xlabel(temp2+' K', fontweight = 'bold', fontsize =size_fonts )
        ax2.xaxis.set_label_position('top')
    elif variable == 'T':
        major_xticks = np.arange(0, 10000, 1000) 
        for axis in [ax,ax2]:
            axis.set_xticks(major_xticks)
            axis.xaxis.set_ticks_position('both')
            ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax0.set_xlabel(r"Temperature (K)", fontweight = 'bold', 
                       fontsize = size_fonts , labelpad = shift_labelpad)
        ax.set_xlabel(str(round(MN/(Na*float(acell1)**3*10**(-24)),2))+r' (g.cm$^{-3}$)',
                      fontweight = 'bold', fontsize = size_fonts )
        ax.xaxis.set_label_position('top')
        ax2.set_xlabel(str(round(MN/(Na*float(acell2)**3*10**(-24)),2))+r' (g.cm$^{-3}$)', 
                       fontweight = 'bold', fontsize =size_fonts )
        ax2.xaxis.set_label_position('top')
    elif variable == 'P':
        for axis in [ax,ax2]:
            axis.set_xscale('log')
            axis.xaxis.set_ticks_position('both')
            ax.set_xlim(1,275)
            axis.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax0.set_xlabel(r"Pressure (GPa)", fontweight = 'bold', 
                       fontsize = size_fonts , labelpad = shift_labelpad*2)
        ax.set_xlabel(temp1+' K', fontweight = 'bold', fontsize = size_fonts )
        ax.xaxis.set_label_position('top')
        ax2.set_xlabel(temp2+' K', fontweight = 'bold', fontsize =size_fonts )
        ax2.xaxis.set_label_position('top')
    else:
        ax.set_xlabel(r"??? (???)", fontweight = 'bold', 
                      fontsize = size_fonts , labelpad = shift_labelpad)

    return fig, ax, ax2


def creation_plot(variable,plot_parameters):
    """     ********** Creation of the plot  **********    """
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    plt.close()
    fig, ax = plt.subplots(figsize=(12,7))
    ax.set_ylabel(r'Coordination polyhedra proportion', fontweight = 'bold', 
                  fontsize = size_fonts , labelpad = shift_labelpad)
    #Adjustment of ticks
    if variable == 'rho':
        major_xticks = np.arange(0, 7, 0.5) 
        minor_xticks = np.arange(0, 7, 0.1)    
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_xlabel(r"Density (g.cm$^{-3}$)", fontweight = 'bold', 
                      fontsize = size_fonts , labelpad = shift_labelpad)
    elif variable == 'T':
        major_xticks = np.arange(0, 10000, 1000) 
        ax.set_xticks(major_xticks)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_xlabel(r"Temperature (K)", fontweight = 'bold', 
                      fontsize = size_fonts , labelpad = shift_labelpad)
    elif variable == 'P':
        ax.set_xscale('log')
        ax.xaxis.set_ticks_position('both')
        ax.set_xlim(1,275)
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_xlabel(r"Pressure (GPa)", fontweight = 'bold', 
                      fontsize = size_fonts , labelpad = shift_labelpad)
    else:
        ax.set_xlabel(r"??? (???)", fontweight = 'bold', 
                      fontsize = size_fonts , labelpad = shift_labelpad)

    major_yticks = np.arange(0, 1.1, 0.1) 
    minor_yticks = np.arange(0, 1.05, 0.05)    
    ax.set_yticks(major_yticks)
    ax.set_yticks(minor_yticks, minor=True) 
    
    ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
    ax.grid(True, which='major',axis = 'y' , linestyle=':', linewidth=size_lines/2)
    #plt.autoscale(enable=True,axis='x',tight=False)
    ax.set_ylim(0,1)
    return fig,ax


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
    xdata = {'T':[],'rho':[],'P':[]} #dictionnary containing the x data
    data = [] #list containing y data
    xdata2 = {'T':[],'rho':[],'P':[]} #dictionnary containing the x data
    data2 = [] #list containing y data
    TP = {} #dictionnary with P and T in each simufile
    colors_species = {}
    #other dictionnaries and parameters for the figure for article version
    #plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 20}
    plot_parameters = {"size_fonts" : 16,"size_font_ticks":14,"size_figure" : (12,7),
                       "size_markers" : 10,"size_lines" : 2,"shift_labelpad" : 10}
    #dictionnary for fullaverages file in full version
    column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,
                     'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,
                     'Cvm':25,'stdev_Cvm':26,'testCv':31,'stdev_testCv':32}
    #other parameters
    Na=6.022*10**23
    statfile2 = ''
    try:
        options,arg = getopt.getopt(argv,"hf:g:a:v:m:",["file1","gfile2","atom","variable","mineralfile"])
    except getopt.GetoptError:
        print("plot_speciation-r0.py -a <first atom determining the type of cluster to plot> -v <variable (rho,P,T)> -m <mineralfile with elements> -f <stat-concentrate_r0_perc_filename> -g <stat-concentrate_perc_filename n°2 (option)>  ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_speciation-r0.py program to  Plot abundance of coordinations for selected cation as a function of rho,P  or T')
            print("plot_speciation-r0.py -a <first atom determining the type of cluster to plot> -v <variable (rho,P,T)> -m <mineralfile with elements> -f <stat-concentrate_r0_perc_filename> -g <stat-concentrate_perc_filename n°2 (option)> ")
            print('requires the file containing elements and number (in order to compute the densities)')
            print('')
            print('For plots as function of P, make sure the files from fullaverages.py are in the current folder')
            sys.exit()
        elif opt in ("-f", "--file1"):
            statfile = str(arg)
        elif opt in ("-g", "--gfile2"):
            statfile2 = str(arg)
        elif opt in ("-a","--atom"):
            atom = str(arg)
        elif opt in ("-v", "--variable"):
            variable = str(arg)
        elif opt in ("-m","--mineralfile"):
            mineralfile = str(arg)
    #***** Calculation of the molecular mass
    with open(mineralfile,'r') as mf:
        entry = mf.readline()
        elements = entry.split()[1:]
        entry = mf.readline()
        number = entry.split()[1:]
    MN = 0
    for i in range(len(elements)):
        MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
    #***** Creation of the plot
    if statfile2 != '':
        with open(statfile,'r') as f:
            line = f.readline()
            file1=line.split()[1]
        with open(statfile2,'r') as f2:
            line = f2.readline()
            file2=line.split()[1]
        fig, ax, ax2 = creation_plot2(variable, file1, file2, MN,plot_parameters)
        axrange = ax.get_xlim()
        ax2range = ax2.get_xlim()
    else:
        fig, ax = creation_plot(variable,plot_parameters)
        axrange = ax.get_xlim()
    #***** Count of the number of clusters for the selected atom type (needed for the automatic color change)
    #first store every cluster type into the dictionary
    with open(statfile,'r') as f:
        while True:
            line = f.readline()
            if not line: break
            else:
                entry=line.split('\n')[0].split('\t')
                if len(entry) > 1:
                    if entry[0][:len(atom)] == atom:
                        colors_species[entry[0]] = []
    mineralname = statfile.split('/')[-1].split('_')[0]
    #then, in case of a second file, we check if there is additional cluster and add them to the dictionary if applicable
    if statfile2 != '':
        with open(statfile2,'r') as f2:
            while True:
                line = f2.readline()
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    if len(entry) > 1:
                        if entry[0][:len(atom)] == atom:
                            if not entry[0] in colors_species:
                                colors_species[entry[0]] = []
    name = 'batlow'#'lajolla' #'imola'
    if 'Al_1O' in atom or 'Si_1O' in atom:
        #Creation of the color list for 11 coordination
        list_species = []
        cm_data = np.loadtxt("/Users/akobsch/Dropbox/Recherche/PhD-MD-silicates/simulations/scripts/Analysis/ScientificColourMaps6/"+name+"/"+name+".txt")
        #cm_data = cm_data[::-1] #for reverse colors (lajolla)
        new_map = LinearSegmentedColormap.from_list('new', cm_data)#[:-20]) #no cut for imola
        color = iter(new_map(np.linspace(0,1,11))) #Creation of the color list  
        #color = iter(plt.cm.jet(np.linspace(0,1,11)))
        for i in range(1,12,1):
            c = next(color)
            curr_species = atom+'_'+str(i)
            if curr_species in colors_species:
                colors_species[curr_species] = c
                list_species.append(curr_species)
        print(list_species, len(list_species))
    else:
        #Creation of the color list for 19 coordination
        list_species = []
        cm_data = np.loadtxt("/Users/akobsch/Dropbox/Recherche/PhD-MD-silicates/simulations/scripts/Analysis/ScientificColourMaps6/"+name+"/"+name+".txt")
        #cm_data = cm_data[::-1] #for reverse colors
        new_map = LinearSegmentedColormap.from_list('new', cm_data)
        color = iter(new_map(np.linspace(0,1,19))) #Creation of the color list  
        #color = iter(plt.cm.jet(np.linspace(0,1,19)))
        for i in range(1,20,1):
            c = next(color)
            curr_species = atom+'_'+str(i)
            if curr_species in colors_species:
                colors_species[curr_species] = c
                list_species.append(curr_species)
        print(list_species, len(list_species))        
    #*****************************
    #*******************
    #*******
    # Extraction of  P and T for each simufile of each thermo file
    thermofiles = sorted(glob.glob('thermo_'+mineralname+'*.txt'))
    for file in thermofiles:
        TP = extract_TP(file, column_number, TP,'')
    #*****************************
    #*******************
    #*******
    #Read & plot for the first file
    #read the file1
    with open(statfile,'r') as f:
        line = f.readline()
        entry=line.split()
        files = entry[1:]
        #creation of the x data list
        for file in files:
            temperature, acell = split_name(file)
            xdata['rho'].append(MN/(Na*float(acell)**3*10**(-24)) )    #calculation density
            xdata['T'].append(int(temperature))
            xdata['P'].append(TP[file.split('/')[-1]][1])
            #print(file.split('/')[-1],MN/(Na*float(acell)**3*10**(-24)),int(temperature),TP[file.split('/')[-1]][1])
    with open(statfile,'r') as f:
        #creation of the y data list and plot for each cluster
        while True:
            line = f.readline()
            if not line: break
            else:
                entry=line.split('\n')[0].split('\t')
                if len(entry) > 1:
                    if entry[0][:len(atom)] == atom:
                        for i in range(1,len(entry)): #we loop over all the files
                            if entry[i] == '':
                                entry[i] = 0
                            data.append(float(entry[i])) 
                        #now the data array is filled, we plot the line for this species
                        x, y = zip(*sorted(zip( xdata[variable], data)))
                        line, = ax.plot(x,y, '.--', color = colors_species[entry[0]], 
                                        markersize = plot_parameters["size_markers"], 
                                        linewidth = plot_parameters["size_lines"])
                        #if variable == 'T':
                        #    label_line(ax, line, entry[0], halign='right')   
                        #elif list_species.index(entry[0]) < int(len(list_species)/2):
                        #    label_line(ax, line, entry[0], halign='left')   
                        #else:
                        #    label_line(ax, line, entry[0], halign='right')
                        #and we plot the name of the species near the max of the curve
                        if max(data) > 0.05:
                            x_max = xdata[variable][data.index(max(data))]
                            if x_max > axrange[0] and x_max < axrange[1]:
                                y_max = max(data)
                                label = entry[0]
                                #**** formatage of label with removing _1 
                                label = format1label(label)
                                ax.text(x_max, y_max, label, color='0.35', 
                                        fontsize = plot_parameters['size_fonts'])
                        
                        
                        data = [] #reinitialization of data array for the next species
    figurename = statfile.split('.dat')[0]+'_'+atom+'-clusters_'+variable    
    #*****************************
    #*******************
    #*******
    #Same as before but for the second file if it exists
    if statfile2 != '':
        with open(statfile2,'r') as f2:
            line = f2.readline()
            entry=line.split()
            files = entry[1:]
            #creation of the x data list
            for file in files:
                temperature, acell = split_name(file)
                xdata2['rho'].append(MN/(Na*float(acell)**3*10**(-24)) )    #calculation density
                xdata2['T'].append(int(temperature))
                xdata2['P'].append(TP[file.split('/')[-1]][1])
        with open(statfile2,'r') as f2:
            #creation of the y data list and plot for each cluster
            while True:
                line = f2.readline()
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    if len(entry) > 1:
                        if entry[0][:len(atom)] == atom:
                            for i in range(1,len(entry)): #we loop over all the files
                                if entry[i] == '':
                                    entry[i] = 0
                                data2.append(float(entry[i]))
                            #now the data array is filled, we plot the line for this species
                            x2, y2 = zip(*sorted(zip( xdata2[variable], data2)))
                            line, = ax2.plot(x2,y2, '.--', color = colors_species[entry[0]],
                                             markersize = plot_parameters["size_markers"], 
                                             linewidth =plot_parameters["size_lines"])
                            #if variable == 'T':
                            #    label_line(ax2, line, entry[0], halign='right')  
                            #elif list_species.index(entry[0]) < int(len(list_species)/2):
                            #    label_line(ax2, line, entry[0], halign='left')  
                            #else:
                            #    label_line(ax2, line, entry[0], halign='right') 
                            #and we plot the name of the species near the max of the curve
                            if max(data2) > 0.05:
                                x2_max = xdata2[variable][data2.index(max(data2))]
                                if x2_max > ax2range[0] and x2_max < ax2range[1]:
                                    y2_max = max(data2)
                                    label2 = entry[0]
                                    #**** formatage of label with removing _1 
                                    label2 = format1label(label2)
                                    ax2.text(x2_max, y2_max, label2, color='0.35', 
                                             fontsize = plot_parameters['size_fonts'])
                            
                            data2 = [] #reinitialization of data array for the next species
                            
        figurename = statfile.split('/')[-1].split('.dat')[0]+'+'+statfile2.split('/')[-1].split('.dat')[0]+'_'+atom+'-clusters_'+variable
    #********* Legend
    #* formatage of label with removing _1 
    for j in range(len(list_species)):
        if list_species[j] != '':
            list_species[j] = format1label(list_species[j])                
    custom_lines = [Line2D([0],[0],color = colors_species[key], ls = '--', 
                           marker = '.', markersize = plot_parameters["size_markers"], 
                           linewidth = plot_parameters["size_lines"]) for key in natsort.natsorted(colors_species)]
    legend = plt.legend([line for line in custom_lines],[label for label in natsort.natsorted(list_species)],
                        title = '$\\bf{Polyhedra}$', bbox_to_anchor=(1.05, 1), loc=2, 
                        fontsize = plot_parameters["size_fonts"],  borderaxespad=0.)
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_fonts"])

    figurename = figurename+'.pdf'
    print(figurename, 'is created')
    plt.savefig(figurename, bbox_inches='tight', dpi=150)
#    plt.show() # à mettre après savefig sinon on ne sauvegarde rien !!  

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



