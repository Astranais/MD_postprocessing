#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3


        Plot species % in the gas phase for any type of clusters as a function of rho or T
                *************  ARTICLE VERSION  *************

"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D #useful to create a custom legend
from matplotlib.colors import LinearSegmentedColormap
import crystallography as cr
import natsort
import re


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

def format_label(list_species):
    """ formatage of all labels of a list """
    for j in range(len(list_species)):
        if list_species[j] != '':
            list_species[j] = format1label(list_species[j])
    return list_species

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
    label = format1label(label)


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

def color_change(color_scheme, color_list, color_dict):
    """apply colors to values of a dictionnary by folowing the color_list and color_scheme """
    for key in natsort.natsorted(color_list):
        c = next(color_scheme)
        color_dict[key] = c
    return color_dict

def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.split('/')[-1]
    temperature = str(int(float(filename.split('_')[1].strip('T'))*1000))
    acell = filename.split('_')[3].strip('a')[:-1]
    #print(acell)
    return temperature, acell


def creation_plot2(variable, file1, file2, MN, plot_parameters,max_den,letter):
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
    ax.set_ylabel(r'Species proportion', fontweight = 'bold', fontsize = size_fonts )

    for axis in [ax,ax2]:
        #Adjustment of ticks
        major_yticks = np.arange(0, 1.1, 0.1) 
        minor_yticks = np.arange(0, 1.05, 0.05)    
        axis.set_yticks(major_yticks)
        axis.set_yticks(minor_yticks, minor=True) 
        axis.yaxis.set_ticks_position('both')
        axis.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        axis.grid(True, which='major',axis = 'y', linestyle=':', linewidth=size_lines/2 )
        axis.set_facecolor((1,1,1,0))
                
    ax.set_ylim(0,1)
    #Fine-tune figure and addition of a title per graph
    temp1, acell1 = split_name(file1)
    temp2, acell2 = split_name(file2)
    ax0 = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    if variable == 'rho':
        major_xticks = np.arange(0, 7, 0.1) 
        minor_xticks = np.arange(0, 7, 0.05)    
        for axis in [ax,ax2]:
            axis.set_xticks(major_xticks)
            axis.set_xticks(minor_xticks, minor=True)
            axis.xaxis.set_ticks_position('both')
            ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
            axis.set_xlim(1,max_den)
        ax0.set_xlabel(r"Density (g.cm$^{-3}$)", fontweight = 'bold', 
                       fontsize = size_fonts , labelpad = shift_labelpad)
        ax.set_xlabel(temp1+' K', fontweight = 'bold', fontsize = size_fonts )
        ax.xaxis.set_label_position('top')
        ax2.set_xlabel(temp2+' K', fontweight = 'bold', fontsize =size_fonts )
        ax2.xaxis.set_label_position('top')
    if variable == 'T':
        major_xticks = np.arange(0, 10000, 1000) 
        for axis in [ax,ax2]:
            axis.set_xticks(major_xticks)
            axis.xaxis.set_ticks_position('both')
            ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
            plt.autoscale(enable=True,axis='x',tight=False)
        ax0.set_xlabel(r"Temperature (K)", fontweight = 'bold', 
                       fontsize = size_fonts , labelpad = shift_labelpad)
        ax.set_xlabel(str(round(MN/(Na*float(acell1)**3*10**(-24)),2))+r' (g.cm$^{-3}$)', 
                      fontweight = 'bold', fontsize = size_fonts )
        ax.xaxis.set_label_position('top')
        ax2.set_xlabel(str(round(MN/(Na*float(acell2)**3*10**(-24)),2))+r' (g.cm$^{-3}$)', 
                       fontweight = 'bold', fontsize =size_fonts )
        ax2.xaxis.set_label_position('top')
        
    if letter != '':
        ax.text(-0.18,0.99, letter , transform=ax.transAxes, horizontalalignment = 'left', 
                fontweight = 'bold', fontsize = plot_parameters["size_fonts"],
                bbox=dict(facecolor='none', edgecolor='k', pad=3.0))  
    
    return fig, ax, ax2


def creation_plot(variable,plot_parameters,max_den,letter):
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
    ax.set_ylabel(r'Species proportion', fontweight = 'bold', 
                  fontsize = size_fonts , labelpad = shift_labelpad)
    #Adjustment of ticks
    if variable == 'rho':
        major_xticks = np.arange(0, 7, 0.1) 
        minor_xticks = np.arange(0, 7, 0.05)    
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_xlabel(r"Density (g.cm$^{-3}$)", fontweight = 'bold', 
                      fontsize = size_fonts , labelpad = shift_labelpad/2)
        ax.set_xlim(1,max_den)
    if variable == 'T':
        major_xticks = np.arange(0, 10000, 1000) 
        ax.set_xticks(major_xticks)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_xlabel(r"Temperature (K)", fontweight = 'bold', 
                      fontsize = size_fonts , labelpad = shift_labelpad/2)
        plt.autoscale(enable=True,axis='x',tight=False)

    major_yticks = np.arange(0, 1.1, 0.1) 
    minor_yticks = np.arange(0, 1.05, 0.05)    
    ax.set_yticks(major_yticks)
    ax.set_yticks(minor_yticks, minor=True) 
    
    ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
    ax.grid(True, which='major',axis = 'y' , linestyle=':', linewidth=size_lines/2)
    ax.set_ylim(0,1)
    
    if letter != '':
        ax.text(-0.08,0.99, letter , transform=ax.transAxes, horizontalalignment = 'left',
                fontweight = 'bold', fontsize = plot_parameters["size_fonts"],
                bbox=dict(facecolor='none', edgecolor='k', pad=3.0))  
        
    return fig,ax




def count_natom(species):
    """ count the number of atoms in a species (works if numbers are indicated by _)"""
    natom = 0
    if species != '':
        i=0
        while True:
            if i <= len(species)-1:
                #print('species[',i,'] = ', species[i])
                if re.match('_',species[i]):
                    num=0
                    value = ''
                    try:
                        while re.match('[0-9]',species[i+1+num]):
                            value = value + species[i+1+num]
                            num +=1
                        natom = natom+int(value)
                    except IndexError: #index error when we arrive at the end of the cluster name
                        #print('end of the cluster')
                        natom = natom+int(value)
                        pass
                i = i+1
            else:break
    return natom


def main(argv):
    """     ********* Main program *********     """
    nmeltatoms = 100 #limit of natoms to consider a species being part of the melt 
    #other dictionnaries and parameters for the figure for article version    
    letter = '' 
    statfile2 = ''
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),
                       "size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 20}
    #plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),"size_markers" : 10,"size_lines" : 2,"shift_labelpad" : 20}
    #other parameters
    Na=6.022*10**23
    try:
        options,arg = getopt.getopt(argv,"hf:g:v:m:d:l:",["file1","gfile2","variable","mineralfile","density_max","letter"])
    except getopt.GetoptError:
        print("plot_speciation-r1-species.py  -v <variable (rho,T)> -m <mineralfile with elements> -f <stat-concentrate_r1_abso_filename>  -g <stat-concentrate_r1_abso_filename2>  -d <maximum density to plot in g/cm3> -l <letter for article subplot, default = ''>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_speciation-r1-species.py program to  Plot abundance of chemical species for selected cation as a function of rho or T')
            print("plot_speciation-r1-species.py -v <variable (rho,T)> -m <mineralfile with elements> -f <stat-concentrate_r1_abso_filename>  -g <stat-concentrate_r1_abso_filename2>  -d <maximum density to plot in g/cm3> -l <letter for article subplot, default = ''>")
            print('requires the file containing elements and number (in order to compute the densities)')
            print('')
            sys.exit()
        elif opt in ("-f", "--file1"):
            statfile1 = str(arg)
        elif opt in ("-g", "--gfile2"):
            statfile2 = str(arg)
        elif opt in ("-v", "--variable"):
            variable = str(arg)
        elif opt in ("-m","--mineralfile"):
            mineralfile = str(arg)
        elif opt in ("-d","--density_max"):
            max_den = float(arg)
        elif opt in ("-l","--letter"):
            letter = str(arg)
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
        allstatfiles = [statfile1,statfile2]
        with open(statfile1,'r') as f:
            line = f.readline()
            file1=line.split()[1]
        with open(statfile2,'r') as f2:
            line = f2.readline()
            file2=line.split()[1]
        fig, ax1, ax2 = creation_plot2(variable, file1, file2, MN,plot_parameters, max_den,letter)
        figurename = statfile1.split('/')[-1].split('.dat')[0]+'+'+statfile2.split('/')[-1].split('.dat')[0]+'_'+variable+'-species'
    else:
        allstatfiles=[statfile1]
        fig, ax1 = creation_plot(variable,plot_parameters,max_den,letter)
        figurename = statfile1.split('.dat')[0]+'_'+variable  + '-species'
    
    for ii in range(len(allstatfiles)):
        statfile = allstatfiles[ii]
        #selection of the plot
        if ii == 1:
            ax = ax2
        else:
            ax = ax1
            
        #initialisation
        xdata = {'T':[],'rho':[]} #dictionnary containing the x data
        data = [] #list containing y data
        melt = [] #list containing y data for species with 100 atoms or more
        selected_files = [] #list containing the files we use for the plot (depends on the density limit)
        species1 = []
        list_colored_species = []
        list_colored_red = []
        list_colored_blue = []
        list_colored_purple = []
        list_colored_green = []
        list_grey_species = []
        list_black_species = ["melt-like"]
        colors_species = {}
        
        #***** Count of the number of clusters for the selected atom type (needed for the automatic color change)
        #first store every cluster type into the dictionary
        if re.search("L208",statfile):
            #print("we group all the species of more than nmeltatoms atoms under the denomination 'melt-like'")
            with open(statfile,'r') as f:
                line = f.readline()
                while True:
                    line = f.readline()
                    if not line: break
                    else:
                        entry=line.split('\n')[0].split('\t')
                        if len(entry) > 1:
                            natoms = count_natom(entry[0])
                            if natoms < nmeltatoms:
                                species1.append(entry[0])
                                colors_species[entry[0]] = []
                species1.append("melt-like")
                colors_species["melt-like"] = []
        else:
            #print("we use the standard script")
            with open(statfile,'r') as f:
                line = f.readline()
                while True:
                    line = f.readline()
                    if not line: break
                    else:
                        entry=line.split('\n')[0].split('\t')
                        if len(entry) > 1:
                            species1.append(entry[0])
                            colors_species[entry[0]] = []
        #*****************************
        #*******************
        #*******
        #Read and compute percentages
        #in case of speciation r1 L208, we group all the species of more than 100 atoms under the denomination "melt-like"
        with open(statfile,'r') as f:
            line = f.readline()
            entry=line.split()
            files = entry[1:]
            #creation of the x data list
            for file in files:
                temperature, acell = split_name(file)
                density = MN/(Na*float(acell)**3*10**(-24))
                if density <= max_den:    
                    xdata['rho'].append(density )    #calculation density
                    xdata['T'].append(int(temperature))   
                    selected_files.append(file)
            if re.search("L208",statfile):
                print("we group all the species of more than", nmeltatoms,"atoms under the denomination 'melt-like'")
                #creation of the y data matrix 
                while True:
                    line = f.readline()
                    if not line: break
                    else:
                        entry=line.split('\n')[0].split('\t')[:-1]
                        linedata = [] #we initialize the list of data per line
                        if len(entry) > 1:   
                            for i in range(1,len(entry)):
                                if files[i-1] in selected_files: #if the column which corresponding file is in the selected files, then we use the data
                                    if entry[i] == '':
                                        entry[i] = 0
                                    else:
                                        entry[i] = float(entry[i])
                                    linedata.append(entry[i]) #we create a list of the data for this line (cluster)
                            natoms = count_natom(entry[0])
                            if natoms < nmeltatoms:
                                data.append(linedata[:]) #we add this list to the big list of data in order to have nested lists
                            else:
                                melt.append(linedata[:])
                print('for files', selected_files)
                print('all species are', species1, len(species1))
                print('and density', xdata['rho'])
                #calculation of the percentages per density
                totalsmelt = np.ndarray.tolist(np.sum(melt,0)) #we sum all the melt species together
                data.append(totalsmelt[:]) #and we add this total to the data matrix
                print('initial data of len',np.size(data,0),np.size(data,1))
                totals = np.sum(data,0)
                #print('totals', totals)
                for i in range(np.size(data,0)):
                    for j in range(np.size(data,1)):
                        if  totals[j] == 0:
                            data[i][j] = 0
                        else:
                            data[i][j] = data[i][j]/totals[j]                    
                #print('calculated percentages',data)
                newtotals= np.sum(data,0)
                print('total should be 1',newtotals)
            else:
                print("we use the standard script")
                #creation of the y data matrix 
                while True:
                    line = f.readline()
                    if not line: break
                    else:
                        entry=line.split('\n')[0].split('\t')[:-1]
                        linedata = [] #we initialize the list of data per line
                        if len(entry) > 1:   
                            for i in range(1,len(entry)):
                                if files[i-1] in selected_files: #if the column which corresponding file is in the selected files, then we use the data
                                    if entry[i] == '':
                                        entry[i] = 0
                                    else:
                                        entry[i] = float(entry[i])
                                    linedata.append(entry[i]) #we create a list of the data for this line (cluster)
                            data.append(linedata[:]) #we add this list to the big list of data in order to have nested lists
                print('for files', selected_files)
                print('all species are', species1)
                #print('initial data',data)
                #print('and density', xdata['rho'])
                #calculation of the percentages per density
                totals = np.sum(data,0)
                #print('totals', totals)
                for i in range(np.size(data,0)):
                    for j in range(np.size(data,1)):
                        if  totals[j] == 0:
                            data[i][j] = 0
                        else:
                            data[i][j] = data[i][j]/totals[j]                    
                #print('calculated percentages',data)
                newtotals= np.sum(data,0)
                print('total should be 1',newtotals)
        #*****************************
        #*******************
        #*******
    #    #now the dictionary has every possible cluster, we attribute the colors to each cluster
    #    list_species = []
    #    color = iter(plt.cm.jet(np.linspace(0,1,len(colors_species)))) #Creation of the color list
    #    for key in natsort.natsorted(colors_species):
    #        c = next(color)
    #        colors_species[key] = c
    #        list_species.append(key)
    #    print(list_species, len(list_species))
        #*** check the abundance to put in grey the species less abundant than 1% and put a label only for species more abundante than 5% at least one time
        for i in range(np.size(data,axis=0)): #loop on each cluster
            if max(data[i]) < 0.01:
                colors_species[species1[i]]='0.5'
                list_grey_species.append(species1[i])
            else:
                #this is for basic automatic color change
                list_colored_species.append(species1[i])
                #this is for custom color change relative to the feldspars
                if (species1[i][0] == 'N' or species1[i][0] == 'K' or species1[i][0] == 'C') and re.search('Si_1O',species1[i]):
                    list_colored_blue.append(species1[i])
                elif re.match('Si_1O',species1[i]):
                    list_colored_red.append(species1[i])
                elif (species1[i][0] == 'N' or species1[i][0] == 'K' or species1[i][0] == 'C') and re.search('O',species1[i]):
                    list_colored_purple.append(species1[i])
                elif species1[i][0:5] == 'Al_1O':
                    list_colored_green.append(species1[i])
                #We add the label
                if max(data[i]) > 0.05:
                    x_max = xdata[variable][data[i].index(max(data[i]))]
                    y_max = max(data[i])
                    label = species1[i]
                    #formatage of label
                    label = format1label(label)
                    ax.text(x_max, y_max, label, color='0.35')            
        print('list red', list_colored_red)
        print('list blue', list_colored_blue)
        print('list purple', list_colored_purple)
        print('list green', list_colored_green)
        #**** we attribute the colors to each cluster with abundance > 0.01
        #this is for basic automatic color change
        #color = iter(plt.cm.jet(np.linspace(0,1,len(list_colored_species)))) #Creation of the color list
        #for key in natsort.natsorted(list_colored_species):
        #    c = next(color)
        #    colors_species[key] = c
        #and this is for custom color change relative to the feldspars
        dict_r = {'red':   ((0.0, 1.0, 1.0),  # <- at 0.0 (first value), the red component is 1 (second value) (the third value may be alpha)
                        (1.0, 0.3, 1.0)),     # <- at 1.0, the red component is 0.3
             'green': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
             'blue':  ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0))
             }
        cmreds = LinearSegmentedColormap('reds', dict_r)
        color_r = iter(cmreds(np.linspace(0,1,len(list_colored_red)))) #Creation of the color list for reds
        dict_b = {'red':   ((0.0, 0.0, 0.0),  
                        (1.0, 0.4, 1.0)),     
             'green': ((0.0, 0.0, 0.0),
                       (1.0, 1.0, 1.0)),
             'blue':  ((0.0, 1.0, 1.0),  # <- at 0.0 (first value), the blue component is 1 (second value) (the third value may be alpha)
                       (1.0, 1.0, 1.0))  # <- at 1.0, the blue component is 0.3
             }
        cmblues = LinearSegmentedColormap('blues', dict_b)
        color_b = iter(cmblues(np.linspace(0,1,len(list_colored_blue)))) #Creation of the color list for blues
        dict_p = {'red':   ((0.0, 0.6, 1.0),  # <- at 0.0 (first value), the red component is 0.6 (second value) (the third value may be alpha)
                        (1.0, 0.3, 1.0)),     # <- at 1.0, the red component is 0.3
             'green': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
             'blue':  ((0.0, 0.6, 1.0),  # <- at 0.0 (first value), the blue component is 0.6 (second value) (the third value may be alpha)
                       (1.0, 0.3, 1.0))  # <- at 1.0, the blue component is 0.3
             }
        cmpurples = LinearSegmentedColormap('purples', dict_p)
        color_p = iter(cmpurples(np.linspace(0,1,len(list_colored_purple)))) #Creation of the color list for purples
        dict_g = {'red':   ((0.0, 0.0, 0.0),  
                        (1.0, 0.0, 0.0)),     
             'green': ((0.0, 1.0, 1.0),   # <- at 0.0 (first value), the green component is 0.7 (second value) (the third value may be alpha)
                       (1.0, 0.3, 1.0)),  # <- at 1.0, the green component is 0.3
             'blue':  ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0))
             }
        cmgreens = LinearSegmentedColormap('greens', dict_g)
        color_g = iter(cmgreens(np.linspace(0,1,len(list_colored_green)))) #Creation of the color list for greens
        for key in natsort.natsorted(list_colored_species):
            if key == 'O_1':
                colors_species[key] = 'gold'
            elif key == 'O_2':
                colors_species[key] = 'darkorange'
            elif key == 'O_3':
                colors_species[key] = 'tomato'
            elif key == 'Na_1' or key == 'K_1' or key == 'Ca_1' :
                colors_species[key] = (0.5,0,1.0)
            elif key == "melt-like":
                colors_species[key] = '0'
            else:
                colors_species[key] = '0.5'
        colors_species = color_change(color_r, list_colored_red, colors_species)
        colors_species = color_change(color_b, list_colored_blue, colors_species)
        colors_species = color_change(color_p, list_colored_purple, colors_species)
        colors_species = color_change(color_g, list_colored_green, colors_species)
        #*****************************
        #*******************
        #*******
        #***plot for each cluster and write file with percentages
        newfilename = statfile.split('.dat')[0]+'_'+variable  + '-species' + '.txt'
        nf = open(newfilename,'w')
        if variable == 'rho':
            nf.write('Density(g/cm3)\t' + '\t'.join(str(round(density,2)) for density in xdata[variable]) + '\n')
        else:
            nf.write('Temperature(K)\t' + '\t'.join(str(temp) for temp in xdata[variable]) + '\n')
        #plot stacked area
        #    plt.stackplot(xdata[variable],data)
        #plot lines
        if re.search("L208",statfile):
            #print("we group all the species of more than nmeltatoms atoms under the denomination 'melt-like'")
            i=-1
            with open(statfile,'r') as f:
                line = f.readline()
                while True:
                    line = f.readline()
                    if not line: break
                    else:
                        entry=line.split('\n')[0].split('\t')[:-1]
                        if len(entry) > 1:
                            natoms = count_natom(entry[0])
                            if natoms < nmeltatoms:
                                i+=1
                                x, y = zip(*sorted(zip( xdata[variable], data[i])))
                                if entry[0] in list_colored_species: #we plot only the colored lines. delete this line to plot also the grey species
                                    line, = ax.plot(x,y, '.--', color = colors_species[entry[0]], 
                                                    markersize = plot_parameters["size_markers"],
                                                    linewidth = plot_parameters["size_lines"])
                                nf.write(format1label(entry[0]) + '\t' + '\t'.join(str(round(perc,4)) for perc in data[i]) + '\n')
                #we plot the last line, for the melt-like
                i+=1
                x, y = zip(*sorted(zip( xdata[variable], data[i])))
                line, = ax.plot(x,y, '.--', color = colors_species['melt-like'], 
                                markersize = plot_parameters["size_markers"], 
                                linewidth = plot_parameters["size_lines"])
                nf.write(format1label('melt-like') + '\t' + '\t'.join(str(round(perc,4)) for perc in data[i]) + '\n')
        else:
            i=-1
            with open(statfile,'r') as f:
                line = f.readline()
                while True:
                    line = f.readline()
                    if not line: break
                    else:
                        entry=line.split('\n')[0].split('\t')[:-1]
                        if len(entry) > 1:
                            i+=1
                            x, y = zip(*sorted(zip( xdata[variable], data[i])))
                            if entry[0] in list_colored_species: #we plot only the colored lines. delete this line to plot also the grey species
                                line, = ax.plot(x,y, '.--', color = colors_species[entry[0]], 
                                                markersize = plot_parameters["size_markers"], 
                                                linewidth = plot_parameters["size_lines"])
                            nf.write(format1label(entry[0]) + '\t' + '\t'.join(str(round(perc,4)) for perc in data[i]) + '\n')
        #                    if variable == 'T':
        #                        label_line(ax, line, entry[0], halign='right')   
        #                    elif list_species.index(entry[0]) < int(len(list_species)/2):
        #                        label_line(ax, line, entry[0], halign='left')   
        #                    else:
        #                        label_line(ax, line, entry[0], halign='right')    
    #    ax.hlines(y=0.05, xmin = 1, xmax=2.1, color = 'r')
    #    ax.hlines(y=0.01, xmin = 1, xmax=2.1, color = 'm')
        #*****************************
        #*******************
        #*******
        #legend           
        #custom_lines_grey = [Line2D([0],[0],color = colors_species[key], ls = '--', marker = '.', markersize = plot_parameters["size_markers"], linewidth = plot_parameters["size_lines"]) for key in natsort.natsorted(list_grey_species)]
        list_grey_species = format_label(natsort.natsorted(list_grey_species))       
        #legend_grey = plt.legend([line for line in custom_lines_grey],[label for label in list_grey_species],title = '$\\bf{Clusters <1\%}$', bbox_to_anchor=(1.25, 1), loc='upper left', fontsize = plot_parameters["size_fonts"],  borderaxespad=0.)
        #plt.gca().add_artist(legend_grey)
        #plt.setp(legend_grey.get_title(),fontsize= plot_parameters["size_fonts"])
        
        custom_lines = [Line2D([0],[0],color = colors_species[key], ls = '--', marker = '.', 
                               markersize = plot_parameters["size_markers"], 
                               linewidth = plot_parameters["size_lines"]) for key in natsort.natsorted(list_colored_species)]
        list_colored_species = format_label(natsort.natsorted(list_colored_species))
        if statfile2 != '':
            legend = ax.legend([line for line in custom_lines],[label for label in list_colored_species], 
                               bbox_to_anchor=(0.5, 0.99), loc='upper center', 
                               fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0.,ncol = 3)
        #if len(list_colored_species) > 16:
    #        legend = plt.legend([line for line in custom_lines],[label for label in list_colored_species],title = '$\\bf{Species}$', bbox_to_anchor=(1.01, 1), loc='upper left', fontsize = plot_parameters["size_fonts"],  borderaxespad=0., ncol = 2)
         #   legend = plt.legend([line for line in custom_lines],[label for label in list_colored_species],title = '$\\bf{Species}$', bbox_to_anchor=(0.5, 0.99), loc='upper center', fontsize = plot_parameters["size_fonts"],  borderaxespad=0.,ncol = 8)
        else:
    #        legend = plt.legend([line for line in custom_lines],[label for label in list_colored_species],title = '$\\bf{Species}$', bbox_to_anchor=(1.01, 1), loc='upper left', fontsize = plot_parameters["size_fonts"],  borderaxespad=0.)
            legend = ax.legend([line for line in custom_lines],[label for label in list_colored_species],
                               title = '$\\bf{Species}$', bbox_to_anchor=(0.5, 0.99), 
                               loc='upper center', fontsize = plot_parameters["size_fonts"],  borderaxespad=0.,ncol = 8)
            plt.setp(legend.get_title(),fontsize= plot_parameters["size_fonts"])
        
        nf.write('\n')
        nf.write('Species<1%\t'+ '\t'.join(str(species) for species in list_grey_species) + '\n')
        print(newfilename, 'is created')

    
    
    figurename = figurename +'.pdf'
    plt.savefig(figurename, bbox_inches='tight', dpi=300)
    print(figurename, 'is created')
    
    
#    plt.show() # à mettre après savefig sinon on ne sauvegarde rien !!  

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



