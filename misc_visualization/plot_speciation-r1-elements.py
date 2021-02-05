#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3


        Plot element % in the gas phase as a function of rho or T
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
    ax.set_ylabel(r'Element proportion', fontweight = 'bold', fontsize = size_fonts )

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
    ax.set_ylabel(r'Element proportion', fontweight = 'bold', fontsize = size_fonts , labelpad = shift_labelpad)
    #Adjustment of ticks
    if variable == 'rho':
        major_xticks = np.arange(0, 7, 0.1) 
        minor_xticks = np.arange(0, 7, 0.05)    
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_xlabel(r"Density (g.cm$^{-3}$)", fontweight = 'bold', fontsize = size_fonts , labelpad = shift_labelpad/2)
        ax.set_xlim(1,max_den)
    if variable == 'T':
        major_xticks = np.arange(0, 10000, 1000) 
        ax.set_xticks(major_xticks)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_xlabel(r"Temperature (K)", fontweight = 'bold', fontsize = size_fonts , labelpad = shift_labelpad/2)
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



def main(argv):
    """     ********* Main program *********     """
    #param for the plot
    colors_elem={'Al':'pink','C':'0.25','Ca':'c','H':'w','K':'m','Na':'b','O':'r','Si':'y'}  #other dictionnaries and parameters for the figure for article version
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
        print("plot_speciation-r1-elements.py  -v <variable (rho,T)> -m <mineralfile with elements> -f <stat-concentrate_r1_abso_filename1> -g <stat-concentrate_r1_abso_filename2> -d <maximum density to plot in g/cm3> -l <letter for article subplot, default = ''>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_speciation-r1-elements.py program to  Plot abundance of elements as a function of rho or T')
            print("plot_speciation-r1-elements.py -v <variable (rho,T)> -m <mineralfile with elements> -f <stat-concentrate_r1_abso_filename1> -g <stat-concentrate_r1_abso_filename2>  -d <maximum density to plot in g/cm3> -l <letter for article subplot, default = ''>")
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
        number = list(map(int,entry.split()[1:]))
    MN = 0
    for i in range(len(elements)):
        MN = MN + number[i] * cr.Elements2rest(elements[i])[3]
    #**** Calculation of congruent proportions
    congruent = {}
    for ii in range(len(elements)):
        congruent[elements[ii]] = number[ii]/sum(number)
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
        figurename = statfile1.split('/')[-1].split('.dat')[0]+'+'+statfile2.split('/')[-1].split('.dat')[0]+'_'+variable+'-elements'
    else:
        allstatfiles=[statfile1]
        fig, ax1 = creation_plot(variable,plot_parameters,max_den,letter)
        figurename = statfile1.split('.dat')[0]+'_'+variable  + '-elements'
    
    for ii in range(len(allstatfiles)):
        statfile = allstatfiles[ii]
        #selection of the plot
        if ii == 1:
            ax = ax2
        else:
            ax = ax1
        #initialisation
        atomsnumbers = {}  #create dictionnary which will contain the number of atom for each element in each species
        filescolumns = {}   #create dictionnary containing the column number for the different files
        xdata = {'T':[],'rho':[]} #dictionnary containing the x data
        lifetime = {} #dictionnary containing lifetime data for each element and file
        perc = {}       #idem for percentages
        selected_files = [] #list containing the files we use for the plot (depends on the density limit)            
        #***** Extraction of all files and cluster and count of their atoms 
        with open(statfile,'r') as f:
            line = f.readline()
            entry = line.split('\n')[0].split('\t')[1:-1]
            for file in entry:
                filescolumns[file] = entry.index(file)
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    if len(entry) > 1:
                        atomsnumbers[entry[0]] = {} #create dictionnary which will contain the number of atom for each element in each species
        for elem in elements:
            #print("for elem",elem)
            for species in atomsnumbers:
                #print("for species",species)
                m = re.search(elem, species)
                if m:                   #if there is the current element in the species
                    ii = m.end()        #then we use the index of the end of the matching pattern m
                    try:                #to try to see if there is several iteration of this element using '_'
                        if species[ii] == '_': 
                            num=''
                            try:
                                while re.match('[0-9]',species[ii+1]):
                                    num = num + species[ii+1]
                                    ii +=1
                            except IndexError: #index error when we arrive at the end of the cluster name
                                pass
                                #print('end of the cluster')
                            #print("num after '_':", num)
                            atomsnumbers[species][elem] = int(num)
                        else:           #if there is no '_' after the matching pattern it means there is only one atom of this element
                            atomsnumbers[species][elem] = 1 
                    except IndexError:  #if the matching pattern is actually at the end of the string we raise an indexerror --> there is only one instance of this element
                        atomsnumbers[species][elem] = 1
                else:
                    atomsnumbers[species][elem] = 0
        #print(atomsnumbers)
        #**** Extraction of selected files
        for file in filescolumns:
            temperature, acell = split_name(file)
            density = MN/(Na*float(acell)**3*10**(-24))
            if density <= max_den:    
                xdata['rho'].append(density )    #calculation density
                xdata['T'].append(int(temperature))   
                selected_files.append(file)
        print("selected files are",selected_files)
        #**** Initialization of perc and lifetime dictionnary for each element and file
        for elem in elements:
            perc[elem] = {}
            lifetime[elem] = {}
            for file in selected_files:
                perc[elem][file] = []
                lifetime[elem][file] = []
        #*****************************
        #*******************
        #*******
        #**** Extract the lifetimes
        with open(statfile,'r') as f:
            f.readline()
            while True:
                line = f.readline()  #for each line of the file we extract the data
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')[:-1]
                    #print("entry is",entry)
                    for elem in atomsnumbers[entry[0]]:   #we skip the analysis for element that are not present in the species
                        if atomsnumbers[entry[0]][elem] == 0:
                            continue
                        else:                           #and we store the lifetime for each element and file
                            for file in selected_files:
                                try:
                                    lifetime[elem][file].append(float(entry[filescolumns[file]+1]) * atomsnumbers[entry[0]][elem])
                                except ValueError: #if there is '' instead of a numerical value
                                    lifetime[elem][file].append(0)
        #**** Compute percentages
        #1st: sum lifetime per element
        for elem in elements:
            perc[elem] = []
            for file in selected_files:
                perc[elem].append(sum(lifetime[elem][file]))
        #2nd: sum all tot lifetime
        totlifetime = []
        for ii in range(len(selected_files)):
            sumlifetime = 0
            for elem in perc:
                sumlifetime = sumlifetime + perc[elem][ii]
            totlifetime.append( sumlifetime )
        #3rd: compute percentages
        for elem in perc:
            for ii in range(len(selected_files)):
                perc[elem][ii] = perc[elem][ii] / totlifetime[ii]
        #*****************************
        #*******************
        #*******
        #**** Plot the percentages and write file with percentages
        newfilename = statfile.split('.dat')[0]+'_'+variable  + '-elements' + '.txt'
        nf = open(newfilename,'w')
        if variable == 'rho':
            nf.write('Density(g/cm3)\t' + '\t'.join(str(round(density,2)) for density in xdata[variable]) + '\n')
        else:
            nf.write('Temperature(K)\t' + '\t'.join(str(temp) for temp in xdata[variable]) + '\n')
        ydata = []
        ycolors = []
        labels = []
        for elem in perc:
            #sort the data by the x value
            x, y = zip(*sorted(zip( xdata[variable], perc[elem])))
            #plot lines
            line, = ax.plot(x,y, '.--', color = colors_elem[elem], markersize = plot_parameters["size_markers"], 
                            linewidth = plot_parameters["size_lines"], label = elem)        
            label_line(ax, line, elem, halign='center')    
            #write data in file
            nf.write(elem + '\t' + '\t'.join(str(round(data,4)) for data in perc[elem]) + '\n')
            #for stackplot only
            ydata.append(y)
            ycolors.append(colors_elem[elem])
            labels.append(elem)
        #plot stacked area
        #plt.stackplot(x,ydata[0],ydata[1],ydata[2],ydata[3], colors = ycolors, labels = labels )
        
        #plot lines of congruent gas
        for elem in congruent:
            ax.axhline(y=congruent[elem], color = colors_elem[elem],  
                       linewidth = plot_parameters["size_lines"]/1.5, linestyle = ':' )
        
        
        #legend = plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize = plot_parameters["size_fonts"], title = '$\\bf{Elements}$', borderaxespad=0., ncol = 1)
        #plt.setp(legend.get_title(),fontsize= plot_parameters["size_fonts"])
        print(newfilename, 'is created')
        
    figurename = figurename +'.pdf'
    plt.savefig(figurename, bbox_inches='tight', dpi=300)
    print(figurename, 'is created')
#    plt.show() # à mettre après savefig sinon on ne sauvegarde rien !!  

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



