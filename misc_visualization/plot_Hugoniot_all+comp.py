#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the Hugoniot diagrams      ****
               *************  ARTICLE VERSION  *************

"""


#     ********* Importation of the packages and modules used here *********
import glob
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import natsort
from matplotlib.lines import Line2D
from matplotlib.colors import to_hex
from matplotlib.ticker import FormatStrFormatter
import re


def creation_plot(plot_parameters, xvariable, yvariable):
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
    
    #tick definition
    if xvariable =='rho':
        major_xticks = np.arange(2.25, 7, 1) 
        minor_xticks = np.arange(2.25, 7, 0.5)    
        (xxmin,xxmax) = (2.25,6.5)
#        (xxmin,xxmax) = (3.25,5.25)
        xlabel = r'Density (g.cm$^{-3}$)'
    elif  xvariable =='T':
        major_xticks = np.arange(0, 25000, 5000) 
        minor_xticks = np.arange(0, 21000, 1000)
        (xxmin, xxmax) = (0, 20500)
#        (xxmin, xxmax) = (0, 6000)
        xlabel = r'Temperature (K)'
    else:
        major_xticks = np.arange(0,1,1)
        minor_xticks = np.arange(0,1, 0.5)    
        (xxmin,xxmax) = (0,1)
        xlabel = 'Not defined yet'
        
    if yvariable == 'T':
        major_yticks = np.arange(0, 25000, 5000) 
        minor_yticks = np.arange(0, 21000, 1000)
        (yymin, yymax) = (0, 20500)
#        (yymin, yymax) = (0, 3500)
        ylabel = r'Temperature (K)'
    elif yvariable == 'P':
        major_yticks = np.arange(0, 500, 50) 
        minor_yticks = np.arange(0, 500, 25)
        (yymin, yymax) = (0,450)
#        (yymin, yymax) = (0,130)
        ylabel = r'Pressure (GPa)'
    else:
        major_yticks = np.arange(0, 1, 1) 
        minor_yticks = np.arange(0, 1, 0.5) 
        (yymin, yymax) = (0, 1)
        ylabel = 'Not defined yet'
    
    
    #apply setup
    ax.set_xticks(major_xticks)
    ax.set_xticks(minor_xticks, minor=True)
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)                                            
    ax.set_yticks(major_yticks)
    ax.set_yticks(minor_yticks, minor=True)
    ax.yaxis.set_ticks_position('both') 
    #ax.grid(True, axis = 'x')
    ax.set_xlim(xxmin,xxmax)
    ax.set_ylim(yymin,yymax)
#    plt.autoscale()
    
    ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    ax.set_ylabel(ylabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    return fig, ax


def add_insert(fig,plot_parameters, xvariable, yvariable):
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]

    # tick and insert definition
    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    if yvariable =='P' and xvariable =='rho':
        left1, bottom1, width1, height1 = [0.2, 0.44, 0.45, 0.42]     
        major_xticks = np.arange(2.25, 7, 1) 
        minor_xticks = np.arange(2.25, 7, 0.5)    
        (xxmin,xxmax) = (2.8,5.5)
        xlabel = r'Density (g.cm$^{-3}$)'     
        major_yticks = np.arange(0, 200, 50) 
        minor_yticks = np.arange(0, 200, 25)
        (yymin, yymax) = (0,150)
        ylabel = r'Pressure (GPa)'
    elif yvariable =='T' and xvariable =='rho':
        left1, bottom1, width1, height1 = [0.215, 0.54, 0.35, 0.32] 
        major_xticks = np.arange(2.25, 7, 1) 
        minor_xticks = np.arange(2.25, 7, 0.5)    
        (xxmin,xxmax) = (2.8,5.5)
        xlabel = r'Density (g.cm$^{-3}$)'
        major_yticks = np.arange(0, 5000, 1000) 
        minor_yticks = np.arange(0, 5000, 500)
        (yymin, yymax) = (0, 3500)
        ylabel = r'Temperature (K)'
    elif yvariable =='P' and xvariable =='T':
        left1, bottom1, width1, height1 = [0.2, 0.54, 0.35, 0.32] 
        major_xticks = np.arange(0, 6000, 2000) 
        minor_xticks = np.arange(0, 6000, 500)
        (xxmin, xxmax) = (0, 6000)
        
        xlabel = r'Temperature (K)'
        major_yticks = np.arange(0, 200, 50) 
        minor_yticks = np.arange(0, 200, 25)
        (yymin, yymax) = (0,150)
        ylabel = r'Pressure (GPa)'
    else:
        left1, bottom1, width1, height1 = [0.215, 0.54, 0.35, 0.32]     
        major_xticks = np.arange(0,1,1)
        minor_xticks = np.arange(0,1, 0.5)    
        (xxmin,xxmax) = (0,1)
        xlabel = 'Not defined yet'     
        major_yticks = np.arange(0, 1, 1) 
        minor_yticks = np.arange(0, 1, 0.5) 
        (yymin, yymax) = (0, 1)
        ylabel = 'Not defined yet'
    
    #apply setup
    ax = fig.add_axes([left1, bottom1, width1, height1])    
    ax.set_xticks(major_xticks)
    ax.set_xticks(minor_xticks, minor=True)
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(which = 'both', labelsize = size_font_ticks-3, width = size_lines/2)
    ax.tick_params(axis = 'x', rotation = 45 )
    ax.set_yticks(major_yticks)
    ax.set_yticks(minor_yticks, minor=True)
    ax.yaxis.set_ticks_position('both') 
    ax.set_xlim(xxmin,xxmax)
    ax.set_ylim(yymin,yymax)
    #ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  #set one decimal to ticks label      
        
    #ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    #ax.set_ylabel(ylabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    return ax



def format_1label(title):
    """ formatage of mineralname"""
    i=0
    while True:
        if i <= len(title)-1:
            if re.match('[0-9]',title[i]):
                num=0
                try:
                    while re.match('[0-9]',title[i+num]):
                        num +=1
                except IndexError: #index error when we arrive at the end of the cluster name
                    pass
                    #print('end of the cluster')
                title = title[:i]+'$_{'+title[i:i+num]+'}$' + title[i+num:]
                i = i+5
            i = i+1
        else:break
    return title


def main(argv):
    """     ********* Main program *********     """
    column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,'Cvm':25,'stdev_Cvm':26}
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    markertype = ['o','s','X','P','d','p'] #one marker per source
    plot_params = {}
    #other dictionnaries and parameters for the figure
    mineralname = 'all'
    minerals = []
    letter = ''
    filename2 = ''
    comparison = 0
    try:
        options,arg = getopt.getopt(argv,"hf:g:t:l:",["fmineralname","gcomparisonfilename",'type', 'letter'])
    except getopt.GetoptError:
        print("plot_Hugoniot_all.py -f <mineralname>(default = all files of every compound) -g <comparison-file> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_Hugoniot_all.py program to plot the Hugoniot ')
            print("plot_Hugoniot_all.py -f <mineralname>(default = all files of every compound) -g <comparison-file> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> -c < =1 for comparison with thermo data points, default = 0>")
            print("plot_Hugoniot_all.py requires to be lauched from the folder containing every Hugoniot .txt file created by analyze_Hugoniot.py")
            print('')
            sys.exit()
        if opt in ('-f','--fmineralname'):
            mineralname = str(arg)
        elif opt in ('-g','--gfilename2'):
            filename2 = str(arg)
        elif opt in ('-t','--type'):
            type_plot = str(arg)
            xvariable = type_plot.split('-')[1]
            yvariable = type_plot.split('-')[0]
        elif opt in ('-l','--letter'):
            letter = str(arg)
    if mineralname == 'all':
        files = sorted(glob.glob('Hugoniot_*_ground-state_*.txt')) #I list every  files
        removefiles = glob.glob('Hugoniot_*_ground-state_3000-0.txt')
        figurename = 'Hugoniot_all+comp'
        for file in sorted(files):
            mineral = file.split('_')[1]
            if mineral in minerals:
                continue
            else:
                minerals.append(mineral)
    else:
        files = sorted(glob.glob('Hugoniot_'+mineralname+'_ground-state_*.txt')) #I list only the files we want
        removefiles = ['Hugoniot_'+mineralname+'_ground-state_3000-0.txt']
        figurename = 'Hugoniot_'+mineralname+'+comp'
        minerals = [mineralname]
    for file in removefiles:
        files.remove(file)
    #**** Creation of the plot 
    fig, ax1 = creation_plot(plot_parameters, xvariable, yvariable)
    #**** Create insert
    if filename2 != '':
        pass
        #ax11 = add_insert(fig,plot_parameters, xvariable, yvariable)
    #creation of color changes: 1 color per ground state
    colors = {'2260-3000':'#ff0000','2500-0':'#b1dbff','2585-1932':'#ff8a00','2600-0':'#0080ed','2700-0':'#004f92'}
    #definition of line and marker style
    for mineral in minerals:
        if mineral == 'CaAl2Si2O8':
            plot_params[mineral] = 'o-'  
        elif mineral == 'KAlSi3O8':
            plot_params[mineral] = 'o--'   
        else:
            plot_params[mineral] = 'o:'  
    #************ creation of arrays and plot at the same time our data
    for file in files:
        #**initialisation of data dicitonnaries for each column (marker type)  
        Hugoniot = {'rho':[],'T':[],'P':[]}
        thermo = {'rho':[],'T':[],'P':[]}
        mineral = file.split('_')[1]
        gs =  file.split('_')[-1].split('.txt')[0]
        #********* fill the dictionnaries with data
        with open(file,'r') as f:
            f.readline()
            f.readline()
            while True:
                line = f.readline()
                if not line: break
                else:    
                    entry=line.split('\t')
                    try:
                        Hugoniot['P'].append(float(entry[2])) #pressure in GPa
                        Hugoniot['rho'].append(float(entry[1])/1000) #density in g/cm3
                        Hugoniot['T'].append(float(entry[0])) #temperature in K
                    except ValueError:
                        print("Missing data in your file", file)
                        sys.exit()
        #******* change of markers and colors style and selection of thermofile
        if mineral == 'NaAlSi3O8':
            thermofile = 'thermo_NaAlSi3O8_hard_all.txt'    
        else:
            thermofile = 'thermo_'+mineral+'_all.txt'
        #******* extraction of thermo data for comparison
        if comparison == 1:
            print("I compare with thermo data")
            try:
                with open(thermofile,'r') as tf:
                    [tf.readline() for i in range(3)]
                    while True:
                        line = tf.readline()
                        if not line: break
                        else:
                            entry=line.split()
                            thermo['rho'].append(float(entry[column_number['rho']]))
                            thermo['T'].append(float(entry[column_number['T']]))
                            thermo['P'].append(float(entry[column_number['P']]))
                ax1.plot(thermo[xvariable],thermo[yvariable],'k+',label='data')
            except FileNotFoundError:
                print("file",thermofile,"does not exist")
        #******* plot the hugoniot lines     
        ax1.plot(Hugoniot[xvariable],Hugoniot[yvariable], plot_params[mineral], marker = None, linewidth = plot_parameters["size_lines"], color =colors[gs])
        if filename2 != '':
            pass
            #ax11.plot(Hugoniot[xvariable],Hugoniot[yvariable], plot_params[mineral], marker = None, linewidth = plot_parameters["size_lines"], color =colors[gs])
    
    #************ creation of arrays and plot at the same time the litterature data
    compdata = {}
    sources = []
    try:
        #**** Extract data
        with open(filename2,'r') as f:
            f.readline()
            while True:
                line = f.readline()
                if not line:break
                else:   
                    entry = line.split('\n')[0].split('\t')
                    try:
                        compdata[entry[0]]['mineral'].append(entry[2])
                        if entry[4] == '0':
                            compdata[entry[0]]['rho'].append(float('nan'))
                        else:
                            compdata[entry[0]]['rho'].append(float(entry[4]))
                        compdata[entry[0]]['P'].append(float(entry[5])/10)
                        if entry[6] == '0':
                            compdata[entry[0]]['T'].append(float('nan'))
                        else:
                            compdata[entry[0]]['T'].append(float(entry[6]))
                    except KeyError:
                        compdata[entry[0]] = {}
                        compdata[entry[0]]['mineral'] = [entry[2]]
                        if entry[4] == '0':
                            compdata[entry[0]]['rho'] = [float('nan')]
                        else:
                            compdata[entry[0]]['rho'] = [float(entry[4])]
                        compdata[entry[0]]['P'] = [float(entry[5])/10]
                        if entry[6] == '0':
                            compdata[entry[0]]['T'] = [float('nan')]
                        else:
                            compdata[entry[0]]['T'] = [float(entry[6])]
                        sources.append(entry[0])
        #**** attribution markers and colors to each source
        markersources = {}
        colorsources = {}
        for ii in range(len(sources)):
            markersources[sources[ii]]=markertype[ii]
            if sources[ii][0] == 'M':
                colorsources[sources[ii]] = '#44aa00' #green
            elif sources[ii][0] == 'B':
                colorsources[sources[ii]] = '#7137c8' #purple
            elif sources[ii][1] == 'h':
                colorsources[sources[ii]] = '#000000' 
            elif sources[ii][1] == 's':
                colorsources[sources[ii]] = '#ff00cc' #pink
        #**** plot
        for source in sources:
            for ii in range(len(compdata[source][xvariable])):
                if compdata[source]['mineral'][ii] == 'NaAlSi3O8':
                    fillcolor = 'w'
                elif compdata[source]['mineral'][ii] == 'CaAl2Si2O8':
                    fillcolor = colorsources[source]+'ff'
                elif compdata[source]['mineral'][ii] == 'KAlSi3O8':
                    fillcolor = colorsources[source]+'7f'
                #ax11.scatter(compdata[source][xvariable][ii],compdata[source][yvariable][ii], marker = markersources[source], edgecolor = colorsources[source], facecolor = fillcolor, s = plot_parameters["size_markers"]*5)
                ax1.scatter(compdata[source][xvariable][ii],compdata[source][yvariable][ii], marker = markersources[source], edgecolor = colorsources[source], facecolor = fillcolor, s = plot_parameters["size_markers"]*5)
    except FileNotFoundError:
        print("Litterature file '",filename2,"' not found")
        
    #****** Create legend from custom artist/label lists  
    custom_markers = [Line2D([0], [0], marker = markersources[key], linestyle = None, markeredgecolor = colorsources[key], markerfacecolor = colorsources[key],  markersize = plot_parameters["size_markers"], markeredgewidth = 0.5, linewidth = 0) for key in sources]
    legend3 = ax1.legend([marker for marker in custom_markers],[label for label in sources], bbox_to_anchor=(1.05, 0.3), loc=2, fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0.)

    #definition of line and marker style
    fillcolor = {}
    for mineral in minerals:
        if mineral == 'CaAl2Si2O8':
            plot_params[mineral] = 'o-'  
            fillcolor[mineral] = '#000000ff'
        elif mineral == 'KAlSi3O8':
            plot_params[mineral] = 'o--'   
            fillcolor[mineral] = '#0000007f'
        else:
            plot_params[mineral] = 'o:'  
            fillcolor[mineral] = 'w'
    #format mineral names
    mineralnames = []
    for label in natsort.natsorted(plot_params):
        mineralnames.append(format_1label(label))
    custom_lines = [Line2D([0], [0], marker = plot_params[key][0], linestyle = plot_params[key][1:], markeredgecolor = 'k', markerfacecolor = fillcolor[key],  color = 'k', markersize = plot_parameters["size_markers"], markeredgewidth = 0.5, linewidth = plot_parameters["size_lines"]) for key in natsort.natsorted(plot_params)]
    legend2 = ax1.legend([line for line in custom_lines],[label for label in mineralnames], bbox_to_anchor=(1.05, 0.55), loc=2, fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0.)

    orderedkeys = ['2500-0','2600-0','2700-0','2585-1932','2260-3000']
    custom_patch = [mpatches.Patch(color=colors[key]) for key in orderedkeys]
    legend = ax1.legend([col for col in custom_patch],[str(round(float(label.split('-')[0])/1000,2)) + '  -  ' + label.split('-')[1] for label in orderedkeys],title = '$\\rho_0$ (g$.cm^{-3}$) - T$_0$ (K)', bbox_to_anchor=(1.05, 1), loc=2, fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0.)
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])

    ax1.add_artist(legend2)
    
    ax1.add_artist(legend3)
     
    #text on plot
    #ax1.text(0.015,0.945, letter , transform=ax1.transAxes, horizontalalignment = 'left', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='none', edgecolor='k', pad=3.0)) 
    ax1.text(0.985,0.945, letter , transform=ax1.transAxes, horizontalalignment = 'right', verticalalignment = 'baseline', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='w', edgecolor='k', pad=3.0)) 
    
    
    #***** Save figure
    figurename = figurename+'_'+yvariable+'-'+xvariable + '.svg'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












