#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the spinodal diagrams      ****
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
        major_xticks = np.arange(0.5, 2.6, 0.5) 
        minor_xticks = np.arange(0.5, 2.1, 0.1)    
        (xxmin,xxmax) = (0.5,1.8)
        xlabel = r'Density (g.cm$^{-3}$)'
    elif  xvariable =='T':
        major_xticks = np.arange(1000, 8000, 1000) 
        minor_xticks = np.arange(1000, 8000, 500)
        (xxmin, xxmax) = (1500, 7500)
        xlabel = r'Temperature (K)'
    else:
        major_xticks = np.arange(0,1,1)
        minor_xticks = np.arange(0,1, 0.5)    
        (xxmin,xxmax) = (0,1)
        xlabel = 'Not defined yet'
        
    if yvariable == 'T':
        major_yticks = np.arange(1000, 8000, 1000) 
        minor_yticks = np.arange(1000, 8000, 500)
        (yymin, yymax) = (1500, 7500)
        ylabel = r'Temperature (K)'
    elif yvariable == 'P':
        major_yticks = np.arange(-3, 1, 0.5) 
        minor_yticks = np.arange(-3, 1, 0.1)
        (yymin, yymax) = (-2.5,0.6)
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
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),"size_markers" : 6,"size_lines" : 1,"shift_labelpad" : 10}
    letter = ''
    try:
        options,arg = getopt.getopt(argv,"hf:t:l:",["filename1",'type', 'letter'])
    except getopt.GetoptError:
        print("plot_Hugoniot_all.py -f <mineralname> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_spinodal_all.py program to plot the spinodal ')
            print("plot_spinodal_all.py -f <mineralname> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> ")
            print('')
            sys.exit()
        if opt in ('-f','--filename1'):
            mineralname = str(arg)
        elif opt in ('-t','--type'):
            type_plot = str(arg)
            xvariable = type_plot.split('-')[1]
            yvariable = type_plot.split('-')[0]
        elif opt in ('-l','--letter'):
            letter = str(arg)
    files = sorted(glob.glob('spinodal_*_'+mineralname+'_*.txt')) #I list only the files we want
    figurename = 'spinodal_'+mineralname
    #creation plot
    fig, ax = creation_plot(plot_parameters, xvariable, yvariable) 
    #************ creation of arrays and plot at the same time
    for file in files:
        #**initialisation of data dicitonnaries for each column (marker type)  
        data = {'T' : [], 'rho':[],'P':[]}
        #********* fill the dictionnaries with data
        with open(file,'r') as f:
            f.readline()
            while True:
                line = f.readline()
                if not line: break
                else:    
                    entry=line.split('\n')[0].split('\t')
                    data['T'].append(float(entry[0].strip('T'))*1000)
                    data['rho'].append(float(entry[1]))
                    data['P'].append(float(entry[2]))
        print(data['T'])
        print(data['rho'])
        #******* change of markers and colors style
        filename = file.split('_')
        if filename[1] == 'constrained':
            bordercolor = '#ffc900'
        else:
            bordercolor = 'k'
        if filename[-2] == 'soft':
            style = 'o:'
            fillcolor = '0.90'
        else:
            style = 'o--'
            fillcolor = bordercolor
        #******* plot the spinodal lines     
        ax.plot(data[xvariable],data[yvariable], style, markersize = plot_parameters["size_markers"]+3, markeredgecolor = bordercolor, markeredgewidth = 0.5, markerfacecolor = fillcolor, linewidth = plot_parameters["size_lines"], color = bordercolor)
    #text on plot
    if letter != '':
        ax.text(0.02,0.94, letter , transform=ax.transAxes, horizontalalignment = 'left', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
    
    figurename = figurename+'_'+yvariable+'-'+xvariable + '.png'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












