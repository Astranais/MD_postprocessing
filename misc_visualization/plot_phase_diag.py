#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the phase diagram in any projection     ****
               *************  ARTICLE VERSION  *************

"""


#     ********* Importation of the packages and modules used here *********
import glob
import sys
import getopt
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import re, os





def update_dictionnary(fillcolor):
    "update dictionnary of plot params from fillcolor"
    markercolor= {'NaAlSi3O8':fillcolor,'KAlSi3O8':fillcolor,'CaAl2Si2O8':'#ffffff'}
    return  markercolor


def format_1label(compound):
    #and formatage of label
    i=0
    while True:
        if i <= len(compound)-1:
            if re.match('[0-9]',compound[i]):
                num=0
                try:
                    while re.match('[0-9]',compound[i+num]):
                        num +=1
                except IndexError: #index error when we arrive at the end of the cluster name
                    pass
                    #print('end of the cluster')
                compound = compound[:i]+'$_{'+compound[i:i+num]+'}$' + compound[i+num:]
                i = i+5
            i = i+1
        else:break
    return compound




def creation_plot_3(plot_parameters,  major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax, yvariable):
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
       figsize = (size_figure[0]*2.8,size_figure[1]))  
    
    
    #apply setup
    for ax in [ax1,ax2,ax3]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_yticks(major_yticks)                                        
        ax.set_yticks(minor_yticks, minor=True)
        ax.yaxis.set_ticks_position('both')
    
    #for ax in [ax2,ax3]:
    #    if yvariable != 'T':
    #        plt.setp(ax.get_yticklabels()[-1], visible=False)
    for ax in [ax1,ax2]:
        plt.setp(ax.get_xticklabels()[-1], visible=False)
    
    ax1.set_xlim(xxmin,xxmax)
    ax1.set_ylim(yymin,yymax)
#    plt.autoscale()
    
    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, 
                      hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    ax.set_ylabel(ylabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*4)
    return f, ax1, ax2, ax3, ax


def creation_plot_2_vertical(plot_parameters,  major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax, yvariable):
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
    f, (ax1, ax2) = plt.subplots(2,1, sharex=True, sharey=True, 
       figsize = (size_figure[0],size_figure[1]*2))  
    
    
    #apply setup
    for ax in [ax1,ax2]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_yticks(major_yticks)                                        
        ax.set_yticks(minor_yticks, minor=True)
        ax.yaxis.set_ticks_position('both')
    
    for ax in [ax2]:
        if yvariable != 'T':
            plt.setp(ax.get_yticklabels()[-1], visible=False)
    
    ax1.set_xlim(xxmin,xxmax)
    ax1.set_ylim(yymin,yymax)
#    plt.autoscale()
    
    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    ax.set_ylabel(ylabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*4)
    return f, ax1, ax2, ax


def creation_plot_2_horizontal(plot_parameters,  major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax, yvariable):
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
    f, (ax1, ax2) = plt.subplots(1,2, sharex=True, sharey=True, 
       figsize = (size_figure[0]*2,size_figure[1]*2))  
     
    #apply setup
    for ax in [ax1,ax2]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_yticks(major_yticks)                                        
        ax.set_yticks(minor_yticks, minor=True)
        ax.yaxis.set_ticks_position('both')
        
    ax1.set_xlim(xxmin,xxmax)
    ax1.set_ylim(yymin,yymax)
#    plt.autoscale()
    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, 
                      hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    ax.set_ylabel(ylabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*4)
    return f, ax1, ax2, ax

def creation_plot(plot_parameters, major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax):
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
    fig, ax = plt.subplots(1,1, figsize = (9,5), linewidth = size_lines /2) 
    
    #apply setup
    if (xxmin,xxmax) == (0.1,275):
        ax.set_xscale('log')
    else:
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
    
    ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad-5)
    ax.set_ylabel(ylabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    return fig, ax



def creation_plot_break(plot_parameters, major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax):
    """     ********** Creation of the plot  **********    """
    print("I use creation plot_break")
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    #plot
    plt.close()
    fig = plt.figure(figsize = (size_figure[0],size_figure[1]))
    
    gs = GridSpec(1,2,width_ratios=[1,2],top = 0.95, bottom = 0.2, 
                  right = 0.95, left = 0.18, wspace = 0, hspace = 0)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])    
     
    ax1.set_xticks(np.arange(-5,1.5,1))
    ax1.set_xticks(np.arange(-5,1.5,0.2), minor=True)
    
    #apply setup
    for ax in [ax1,ax2]:
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks,
                       width = size_lines/2)
        ax.set_yticks(major_yticks)                                        
        ax.set_yticks(minor_yticks, minor=True)
        ax.yaxis.set_ticks_position('both')
        ax.set_ylim(yymin,yymax)
        
    ax1.set_xlim(-3,1)
    ax2.set_xlim(xxmin,xxmax)
    
    plt.setp(ax2.get_xticklabels()[1], visible=False) 
    ax2.tick_params(labelleft=False)
    ax2.set_xscale('log')
    
    ax = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, 
                    bottom=False, labelleft=False, left=False, 
                    labelright = False, right=False)
    ax.set_xlabel(r'Pressure (GPa)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = 0)#shift_labelpad)
    ax.set_ylabel(r'Temperature (K)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad/2)
    return fig, ax1, ax2, ax



def append_value(data,entry,i):
    """      ******** append a float to the list, replacing empty str by NaN ****** """
    if entry[i] == '':
        data.append(float('nan'))
    elif entry[i] == '\n':
        data.append(float('nan'))
    else:
        data.append(int(entry[i]))
    return data


def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    #plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (9,5)
    #,"size_markers" : 10,"size_lines" : 1,"shift_labelpad" : 10} #for beamer figures
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4)
    ,"size_markers" : 10,"size_lines" : 1,"shift_labelpad" : 10} #for article figures
    markeralpha= {'NaAlSi3O8':'ff','KAlSi3O8':'7f','CaAl2Si2O8':'ff'}
    markertype= {'NaAlSi3O8':'o','KAlSi3O8':'s','CaAl2Si2O8':'d'}
    linetype= {'NaAlSi3O8':'-','KAlSi3O8':'--','CaAl2Si2O8':':'}
    zorder = 0
    #other dictionnaries and parameters for the figure
    filename = ''
    filename2 = ''
    combined = 0
    fillcolor = '#000000'
    init_letter = ''
    letters = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    #other dictionnaries and parameters for the figure
    compounds = []
    #location of the critical point
    critical_area = {'NaAlSi3O8':{'P':(0.03,0.14),'T':(5500,6000),'rho':(0.51,0.83)},
                'KAlSi3O8':{'P':(0.02,0.11),'T':(5000,5500),'rho':(0.56,0.87)},
                'CaAl2Si2O8':{'P':(0.11,0.21),'T':(7000,7500),'rho':(0.57,0.79)}  }
    # with constrained fit for all
    critical_constrained = {'NaAlSi3O8':{'P':(0.1,0.2),'T':(6000,6500),'rho':(0.44,0.62)}, 
            'KAlSi3O8':{'P':(0.1,0.15),'T':(5500,6000),'rho':(0.36,0.64)} , 
            'CaAl2Si2O8':{'P':(0.1,0.17),'T':(7000,7500),'rho':(0.34,0.77)}  }
    try:
        options,arg = getopt.getopt(argv,"hf:g:t:l:v:c:",["filename","gfilename",'type', 'letter', "view", "combined"])
    except getopt.GetoptError:
        print("plot_phase_diag.py -f <filename> -g <filename 2 (if we want to plot 2 files)> -t <type of plot: y-x variables> -l <letter of first plot, default = ''> -v <view = zoom,classic,full> -c <=1 to combine all data on 1 plot, default = 0>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_phase_diag.py program to plot the phase diagram ')
            print("plot_phase_diag.py -f <filename (or 'all')> -g <filename 2 (if we want to plot 2 files)> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> -v <view = zoom,classic,full> -c <=1 to combine all data on 1 plot, default = 0>")
            print("plot_phase_diag.py requires to be lauched from the folder containing every phase_diag_compound.txt file")
            print('')
            sys.exit()
        if opt in ('-f','--filename'):
            filename = str(arg)
        elif opt in ('-g','--gfilename'):
            filename2 = str(arg)
        elif opt in ('-t','--type'):
            type_plot = str(arg)
            xvariable = type_plot.split('-')[1]
            yvariable = type_plot.split('-')[0]
        elif opt in ('-l','--letter'):
            init_letter = str(arg)
        elif opt in ('-v','--view'):
            view = str(arg)
        elif opt in ('-c','--combined'):
            combined = int(arg)
    
    #************ tick definition
    tickdefinitions = {'P':{'zoom':(-1,1.6,0.5,0.1),'classic':(0,300,50,10),'full':(0,550,10,5)} ,
                     'rho':{'zoom':(0,5,0.5,0.1),'classic':(0,7,0.5,0.1),'full':(0,9,1,0.5)}, 
                       'T':{'zoom':(0,9000,1000,500),'classic':(0,8000,1000,500),'full':(0,21000,1000,500)} }   #begin/end/stepmajor/stepminor
    limitaxis =  {'P':{'zoom':(-1,1),'classic':(1,275),'full':(0,25)} , 
                'rho':{'zoom':(0.3,2),'classic':(0.2,6),'full':(0.3,6.5)}, 
                  'T':{'zoom':(1500,7700),'classic':(273,7700),'full':(0,7200)} }  #begin/end
    axislabel =  {'P':r'Pressure (GPa)', 'rho':r'Density (g.cm$^{-3}$)', 'T':r'Temperature (K)' } 
    major_xticks = np.arange(tickdefinitions[xvariable][view][0], tickdefinitions[xvariable][view][1], 
                             tickdefinitions[xvariable][view][2]) 
    minor_xticks = np.arange(tickdefinitions[xvariable][view][0], tickdefinitions[xvariable][view][1], 
                             tickdefinitions[xvariable][view][3]) 
    (xxmin,xxmax) = (limitaxis[xvariable][view][0], limitaxis[xvariable][view][1])
    xlabel = axislabel[xvariable]
    major_yticks = np.arange(tickdefinitions[yvariable][view][0], tickdefinitions[yvariable][view][1], 
                             tickdefinitions[yvariable][view][2]) 
    minor_yticks = np.arange(tickdefinitions[yvariable][view][0], tickdefinitions[yvariable][view][1], 
                             tickdefinitions[yvariable][view][3]) 
    (yymin, yymax) = (limitaxis[yvariable][view][0], limitaxis[yvariable][view][1])
    ylabel = axislabel[yvariable]
    
    #************ initialization of the plot
    if filename == 'all':
        files = sorted(glob.glob('phase_diag_*O8.txt'),reverse=False) #I list every phase diag files
        if combined == 1:
            if type_plot == 'T-P':
                fig, ax1,axbis, ax0 = creation_plot_break(plot_parameters, 
                                                          major_xticks,  minor_xticks,
                                                          major_yticks,  minor_yticks, 
                                                          xlabel, ylabel, 
                                                          xxmin, xxmax, yymin, yymax)
                axis = [ax1,axbis]
            else:
                fig, ax1 = creation_plot(plot_parameters, major_xticks,  minor_xticks, 
                                         major_yticks,  minor_yticks, xlabel, ylabel, 
                                         xxmin, xxmax, yymin, yymax)
                axis = [ax1]
            figurename = 'phase_diag_combined_all'
        else:
            fig,ax1,ax2,ax3, ax0 = creation_plot_3(plot_parameters, major_xticks,  minor_xticks, 
                                                   major_yticks,  minor_yticks, xlabel, ylabel, 
                                                   xxmin, xxmax, yymin, yymax, yvariable)
            #print("axis are ax1",ax1,"ax2",ax2,"ax3",ax3,"ax0",ax0)
            figurename = 'phase_diag_all'
    elif filename2 != '':
        files = [filename, filename2] #we will plot 2 files
        if combined == 1:
            fig, ax1 = creation_plot(plot_parameters, major_xticks,  minor_xticks, 
                                     major_yticks,  minor_yticks, xlabel, ylabel,
                                     xxmin, xxmax, yymin, yymax)
            axis=[ax1]
            figurename = figurename =  'phase_diag_combined_' + filename.split('.txt')[0].split('_')[2] +'_'+ filename2.split('.txt')[0].split('_')[2]
        else:
            #fig,ax1,ax2, ax0 = creation_plot_2_vertical(plot_parameters, major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax, yvariable)
            fig,ax1,ax2, ax0 = creation_plot_2_horizontal(plot_parameters, major_xticks,  minor_xticks,
                                                          major_yticks,  minor_yticks, 
                                                          xlabel, ylabel, xxmin, xxmax,
                                                          yymin, yymax, yvariable)
            print("axis are ax1",ax1,"ax2",ax2,"ax0",ax0)
            figurename = filename.split('.txt')[0] +'_'+ filename2.split('.txt')[0].split('_')[2]
    else:
        fig, ax1 = creation_plot(plot_parameters, major_xticks,  minor_xticks, 
                                 major_yticks,  minor_yticks, xlabel, ylabel,
                                 xxmin, xxmax, yymin, yymax)
        axis=[ax1]
        files = [filename] #I take only the file we want
        figurename = filename.split('.txt')[0]
    
    #************ creation of arrays and plot at the same time
    for file in files:
        print("****** For file",file)
        #******** Initialisations
        #**change of subplot
        if combined == 0:
            if filename == 'all' or filename2 !='':
                if files.index(file) == 0:
                    axis=[ax1]
                    letter = init_letter                    
                if files.index(file) == 1:
                    axis=[ax2]
                    try:
                        #letter = letters[letters.index(init_letter)+3] #for all 3 projections
                        #letter = letters[letters.index(init_letter)+2] #for only 2 projections
                        letter = letters[letters.index(init_letter)+1] #for naming  from up to down
                    except ValueError:
                        print("You haven't indicated any letter")
                        letter = ''
                if files.index(file) == 2:
                    axis=[ax3] 
                    try:
                        #letter = letters[letters.index(init_letter)+6] #for all 3 projections
                        #letter = letters[letters.index(init_letter)+4] #for only 2 projections
                        letter = letters[letters.index(init_letter)+2] #for naming  from up to down
                    except ValueError:
                        print("You haven't indicated any letter")
        else:
            letter = init_letter
            #print("I plot on axis",ax)
        #**extraction compound
        mineral = file.split('_')[2].split('.txt')[0]  
        compounds.append(mineral)
         
        #************ Extraction of data from phase_diag.txt file
        #**initialisation of data dictionnaries for each column (marker type)
        not_started = {'rho':[],'T':[],'P':[]}
        ongoing = {'rho':[],'T':[],'P':[]}
        viscous_fluid = {'rho':[],'T':[],'P':[]}
        fluid = {'rho':[],'T':[],'P':[]}
        mix_liq_gaz = {'rho':[],'T':[],'P':[]}
        O2 = {'rho':[],'T':[],'P':[]}        
        liquid_spinodal = {'rho':[],'T':[],'P':[]}
        gas_spinodal = {'rho':[],'T':[],'P':[]}
        #********* fill the dictionnaries with data
        with open(file,'r') as f:
            line = f.readline()
            while True:
                line = f.readline()
                if not line: break
                else:    
                    entry=line.split('\t')
                    not_started['P'].append(float(entry[0]))
                    not_started['rho'].append(float(entry[1]))
                    not_started['T'] = append_value(not_started['T'],entry,2)
                    ongoing['P'].append(float(entry[0]))
                    ongoing['rho'].append(float(entry[1]))
                    ongoing['T'] = append_value(ongoing['T'],entry,3)
                    viscous_fluid['P'].append(float(entry[0]))
                    viscous_fluid['rho'].append(float(entry[1]))
                    viscous_fluid['T'] = append_value(viscous_fluid['T'],entry,4)
                    fluid['P'].append(float(entry[0]))
                    fluid['rho'].append(float(entry[1]))
                    fluid['T'] = append_value(fluid['T'],entry,5)
                    mix_liq_gaz['P'].append(float(entry[0]))
                    mix_liq_gaz['rho'].append(float(entry[1]))
                    mix_liq_gaz['T'] = append_value(mix_liq_gaz['T'],entry,6)
                    O2['P'].append(float(entry[0]))
                    O2['rho'].append(float(entry[1]))
                    O2['T'] = append_value(O2['T'],entry,7)
                    liquid_spinodal['P'].append(float(entry[0]))
                    liquid_spinodal['rho'].append(float(entry[1]))
                    liquid_spinodal['T'] = append_value(liquid_spinodal['T'],entry,8)
                    gas_spinodal['P'].append(float(entry[0]))
                    gas_spinodal['rho'].append(float(entry[1]))
                    gas_spinodal['T'] = append_value(gas_spinodal['T'],entry,9)
        #print('not_started',not_started)
        #print('fluid',fluid)
        #print('O2',O2)
        #print('liquid_spinodal',liquid_spinodal)
        #*********we remove data in case of P-rho plot
        if yvariable == 'P':
            dictionaries = [not_started,ongoing,viscous_fluid,fluid,mix_liq_gaz,O2,liquid_spinodal,gas_spinodal]
            for dicname in dictionaries:
                #print("*****************Change dicname")
                for ii in range(len(dicname['P'])):
                    #print(dicname['T'][ii])
                    if math.isnan(dicname['T'][ii]):
                        #print("NaN detected")
                        dicname['P'][ii] = float('nan')
        #create list for spinodal curve
        spinodal = {'P':[],'rho':[],'T':[]}
        #for ii in range(len(gas_spinodal['T'])):
        #    if math.isnan(gas_spinodal['T'][ii]):
        #        continue
        #    else:
        #        spinodal[xvariable].append(gas_spinodal[xvariable][ii])
        #        spinodal[yvariable].append(gas_spinodal[yvariable][ii])
        for ii in range(len(liquid_spinodal['T'])):
            if math.isnan(liquid_spinodal['T'][ii]):
                continue
            else:
                spinodal[xvariable].append(liquid_spinodal[xvariable][ii])
                spinodal[yvariable].append(liquid_spinodal[yvariable][ii])
    
        #************ Extraction of litterature melting curves and phase domains
        #data are stored in one .txt file
        litteraturefile = mineral+'/litterature-data_'+mineral+'.txt'
        if (os.path.isfile(litteraturefile)):
            print("data from litterature in file",litteraturefile)
            with open(litteraturefile,'r') as f:
                line = f.readline()
            entry = line.split('\n')[0].split('\t')
        else:
            print('no litterature data for melting curves and phase stability regions')
            entry = ''
        curves = []
        for i in range(0,len(entry),2):
            Pcurve,Tcurve = np.genfromtxt(litteraturefile, delimiter = '\t', usecols = (i,i+1), skip_header = 2, unpack=True)
            curves.append((Pcurve,Tcurve,entry[i+1]))
        
        for ax in axis: 
            #************ Plot points
            zorder +=1
            fillcolor = '#ff6600'
            markercolor = update_dictionnary(fillcolor)
            rect = mpatches.Rectangle((critical_area[mineral][xvariable][0],
                                       critical_area[mineral][yvariable][0]),
        width = critical_area[mineral][xvariable][1]-critical_area[mineral][xvariable][0], 
        height = critical_area[mineral][yvariable][1]-critical_area[mineral][yvariable][0], 
                                        linewidth = 1, fill = True,
                                        edgecolor = fillcolor,
                                        facecolor = markercolor[mineral]+markeralpha[mineral],
                                        joinstyle = 'round', zorder = -99 )
            ax.add_patch(rect)
            #ax.scatter(not_started[xvariable],not_started[yvariable],s=plot_parameters['size_markers'], facecolor = 'r', edgecolor='r', label = 'not started')
            #ax.scatter(ongoing[xvariable],ongoing[yvariable],s=plot_parameters['size_markers'], facecolor = 'orange', edgecolor='orange', label = 'ongoing')
            fillcolor = '#0000ff'
            markercolor = update_dictionnary(fillcolor)
            ax.scatter(fluid[xvariable],fluid[yvariable],
                       s=plot_parameters['size_markers']*3, marker = markertype[mineral],
                       facecolor = markercolor[mineral]+markeralpha[mineral], 
                       edgecolor=fillcolor, label = 'fluid',
                       linewidths =1, zorder = zorder)
            fillcolor = '#008000'
            markercolor = update_dictionnary(fillcolor)
            ax.scatter(mix_liq_gaz[xvariable],mix_liq_gaz[yvariable],
                       s=plot_parameters['size_markers']*3,marker = markertype[mineral],
                       facecolor = markercolor[mineral]+markeralpha[mineral], 
                       edgecolor=fillcolor, label = 'mix liquid + gaz',
                       linewidths =1,zorder = zorder)
            ax.scatter(viscous_fluid[xvariable],viscous_fluid[yvariable],
                       marker = 'x', s=plot_parameters['size_markers']*8, 
                       facecolor = 'grey', edgecolor='grey', 
                       label = 'viscous fluid',zorder = -90)
            #ax.scatter(O2[xvariable],O2[yvariable],s=plot_parameters['size_markers']*10, facecolor = 'none', edgecolor='r', label = 'O$_2$')
    #        ax.scatter(liquid_spinodal[xvariable],liquid_spinodal[yvariable],s=plot_parameters['size_markers']*2, marker = 'P',  facecolor = 'k', edgecolor='k', label = 'liquid spinodal')
            fillcolor = '#000000'
            markercolor = update_dictionnary(fillcolor)
            ax.plot(spinodal[xvariable],spinodal[yvariable],linestyle = linetype[mineral], 
                    linewidth =plot_parameters['size_lines'], 
                    color = fillcolor)
            ax.scatter(liquid_spinodal[xvariable],liquid_spinodal[yvariable],
                       s=plot_parameters['size_markers']*8, marker = 'P',
                       facecolor = markercolor[mineral]+markeralpha[mineral], 
                       edgecolor=fillcolor, label = 'liquid spinodal',
                       linewidths =1,zorder = zorder +10)
    #        ax.plot(gas_spinodal[xvariable],gas_spinodal[yvariable],linestyle = '', linewidth =plot_parameters['size_lines']*2, color='k' , markersize=plot_parameters['size_markers']*0.5, marker = 'P',  markerfacecolor = '0.90', markeredgecolor='k', label = 'liquid spinodal')
             
            #************ Plot literature curves
            if mineral == 'NaAlSi3O8':
                colorstyles = {'Zhang1996':'0.5','Zhang1996s':'g','Tsuchiya2011':'g',
                               'Litvin1993':'g', 'Bell1969s':'g', 'Bell1969':'0.5',
                               'Tutti2000':'0.5','Newton1967':'0.5'}
                linestyles =   {'Zhang1996':'-.','Zhang1996s':'-.','Tsuchiya2011':'-.',
                                'Litvin1993':'-', 'Bell1969s':'-', 'Bell1969':'-',
                                'Tutti2000':'-','Newton1967':'-'} 
            elif mineral == 'KAlSi3O8':
                colorstyles = {'Zhang1996':'0.5','Zhang1996s':'g','Tsuchiya2011':'g',
                               'Urakawa1994':'g', 'Lindsley1966':'g','Akaogi2004':'0.5'}
                linestyles = {'Zhang1996':'-.','Zhang1996s':'-.','Tsuchiya2011':'-.',
                              'Urakawa1994':'--', 'Lindsley1966':'--','Akaogi2004':'--'}
            elif mineral == 'CaAl2Si2O8':
                colorstyles = {'Zhang1996':'0.5','Zhang1996s':'g','Tsuchiya2011':'g',
                               'Hariya1968':'0.5', 'Hariya1968s':'g'}
                linestyles = {'Zhang1996':'-.','Zhang1996s':'-.','Tsuchiya2011':'-.',
                              'Hariya1968':':', 'Hariya1968s':':'}
            
            if yvariable == 'T' and xvariable == 'P':
                for i in range(len(curves)):
                    ref = curves[i][2]
                    print(ref)
                    ax.plot(curves[i][0],curves[i][1],linestyle = linestyles[ref], 
                            color = colorstyles[ref], zorder = -80,
                            marker = None, markersize = 0 , linewidth = plot_parameters['size_lines'])
    
            #************ Text on plot        
            if yvariable == 'T':
                positiontext = {'zoom':{'supercrit':(0.98,0.91),'liquid':(0.98,0.42),
                                        'liquidg':(0.9,3500),'+gas':(0.9,3200)} ,
                             'classic':{'supercrit':(0.52,0.95),'liquid':(0.47,0.52),
                                        'liquidg':(0.7,3500),'+gas':(0.7,3100)}, 
                                'full':{'supercrit':(0.4,0.91),'liquid':(0.5,0.21)}   }
                ax.text(positiontext[view]['supercrit'][0],positiontext[view]['supercrit'][1],
                        "Supercritical" , transform=ax.transAxes, horizontalalignment = 'right', 
                        fontsize = plot_parameters["size_fonts"])
                ax.text(positiontext[view]['liquid'][0],positiontext[view]['liquid'][1], "Liquid" ,
                        transform=ax.transAxes, horizontalalignment = 'right', 
                        fontsize = plot_parameters["size_fonts"])
                if xvariable == 'rho' and view != 'full':
                    ax.text(positiontext[view]['liquidg'][0],positiontext[view]['liquidg'][1], "Liquid" ,
                            horizontalalignment = 'left', fontsize = plot_parameters["size_fonts"])
                    ax.text(positiontext[view]['+gas'][0],positiontext[view]['+gas'][1], "+ Gas" , 
                            horizontalalignment = 'left', fontsize = plot_parameters["size_fonts"])                
            elif yvariable == 'P' and xvariable == 'rho':
                positiontext = {'zoom':{'supercrit':(0.4,0.85),'liquidg':(0.4,0.25),'+gas':(0.4,0.2)} ,
                             'classic':{'supercrit':(0.4,0.85),'liquidg':(0.6,0.25),'+gas':(0.6,0.2)},
                                'full':{'supercrit':(0.4,0.55)}   }
                ax.text(positiontext[view]['supercrit'][0],positiontext[view]['supercrit'][1], "Supercritical" ,
                        transform=ax.transAxes, horizontalalignment = 'right', fontsize = plot_parameters["size_fonts"])
                if view != 'full':
                    ax.text(positiontext[view]['liquidg'][0],positiontext[view]['liquidg'][1], "Liquid" , 
                            transform=ax.transAxes, horizontalalignment = 'left', fontsize = plot_parameters["size_fonts"])
                    ax.text(positiontext[view]['+gas'][0],positiontext[view]['+gas'][1], "+ Gas" , 
                            transform=ax.transAxes, horizontalalignment = 'left', fontsize = plot_parameters["size_fonts"])        
            elif yvariable == 'P' and xvariable == 'T':
                positiontext = {'zoom':{'supercrit':(7150,0.4),'liquid':(4000,0.6)} ,
                             'classic':{'supercrit':(7150,1),'liquid':(4000,1)}, 
                                'full':{'supercrit':(8000,10)}   }
                ax.text(positiontext[view]['supercrit'][0],positiontext[view]['supercrit'][1], "Supercritical" ,
                        horizontalalignment = 'left', rotation = 90, fontsize = plot_parameters["size_fonts"])
                if view != 'full':
                    ax.text(positiontext[view]['liquid'][0],positiontext[view]['liquid'][1], "Liquid" , 
                            horizontalalignment = 'left', fontsize = plot_parameters["size_fonts"])
            if (filename == 'all' or filename2 != '') and init_letter != '':
                ax.text(0.02,0.95, letter , transform=ax.transAxes, horizontalalignment = 'left', 
                        fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                        bbox=dict(facecolor='w', edgecolor='k', pad=3.0))
            elif init_letter != '':
                ax.text(0.02,0.94, init_letter , transform=ax.transAxes, horizontalalignment = 'left', 
                        fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                        bbox=dict(facecolor='w', edgecolor='k', pad=3.0))
        
    #***** Create legend from custom artist/label lists  
    custom_lines = [Line2D([0], [0], marker = markertype[key], 
                           linestyle = linetype[key], markeredgecolor = fillcolor,
                           markerfacecolor = markercolor[key]+markeralpha[key], 
                           color = fillcolor, markersize = 5, 
                           markeredgewidth = 1, linewidth = plot_parameters["size_lines"]) 
                    for key in compounds]
    #format mineral names
    mineralnames = []
    for label in compounds:
        mineralnames.append(format_1label(label))
    legend2 = plt.legend([line for line in custom_lines],[label for label in mineralnames],
                         bbox_to_anchor=(0.5, 1.08), loc="upper center", 
                         fontsize = plot_parameters["size_font_ticks"], 
                         borderaxespad=0., ncol=len(mineralnames))
    
    labelcolor = ['Fluid','Liquid + Gas','Critical point']
    colors = ['#0000ff','#008000','#ff6600']
    custom_patch = [mpatches.Patch(color=colors[ii]) for ii in range(len(colors))]
    legend = plt.legend([col for col in custom_patch],[label for label in labelcolor],
                        bbox_to_anchor=(0.5, 1.08), loc="lower center", 
                        fontsize = plot_parameters["size_font_ticks"],  
                        borderaxespad=0.,ncol=len(labelcolor))
    plt.gca().add_artist(legend2)
    
    #***** Save
    figurename = figurename+'_'+yvariable+'-'+xvariable +'_'+ view + '.pdf'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, dpi = 300, bbox_inches = 'tight')
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












