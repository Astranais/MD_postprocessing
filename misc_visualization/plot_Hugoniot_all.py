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
import re
from scipy.interpolate import interp1d


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
        xlabel = r'Density (g.cm$^{-3}$)'
    elif  xvariable =='T':
        major_xticks = np.arange(0, 25000, 5000) 
        minor_xticks = np.arange(0, 21000, 1000)
        (xxmin, xxmax) = (0, 20500)
        xlabel = r'Temperature (K)'
    elif xvariable == 'P':
        major_xticks = np.arange(0, 500, 50) 
        minor_xticks = np.arange(0, 500, 25)
        (xxmin, xxmax) = (0.1,500)
        #(xxmin, xxmax) = (0,275)
        xlabel = r'Pressure (GPa)'
        ax.set_xscale('log')
    else:
        major_xticks = np.arange(0,1,1)
        minor_xticks = np.arange(0,1, 0.5)    
        (xxmin,xxmax) = (0,1)
        xlabel = 'Not defined yet'
        
    if yvariable == 'T':
        major_yticks = np.arange(0, 25000, 5000) 
        minor_yticks = np.arange(0, 21000, 1000)
        (yymin, yymax) = (0, 20500)
#        major_yticks = np.arange(0, 77000, 1000) 
#        minor_yticks = np.arange(0, 77000, 500)
#        (yymin, yymax) = (273, 7700)
        ylabel = r'Temperature (K)'
    elif yvariable == 'P':
        major_yticks = np.arange(0, 500, 50) 
        minor_yticks = np.arange(0, 500, 25)
        (yymin, yymax) = (0,500)
        ylabel = r'Pressure (GPa)'
    else:
        major_yticks = np.arange(0, 1, 1) 
        minor_yticks = np.arange(0, 1, 0.5) 
        (yymin, yymax) = (0, 1)
        ylabel = 'Not defined yet'
    
    
    #apply setup
    if xvariable != 'P':
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
    column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,
                     'err_T':10,'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,
                     'stdev_Cvm_Nkb':24,'Cvm':25,'stdev_Cvm':26}
    #parameters for the figure for article
    #plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),
                       "size_markers" : 4,"size_lines" :2,"shift_labelpad" : 10}
    markertype = ['>','^','D','v','<'] #one marker per impactor velocity
    impactor_velocities = [12.9,15.2,18.1,8.3,11.5] #in km/s, see bibliography excel file
    plot_params = {}
    fillcolor = {}
    #other dictionnaries and parameters for the figure
    mineralname = 'all'
    mineralname2 = ''
    minerals = []
    letter = ''
    comparison = 0
    try:
        options,arg = getopt.getopt(argv,"hf:g:t:l:c:",["fmineralname","gmineralname",'type', 'letter','comparison' ])
    except getopt.GetoptError:
        print("plot_Hugoniot_all.py -f <mineralname>(default = all files of every compound) -g <mineralname 2 (if we want to plot 2 files)> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> -c < =1 for comparison with thermo data points, default = 0>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_Hugoniot_all.py program to plot the Hugoniot ')
            print("plot_Hugoniot_all.py -f <mineralname>(default = all files of every compound) -g <mineralname 2 (if we want to plot 2 files)> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> -c < =1 for comparison with thermo data points, default = 0>")
            print("plot_Hugoniot_all.py requires to be lauched from the folder containing every Hugoniot .txt file created by analyze_Hugoniot.py")
            print('')
            sys.exit()
        if opt in ('-f','--fmineralname'):
            mineralname = str(arg)
        elif opt in ('-g','--gmineralname'):
            mineralname2 = str(arg)
        elif opt in ('-t','--type'):
            type_plot = str(arg)
            xvariable = type_plot.split('-')[1]
            yvariable = type_plot.split('-')[0]
        elif opt in ('-l','--letter'):
            letter = str(arg)
        elif opt in ('-c','--comparison'):
            comparison = int(arg)
    if mineralname == 'all':
        files = sorted(glob.glob('Hugoniot_*_ground-state_*.txt')) #I list every  files
        removefiles = glob.glob('Hugoniot_*_ground-state_3000-0.txt')
        shockfiles = sorted(glob.glob('Shock-state_*_ground-state_*.txt')) #I list every  files
        figurename = 'Hugoniot_all'
        for file in sorted(files):
            mineral = file.split('_')[1]
            if mineral in minerals:
                continue
            else:
                minerals.append(mineral)
    elif mineralname2 != '':
        files = sorted(glob.glob('Hugoniot_'+mineralname+'_ground-state_*.txt'))
        files.extend(sorted(glob.glob('Hugoniot_'+mineralname2+'_ground-state_*.txt')))
        removefiles = glob.glob('Hugoniot_'+mineralname+'_ground-state_3000-0.txt')
        removefiles.extend(sorted(glob.glob('Hugoniot_'+mineralname2+'_ground-state_3000-0.txt')))
        shockfiles = glob.glob('Shock-state_'+mineralname+'_ground-state_*.txt')
        shockfiles.extend(glob.glob('Shock-state_'+mineralname2+'_ground-state_*.txt'))
        figurename = 'Hugoniot_'+mineralname+'_'+mineralname2
        minerals = [mineralname,mineralname2]
    else:
        files = sorted(glob.glob('Hugoniot_'+mineralname+'_ground-state_*.txt')) #I list only the files we want
        removefiles = ['Hugoniot_'+mineralname+'_ground-state_3000-0.txt']
        shockfiles = sorted(glob.glob('Shock-state_'+mineralname+'_ground-state_*.txt'))
        figurename = 'Hugoniot_'+mineralname
        minerals = [mineralname]
    for file in removefiles:
        files.remove(file)
    fig, ax = creation_plot(plot_parameters, xvariable, yvariable)
    #creation of color changes: 1 color per ground state
    colors = {'2260-3000':'#ff0000','2500-0':'#b1dbff','2585-1932':'#ff8a00',
              '2600-0':'#0080ed','2700-0':'#004f92'}
    #colors = {}
    #for file in files:
    #    gs = file.split('_')[-1].split('.txt')[0]
    #    colors[gs] = ''
    #color = iter(plt.cm.nipy_spectral(np.linspace(0,1,len(colors)))) #Creation of the color list
    #for key in natsort.natsorted(colors):
    #    c = next(color)
    #    colors[key] = c
    #    print(key, to_hex(c) )
    #definition of line and marker style
    for mineral in minerals:
        if mineral == 'NaAlSi3O8':
            plot_params[mineral] = 'o-'  
            fillcolor[mineral] = '#000000ff'
        elif mineral == 'KAlSi3O8':
            plot_params[mineral] = 'o--'   
            fillcolor[mineral] = '#0000007f'
        else:
            plot_params[mineral] = 'o:'  
            fillcolor[mineral] = 'w'
    #************ creation of arrays and plot at the same time
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
            colorfill = colors[gs]
            thermofile = 'thermo_NaAlSi3O8_hard_all.txt'    
        else:
            colorfill = 'w'
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
                ax.plot(thermo[xvariable],thermo[yvariable],'k+',label='data')
            except FileNotFoundError:
                print("file",thermofile,"does not exist")
        #******* interpolation
        xnew = np.linspace(Hugoniot[xvariable][1], Hugoniot[xvariable][-1], num=500, endpoint=True)
        f2 = interp1d(Hugoniot[xvariable][1:],Hugoniot[yvariable][1:], kind='cubic')
        A = (Hugoniot[yvariable][0]-Hugoniot[yvariable][1])/(Hugoniot[xvariable][0]**2 - Hugoniot[xvariable][1]**2 -2*Hugoniot[xvariable][0]*(Hugoniot[xvariable][0]-Hugoniot[xvariable][1]) )
        B = -2*A*Hugoniot[xvariable][0]
        C = Hugoniot[yvariable][0]-A*Hugoniot[xvariable][0]**2-B*Hugoniot[xvariable][0]
        xnewsmall = np.linspace(Hugoniot[xvariable][0], Hugoniot[xvariable][1], num=100, endpoint=True)
        ynewsmall = A*xnewsmall**2+B*xnewsmall+C
        #******* plot the hugoniot lines     
        ax.plot(xnewsmall,ynewsmall, plot_params[mineral], markersize = 0, 
                linewidth = plot_parameters["size_lines"], color =colors[gs])
        ax.plot(xnew,f2(xnew), plot_params[mineral], markersize = 0, 
                linewidth = plot_parameters["size_lines"], color =colors[gs])
        ax.plot(Hugoniot[xvariable][1:],Hugoniot[yvariable][1:], 
                plot_params[mineral], marker = None, linewidth = plot_parameters["size_lines"], color =colors[gs])
        
    for file in shockfiles:
        Shock = {'rho':[],'T':[],'P':[]}
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
                        Shock['P'].append(float(entry[2])) #pressure in GPa
                        Shock['rho'].append(float(entry[0])/1000) #density in g/cm3
                        Shock['T'].append(float(entry[1])) #temperature in K
                    except ValueError:
                        print("Missing data in your file", file)
                        sys.exit()
        #******* change of markers and colors style
        if mineral == 'NaAlSi3O8':
            colorfill = colors[gs]
            alphafill = 1
        elif mineral == 'KAlSi3O8':
            colorfill = colors[gs]
            alphafill = 0.5
        else:
            colorfill = 'w'
            alphafill = 1
        #******* plot the shock condition for each impactor velocity 
        #for ii in range(5):
        #    ax.plot(Shock[xvariable][ii],Shock[yvariable][ii], marker = markertype[ii] , markersize = plot_parameters["size_markers"]+2,markeredgecolor = colors[gs], color = 'none',  markeredgewidth = 0.5,  markerfacecolor = colorfill, alpha=alphafill)
    
    #text on plot
    if (mineralname == 'all' or mineralname2 != '') and letter != '':
        ax.text(0.015,0.945, letter , transform=ax.transAxes, horizontalalignment = 'left', 
                fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
    elif letter != '':
        ax.text(0.06,0.94, letter , transform=ax.transAxes, horizontalalignment = 'left', 
                fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
    
    #Create legend from custom artist/label lists  
    custom_lines = [Line2D([0], [0], marker = plot_params[key][0], linestyle = plot_params[key][1:],
                           markeredgecolor = 'k', markerfacecolor = fillcolor[key], 
                           color = 'k', markersize = plot_parameters["size_markers"], 
                           markeredgewidth = 0.5, linewidth = plot_parameters["size_lines"]) for key in natsort.natsorted(plot_params)]
    #custom_lines = [Line2D([0], [0], linestyle = plot_params[key][1:], color = 'k', linewidth = plot_parameters["size_lines"]) for key in natsort.natsorted(plot_params)]
    #format mineral names
    mineralnames = []
    for label in natsort.natsorted(plot_params):
        mineralnames.append(format_1label(label))
    print(mineralnames)
    legend2 = plt.legend([line for line in custom_lines],[label for label in mineralnames], 
                         bbox_to_anchor=(1.05, 0.55), loc=2, fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0.)

    orderedmarkertype = ['v','<','>','^','D'] #one marker per impactor velocity
    custom_markers = [Line2D([0], [0], marker = orderedmarkertype[ii], linestyle = '', 
                             markeredgecolor = 'k', markerfacecolor = 'k',  color = 'k', 
                             markersize = plot_parameters["size_markers"], markeredgewidth = 0.5, 
                             linewidth = plot_parameters["size_lines"]) for ii in range(len(orderedmarkertype))]
    legend3 = plt.legend([line for line in custom_markers],[str(label) + ' km/s' for label in natsort.natsorted(impactor_velocities)],
                         title = 'Impactor velocity', bbox_to_anchor=(1.05, 0.31), loc=2, 
                         fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0.)
    plt.setp(legend3.get_title(),fontsize= plot_parameters["size_font_ticks"])
    
    orderedkeys = ['2500-0','2600-0','2700-0','2585-1932','2260-3000']
    #custom_patch = [mpatches.Patch(color=colors[key]) for key in natsort.natsorted(colors)]
    custom_patch = [mpatches.Patch(color=colors[key]) for key in orderedkeys]
    #legend = plt.legend([col for col in custom_patch],[str(round(float(label.split('-')[0])/1000,2)) + '  -  ' + label.split('-')[1] for label in natsort.natsorted(colors)],title = '$\\rho_0$ (g$.cm^{-3}$) - T$_0$ (K)', bbox_to_anchor=(1.05, 1), loc=2, fontsize = plot_parameters["size_fonts"],  borderaxespad=0.)
    legend = plt.legend([col for col in custom_patch],[str(round(float(label.split('-')[0])/1000,2)) + '  -  ' + label.split('-')[1] for label in orderedkeys],title = '$\\rho_0$ (g$.cm^{-3}$) - T$_0$ (K)', bbox_to_anchor=(1.05, 1), loc=2, fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0.)
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])

    plt.gca().add_artist(legend2)
    
    plt.gca().add_artist(legend3)
    
    figurename = figurename+'_'+yvariable+'-'+xvariable + '.pdf'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












