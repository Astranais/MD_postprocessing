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
        (xxmin,xxmax) = (0.1,500)
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


def interpolate_gaps(yvalues, xvalues, limit=None):
    """
    Fill gaps using linear interpolation, optionally only fill gaps up to a
    size of `limit`.
    """
    xvalues = np.asarray(xvalues)
    yvalues = np.asarray(yvalues)
    valid = np.isfinite(yvalues)
    filled = np.interp(xvalues, xvalues[valid], yvalues[valid])
    if limit is not None:
        invalid = ~valid
        for n in range(1, limit+1):
            invalid[:-n] &= invalid[n:]
        filled[invalid] = np.nan    
    return filled


def extractdata(file,Hugoniot_all,gsT,Hugoniot):
    with open(file,'r') as f:
        f.readline()
        f.readline()
        while True:
            line = f.readline()
            if not line: break
            else:    
                entry=line.split('\t')
                try:
                    Hugoniot_all[gsT]['P'].append(float(entry[2])) #pressure in GPa
                    Hugoniot_all[gsT]['rho'].append(float(entry[1])/1000) #density in g/cm3
                    Hugoniot_all[gsT]['T'].append(float(entry[0])) #temperature in K
                    Hugoniot['P'].append(float(entry[2])) #pressure in GPa
                    Hugoniot['rho'].append(float(entry[1])/1000) #density in g/cm3
                    Hugoniot['T'].append(float(entry[0])) #temperature in K
                except ValueError:
                    print("Missing data in your file", file)
                    sys.exit() 

def main(argv):
    """     ********* Main program *********     """
    column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,
                     'err_T':10,'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,
                     'stdev_Cvm_Nkb':24,'Cvm':25,'stdev_Cvm':26}
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),
                       "size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    #plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),"size_markers" : 4,"size_lines" : 2,"shift_labelpad" : 10}
    
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
        print("plot_Hugoniot_all_enveloppes.py -f <mineralname>(default = all files of every compound) -g <mineralname 2 (if we want to plot 2 files)> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> -c < =1 for comparison with thermo data points, default = 0>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_Hugoniot_all_enveloppes.py program to plot the Hugoniot lines as enveloppes (highest and lowest only) ')
            print("plot_Hugoniot_all_enveloppes.py -f <mineralname>(default = all files of every compound) -g <mineralname 2 (if we want to plot 2 files)> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> -c < =1 for comparison with thermo data points, default = 0>")
            print("plot_Hugoniot_all_enveloppes.py requires to be lauched from the folder containing every Hugoniot .txt file created by analyze_Hugoniot.py")
            print("plot_Hugoniot_all_enveloppes.py requires to be lauched from the folder containing also every Shock-state .txt file created by plot_Hugoniot_impedance-match.py")
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
        figurename = 'Hugoniot_all_enveloppe'
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
        figurename = 'Hugoniot_'+mineralname+'_'+mineralname2 + '_enveloppe'
        minerals = [mineralname,mineralname2]
    else:
        files = sorted(glob.glob('Hugoniot_'+mineralname+'_ground-state_*.txt')) #I list only the files we want
        removefiles = ['Hugoniot_'+mineralname+'_ground-state_3000-0.txt']
        shockfiles = sorted(glob.glob('Shock-state_'+mineralname+'_ground-state_*.txt'))
        figurename = 'Hugoniot_'+mineralname + '_enveloppe'
        minerals = [mineralname]
    for file in removefiles:
        files.remove(file)
    fig, ax = creation_plot(plot_parameters, xvariable, yvariable)
    #creation of color changes: 1 color per ground state
    colors = {'2260-3000':'#ff0000','2500-0':'#b1dbff','2585-1932':'#ff8a00','2600-0':'#0080ed','2700-0':'#004f92'}
    colors_enveloppes = {'3000':'#ff0000','1932':'#ff8a00','0':'#0080ed'}
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
    #************ creation of enveloppes for each hugoniot thermal state (cold, warm and hot)                                             
    Hugoniot_all = {}
    enveloppeUp = {}
    enveloppeDown = {}
    Xenveloppe = {}                
    for file in files:
        #**initialisation of data dicitonnaries
        gsT =  file.split('_')[-1].split('.txt')[0].split('-')[1]
        Hugoniot_all[gsT] = {'rho':[],'T':[],'P':[]}
        enveloppeUp[gsT] = {'x':[],'y':[]}
        enveloppeDown[gsT] = {'x':[],'y':[]}
        Xenveloppe[gsT] = []
    for file in files:
        mineral = file.split('_')[1]
        gs =  file.split('_')[-1].split('.txt')[0]
        gsT =  file.split('_')[-1].split('.txt')[0].split('-')[1]
        Hugoniot = {'rho':[],'T':[],'P':[]}
        #********* fill the dictionnaries with data only for the extreme feldspars (those making the enveloppe)
        f2, xnew = '',''
        if mineralname == 'all':
            if mineral != 'KAlSi3O8' and gs != '2600-0':
                extractdata(file,Hugoniot_all,gsT,Hugoniot)
                #******* interpolation
                xnew = np.linspace(Hugoniot[xvariable][1], Hugoniot[xvariable][-1], num=500, endpoint=True)
                f2 = interp1d(Hugoniot[xvariable][1:],Hugoniot[yvariable][1:], kind='cubic')
                A = (Hugoniot[yvariable][0]-Hugoniot[yvariable][1])/(Hugoniot[xvariable][0]**2 - Hugoniot[xvariable][1]**2 -2*Hugoniot[xvariable][0]*(Hugoniot[xvariable][0]-Hugoniot[xvariable][1]) )
                B = -2*A*Hugoniot[xvariable][0]
                C = Hugoniot[yvariable][0]-A*Hugoniot[xvariable][0]**2-B*Hugoniot[xvariable][0]
                xnewsmall = np.linspace(Hugoniot[xvariable][0], Hugoniot[xvariable][1], num=100, endpoint=True)
                ynewsmall = A*xnewsmall**2+B*xnewsmall+C
                #******* plot the hugoniot lines     
                #ax.plot(xnewsmall,ynewsmall, plot_params[mineral], markersize = 0, linewidth = plot_parameters["size_lines"], color =colors[gs])
                #ax.plot(xnew,f2(xnew), plot_params[mineral], markersize = 0, linewidth = plot_parameters["size_lines"], color =colors[gs])
                #ax.plot(Hugoniot[xvariable][1:],Hugoniot[yvariable][1:], plot_params[mineral], marker = None, linewidth = plot_parameters["size_lines"], color =colors[gs])
                
                #******* interpolation
                #f2 = interp1d(Hugoniot[xvariable],Hugoniot[yvariable], kind='cubic')
                #xnew = np.linspace(Hugoniot[xvariable][0], Hugoniot[xvariable][-1], num=500, endpoint=True)    
                #print('1',mineral,gs)
            if mineral == 'NaAlSi3O8' and gs != '2500-0':
                if f2 != '':
                    enveloppeDown[gsT]['x'].extend(xnewsmall)
                    enveloppeDown[gsT]['x'].extend(xnew)
                    enveloppeDown[gsT]['y'].extend(ynewsmall)
                    enveloppeDown[gsT]['y'].extend(f2(xnew))
                    Xenveloppe[gsT].extend(xnewsmall)
                    Xenveloppe[gsT].extend(xnew)
                else:
                    enveloppeDown[gsT]['x'].extend(Hugoniot[xvariable])
                    enveloppeDown[gsT]['y'].extend(Hugoniot[yvariable])
                    Xenveloppe[gsT].extend(Hugoniot[xvariable])
                #print('2',mineral,gs)
                #print("mineral ", mineral, len(enveloppeDown[gsT]['x']),len(enveloppeDown[gsT]['y']))
            elif mineral == 'CaAl2Si2O8' and gs != '2700-0':
                if f2 != '':
                    enveloppeUp[gsT]['x'].extend(xnewsmall)
                    enveloppeUp[gsT]['x'].extend(xnew)
                    enveloppeUp[gsT]['y'].extend(ynewsmall)
                    enveloppeUp[gsT]['y'].extend(f2(xnew))
                    Xenveloppe[gsT].extend(xnewsmall)
                    Xenveloppe[gsT].extend(xnew)
                else:
                    enveloppeUp[gsT]['x'].extend(Hugoniot[xvariable])
                    enveloppeUp[gsT]['y'].extend(Hugoniot[yvariable])
                    Xenveloppe[gsT].extend(Hugoniot[xvariable])
                #print('3',mineral,gs)                
        else: #for comparison between Na and K only
            if gs != '2600-0':
                extractdata(file,Hugoniot_all,gsT,Hugoniot)
                #print('1',mineral,gs)
            if mineral == 'NaAlSi3O8' and gs != '2500-0':
                enveloppeDown[gsT]['x'].extend(Hugoniot[xvariable])
                enveloppeDown[gsT]['y'].extend(Hugoniot[yvariable])
                Xenveloppe[gsT].extend(Hugoniot[xvariable])
                #print('2',mineral,gs)
                #print("mineral ", mineral, len(enveloppeDown[gsT]['x']),len(enveloppeDown[gsT]['y']))
            elif mineral == 'KAlSi3O8' and gs != '2700-0':
                enveloppeUp[gsT]['x'].extend(Hugoniot[xvariable])
                enveloppeUp[gsT]['y'].extend(Hugoniot[yvariable])
                Xenveloppe[gsT].extend(Hugoniot[xvariable])
                #print('3',mineral,gs)           
            #print("mineral ", mineral, len(enveloppeUp[gsT]['x']),len(enveloppeUp[gsT]['y']))
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
            thermo = {'rho':[],'T':[],'P':[]}         
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
        #******* plot the hugoniot lines
        if gsT !='0':
            if f2 != '':
                if yvariable == 'P' and xvariable == 'T':
                    ax.plot(Hugoniot[xvariable][:2],Hugoniot[yvariable][:2], plot_params[mineral],
                            markersize = 0, linewidth = plot_parameters["size_lines"], color =colors[gs])
                else:
                    ax.plot(xnewsmall,ynewsmall, plot_params[mineral], markersize = 0, 
                            linewidth = plot_parameters["size_lines"], color =colors[gs])
                ax.plot(xnew,f2(xnew), plot_params[mineral], markersize = 0, 
                        linewidth = plot_parameters["size_lines"], color =colors[gs])        
            else:
                ax.plot(Hugoniot[xvariable],Hugoniot[yvariable], plot_params[mineral],
                        markersize = 0, linewidth = plot_parameters["size_lines"], color =colors[gs])
        else:
            if mineralname == 'all':
                if mineral == 'CaAl2Si2O8' and gs == '2500-0'  or mineral == 'NaAlSi3O8' and gs == '2700-0':
                    if f2 != '':
                        if yvariable == 'P' and xvariable == 'T':
                            ax.plot(Hugoniot[xvariable][:2],Hugoniot[yvariable][:2], 
                                    plot_params[mineral], markersize = 0, 
                                    linewidth = plot_parameters["size_lines"], color =colors[gs])
                        else:
                            ax.plot(xnewsmall,ynewsmall, plot_params[mineral], 
                                    markersize = 0, linewidth = plot_parameters["size_lines"], color =colors[gs])
                        ax.plot(xnew,f2(xnew), plot_params[mineral], 
                                markersize = 0, linewidth = plot_parameters["size_lines"], color =colors[gs])        
                    else:
                        ax.plot(Hugoniot[xvariable],Hugoniot[yvariable], plot_params[mineral], 
                                markersize = 0, linewidth = plot_parameters["size_lines"], color =colors[gs])                
            else: #for comparison with only Na and K
                if mineral == 'KAlSi3O8' and gs == '2500-0'  or mineral == 'NaAlSi3O8' and gs == '2700-0':
                    ax.plot(Hugoniot[xvariable],Hugoniot[yvariable], plot_params[mineral], 
                            markersize = 0, linewidth = plot_parameters["size_lines"], color =colors[gs])                
    
    if yvariable == 'P' and xvariable == 'T':
        enveloppeUp = {}
        enveloppeDown = {}
        Xenveloppe = {}
        for file in files:
            #**initialisation of data dicitonnaries
            gsT =  file.split('_')[-1].split('.txt')[0].split('-')[1]
            enveloppeUp[gsT] = []
            enveloppeDown[gsT] = []
            Xenveloppe[gsT] = []
            #** sort all data relative to the x variable
            Hugoniot_all[gsT][xvariable],Hugoniot_all[gsT][yvariable] = zip(*sorted(zip( Hugoniot_all[gsT][xvariable],Hugoniot_all[gsT][yvariable]   )))
            #** create the enveloppes
            for ii in range(len(Hugoniot_all[gsT][xvariable])):
                enveloppeUp[gsT] = [Hugoniot_all[gsT][yvariable][0]]
                enveloppeDown[gsT]  = [Hugoniot_all[gsT][yvariable][0]]
                Xenveloppe[gsT] = [Hugoniot_all[gsT][xvariable][0]]
                j = 0
                for i in range(1,len(Hugoniot_all[gsT][xvariable])):
                    if Hugoniot_all[gsT][xvariable][i] == Hugoniot_all[gsT][xvariable][i-1]:
                        if Hugoniot_all[gsT][yvariable][i] > enveloppeUp[gsT][j]:
                            enveloppeUp[gsT][j] = Hugoniot_all[gsT][yvariable][i]
                        if Hugoniot_all[gsT][yvariable][i] < enveloppeDown[gsT][j]:
                            enveloppeDown[gsT][j] = Hugoniot_all[gsT][yvariable][i]
                    else:
                        j+=1
                        Xenveloppe[gsT].append(Hugoniot_all[gsT][xvariable][i])
                        enveloppeUp[gsT].append(Hugoniot_all[gsT][yvariable][i])
                        enveloppeDown[gsT].append(Hugoniot_all[gsT][yvariable][i])
                print('*************** gsT:', gsT)
                print(Xenveloppe[gsT])
                print(enveloppeUp[gsT])
                print(enveloppeDown[gsT])
        #** plot
        for gsT in Hugoniot_all:
            ax.fill_between(Xenveloppe[gsT],enveloppeUp[gsT],enveloppeDown[gsT], 
                            facecolor = colors_enveloppes[gsT], linewidth=plot_parameters['size_lines'], alpha = 0.1)
    else:
        #**** create the enveloppes
        for gsT in Xenveloppe:
            print(gsT)
            #complete the y data for each x (need to have the same x axis)
            #print(len(enveloppeDown[gsT]['x']), len(enveloppeDown[gsT]['y']), len(enveloppeUp[gsT]['x']), len(enveloppeUp[gsT]['y']))
            #print(enveloppeDown[gsT]['x'],enveloppeDown[gsT]['y'])
            #print(enveloppeUp[gsT]['x'],enveloppeUp[gsT]['y'])
            for x in Xenveloppe[gsT]:
                if Xenveloppe[gsT].count(x) < 2:
                    Xenveloppe[gsT].append(x)
                    if x in enveloppeDown[gsT]['x']:
                        enveloppeUp[gsT]['x'].append(x)
                        enveloppeUp[gsT]['y'].append(float('nan'))
                    else:
                        enveloppeDown[gsT]['x'].append(x)
                        enveloppeDown[gsT]['y'].append(float('nan'))
            #print(len(enveloppeDown[gsT]['x']), len(enveloppeDown[gsT]['y']), len(enveloppeUp[gsT]['x']), len(enveloppeUp[gsT]['y']))
            #sort all data relative to the x variable
            enveloppeDown[gsT]['x'],enveloppeDown[gsT]['y'] = zip(*sorted(zip( enveloppeDown[gsT]['x'],enveloppeDown[gsT]['y']  )))
            enveloppeUp[gsT]['x'],enveloppeUp[gsT]['y'] = zip(*sorted(zip( enveloppeUp[gsT]['x'],enveloppeUp[gsT]['y']  )))
            #print(enveloppeDown[gsT]['x'],'\n', enveloppeDown[gsT]['y'], '\n', enveloppeUp[gsT]['x'], '\n', enveloppeUp[gsT]['y'])
            #fill the gaps with interpolated values
            enveloppeDown[gsT]['y'] = interpolate_gaps(enveloppeDown[gsT]['y'],enveloppeDown[gsT]['x'], limit=None)
            enveloppeUp[gsT]['y']  = interpolate_gaps(enveloppeUp[gsT]['y'],enveloppeUp[gsT]['x'], limit=None)
            #special adjustements for P-rho projection
            if yvariable == 'P':
                #replace last interpolated value by the other one
                if enveloppeUp[gsT]['y'][-1]<enveloppeDown[gsT]['y'][-1]:
                    enveloppeUp[gsT]['y'][-1] = enveloppeDown[gsT]['y'][-1]
                if enveloppeUp[gsT]['y'][-2]<enveloppeDown[gsT]['y'][-2]:
                    enveloppeUp[gsT]['y'][-2] = enveloppeDown[gsT]['y'][-1]
            #print(enveloppeDown[gsT]['x'],enveloppeDown[gsT]['y'])
            #print(enveloppeUp[gsT]['x'],enveloppeUp[gsT]['y'])
        #**** plot the enveloppes
        for gsT in Hugoniot_all:
            #ax.plot(enveloppeDown[gsT]['x'], enveloppeDown[gsT]['y'],color = colors_enveloppes[gsT], linewidth=plot_parameters['size_lines'])
            #ax.plot(enveloppeDown[gsT]['x'], enveloppeUp[gsT]['y'],color = colors_enveloppes[gsT], linewidth=plot_parameters['size_lines'])
            ax.fill_between(enveloppeDown[gsT]['x'],enveloppeUp[gsT]['y'],
                            enveloppeDown[gsT]['y'], interpolate = True, 
                            facecolor = colors_enveloppes[gsT], linewidth=plot_parameters['size_lines'], alpha = 0.3)
    
    
    
    
    #************ creation of arrays and plot at the same time shock data
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
        for ii in range(5):
            ax.plot(Shock[xvariable][ii],Shock[yvariable][ii], marker = markertype[ii], 
                    markersize = plot_parameters["size_markers"]+2,markeredgecolor = colors[gs], 
                    color = 'none',  markeredgewidth = 0.5,  markerfacecolor = colorfill, alpha=alphafill)
    
    #text on plot
    if (mineralname == 'all' or mineralname2 != '') and letter != '':
        ax.text(0.015,0.945, letter , transform=ax.transAxes, horizontalalignment = 'left', 
                fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
    elif letter != '':
        ax.text(0.06,0.94, letter , transform=ax.transAxes, horizontalalignment = 'left', 
                fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
    
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
    legend2 = plt.legend([line for line in custom_lines],[label for label in mineralnames], bbox_to_anchor=(1.05, 0.55), loc=2, fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0.)

    orderedmarkertype = ['v','<','>','^','D'] #one marker per impactor velocity
    custom_markers = [Line2D([0], [0], marker = orderedmarkertype[ii], linestyle = '', 
                             markeredgecolor = 'k', markerfacecolor = 'k',  color = 'k', 
                             markersize = plot_parameters["size_markers"], markeredgewidth = 0.5, 
                             linewidth = plot_parameters["size_lines"]) for ii in range(len(orderedmarkertype))]
    legend3 = plt.legend([line for line in custom_markers],[str(label) + ' km/s' for label in natsort.natsorted(impactor_velocities)], 
                         title = 'Impactor velocity', bbox_to_anchor=(1.05, 0.31), 
                         loc=2, fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0.)
    plt.setp(legend3.get_title(),fontsize= plot_parameters["size_font_ticks"])
    
    orderedkeys = ['2500-0','2600-0','2700-0','2585-1932','2260-3000']
    #custom_patch = [mpatches.Patch(color=colors[key]) for key in natsort.natsorted(colors)]
    custom_patch = [mpatches.Patch(color=colors[key]) for key in orderedkeys]
    #legend = plt.legend([col for col in custom_patch],[str(round(float(label.split('-')[0])/1000,2)) + '  -  ' + label.split('-')[1] for label in natsort.natsorted(colors)],title = '$\\rho_0$ (g$.cm^{-3}$) - T$_0$ (K)', bbox_to_anchor=(1.05, 1), loc=2, fontsize = plot_parameters["size_fonts"],  borderaxespad=0.)
    legend = plt.legend([col for col in custom_patch],[str(round(float(label.split('-')[0])/1000,2)) + '  -  ' + label.split('-')[1] for label in orderedkeys],
                        title = '$\\rho_0$ (g$.cm^{-3}$) - T$_0$ (K)', 
                        bbox_to_anchor=(1.05, 1), loc=2, fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0.)
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







