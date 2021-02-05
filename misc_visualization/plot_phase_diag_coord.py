#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the phase diagram in T-P projection   ****
                            with symbols of major coordinence
                                and melting curves

"""


#     ********* Importation of the packages and modules used here *********
import glob
import sys
import getopt
import os
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap








def creation_plot_2_horizontal(plot_parameters):
    """     ********** Creation of the plot  **********    """
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
     
    ax1.set_xticks(np.arange(-1,1.5,0.5))
    ax1.set_xticks(np.arange(-1,1.5,0.1), minor=True)

    major_yticks = np.arange(0,8000,1000)
    minor_yticks = np.arange(0,8000,500)
    
    #apply setup
    for ax in [ax1,ax2]:
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks,
                       width = size_lines/2)
        ax.set_yticks(major_yticks)                                        
        ax.set_yticks(minor_yticks, minor=True)
        ax.yaxis.set_ticks_position('both')
        ax.set_ylim(273,7700)
        
    ax1.set_xlim(-1,1)
    ax2.set_xlim(1,275)
    
    plt.setp(ax2.get_xticklabels()[1], visible=False) 
    ax2.tick_params(labelleft=False)
    ax2.set_xscale('log')
    
    ax = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, 
                    bottom=False, labelleft=False, left=False, 
                    labelright = False, right=False)
    ax.set_xlabel(r'Pressure (GPa)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad)
    ax.set_ylabel(r'Temperature (K)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad)
    return fig, ax1, ax2, ax


def creation_plot_horizontal(plot_parameters):
    """     ********** Creation of the plot  **********    """
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    #plot
    plt.close()
    fig, ax1 = plt.subplots(1,1, figsize = (size_figure[0],size_figure[1]))
    
    major_yticks = np.arange(0,8000,1000)
    minor_yticks = np.arange(0,8000,500)
    
    #apply setup
    for ax in [ax1]:
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks,
                       width = size_lines/2)
        ax.set_yticks(major_yticks)                                        
        ax.set_yticks(minor_yticks, minor=True)
        ax.yaxis.set_ticks_position('both')
        ax.set_ylim(273,7700)
        
    #ax1.set_xlim(0.1,375)
    ax1.set_xlim(0.1,275)
    ax1.set_xscale('log')
    
    ax = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, 
                    bottom=False, labelleft=False, left=False, 
                    labelright = False, right=False)
    ax.set_xlabel(r'Pressure (GPa)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad*20)
    ax.set_ylabel(r'Temperature (K)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad*40)
    return fig, ax1, ax

def append_value(data,entry,i):
    """      ******** append a float to the list, replacing empty str by NaN ****** """
    if entry[i] == '':
        data.append(float('nan'))
    elif entry[i] == '\n':
        data.append(float('nan'))
    else:
        data.append(int(entry[i]))
    return data





def extract_major_coordinence(percfile, atoms, major_coord, max_coord,addtxt):
    """extract the major coordination polyhedra for each file and the selected type of polehydra """
    data_coord = {}
    coord_polyhedra = []
    with open(percfile,'r') as f:        
        line = f.readline()
        entry=line.split()
        files = entry[1:]
        #creation of the x data list
        for file in files:
            data_coord[file.split('/')[-1]+addtxt] = []
    with open(percfile,'r') as f:
        #creation of the y data list and plot for each cluster
        while True:
            line = f.readline()
            if not line: break
            else:
                entry=line.split('\n')[0].split('\t')
                if len(entry) > 1:
                    if entry[0][:len(atoms)] == atoms:
                        coord_polyhedra.append(entry[0])
                        for i in range(1,len(entry)): #we loop over all the files
                            if entry[i] == '':
                                entry[i] = 0
                            data_coord[files[i-1].split('/')[-1]+addtxt].append(float(entry[i]))
    for file in data_coord:
        major_coord[file] = int(coord_polyhedra[np.argmax(np.array(data_coord[file]))].split('_')[-1])
        if major_coord[file] > max_coord:
            max_coord = major_coord[file]
    return major_coord, max_coord



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
    filename = ''
    major_coord = {} #dictionnary with the number of the major coordination polyhedra in each simufile
    max_coord= 0
    TP = {} #dictionnary with P and T in each simufile
    column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,
                     'err_T':10,'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,
                     'stdev_Cvm_Nkb':24,'Cvm':25,'stdev_Cvm':26,'testCv':31,'stdev_testCv':32}
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (9,4)
    ,"size_markers" : 50,"size_lines" : 2,"shift_labelpad" : 1} 
    letter = ''
    coord_marker = ['','o','X','^','s','p','h',(7,0,0),(8,0,0),(9,0,0),(10,0,0),(11,0,0),(12,0,0)]
    #location of the critical point
    critical_area = {'NaAlSi3O8':{'P':(0.05,0.15),'T':(5500,6000),'rho':(0.51,0.83)} ,
                                  'KAlSi3O8':{'P':(0,0.11),'T':(5000,5500),'rho':(0.56,0.87)} , 
                                  'CaAl2Si2O8':{'P':(0.1,0.2),'T':(7000,7500),'rho':(0.57,0.79)}  }
    # with constrained fit {'NaAlSi3O8':{'P':(0.1,0.2),'T':(6000,6500),'rho':(0.44,0.62)} 
    try:
        options,arg = getopt.getopt(argv,"hf:l:a:t:",["filename",'letter','atoms','type'])
    except getopt.GetoptError:
        print("plot_phase_diag_coord.py -f <phase diag filename> -l <letter of first plot, default = ''> -a <Al_1O or Si_1O (for major coordination)> -t <plot type: 1=one log plot P>0.1, 2=two side by side plots (linear+log)>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_phase_diag_coord.py program to plot the T-P projection of the phase diagram of one feldspar with indication of major coordination')
            print("plot_phase_diag_coord.py -f <phase diag filename> -l <letter of first plot, default = ''> -a <Al_1O or Si_1O (for major coordination)> -t <plot type: 1=one log plot P>0.1, 2=two side by side plots (linear+log)>")
            print("plot_phase_diag_coord.py requires to be launched from the folder with every file needed (perc files from stat-concentrate, melting curves from litterature, fullthermo.txt)")
            print('')
            sys.exit()
        if opt in ('-f','--filename'):
            filename = str(arg)
        elif opt in ('-l','--letter'):
            letter = str(arg)
        elif opt in ('-a','--atoms'):
            atoms = str(arg)
        elif opt in ('-t','--type'):
            plottype=int(arg)
            
    #************ initialization of the plot
    if plottype == 2:
        fig,ax1,ax2, ax0 = creation_plot_2_horizontal(plot_parameters)
        axis = [ax1,ax2]
        if letter != '':
            ax0.text(-0.03,1.03, letter , transform=ax0.transAxes, horizontalalignment = 'left',
                     fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                     bbox=dict(facecolor='w', edgecolor='k', pad=3.0))
    else:
        fig, ax1, ax0 = creation_plot_horizontal(plot_parameters)
        axis = [ax1]
        ax2 = ax1
        if letter != '':
            ax0.text(-0.1,0.95, letter , transform=ax0.transAxes, horizontalalignment = 'left', 
                     fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                     bbox=dict(facecolor='w', edgecolor='k', pad=3.0))
    figurename = filename.split('/')[-1].split('.txt')[0]+'_T-P_'+atoms.split('_')[0]+atoms.split('_')[1]+'_-t'+str(plottype)
    colors_coord = []
    cm_data = np.loadtxt("/Users/akobsch/Dropbox/Recherche/PhD-MD-silicates/simulations/scripts/Analysis/ScientificColourMaps6/romaO/romaO.txt")
    cm_data = cm_data[::-1] #for reverse colors
    new_map = LinearSegmentedColormap.from_list('new', cm_data)
    color = iter(new_map(np.linspace(0,1,8)))
    for i in range(8):
        c = next(color)
        colors_coord.append(c)
    #print(colors_coord,len(colors_coord))
    
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
    #** fill the dictionnaries with data
    with open(filename,'r') as f:
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
    xvariable = 'P'
    yvariable = 'T'
    #create list for spinodal curve
    spinodal = {'P':[],'rho':[],'T':[]}
    for ii in range(len(gas_spinodal['T'])):
        if math.isnan(gas_spinodal['T'][ii]):
            continue
        else:
            spinodal[xvariable].append(gas_spinodal[xvariable][ii])
            spinodal[yvariable].append(gas_spinodal[yvariable][ii])
    for ii in range(len(liquid_spinodal['T'])):
        if math.isnan(liquid_spinodal['T'][ii]):
            continue
        else:
            spinodal[xvariable].append(liquid_spinodal[xvariable][ii])
            spinodal[yvariable].append(liquid_spinodal[yvariable][ii])
    
    #************ Extraction of  major coordination for each simufile of each perc file
    mineralname = filename.split('_')[-1].split('.txt')[0]
    #** we identify all the perc files    
    coordfiles = sorted(glob.glob(mineralname+'*perc-T0-L208-R0.dat'))
    soft_coordfiles = sorted(glob.glob('soft_pseudo/'+mineralname+'*perc-T0-L208-R0.dat'))
    #** we extract the major coordination for each simufile in each coordfile
    for file in coordfiles:
        major_coord, max_coord = extract_major_coordinence(file, atoms, major_coord,max_coord,'')
    for file in soft_coordfiles:
        major_coord, max_coord = extract_major_coordinence(file, atoms, major_coord,max_coord,'soft')    
    print('the max coordination is', max_coord)
    #************ Extraction of  P and T for each simufile of each thermo file
    thermofiles = sorted(glob.glob('thermo_'+mineralname+'*.txt'))
    soft_thermofiles = sorted(glob.glob('soft_pseudo/thermo_'+mineralname+'*.txt'))
    for file in thermofiles:
        TP = extract_TP(file, column_number, TP,'')
    for file in soft_thermofiles:
        TP = extract_TP(file, column_number, TP,'soft')
    #print(TP)
    
    #************ Extraction of litterature melting curves and phase domains
    #data are stored in one .txt file
    litteraturefile = 'litterature-data_'+mineralname+'.txt'
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
    
    #************ Extraction of phases or coordination number
    phasefile = 'phases_'+mineralname+'.txt'
    if (os.path.isfile(phasefile)):
        print("data from litterature in file",phasefile)
        phases = np.loadtxt(phasefile,delimiter = '\t', skiprows=1,usecols=(2,3,4,5))
    else:
        print('no phase data')
        phases = [[1,1,'0','0']]
    
    #************ Extraction of additional solid points
    solidpointsfile = 'solid-points_'+mineralname+'.txt'
    if (os.path.isfile(solidpointsfile)):
        print("additional solid data points from litterature in file",solidpointsfile)
        solidpoints_P, solidpoints_T = np.loadtxt(solidpointsfile,delimiter = '\t', skiprows=2,usecols=(0,1), unpack = True)
    else:
        print('no additional solid data points')
        solidpoints_P, solidpoints_T = [np.nan],[np.nan]
    
    #************ Plot
    #** first plot = zoom on low pressures      
    rect = mpatches.Rectangle((critical_area[mineralname][xvariable][0],
                               critical_area[mineralname][yvariable][0]), 
        width = critical_area[mineralname][xvariable][1]-critical_area[mineralname][xvariable][0], 
        height = critical_area[mineralname][yvariable][1]-critical_area[mineralname][yvariable][0], 
        linewidth = None, fill = True, color = 'orange', joinstyle = 'round', zorder = -99 )
    #ax1.add_patch(rect)
       
    #ax1.plot(spinodal[xvariable],spinodal[yvariable],linestyle = '--', 
    #         linewidth =plot_parameters['size_lines']*2, color='k' )
    #ax1.plot(liquid_spinodal[xvariable],liquid_spinodal[yvariable],
    #         linestyle = '', linewidth =plot_parameters['size_lines']*2, 
    #         color='k' , markersize=plot_parameters['size_markers']*0.25, 
    #         marker = 'P',  markerfacecolor = 'k', markeredgecolor='k', 
    #         label = 'liquid spinodal')
    
    #ax1.plot(gas_spinodal[xvariable],gas_spinodal[yvariable],linestyle = '',
    #         linewidth =plot_parameters['size_lines']*2, color='k' , 
    #         markersize=plot_parameters['size_markers']*0.5, marker = 'P', 
    #         markerfacecolor = '0.90', markeredgecolor='k', label = 'gas spinodal')
    
    #for file in major_coord:
    #    if TP[file][1] < 1:
    #        ax1.scatter(TP[file][1],TP[file][0],marker = coord_marker[major_coord[file]], 
    #                s=plot_parameters['size_markers'], edgecolors='k',
    #                facecolor = colors_coord[major_coord[file]-1],
    #                linewidths =0.5 , zorder = -major_coord[file] )
    if atoms[0] == 'A':
        idx = 2
    else:
        idx = 3
    #ax1.text(phases[0][0],phases[0][1],str(int(phases[0][idx])),fontsize=plot_parameters['size_font_ticks'])
    
    #** second plot = coordination change
    for file in major_coord:
        ax2.scatter(TP[file][1],TP[file][0],marker = coord_marker[major_coord[file]], 
                    s=plot_parameters['size_markers'], edgecolors='k',
                    facecolor = colors_coord[major_coord[file]-1],
                    linewidths =0.5 ,zorder = -major_coord[file] )    
    #************ Plot literature curves
    if mineralname == 'NaAlSi3O8':
        colorstyles = {'Zhang1996':'0.5','Zhang1996s':'g','Tsuchiya2011':'g','Litvin1993':'g', 
                       'Bell1969s':'g', 'Bell1969':'0.5','Tutti2000':'0.5','Newton1967':'0.5'}
        linestyles =   {'Zhang1996':'-.','Zhang1996s':'-.','Tsuchiya2011':'-.','Litvin1993':'-', 
                        'Bell1969s':'-', 'Bell1969':'-','Tutti2000':'-','Newton1967':'-'} 
    elif mineralname == 'KAlSi3O8':
        colorstyles = {'Zhang1996':'0.5','Zhang1996s':'g','Tsuchiya2011':'g','Urakawa1994':'g', 
                       'Lindsley1966':'g','Akaogi2004':'0.5'}
        linestyles = {'Zhang1996':'-.','Zhang1996s':'-.','Tsuchiya2011':'-.','Urakawa1994':'--',
                      'Lindsley1966':'--','Akaogi2004':'--'}
    elif mineralname == 'CaAl2Si2O8':
        colorstyles = {'Zhang1996':'0.5','Zhang1996s':'g','Tsuchiya2011':'g',
                       'Hariya1968':'0.5', 'Hariya1968s':'g'}
        linestyles = {'Zhang1996':'-.','Zhang1996s':'-.','Tsuchiya2011':'-.',
                      'Hariya1968':':', 'Hariya1968s':':'}
    #ax2.plot(solidpoints_P, solidpoints_T, linestyle = None, linewidth = 0, marker = 'o', markersize = plot_parameters['size_markers']*0.1, color = 'k')
    for i in range(1,np.size(phases,0)):
        ax2.text(phases[i][0],phases[i][1],str(int(phases[i][idx])),fontsize=plot_parameters['size_font_ticks'])
    #** on both plots
    for ax in axis:
    #    ax.scatter(viscous_fluid[xvariable],viscous_fluid[yvariable],marker = 'x',
    #            s=plot_parameters['size_markers']*3, facecolor = 'grey',
    #            edgecolor='grey', label = 'viscous fluid',zorder = -90)
        for i in range(len(curves)):
            ref = curves[i][2]
            print(ref)
            ax.plot(curves[i][0],curves[i][1],linestyle = linestyles[ref], 
                    color = colorstyles[ref], zorder = -80,
                    marker = None, markersize = 0 , linewidth = plot_parameters['size_lines'])

    #** save the figure
    figurename = figurename+'.pdf'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












