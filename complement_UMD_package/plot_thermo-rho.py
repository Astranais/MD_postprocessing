#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3


        Plot any properties of the thermo file as a function of density, and fit 3rd order polynome on P-rho
        Produces the spinodal.txt file with the values of the fit
         *************  ARTICLE VERSION  *************
         
         
         For use after the fullaverages.py script !!!!!!!!!!
         
         For use after the fullaverages+Cvanalysis change lines 737-738 (search of every fullaverage filename based on name pattern)

"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D #useful to create a custom legend
import natsort
from matplotlib.ticker import AutoLocator,AutoMinorLocator
import crystallography as cr
import glob
import re
from scipy.interpolate import make_interp_spline, BSpline

def plot_deKoker2010(ax,plot_parameters, prop):
    """ Add  Cvm_Nkb values computed by deKoker2010 """
    try:
        if prop == 'Cvm_Nkb':
            rho, Cv, stdev = np.loadtxt('Cv_deKoker2010.txt',usecols = (0,1,2), skiprows = 1, unpack =True)
        else:
            rho, Cv, stdev = np.loadtxt('Cv_deKoker2010.txt',usecols = (0,3,4), skiprows = 1, unpack =True)
        ax.errorbar(rho,Cv, yerr=stdev, fmt='none',  ecolor = 'k', elinewidth = 1, capsize=2)
        ax.plot(rho,Cv, linestyle = '', marker='P', markeredgecolor = 'k', markeredgewidth = 0.5, markersize = plot_parameters["size_markers"]-1, markerfacecolor = 'k', label = 'de Koker (2010)')
        line1 = Line2D([0],[0],color = 'k', ls = '', marker = 'P',markeredgewidth = 0.5, markersize = plot_parameters["size_markers"]-1, markerfacecolor = 'k' )
        ax.legend([line1],['de Koker (2010)'])
    except FileNotFoundError:
        print('file Cv_deKoker2010.txt not found')
    
    
    
    
    
def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.outcar.umd.dat
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.strip('./')
    temperature = filename.split('_')[1]
    acell = filename.split('.outcar.umd.dat')[0].split('_')[3].strip('a')
    return temperature, acell

def spinodal(thermofile,addtxt):
    """creation of the newfile for spinodal values"""
    print(thermofile)
    newfilename = 'spinodal_'+addtxt+ thermofile.split('.txt')[0] +'.txt'
    print('The file ', newfilename, 'is also created')
    sp = open(newfilename, 'w')
    sp.write('T\trho(g/cm3)\tP(GPa)\tChi2_red\trho(gas_spinodal)\tP(gas_spinodal)\trho_crit\tP_crit\n')
    return sp

def coefficients(thermofile,addtxt):
    """creation of the newfile for coefficient values"""
    print(thermofile)
    newfilename = 'poly3_fit_'+addtxt+ thermofile.split('.txt')[0] +'.txt'
    print('The file ', newfilename, 'is also created')
    coeff = open(newfilename, 'w')
    coeff.write('T\ta\tb\tc\td\tstdev_a\tstdev_b\tstdev_c\tstdev_d\n')
    return coeff

def poly3_constrained(x,a,b,c):
    """ 3rd order polynomial """
    y = a * x**3 + b * x**2 + c * x
    return y

def poly3(x,a,b,c,d):
    """ 3rd order polynomial """
    y = a * x**3 + b * x**2 + c * x + d
    return y

def VdW(x,a,b,R,T):
    """ Van der Waals"""
    y = (R * T) / (1/x - b) - a * x**2
    return y


def chi2red(rho,data,stdev_data,popt,functiontofit):
    """ function chi2 which has to be minimized"""
    y = functiontofit(rho,*popt)
    chi2 = 0.0
    for i in range(len(rho)):
        chi2 = chi2 + ( ( data[i] - y[i]  ) / stdev_data[i] )**2        
    chi2red = chi2 / (len(rho) - len(popt))
    return chi2red

def calculation_min_max_constrained(popt,pcov):
    """ calculation of the x,y coordinates of the minimum of the 3rd order polynomial"""
    a = popt[0]
    b = popt[1]
    c = popt[2]
    xmin = (-2*b + np.sqrt(4*b**2 - 12*a*c)) /(6*a)
    ymin = a*xmin**3 + b*xmin**2 + c*xmin
    xmax = (-2*b - np.sqrt(4*b**2 - 12*a*c)) /(6*a)
    ymax = a*xmax**3 + b*xmax**2 + c*xmax
    xcrit = -b/(3*a)
    ycrit = a*xcrit**3 + b*xcrit**2 + c*xcrit
    return xmin, ymin, xmax, ymax, xcrit, ycrit

def calculation_min_max(popt,pcov):
    """ calculation of the x,y coordinates of the minimum and max of the 3rd order polynomial"""
    a = popt[0]
    b = popt[1]
    c = popt[2]
    d = popt[3] 
    xmin = (-2*b + np.sqrt(4*b**2 - 12*a*c)) /(6*a)
    ymin = a*xmin**3 + b*xmin**2 + c*xmin + d
    xmax = (-2*b - np.sqrt(4*b**2 - 12*a*c)) /(6*a)
    ymax = a*xmax**3 + b*xmax**2 + c*xmax + d
    xcrit = -b/(3*a)
    ycrit = a*xcrit**3 + b*xcrit**2 + c*xcrit + d
    return xmin, ymin, xmax, ymax, xcrit, ycrit
 
    
def convert_to_float(string):
    if string[0] == '>':
        value = float(string[1:])
    else:
        value = float(string[:])
    return value 


def creation_plot_3_vertical(prop, plot_parameters,major_xticks,minor_xticks,major_yticks,minor_yticks,xxmin, xxmax,yymin, yymax):
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
    global ax1, ax2, ax3
    f, (ax1, ax2, ax3) = plt.subplots(3,1, sharex=True, sharey=True, figsize = (size_figure[1],size_figure[0]*2.8))
    #Adjustment of ticks
    for ax in [ax1,ax2,ax3]:
        if prop == 'P':
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True) 
            ax.xaxis.set_ticks_position('both')
            ax.set_yticks(major_yticks)
            ax.set_yticks(minor_yticks, minor=True) 
            ax.yaxis.set_ticks_position('both') 
            ax.set_ylim(yymin,yymax)
        else:
            majorLocator = AutoLocator()
            minorLocator = AutoMinorLocator()
            ax.yaxis.set_major_locator(majorLocator)
            ax.yaxis.set_minor_locator(minorLocator) 
            ax.yaxis.set_ticks_position('both') 
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True)
            ax.xaxis.set_ticks_position('both')
            plt.autoscale(enable=True,axis='y',tight=False)
            #ax.grid(True, which='both',axis = 'x', )
        ax.set_xlim(xxmin,xxmax)
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)    
        if prop == 'testCv':
            ax.autoscale(enable=True,axis='both',tight=False)
    for ax in [ax1,ax2]:
        plt.setp(ax.get_yticklabels()[0], visible=False) 
        #ax.grid(True, axis = 'x')

    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    
    #labels
    ax.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    if prop == 'P':
        ax.set_ylabel(r"Pressure (GPa)", fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*3)
    if prop == 'stdev_P':
        ax.set_ylabel(r"Standard deviation of pressure (GPa)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'E':
        ax.set_ylabel(r"Internal Energy (eV)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'stdev_E':
        ax.set_ylabel(r"Standard deviation of internal energy (eV)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'T':
        ax.set_ylabel(r"Temperature (K)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'stdev_T' or prop == 'testCv':
        ax.set_ylabel(r"$\sigma^{2}_{T}/T^{2}$", fontweight = 'bold', fontsize = size_fonts,labelpad = shift_labelpad*5)
    if prop == 'Cv':
        ax.set_ylabel(r"Cv (J/K)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'Cvm':
        ax.set_ylabel(r"Cv (J/K/mol)", fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*3)    
    if prop == 'Cvm_Nkb':
        ax.set_ylabel(r"Cvm ($\mathcal{N}_{A}k_{B}$)", fontweight = 'bold', fontsize = size_fonts,labelpad = shift_labelpad*3)
    return f, ax1, ax2, ax3, ax


def creation_plot_3_horizontal(prop, plot_parameters,major_xticks,minor_xticks,major_yticks,minor_yticks,xxmin, xxmax,yymin, yymax):
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
    global ax1, ax2, ax3
    f, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=False, sharey=True, figsize = (size_figure[0]*2.8,size_figure[1]))
    #Adjustment of ticks
    for ax in [ax1,ax2,ax3]:
        if prop == 'P':
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True) 
            ax.xaxis.set_ticks_position('both')
            ax.set_yticks(major_yticks)
            ax.set_yticks(minor_yticks, minor=True) 
            ax.yaxis.set_ticks_position('both') 
            ax.set_ylim(yymin,yymax)
        else:
            majorLocator = AutoLocator()
            minorLocator = AutoMinorLocator()
            ax.yaxis.set_major_locator(majorLocator)
            ax.yaxis.set_minor_locator(minorLocator)   
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True) 
            plt.autoscale(enable=True,axis='both',tight=False)
            #ax.grid(True, which='both',axis = 'x', )
        ax.set_xlim(xxmin,xxmax)
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)    
        if prop == 'testCv':
            ax.autoscale(enable=True,axis='both',tight=False)
    for ax in [ax1,ax2]:
        plt.setp(ax.get_xticklabels()[-1], visible=False) 
        #ax.grid(True, axis = 'x')

    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    
    #labels
    ax.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    if prop == 'P':
        ax.set_ylabel(r"Pressure (GPa)", fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*3)
    if prop == 'stdev_P':
        ax.set_ylabel(r"Standard deviation of pressure (GPa)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'E':
        ax.set_ylabel(r"Internal Energy (eV)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'stdev_E':
        ax.set_ylabel(r"Standard deviation of internal energy (eV)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'T':
        ax.set_ylabel(r"Temperature (K)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'stdev_T' or prop == 'testCv':
        ax.set_ylabel(r"$\sigma^{2}_{T}/T^{2}$", fontweight = 'bold', fontsize = size_fonts,labelpad = shift_labelpad*5)
    if prop == 'Cv':
        ax.set_ylabel(r"Cv (J/K)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'Cvm':
        ax.set_ylabel(r"Cv (J/K/mol)", fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*3)    
    if prop == 'Cvm_Nkb':
        ax.set_ylabel(r"Cvm ($\mathcal{N}_{A}k_{B}$)", fontweight = 'bold', fontsize = size_fonts,labelpad = shift_labelpad*3)
    return f, ax1, ax2, ax3, ax

def creation_plot_2(prop, plot_parameters,major_xticks,minor_xticks,major_yticks,minor_yticks,xxmin, xxmax,yymin, yymax):
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
    global ax1, ax2
    f, (ax1,ax2) = plt.subplots(2,1, sharex=True, sharey=True, figsize = (size_figure[0],size_figure[1]*2))
    #Adjustment of ticks
    for ax in [ax1,ax2]: 
        if prop == 'P':
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True) 
            ax.xaxis.set_ticks_position('both')
            ax.set_yticks(major_yticks)
            ax.set_yticks(minor_yticks, minor=True) 
            ax.yaxis.set_ticks_position('both') 
            ax.set_ylim(yymin,yymax)
        else:
            majorLocator = AutoLocator()
            minorLocator = AutoMinorLocator()
            ax.yaxis.set_major_locator(majorLocator)
            ax.yaxis.set_minor_locator(minorLocator)   
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True) 
            plt.autoscale(enable=True,axis='y',tight=False)
            #ax.grid(True, which='both',axis = 'x', )
        ax.set_xlim(xxmin,xxmax)
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)    
    for ax in [ax2]:
        plt.setp(ax.get_yticklabels()[-1], visible=False) 
        #ax.grid(True, axis = 'x')

    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    
    #labels
    ax.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    if prop == 'P':
        ax.set_ylabel(r"Pressure (GPa)", fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*4)
    if prop == 'stdev_P':
        ax.set_ylabel(r"Standard deviation of pressure (GPa)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'E':
        ax.set_ylabel(r"Internal Energy (eV)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'stdev_E':
        ax.set_ylabel(r"Standard deviation of internal energy (eV)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'T':
        ax.set_ylabel(r"Temperature (K)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'stdev_T' or prop == 'testCv':
        ax.set_ylabel(r"$\sigma^{2}_{T}/T^{2}$", fontweight = 'bold', fontsize = size_fonts,labelpad = shift_labelpad)
    if prop == 'Cv':
        ax.set_ylabel(r"Cv (J/K)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'Cvm':
        ax.set_ylabel(r"Cv (J/K/mol)", fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)    
    if prop == 'Cvm_Nkb':
        ax.set_ylabel(r"Cvm ($\mathcal{N}_{A}k_{B}$)", fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    return f, ax1, ax2, ax


def creation_plot(prop, plot_parameters,major_xticks,minor_xticks,major_yticks,minor_yticks,xxmin, xxmax,yymin, yymax):
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
    global ax1
    fig, ax1 = plt.subplots(1,1, figsize = size_figure, linewidth = size_lines /2)  
    #Adjustment of ticks
    if prop == 'P':
        ax1.set_xticks(major_xticks)
        ax1.set_xticks(minor_xticks, minor=True) 
        ax1.xaxis.set_ticks_position('both')
        ax1.set_yticks(major_yticks)
        ax1.set_yticks(minor_yticks, minor=True) 
        ax1.yaxis.set_ticks_position('both') 
        ax1.set_ylim(yymin,yymax)
    else:
        majorLocator = AutoLocator()
        minorLocator = AutoMinorLocator()
        ax1.yaxis.set_major_locator(majorLocator)
        ax1.yaxis.set_minor_locator(minorLocator)   
        ax1.set_xticks(major_xticks)
        ax1.set_xticks(minor_xticks, minor=True) 
        plt.autoscale(enable=True,axis='both',tight=False)
        #ax.grid(True, which='both',axis = 'x', )
    ax1.set_xlim(xxmin,xxmax)
    ax1.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
    #plt.autoscale()
    #labels
    ax1.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    if prop == 'P':
        ax1.set_ylabel(r"Pressure (GPa)", fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    if prop == 'stdev_P':
        ax1.set_ylabel(r"Standard deviation of pressure (GPa)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'E':
        ax1.set_ylabel(r"Internal Energy (eV/atom)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'stdev_E':
        ax1.set_ylabel(r"Standard deviation of internal energy (eV)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'T':
        ax1.set_ylabel(r"Temperature (K)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'stdev_T' or prop == 'testCv':
        ax1.set_ylabel(r"$\sigma^{2}_{T}/T^{2}$", fontweight = 'bold', fontsize = size_fonts,labelpad = shift_labelpad)
    if prop == 'Cv':
        ax1.set_ylabel(r"Cv (J/K)", fontweight = 'bold', fontsize = size_fonts)
    if prop == 'Cvm':
        ax1.set_ylabel(r"Cv (J/K/mol)", fontweight = 'bold', fontsize = size_fonts,labelpad = shift_labelpad)    
    if prop == 'Cvm_Nkb':
        ax1.set_ylabel(r"Cvm ($\mathcal{N}_{A}k_{B}$)", fontweight = 'bold', fontsize = size_fonts)
    return fig, ax1

def initialization(files,file,thermofile,thermofile2,prop,compounds, view,addtxt):
    #********initialisations
    sp = 0
    coeff = 0
    #**extraction compound
    compound = file.split('_')[1]
    compounds.append(compound)
    #**change of subplot
    if (thermofile == 'all'): 
        if files.index(file) == 0:
            ax=ax1
            if view == 'classic':
                letter = 'a'
            else:
                letter = 'd'
        elif files.index(file) == 1:
            ax=ax2
            if view == 'classic':
                letter = 'b'
            else:
                letter ='e'
        elif files.index(file) == 2:
            ax=ax3 
            if view == 'classic':
                letter = 'c'
            else:
                letter = 'f'
        print("I plot on axis",ax)
    elif (thermofile2 != '') :
        if files.index(file) == 0:
            ax=ax1
            if view == 'classic':
                letter = 'a'
            else:
                letter = 'c'
        elif files.index(file) == 1:
            ax=ax2
            if view == 'classic':
                letter = 'b'
            else:
                letter ='d'
        print("I plot on axis",ax)
    else:
        ax = ax1
        letter = compound
    #creation of the spinodal txt file if P selected
    if prop == 'P':
        sp = spinodal(file,addtxt)
        coeff = coefficients(file,addtxt)
    #**creation of elements and number lists and initialization of T
    with open(file,'r') as f:
        [f.readline() for i in range(3)]
        line = f.readline()
        entry=line.split()
        temperature0, acell0 = split_name(entry[0])    
    return ax, letter, sp, coeff, compounds, temperature0

def extractandplot(file, temperature0, Temperatures, prop,ax,plot_parameters,colors_T,colors_Tfill, colors_spfill, sp, coeff,column_number,thermofile,thermofile2,letter,xxmin, xxmax, yymin, yymax,markertype,linetype,constrained,den_lim):
    """ we extract and plot data at the same time"""
    #**************initialisation of data dictionnaries
    rho = {}
    Y = {} #y axis data
    stdev = {}
    err = {}
    data = {}  
    #************* define append and plot functions depending on the prop to plot
    if (prop == "P") or (prop == "T") or (prop == "E"):
        def append_function(data,temperature,entry, column_number, prop):
            data[temperature].append([  float(entry[column_number['rho']]) , float(entry[column_number[prop]]) , float(entry[column_number[prop]+1]) , convert_to_float(entry[column_number[prop]+2])   ])
        def append_function2(data, temperature, ii, prop, Y, rho, stdev, err):
            rho[temperature].append(data[temperature][ii][0] )
            Y[temperature].append(data[temperature][ii][1]  )
            stdev[temperature].append(data[temperature][ii][2] )
            err[temperature].append(data[temperature][ii][3] )
        if prop == "P" :
            def plot_function(ax,sp,coeff,temperature, rho, Y, err, stdev, colors_T, markertype, plot_parameters, colors_Tfill,constrained):
                #transform to arrays
                rho_data = np.asarray(rho[temperature])
                Y_data = np.asarray(Y[temperature])
                stdev_data = np.asarray(stdev[temperature])
                err_data = np.asarray(err[temperature])
                #fit and plot
                ax.errorbar(rho_data,Y_data, yerr=2*err_data, fmt='none',  ecolor = '0.5', elinewidth = 1, capsize=2)
                try:
                    ax.plot(rho_data,Y_data, linestyle = '', marker=markertype[0], markeredgecolor = colors_T[temperature], markeredgewidth = 0.5, markersize = plot_parameters["size_markers"], markerfacecolor = colors_Tfill[temperature], label = temperature)
                except KeyError:
                    pass
                if len(Y[temperature]) > 4: 
                    if constrained == 0:
                        functiontofit = poly3
                        functionmin = calculation_min_max
                    else: 
                        functiontofit = poly3_constrained
                        functionmin = calculation_min_max_constrained
                    #if there is enough data at this T then we fit a 3rd order polynome to the data 
                    popt, pcov = curve_fit(functiontofit, rho_data,Y_data, p0=None, sigma=stdev_data, absolute_sigma=False )
                    if constrained == 1:
                        d = 0
                        stdev_d = 0
                    else:
                        d = popt[3]
                        stdev_d = np.sqrt(np.diag(pcov))[3]
                    #and then we find the mininum of the corresponding curves and associated stdev       
                    xmin, ymin, xmax, ymax, xcrit, ycrit = functionmin(popt,pcov)
                    #plot of the fit and spinodal
                    rho_extended = np.arange(0.2,max(rho_data)+0.05, 0.05)
                    P_fit = functiontofit(rho_extended, *popt)
                    ax.plot(rho_extended, P_fit, linestyle = linetype, color = colors_T[temperature], linewidth = plot_parameters["size_lines"])
                    ax.plot(xmin,ymin, marker = markertype[1],  markersize = plot_parameters["size_markers"]+3, markeredgecolor = colors_T[temperature], markeredgewidth = 0.5, markerfacecolor = colors_spfill)
                    sp.write(str(temperature)+'\t'+str(round(xmin,2))+'\t'+str(round(ymin,2))+'\t'+str(round(chi2red(rho_data,Y_data,stdev_data,popt,functiontofit),3))+'\t'+str(round(xmax,2))+'\t'+str(round(ymax,2))+'\t'+str(round(xcrit,2))+'\t'+str(round(ycrit,2))+'\n')
                    coeff.write(str(temperature)+'\t'+str(round(popt[0],3))+'\t'+str(round(popt[1],3))+'\t'+str(round(popt[2],3))+'\t'+str(round(d,3))+'\t'+str(round(np.sqrt(np.diag(pcov))[0],3))+'\t'+str(round(np.sqrt(np.diag(pcov))[1],3))+'\t'+str(round(np.sqrt(np.diag(pcov))[2],3))+'\t'+str(round(stdev_d,3))+'\n')
        else:
            def plot_function(ax,sp,coeff,temperature, rho, Y, err, stdev, colors_T, markertype, plot_parameters, colors_Tfill,constrained ):
                #transform to arrays
                err_data = np.asarray(err[temperature])
                #plot
                ax.errorbar(rho[temperature],Y[temperature], yerr=2*err_data, fmt='none',  ecolor = '0.5', elinewidth = 1, capsize=2)
                try:
                    ax.plot(rho[temperature],Y[temperature], linestyle = '', marker=markertype[0], markeredgecolor = colors_T[temperature], markeredgewidth = 0.5, markersize = plot_parameters["size_markers"], markerfacecolor = colors_Tfill[temperature], label = temperature)
                except KeyError:
                    pass
    elif (prop == "stdev_P") or (prop == "stdev_E"):
        def append_function(data,temperature,entry, column_number, prop):
            data[temperature].append([  float(entry[column_number['rho']]) , float(entry[column_number[prop]]) , 0 , 0  ])
        def append_function2(data, temperature, ii, prop, Y, rho, stdev, err):
            rho[temperature].append(data[temperature][ii][0] )
            Y[temperature].append(data[temperature][ii][1]  )
        def plot_function(ax,sp,coeff,temperature, rho, Y, err, stdev, colors_T, markertype, plot_parameters, colors_Tfill,constrained ):
            try:
                ax.plot(rho[temperature],Y[temperature], linestyle = '', marker=markertype[0], markeredgecolor = colors_T[temperature], markeredgewidth = 0.5, markersize = plot_parameters["size_markers"], markerfacecolor = colors_Tfill[temperature], label = temperature)
            except KeyError:
                pass
    elif (prop == "stdev_T"): #to check if sigmaT **2 / T**2 = 2/3N  (if no then fluctuations are not ok, and then Cv too)
        def append_function(data,temperature,entry, column_number, prop):
            data[temperature].append([  float(entry[column_number['rho']]) , float(entry[column_number[prop]]) , float(entry[column_number['T']]) , 0  ])
        def append_function2(data, temperature, ii, prop, Y, rho, stdev, err):
            rho[temperature].append(data[temperature][ii][0] )
            Y[temperature].append(data[temperature][ii][1]**2/data[temperature][ii][2]**2 ) 
        def plot_function(ax,sp,coeff,temperature, rho, Y, err, stdev, colors_T, markertype, plot_parameters, colors_Tfill,constrained ):
            try:
                ax.plot(rho[temperature],Y[temperature], linestyle = '', marker=markertype[0], markeredgecolor = colors_T[temperature], markeredgewidth = 0.5, markersize = plot_parameters["size_markers"], markerfacecolor = colors_Tfill[temperature], label = temperature)
                ax.axhline(y=2/(3*208), color = 'k', linewidth = 2)
            except KeyError:
                pass
    elif (prop == "testCv"): #to check if sigmaT **2 / T**2 = 2/3N  (if no then fluctuations are not ok, and then Cv too)
        def append_function(data,temperature,entry, column_number, prop):
            data[temperature].append([  float(entry[column_number['rho']]) , float(entry[column_number[prop]]) , 0 , convert_to_float(entry[column_number[prop]+1])  ])
        def append_function2(data, temperature, ii, prop, Y, rho, stdev, err):
            rho[temperature].append(data[temperature][ii][0] )
            Y[temperature].append(data[temperature][ii][1])
            err[temperature].append(data[temperature][ii][3] )
        def plot_function(ax,sp,coeff,temperature, rho, Y, err, stdev, colors_T, markertype, plot_parameters, colors_Tfill,constrained ):
            err_data = np.asarray(err[temperature])
            try:
                #ax.errorbar(rho[temperature],Y[temperature], yerr=2*err_data, fmt='none',  ecolor = '0.5', elinewidth = 1, capsize=2)
                #ax.plot(rho[temperature],Y[temperature], linestyle = '', marker=markertype[0], markeredgecolor = colors_T[temperature], markeredgewidth = 0.5, markersize = plot_parameters["size_markers"], markerfacecolor = colors_Tfill[temperature], label = temperature)
                val = 2/(3*208)
                ax.axhline(y=val, color = 'k', linewidth = 2)
                transp = '19' #10% #'55'#25% #'7f' #50%
                count = 0
                totcount = 0
                for i in range(len(rho[temperature])):
                    inflim = Y[temperature][i] - 2*err_data[i]
                    suplim = Y[temperature][i] + 2*err_data[i]
                    totcount +=1
                    if val < suplim and val > inflim:
                        ax.errorbar(rho[temperature][i],Y[temperature][i], yerr=2*err_data[i], fmt='none',  ecolor = colors_T[temperature], elinewidth = 1, capsize=2)
                        ax.plot(rho[temperature][i],Y[temperature][i], linestyle = '', marker=markertype[0], markeredgecolor = colors_T[temperature], markeredgewidth = 0.5, markersize = plot_parameters["size_markers"], markerfacecolor = colors_Tfill[temperature], label = temperature)
                        count +=1
                    else:
                        ax.errorbar(rho[temperature][i],Y[temperature][i], yerr=2*err_data[i], fmt='none',  ecolor = colors_T[temperature]+transp, elinewidth = 1, capsize=2)
                        ax.plot(rho[temperature][i],Y[temperature][i], linestyle = '', marker=markertype[0], markeredgecolor = colors_T[temperature]+transp, markeredgewidth = 0.5, markersize = plot_parameters["size_markers"], markerfacecolor = colors_Tfill[temperature]+transp, label = temperature)
                print(count,totcount)
            except KeyError:
                pass
    else: #for Cv, stdev estimated using the 4 test simulations of NaAlSi3O8 at T5 and 16.5 A and err using bootstrap method
        if prop == "Cv":
            stdev_data =  1.029*10**(-21) 
        elif prop == "Cvm":
            stdev_data = 2.98 
        else: #"Cvm_Nkb":
            stdev_data = 0.36
        def append_function(data,temperature,entry, column_number, prop):
            if prop == 'Cvm':
                data[temperature].append([  float(entry[column_number['rho']]) , float(entry[column_number[prop]]) , (float(entry[column_number['testCv']]),float(entry[column_number['stdev_testCv']])) , convert_to_float(entry[column_number[prop]+1]) ])
            else:
                #data[temperature].append([  float(entry[column_number['rho']]) , float(entry[column_number[prop]]) , stdev_data , convert_to_float(entry[column_number[prop]+1])  ])
                data[temperature].append([  float(entry[column_number['rho']]) , float(entry[column_number[prop]]), (float(entry[column_number['testCv']]),float(entry[column_number['stdev_testCv']])) , convert_to_float(entry[column_number[prop]+1])  ])
        def append_function2(data, temperature, ii, prop, Y, rho, stdev, err):
            rho[temperature].append(data[temperature][ii][0] )
            Y[temperature].append(data[temperature][ii][1]  )
            #stdev[temperature].append(data[temperature][ii][2] )
            stdev[temperature].append((data[temperature][ii][2][0],data[temperature][ii][2][1])) #testCv and err
            err[temperature].append(data[temperature][ii][3] )
        def plot_function(ax,sp,coeff,temperature, rho, Y, err, stdev, colors_T, markertype, plot_parameters, colors_Tfill,constrained ):
            #transform to arrays
            err_data = np.asarray(err[temperature])
            try:
                #plot
                #ax.errorbar(rho[temperature],Y[temperature], yerr=2*err_data, fmt='none',  ecolor = '0.5', elinewidth = 1, capsize=2)
                #ax.plot(rho[temperature],Y[temperature], linestyle = '-', marker=markertype[0], markeredgecolor = colors_T[temperature], markeredgewidth = 0.5, markersize = plot_parameters["size_markers"], markerfacecolor = colors_Tfill[temperature], color =  colors_T[temperature], label = temperature)
                transp = '19'
                val = 2/(3*208)
                for i in range(len(rho[temperature])):
                    inflim = stdev[temperature][i][0] - 2*stdev[temperature][i][1]
                    suplim = stdev[temperature][i][0] + 2*stdev[temperature][i][1]
                    if val < suplim and val > inflim:                        
                        ax.errorbar(rho[temperature][i],Y[temperature][i], yerr=2*err_data[i], fmt='none',  ecolor = colors_T[temperature], elinewidth = 1, capsize=2)
                        ax.plot(rho[temperature][i],Y[temperature][i], linestyle = '', marker=markertype[0], markeredgecolor = colors_T[temperature], markeredgewidth = 0.5, markersize = plot_parameters["size_markers"], markerfacecolor = colors_Tfill[temperature], label = temperature)
                    else:
                        ax.errorbar(rho[temperature][i],Y[temperature][i], yerr=2*err_data[i], fmt='none',  ecolor = colors_T[temperature]+transp, elinewidth = 1, capsize=2)
                        ax.plot(rho[temperature][i],Y[temperature][i], linestyle = '', marker=markertype[0], markeredgecolor = colors_T[temperature]+transp, markeredgewidth = 0.5, markersize = plot_parameters["size_markers"], markerfacecolor = colors_Tfill[temperature]+transp, label = temperature)    
            except KeyError:
                pass
    #********* fill the dictionnaries with data
    with open(file,'r') as f:
        [f.readline() for i in range(3)]
        #extract data
        while True:
            line = f.readline()
            if not line: break
            else:    
                entry=line.split('\n')[0].split('\t')
                temperature, acell = split_name(entry[0])
                try:
                    append_function(data,temperature,entry, column_number, prop)
                except KeyError:                      
                    data[temperature] = []
                    append_function(data,temperature,entry, column_number, prop)
    #********* remove the data we do not want to use
    for temperature in data:
        if temperature == 'T10' or temperature == 'T15' or temperature == 'T20':
            continue
        else:
            #first we sort the data by density from smallest to biggest
            data[temperature] = sorted(data[temperature], key=lambda x: x[0], reverse=False) #sort by the first column: rho
            #second we create rho,P,stdev_P dictionnaries
            for ii in range(len(data[temperature])):
                if data[temperature][ii][0] > den_lim   :
                    #if the density is too high then we stop the loop to not take the others
                    break
                else:
                    try:
                        append_function2(data, temperature, ii, prop, Y, rho, stdev, err)
                    except KeyError:
                        rho[temperature] = []
                        Y[temperature] = []
                        stdev[temperature] = [] 
                        err[temperature] = []
                        append_function2(data, temperature, ii, prop, Y, rho, stdev, err)
    #******** Fit data and plot
    for temperature in sorted(rho):
        plot_function(ax,sp,coeff,temperature, rho, Y, err, stdev, colors_T, markertype, plot_parameters, colors_Tfill,constrained )
    #******* addition text on figure    
    if (thermofile == 'all'):
        ax.text(0.015,0.945, letter , transform=ax.transAxes, horizontalalignment = 'left',verticalalignment = 'baseline',  fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='none', edgecolor='k', pad=3.0))    
    elif (thermofile2 != ''):
        if prop == 'Cvm_Nkb':
            position = (2.4,8)
        elif prop == 'P':
            position = (xxmin+0.047*(xxmax-xxmin),yymax-0.05*(yymax-yymin))
        ax.text(position[0],position[1], letter, horizontalalignment = 'right', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
    else:
        yymin,yymax = ax.get_ylim()
        position = (xxmin+0.3*(xxmax-xxmin),yymax-0.05*(yymax-yymin))
        ax.text(position[0],position[1], letter,horizontalalignment = 'right', fontweight = 'bold', fontsize = plot_parameters["size_fonts"])            

    
def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),"size_markers" : 8,"size_lines" : 1,"shift_labelpad" : 10}
    #other dictionnaries and parameters for the figure
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200','T5':'#ffcd01','T5.5':'#ff6e00','T6':'#ff0101','T6.5':'#ff00a2','T7':'#ff01de','T7.5':'#ffa6f4'}#,'T10':'#ffe86e','T15':'#ffbf90','T20':'#ff7788'}
    thermofile = 'all'
    thermofile2 = ''
    filetype = 'all'
    constrained = 0
    doubleplot = 0
    #other dictionnaries and parameters for the figure
    Temperatures = []
    compounds = []
    den_lim = 2.5
    try:
        options,arg = getopt.getopt(argv,"hf:g:p:t:d:c:l:",["filename","gfilename","property","type","double",'constrained','limit'])
    except getopt.GetoptError:
        print("plot_thermo-rho.py -f <thermo_filename (default = 'all')>  -g <thermo_filename 2 (if we want to plot 2 files)> -p <property to plot  (P,stdev_P,T,E,stdev_E,Cv,Cvm,Cvm_Nkb)> -t <type of the file: 'short' or 'all'> -d < =1 if we want to plot soft pseudopot values on classic plot. =0 by default>  -c <1=constrained poly3 fit, default = 0> -l <density limit, default = 2.5>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_thermo-rho.py program to plot any properties of the thermo file as a function of density, and fit 3rd order polynome on P-rho. If -s option is filled, then the stdev on averages will be used instead of arithmetic stdev')
            print("plot_thermo-rho.py -f <thermo_filename (default = 'all' -> 3plot)> -g <thermo_filename 2 (if we want to plot 2 files)>  -p <property to plot (P,stdev_P,T,E,stdev_E,Cv,Cvm,Cvm_Nkb)> -t <type of the file: 'short' or 'all'> -d < =1 if we want to plot soft pseudopot values on classic plot. =0 by default> -c <1=constrained poly3 fit, default = 0> -l <density limit, default = 2.5>")
            print('')
            sys.exit()
        elif opt in ("-p","--property"):
            prop = str(arg)
        elif opt in ('-f','--filename'):
            thermofile = str(arg)
        elif opt in ('-g','--gfilename'):
            thermofile2 = str(arg)
        elif opt in ('-t','--type'):
            filetype = str(arg)
        elif opt in ('-d','--double'):
            doubleplot = int(arg)
        elif opt in ('-c','--constrained'):
            constrained = int(arg)  
        elif opt in ('-l','--limit'):
            den_lim = float(arg)
    if constrained == 1:
        addtxt = 'constrained_'
    else:
        addtxt = ''
    if den_lim <2.5: #for zoom view
        linetype1 = ''  
        major_xticks = np.arange(0, den_lim+0.1, 0.5) #zoom on low density 
        minor_xticks = np.arange(0, den_lim+0.1, 0.1)#zoom on low density
        major_yticks = np.arange(-1, 1.6, 0.5) #zoom on low density 
        minor_yticks = np.arange(-1, 1.6, 0.1)#zoom on low density
        (yymin, yymax) = (-1,1.0) #zoom on low density 
        (xxmin, xxmax) = (0,den_lim) #zoom on low density 
        view = 'zoom'
    else: 
        linetype1 = (0, (5, 8))  #dashed line with 5pt line and 8 pt space
        major_xticks = np.arange(0, den_lim+0.5, 0.5) #zoom on low density 
        minor_xticks = np.arange(0, den_lim+0.1, 0.1)#zoom on low density
        major_yticks = np.arange(-10, 100, 2.5) #classic view
        minor_yticks = np.arange(-10, 100, 0.5)#classic view
        (yymin, yymax) = (-3,5) #classic view
        (xxmin, xxmax) = (1,den_lim) #classic view
        view = 'classic'
    #************ initialization of the plot
    if thermofile == 'all':
        fig,ax1,ax2,ax3, ax0 = creation_plot_3_vertical(prop, plot_parameters,major_xticks,minor_xticks,major_yticks,minor_yticks,xxmin, xxmax,yymin, yymax)
        #print("axis are ax1",ax1,"ax2",ax2,"ax3",ax3,"ax0",ax0)
        #files = sorted(glob.glob('fullthermo+Cvanalysis_*O8_'+filetype+'_clean.txt'),reverse=True) #I list every thermo files
        files = sorted(glob.glob('thermo_*O8_'+filetype+'.txt'),reverse=True) #I list every thermo files
        if doubleplot == 1:
            softfiles = sorted(glob.glob('thermo_*O8_soft_'+filetype+'.txt'),reverse=True) #I list every thermo files
        figurename = 'thermo_all_'+addtxt+prop
    elif thermofile2 != '':
        fig,ax1,ax2, ax0  = creation_plot_2(prop, plot_parameters,major_xticks,minor_xticks,major_yticks,minor_yticks,xxmin, xxmax,yymin, yymax)
        #print("axis are ax1",ax1,"ax2",ax2,"ax0",ax0)
        files = [thermofile, thermofile2] #we will plot 2 files
        if doubleplot == 1:
            softfiles = [thermofile.split(filetype)[0]+'soft_'+filetype+'.txt',thermofile2.split(filetype)[0]+'soft_'+filetype+'.txt']
        figurename = thermofile.split('.txt')[0]+'_'+thermofile2.split('.txt')[0].split('_')[1]+'_'+addtxt+prop
    else:
        fig, ax0 = creation_plot(prop, plot_parameters,major_xticks,minor_xticks,major_yticks,minor_yticks,xxmin, xxmax,yymin, yymax)
        files = [thermofile] #I take only the file we want
        if doubleplot == 1:
            softfiles = [thermofile.split(filetype)[0]+'soft_'+filetype+'.txt']
        figurename = thermofile.split('.txt')[0]+'_'+addtxt+prop
    #************ initialization of the column number
    if filetype == 'all':
        column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,'Cvm':25,'stdev_Cvm':26,'testCv':31,'stdev_testCv':32}
    else:
        column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,'E':11,'stdev_E':12,'err_E':13,'Cvm_Nkb':14,'stdev_Cvm_Nkb':15,'testCv':16,'stdev_testCv':17}
    #************ creation of arrays and plot at the same time
    for file in files:
        #print("I use file ", file)
        #initialization of the subplot letter and temperature
        ax, letter, sp,coeff, compounds, temperature0 = initialization(files,file,thermofile,thermofile2,prop,compounds,view,addtxt)
        markertype = ['o','P']        
        colors_Tfill = colors_T
        colors_spfill = 'k'
        extractandplot(file, temperature0, Temperatures, prop,ax,plot_parameters,colors_T,colors_Tfill,colors_spfill, sp,coeff,column_number,thermofile,thermofile2,letter,xxmin, xxmax, yymin, yymax,markertype,linetype1,constrained,den_lim)
        #** plot comparison data
        if prop == 'Cvm' or prop == 'Cvm_Nkb':
            if compounds[-1] == 'CaAl2Si2O8':
                print('plot deKoker2010')
                plot_deKoker2010(ax,plot_parameters,prop)
    if doubleplot == 1:
        for file in softfiles:
            #print("I use file ", file)
            #initialization of the subplot letter and temperature
            ax, letter, sp,coeff, compounds, temperature0 = initialization(softfiles,file,thermofile,thermofile2,prop,compounds,view,addtxt)
            markertype = ['o','P']       
            linetype2 = ':'
            colors_Tfill = {'T2':'w','T3':'w','T4':'w','T4.5':'w','T5':'w','T5.5':'w','T6':'w','T6.5':'w','T7':'w','T7.5':'w'}
            colors_spfill = '0.90'
            extractandplot(file, temperature0, Temperatures, prop,ax,plot_parameters,colors_T,colors_Tfill,colors_spfill, sp,coeff,column_number,thermofile,thermofile2,letter,xxmin, xxmax, yymin, yymax,markertype,linetype2,constrained,den_lim)           
    
    #************ legend
    custom_patch = [mpatches.Patch(color=colors_T[key]) for key in natsort.natsorted(colors_T,reverse=True)]
    #legend = ax0.legend([col for col in custom_patch],[str(int(float(label.strip('T'))*1000)) for label in natsort.natsorted(colors_T)],title = '$\\bf{Temperature~(K)}$', bbox_to_anchor=(0.5, 1.03), loc="lower center", fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol=10)
    #legend = ax0.legend([col for col in custom_patch],[str(int(float(label.strip('T'))*1000)) for label in natsort.natsorted(colors_T,reverse=True)],title = '$\\bf{Temperature~(K)}$', bbox_to_anchor=(1.01, 1), loc="upper left", fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol=1)
    #plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])
    
    
    figurename = figurename+'_'+view + '.pdf'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
#    plt.show() # à mettre après savefig sinon on ne sauvegarde rien !!  

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



