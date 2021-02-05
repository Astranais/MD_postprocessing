#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3


        Fit and plot BM3 from the thermo file as a function of density
        Produces the BM3.txt file with the values of the fit
         *************  ARTICLE VERSION  *************
         
         
         For use after the fullaverages.py script !!!!!!!!!!

"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import fmin as simplex
import scipy.odr as odr
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import make_interp_spline, BSpline
import crystallography as cr


def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.outcar.umd.dat
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.strip('./')
    temperature = filename.split('_')[1]
    acell = filename.split('.outcar.umd.dat')[0].split('_')[3].strip('a')
    return temperature, acell

def fileBM3(thermofile,denlim):
    """creation of the newfile for BM3 values"""
    newfilename = 'BM3_'+str(denlim)+'_'+ thermofile.split('.txt')[0] +'.txt'
    bm = open(newfilename, 'w')
    bm.write('BM3 fit for pressures above 1 GPa or '+str(denlim)+' g/cm3\n')
    bm.write('T\tChi2_red\trho0(g/cm3)\tK0(GPa)\tKp0\tstdev_rho0(g/cm3)\tstdev_K0(GPa)\tstdev_Kp0\n')
    print('File',newfilename,'created')
    return bm

def poly3(x,a,b,c,d):
    """ 3rd order polynomial """
    xarr = np.asarray(x)
    yarr = a * xarr**3 + b * xarr**2 + c * xarr + d
    y = yarr.tolist()
    return y

def BM3_P_V(parameters,V):
    """definition of the function BM3 : 3rd order Birch-Murnaghan  P=f(V)"""
    V0,K0,Kp0 = parameters
    BM3=3.*K0/2.*((V0/V)**(7./3.)-(V0/V)**(5./3.))*(1.+3./4.*(Kp0-4.)*((V0/V)**(2./3.)-1.))
    return BM3

def BM3_P_rho(parameters,rho):
    """definition of the function BM3 : 3rd order Birch-Murnaghan  P=f(V)"""
    rho0,K0,Kp0 = parameters
    BM3=3.*K0/2.*((rho/rho0)**(7./3.)-(rho/rho0)**(5./3.))*(1.+3./4.*(Kp0-4.)*((rho/rho0)**(2./3.)-1.))
    return BM3

def fonction_BM3(rho,rho0,K0,Kp0):
    """definition of the function BM3 : 3rd order Birch-Murnaghan  P=f(V)"""
    BM3=3.*K0/2.*((rho/rho0)**(7./3.)-(rho/rho0)**(5./3.))*(1.+3./4.*(Kp0-4.)*((rho/rho0)**(2./3.)-1.))
    return BM3

def chi2red(popt,rho,data,stdev_data):
    """ function chi2 which has to be minimized"""
    y = fonction_BM3(rho,*popt)
    # compute chi-square
    chi2 = 0.0
    for ii in range(len(rho)):
        chi2 = chi2 + ( (data[ii] - y[ii]) / stdev_data[ii]  )**2
    chi2red = chi2 / (len(rho) - len(popt))
    return chi2red

def fonction(params,X,P,Err):
    """ function chi2 which has to be minimized"""
    # extract current values of fit parameters from input array
    rho0 = params[0]
    K0 = params[1]
    Kp0 = params[2]
    # compute chi-square
    chi2 = 0.0
    for n in range(len(X)):
        rho = X[n]
        gmP=3.*K0/2.*((rho/rho0)**(7./3.)-(rho/rho0)**(5./3.))*(1.+3./4.*(Kp0-4.)*((rho/rho0)**(2./3.)-1.))
        
        chi2 = chi2 + ( (P[n] - gmP) / Err[n]  )**2
    return chi2
 
    
def convert_to_float(string):
    if string[0] == '>':
        value = float(string[1:])
    else:
        value = float(string[:])
    return value 


def creation_plot_3(plot_parameters):
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
    f, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, sharey=True, 
       figsize = (size_figure[0]*2.5,size_figure[1]))
    #Adjustment of ticks
    major_xticks = np.arange(0, 2.6, 0.5) 
    minor_xticks = np.arange(0, 2.6, 0.1)
    for ax in [ax1,ax2,ax3]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True) 
        ax.xaxis.set_ticks_position('both')
        #ax.set_xlim(xxmin,xxmax)
        major_yticks = np.arange(-10, 100, 5) 
        minor_yticks = np.arange(-10, 100, 1)
        ax.set_yticks(major_yticks)
        ax.set_yticks(minor_yticks, minor=True) 
        ax.yaxis.set_ticks_position('both') 
        #ax.set_ylim(yymin,yymax)
        #ax.grid(True, which='both',axis = 'x', )
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)    
    for ax in [ax1,ax2]:
        plt.setp(ax.get_xticklabels()[-1], visible=False) 
        #ax.grid(True, axis = 'x')

    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, 
                      hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False,
                    labelleft=False, left=False, labelright = False, right=False)
    
    #labels
    ax.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad*2)
    ax.set_ylabel(r"Pressure (GPa)", fontweight = 'bold', fontsize = size_fonts, 
                  labelpad = shift_labelpad*2)
    return f, ax1, ax2, ax3, ax



def creation_plot_2(plot_parameters):
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
    f, (ax1,ax2) = plt.subplots(2,1, sharex=True, sharey=True, 
       figsize = (size_figure[0],size_figure[1]*2))
    #Adjustment of ticks
    major_xticks = np.arange(0, 2.6, 0.5) 
    minor_xticks = np.arange(0, 2.6, 0.1)
    for ax in [ax1,ax2]: 
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True) 
        ax.xaxis.set_ticks_position('both')
        #ax.set_yticks(major_yticks)
        #ax.set_yticks(minor_yticks, minor=True) 
        ax.yaxis.set_ticks_position('both') 
        #ax.set_ylim(yymin,yymax)
        #ax.set_xlim(xxmin,xxmax)
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)    
    for ax in [ax2]:
        plt.setp(ax.get_yticklabels()[-1], visible=False) 
        #ax.grid(True, axis = 'x')

    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, 
                      hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    
    #labels
    ax.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad*2)
    ax.set_ylabel(r"Pressure (GPa)", fontweight = 'bold', fontsize = size_fonts, 
                  labelpad = shift_labelpad*4)
    return f, ax1, ax2, ax


def creation_plot(plot_parameters):
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
    major_xticks = np.arange(0, 4.5, 0.5) 
    minor_xticks = np.arange(0, 4.5, 0.1)
    ax1.set_xticks(major_xticks)
    ax1.set_xticks(minor_xticks, minor=True) 
    ax1.xaxis.set_ticks_position('both')
    #ax1.set_yticks(major_yticks)
    #ax1.set_yticks(minor_yticks, minor=True) 
    ax1.yaxis.set_ticks_position('both') 
    #ax1.set_ylim(yymin,yymax)
    #ax1.set_xlim(xxmin,xxmax)
    ax1.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
    #plt.autoscale()
    #labels
    ax1.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', 
                   fontsize = size_fonts, labelpad = shift_labelpad)
    ax1.set_ylabel(r"Pressure (GPa)", fontweight = 'bold', fontsize = size_fonts, 
                   labelpad = shift_labelpad)
    return fig, ax1


def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),
                       "size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    thermofile = 'all'
    thermofile2 = ''
    filetype = 'all'
    init_letter = ''
    letters = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    #figure dictionnary
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200',
                'T5':'#ffcd01','T5.5':'#ff6e00','T6':'#ff0101','T6.5':'#ff00a2',
                'T7':'#ff01de','T10':'#ffa6f4','T15':'#ffe86e','T20':'#ffbf90'}
    #other parameters
    Na=6.022*10**23
    eVtoJ = 1.6e-19        #1eV = 1.6e-19 J
    kb = 1.38064852e-23    #boltzmann constant
    try:
        options,arg = getopt.getopt(argv,"hf:g:t:l:d:",["filename","gfilename","type",'letter','denlim'])
    except getopt.GetoptError:
        print("plot_thermo-rho_BM3.py -f <thermo_filename (default = 'all')>  -g <thermo_filename 2 (if we want to plot 2 files)>  -t <type of the file: 'short' or 'all'> -d <minimum density of the data point used in fit> -l <letter of first plot, default = ''> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_thermo-rho_BM3.py program to fit and plot BM3 as a function of density')
            print("plot_thermo-rho_BM3.py -f <thermo_filename (default = 'all' -> 3plot)> -g <thermo_filename 2 (if we want to plot 2 files)>  -t <type of the file: 'short' or 'all'> -d <minimum density of the data point used in fit> -l <letter of first plot, default = ''>")
            print('')
            sys.exit()
        elif opt in ('-f','--filename'):
            thermofile = str(arg)
        elif opt in ('-g','--gfilename'):
            thermofile2 = str(arg)
        elif opt in ('-t','--type'):
            filetype = str(arg)
        elif opt in ('-l','--letter'):
            init_letter = str(arg)
        elif opt in ('-d','--denlim'):
            denlim = float(arg)
    #selection of the column depending on the type of the thermo file
    if filetype == 'all':
        column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,
                         'err_T':10,'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,
                         'stdev_Cvm_Nkb':24,'Cvm':25,'stdev_Cvm':26}
    else:
        column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,
                         'err_T':10,'E':11,'stdev_E':12,'err_E':13,'Cvm_Nkb':14,
                         'stdev_Cvm_Nkb':15}
    #************ initialization of the plot
    if thermofile == 'all':
        fig,ax1,ax2,ax3, ax0 = creation_plot_3(plot_parameters)
        print("axis are ax1",ax1,"ax2",ax2,"ax3",ax3,"ax0",ax0)
        files = sorted(glob.glob('thermo_*O8_'+filetype+'.txt'),reverse=True) #I list every thermo files
        figurename = 'Fit_BM3_thermo_all_'+str(denlim)
    elif thermofile2 != '':
        fig,ax1,ax2, ax0  = creation_plot_2(plot_parameters)
        print("axis are ax1",ax1,"ax2",ax2,"ax0",ax0)
        files = [thermofile, thermofile2] #we will plot 2 files
        figurename = 'Fit_BM3_'+str(denlim)+'_'+thermofile.split('.txt')[0]+'_'+thermofile2.split('.txt')[0].split('_')[1]
    else:
        fig, ax = creation_plot(plot_parameters)
        files = [thermofile] #I take only the file we want
        figurename = 'Fit_BM3_'+str(denlim)+'_'+thermofile.split('.txt')[0]
    #************ creation of arrays and plot at the same time
    for file in files:
        print("****** For file",file)
        #********initialisations
        #**change of subplot
        if thermofile == 'all' or thermofile2 !='':
            if files.index(file) == 0:
                ax=ax1
                letter = init_letter                    
            if files.index(file) == 1:
                ax=ax2
                try:
                    letter = letters[letters.index(init_letter)+1]
                except ValueError:
                    print("You haven't indicated any letter")
                    letter = ''
            if files.index(file) == 2:
                ax=ax3 
                try:
                    letter = letters[letters.index(init_letter)+2]
                except ValueError:
                    print("You haven't indicated any letter")
            #print("I plot on axis",ax)
        #**initialisation of data dictionnaries
        rho = {}
        V = {}
        P = {}
        stdev_P = {}
        data = {}
        #********* creation of the BM3 txt file
        bm = fileBM3(file,denlim)
        #********* fill the dictionnaries with data
        with open(file,'r') as f:
            line = f.readline()
            entry=line.split()
            elements = entry[1:]
            line = f.readline()
            entry=line.split()
            number = entry[1:]
            f.readline()
            #calculation of M*N nedded for the calculation of volumes
            MN = 0
            for i in range(len(elements)):
                MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]   
            #extract data
            while True:
                line = f.readline()
                if not line: break
                else:    
                    entry=line.split('\n')[0].split('\t')
                    temperature, acell = split_name(entry[0])
                    try:
                        data[temperature].append([  float(entry[column_number['rho']]) , float(entry[column_number['P']]) , float(entry[column_number['stdev_P']])   ])
                    except KeyError:                      
                        data[temperature] = []
                        data[temperature].append([  float(entry[column_number['rho']]) , float(entry[column_number['P']]) , float(entry[column_number['stdev_P']])   ])
        #********* remove the data we do not want to use
        for temperature in data:
            #first we sort the data by density from biggest to smallest
            data[temperature] = sorted(data[temperature], key=lambda x: x[0], reverse=True) #sort by the first column: rho
            #second we create rho,P,stdev_P dictionnaries with only the data with positive P starting below 1000 GPa
            for ii in range(len(data[temperature])):
                if data[temperature][ii][1] < 1000:
                    if data[temperature][ii][1] < 1 or  data[temperature][ii][0] < denlim   :
                        #if the pressure is below 1GPa or below the limit density, then we take this value and STOP the loop to not take the others
                        try:
                            rho[temperature].append(data[temperature][ii][0] )
                            V[temperature].append( MN / (Na * data[temperature][ii][0]*1E-24))
                            P[temperature].append(data[temperature][ii][1]  )
                            stdev_P[temperature].append(data[temperature][ii][2] )
                        except KeyError:
                            rho[temperature] = []
                            rho[temperature].append(data[temperature][ii][0] )
                            V[temperature] = [ MN / (Na * data[temperature][ii][0]*1E-24)]
                            P[temperature] = [data[temperature][ii][1]  ]
                            stdev_P[temperature] = [data[temperature][ii][2] ]
                        break
                    else:
                        try:
                            rho[temperature].append(data[temperature][ii][0] )
                            V[temperature].append( MN / (Na * data[temperature][ii][0]*1E-24))
                            P[temperature].append(data[temperature][ii][1]  )
                            stdev_P[temperature].append(data[temperature][ii][2] )
                        except KeyError:
                            rho[temperature] = []
                            rho[temperature].append(data[temperature][ii][0] )
                            V[temperature] = []
                            V[temperature].append(MN / (Na * data[temperature][ii][0]*1E-24))
                            P[temperature] = []
                            P[temperature].append(data[temperature][ii][1]  )
                            stdev_P[temperature] = []
                            stdev_P[temperature].append(data[temperature][ii][2] )     
        #******** Fit data and plot
        #This is the theoretical model or initial guesses (rho0 in g/cm3, K0 in GPa, Kp0), based on common values 
        #m0 = {'T2':[2.6,15,4],'T3':[1.8,4,6.8],'T4':[1.6,3,6.7],'T5':[1.1,0.4,11],'T6':[0.8,0,66],'T7':[0.05,0,7],'T10':[2.6,15,4],'T15':[2.6,15,4],'T20':[2.6,15,4]} #from Razvan fit (different software)
        m0 = {'T2':[2.6,15,4],'T3':[2.6,15,4],'T4':[2.6,15,4],'T5':[2.6,15,4],
              'T6':[2.6,15,4],'T7':[2.6,15,4],'T10':[2.6,15,4],'T15':[2.6,15,4],'T20':[2.6,15,4]} #rho0 in g/cm3, K0 in GPa, Kp0
        #m0 = {'T3':[3000,19,4],'T4':[3000,19,4],'T5':[3000,15,4],'T6':[3000,15,6]} #V0 in angstrom3, K0 in GPa, Kp0
        drho = 0.1
        for temperature in sorted(rho):
            if len(rho[temperature]) > 4:
                #***************** fit using P=f(rho)
               # try:  #Least square fit
               #     m, pcov = curve_fit(fonction_BM3, rho[temperature], P[temperature], p0=m0[temperature], sigma=stdev_P[temperature], absolute_sigma=False )
               #     chi2 = chi2red(m,rho[temperature],P[temperature],stdev_P[temperature])
               #     DensityX=np.arange(rho[temperature][-1]-drho,rho[temperature][0]+drho,drho) 
               #     ax.plot(DensityX,fonction_BM3(DensityX,m[0],m[1],m[2]),'-', color=colors_T[temperature], label= temperature + ' BM3 fit')
               #     print(temperature, "Least-square fit of BM3")
               #     perr = np.sqrt(np.diag(pcov)) #stdev on the parameters
               #     bm.write(str(temperature) + '\t' + str('{:1.2e}'.format(chi2)) + '\t' + str('{:1.2e}'.format(m[0])) + '\t' + str('{:1.2e}'.format(m[1])) + '\t' + str('{:1.2e}'.format(m[2])) + '\t' + str('{:1.2e}'.format(perr[0])) + '\t' + str('{:1.2e}'.format(perr[1])) + '\t' + str('{:1.2e}'.format(perr[2])) +'\n')
                #try:   #fmin fit    
                #    m=simplex(fonction, m0[temperature], args=(rho[temperature], P[temperature], stdev_P[temperature]), maxiter = 1000000, full_output=0)
                #    chi2=chi2red(m,rho[temperature],P[temperature],stdev_P[temperature])
                #    DensityX=np.arange(rho[temperature][-1]-drho,rho[temperature][0]+drho,drho) 
                #    ax.plot(DensityX,fonction_BM3(DensityX,m[0],m[1],m[2]),'-', color=colors_T[temperature], label= temperature + ' BM3 fit')
                #    print(temperature, "Least-square fit of BM3")
                #    bm.write(str(temperature) + '\t' + str('{:1.2e}'.format(chi2)) + '\t' + str('{:1.2e}'.format(m[0])) + '\t' + str('{:1.2e}'.format(m[1])) + '\t' + str('{:1.2e}'.format(m[2])) + '\n')
                try: #odr fit
                    eos_model_BM3 = odr.Model(BM3_P_rho)  # model for fitting
                    stdev_rho = [i/1000 for i in rho[temperature]] #creation of stdev for rho
                    data_for_ODR = odr.RealData(rho[temperature], P[temperature], stdev_rho, stdev_P[temperature]) # RealData object using our data
                    odr_process = odr.ODR(data_for_ODR, eos_model_BM3, beta0=m0[temperature], maxit = 10000) # Set up orth-dist-reg with the model and data
                    odr_results = odr_process.run() # run the regression
                    print('***************',temperature, "Fit of BM3")
                    print(odr_results.stopreason)
                    #odr_results.pprint() # use the pprint method to display results
                    DensityX=np.arange(rho[temperature][-1]-drho,rho[temperature][0]+drho,drho) 
                    ax.plot(DensityX,BM3_P_rho(odr_results.beta,DensityX),'-', color=colors_T[temperature], label= temperature + ' BM3 fit')
                    string = odr_results.stopreason[0]
                    for ii in range(1,len(odr_results.stopreason)):
                        string=string+' '+odr_results.stopreason[ii]
                    bm.write(str(temperature) + '\t' + str('{:1.2e}'.format(odr_results.res_var)) + '\t' + str('{:1.2e}'.format(odr_results.beta[0])) + '\t' + str('{:1.2e}'.format(odr_results.beta[1])) + '\t' + str('{:1.2e}'.format(odr_results.beta[2])) + '\t' + str('{:1.2e}'.format(odr_results.sd_beta[0])) + '\t' + str('{:1.2e}'.format(odr_results.sd_beta[1])) + '\t' + str('{:1.2e}'.format(odr_results.sd_beta[2])) + '\t'+string+ '\n')
                    #print('chi2red', chi2red(odr_results.beta,rho[temperature],P[temperature],stdev_P[temperature]))
                    #print('residual variance',odr_results.res_var)
                #***************** fit using P=f(V)
                #try: #odr fit
                #    eos_model_BM3 = odr.Model(BM3_P_V)  # model for fitting
                #    stdev_V = [i/1000 for i in V[temperature]] #creation of stdev for rho
                #    data_for_ODR = odr.RealData(V[temperature], P[temperature], stdev_V, stdev_P[temperature]) # RealData object using our data
                #    odr_process = odr.ODR(data_for_ODR, eos_model_BM3, beta0=m0[temperature], maxit = 10000) # Set up orth-dist-reg with the model and data
                #    odr_results = odr_process.run() # run the regression
                #    print('***************',temperature, "Fit of BM3")
                #    odr_results.pprint() # use the pprint method to display results
                #    DensityX=np.arange(rho[temperature][-1]-drho,rho[temperature][0]+drho,drho) 
                #    VolumeX = [MN/(Na*DensityX[ii]*1E-24) for ii in range(len(DensityX))]
                #    ax.plot(DensityX,BM3_P_V(odr_results.beta,VolumeX),'-', color=colors_T[temperature], label= temperature + ' BM3 fit')
                #    bm.write(str(temperature) + '\t' + str('{:1.2e}'.format(odr_results.sum_square)) + '\t' + str('{:1.2e}'.format(MN/(Na*odr_results.beta[0]*1E-24))) + '\t' + str('{:1.2e}'.format(odr_results.beta[1])) + '\t' + str('{:1.2e}'.format(odr_results.beta[2])) + '\n')
                except RuntimeError:
                    print(temperature, "FAIL of fit of BM3")
                    bm.write(str(temperature) + '\t-\t-\t-\t-\n')
                ax.plot(rho[temperature],P[temperature],'o', color=colors_T[temperature], markersize = plot_parameters['size_markers'], label= temperature + ' data')
    #Create legend from custom artist/label lists  
#    legend_labels['not started'] = plt.Line2D((0,1),(0,0),  markersize = plot_parameters["size_markers"], marker='o', facecolor = 'r', edgecolor='r', linestyle='')
#    s = [(k, legend_labels[k]) for k in sorted(legend_labels.keys(),reverse = False)]
#    if filename == 'all':
#        ax=ax0
#    plt.legend([v for k,v in s],[k for k,v in s],loc='upper center', bbox_to_anchor=(0.46, 1.14),fancybox=True, fontsize = plot_parameters["size_fonts"], ncol=len(style_lines))
    figurename = figurename+'.pdf'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    #plt.show()
    
    



#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



