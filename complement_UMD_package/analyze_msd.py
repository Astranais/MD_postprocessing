#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Anaïs
Langage : Python3

                    ****   compute the self diffusion coefficients       ****  
                             and save all the data into a .txt file
                         

This code requires all the msd.dat files created by the script msd_umd.py

This code produces a .txt file with:
    - a column with filename
    - 6 columns per atom type
        - Diffusivity (m^2/s)
        - stdev on diffusivity (m^2/s)
        - R^2  coefficient associated to the fit
        - slope of the fit
        - intercept (y(x=0))
        - time of change of slopes (ballistic to diffusive regime)

and another .txt file with the filename of simulations corresponding to viscous fluids (msd of Si < 9 angstrom^2)
"""



"""     ********* Importation of the packages and modules used here *********     """
import sys
import getopt
import glob
import os
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def headerfile(firstfile, mineralfile,regimeseparation,SkipTime):
    """creation of the newfile with correct header"""
    firstline = ['atom']  #beginning of the first line of the file gofr
    secondline = ['file']
    # creation of the header from the first line of the first file
    with open(firstfile, 'r') as f:
        line = f.readline()
    atoms = line.strip('time_(fs)').split() #list of atoms
    for atom in atoms:
        for i in range(0,6):
            firstline.append(atom)
        secondline.extend(['D(m2/s)','D_stdev','R_squared','slope','intercept','regime_change'])
    newfilename = 'diffusivities.txt' 
    print('The file ',newfilename,' is created')
    f = open(newfilename,'w')
    if mineralfile != '':
        with open(mineralfile,'r') as hf:
            for line in hf:
                f.write(line)
    f.write("Location of ballistic to diffusif regime using linear regression curves (before "
        + str(regimeseparation[0]) + " fs and after " +  str(regimeseparation[1]) + " fs).\n")
    f.write("Calculation of atomic self diffusivities using linear regression after " + str(SkipTime) + " fs.\n")
    f.write("\t".join(x for x in firstline)+ "\n")
    f.write("\t".join(x for x in secondline)+ "\n")
    return f, atoms      #I return the newly created files f along with the list of element couples


def calculation_diffusivities(t,data,steps):
    """     ********* Calculation of the self diffusivities (using linear regression) *********     """
    i = 0
    while t[i] < steps:
        #print('time is',t[i])
        i = i+1
    NewTime = t[i:]
    NewData = data[i:]
    #least squared linear regression giving the correlation coefficient
    slope, intercept, r_value, p_value, std_err = stats.linregress(NewTime, NewData) #intercept = y(x=0), std_err = error on the slope
    R_squared = r_value**2
    diffusivity = slope * 10**(-5) / 6.
    diff_err = std_err * 10**(-5) / 6.
    return diffusivity , diff_err , R_squared , slope , intercept

def poly2_a(x,a):
    return a*x**2

def poly2_b(x,a,b):
    return a*x**2+b*x


def regime_change(t,data,regimeseparation,b):
    """     ********* Calculation of the ballistic to diffusif regime change 
                        (intersection of linear regressions) *********     """
    #***Ballistic part
    i = 0
    while t[i] < int(regimeseparation[0]):
        #print('time is',t[i])
        i = i+1
    Time1 = t[0:i+1]
    Data1 = data[0:i+1]
    #poly2 curve fit
    if b == 0:
        popt, pcov = curve_fit(poly2_a,Time1,Data1)    
        A_1 = popt[0]
        B_1 = 0
    else:
        popt, pcov = curve_fit(poly2_b,Time1,Data1)
        A_1 = popt[0]
        B_1 = popt[1]
    #***Diffusive part
    j = 0
    while t[j] < int(regimeseparation[1]):
        j = j+1
    Time2 = t[j:]
    Data2 = data[j:]
    #least squared linear regression giving the correlation coefficient
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(Time2, Data2) #intercept = y(x=0), std_err = error on the slope
    #***Intersection of the line with the poly2
    delta = (B_1-slope2)**2+4*A_1*intercept2
    if delta < 0:
        intersection = np.nan
        Yintersection = np.nan
    else:
        intersection = (-(B_1-slope2)+np.sqrt(delta))/(2*A_1)
        Yintersection = slope2*intersection+intercept2
    return A_1, B_1, slope2, intercept2, intersection, Yintersection, i, j


def creation_plot(plot_parameters):
    """     ********* Creation of the plot background and titles *********     """
    plt.close(1)
    fig, ax1 = plt.subplots(figsize=(12,7))
    ax1.set_xlabel("Time (fs)",fontsize=plot_parameters['size_fonts'],fontweight='bold')
    ax1.set_ylabel(r"Mean square displacement ($\mathregular{\AA^2}$)",
                   fontsize=plot_parameters['size_fonts'],fontweight='bold')
    
    #addition of a zoom
    # These are in unitless percentages of the figure size. (0,0 is bottom left)
    left1, bottom1, width1, height1 = [0.18, 0.56, 0.25, 0.32]
    ax11 = fig.add_axes([left1, bottom1, width1, height1])    

    for ax in [ax1,ax11]:
        ax.minorticks_on()
        ax.xaxis.set_tick_params(labelsize = plot_parameters['size_font_ticks'], 
                                 width = plot_parameters['size_lines']/2)
        ax.yaxis.set_tick_params(labelsize = plot_parameters['size_font_ticks'],
                                 width = plot_parameters['size_lines']/2) 
        ax.grid(b=True,which='major',axis='x',color='0.85',linewidth=0.5)
    return fig, ax1, ax11

    
    


def main(argv):
    """     ********* Main program *********     """
    SkipTime = 0
    mineralfile = ''
    regimeseparation = ['60','500']
    figure = 0
    Atom_colors={'Al':'#ff9999','C':'0.25','Ca':'#38bbbe','H':'w','K':'#c802c5',
                 'Na':'#0000fd','O':'#f50000','Si':'#ffff00'} 
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),
                       "size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 20}
    try:
        options,arg = getopt.getopt(argv,"hs:m:r:p:",["sSkipTime","mineralfile","regimeseparation","plot"])
    except getopt.GetoptError:
        print('analyze_msd.py -s <SkipTime>(time in fs, default=0) -m <mineralfile (default none)> -r <range of regime change in fs, default 60-500> -p <=1 to plot msd for each file, default = 0>')
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('analyze_msd.py program to compute self diffusivities from msd.dat file and write them in a .txt file per subfolder')
            print('print also in viscous_simulations.txt all the files for which the diffusivity of Si is not large enough to consider the simulation as real liquid')
            print('analyze_msd.py -s <SkipTime>(time in fs, default=0) -m <mineralfile (default none)> -r <range of regime change in fs, default 60-500> -p <=1 to plot msd for each file, default = 0>')
            print('')
            print('mineralfile is a .txt file containing the element and number (useful to compute the densities)')
            print('')
            print('This code produces a diffusivities.txt file with:')
            print('     - a column with filename')
            print('     - 6 columns per atom type')
            print('         - Diffusivity (m^2/s)')
            print('         - stdev on diffusivity (m^2/s)')
            print('         - R^2  coefficient associated to the fit')     
            print('         - slope of the fit') 
            print('         - intercept (y(x=0))') 
            print('         - time of change of slopes (ballistic to diffusive regime)')
            sys.exit()
        elif opt in ("-s","--sSkipTime"):
            SkipTime = int(arg)
        elif opt in ('-m','--mineralfile'):
            mineralfile = str(arg)
        elif opt in ('-r','--regimeseparation'):
            regimeseparation = arg.split('-')
        elif opt in ('-p','--plot'):
            figure = int(arg)
    files = sorted(glob.glob('*.msd.dat')) #I list every msd files in alphabetic order
    if files != []:
        f, atoms = headerfile(files[0], mineralfile,regimeseparation,SkipTime)                          #I create the first newfile for gofr and save the list of element couples 
        visc_f = open('viscous_simulations.txt','w')
        print("The viscous simulations are:")
        for file in files:
            #print("working on file",file)
            if figure == 1 :
                fig, ax1, ax11 = creation_plot(plot_parameters)
                figurename = file.split('/')[-1].split('.outcar.msd.dat')[0]+'_'+'intersection-regimes'
            results = [file]
            maxdata = 0
            maxshortdata = 0
            for atom in atoms:
                data = np.loadtxt(file,skiprows=1,usecols=atoms.index(atom)+1,unpack=True)
                Time = np.loadtxt(file,skiprows=1,usecols=0,unpack=True)
                try:
                    A_1, B_1, slope2, intercept2, intersection, Yintersection, index1, index2 = regime_change(Time,data,regimeseparation,0) #we compute the intercection of the poly2 fit with a linear regression (change of ballistic to diffusif regime)
                    if intersection < 0 : 
                        #print('intersection < 0 for ', atom, ' in ', file)
                        intersection = np.nan
                    if figure == 1:
                        if data[-1] > maxdata:
                            maxdata = data[-1]
                        if data[index2] > maxshortdata:
                            maxshortdata = data[index2]
                        for ax in [ax1,ax11]:
                            ax.plot(Time[0:],data[0:], 'o', markersize = 5,  markeredgewidth = 0.2, 
                                    markeredgecolor = '0.5', markerfacecolor = Atom_colors[atom]+'7f')
                            #ax.plot(Time[:index2], slope1*Time[:index2] + intercept1 , '-', color = Atom_colors[atom], linewidth = 2)
                            ax.plot(Time[:index2], poly2_b(Time[:index2],A_1,B_1) , '-',
                                    color = Atom_colors[atom], linewidth = 2)
                            ax.plot(Time[index1:], slope2*Time[index1:] + intercept2 , '--', 
                                    color = Atom_colors[atom], linewidth = 2)
                            ax.vlines(x=intersection, ymin = 0 , ymax = Yintersection+Yintersection*0.1,
                                      color = Atom_colors[atom], linewidth = 2, linestyle = ':')
                        if atom == 'K' or atom == 'Na' or atom == 'Ca':
                            #ax11.set_ylim([0,slope1*int(regimeseparation[1]) + intercept1 ]) 
                            ax11.set_ylim([0,poly2_b(int(regimeseparation[1]),A_1,B_1)]) 
                    D, D_stdev, R, slope, intercept = calculation_diffusivities(Time,data,SkipTime) #we compute the diffusivities and stdev using linear fit
                    results.extend([str(D), str(D_stdev), str(R), str(slope), str(intercept), str(intersection)])
                    if atom == 'Si':
                        #print(data[-1])
                        if data[-1] <= 9.0: #if during the simu, the Si atoms do not move to the next Si site, then the simu is considered as viscous
                            visc_f.write(file.split('.outcar.msd.dat')[0]+'\t'+str(int(Time[-1]))+'\t'+str(round(data[-1],1))+'\n')                        #we write in the visc_file the filename
                            print(file.split('.outcar.msd.dat')[0]+'\t'+str(int(Time[-1]))+'\t'+str(round(data[-1],1)))
                except IndexError: #for too short simulation we just print the filename without any results
                    results.extend(['','','','',''])
            if figure == 1 :                
                ax1.set_xlim(0,9999)
                ax1.set_ylim(0,maxdata)
                ax11.set_xlim([0,int(regimeseparation[1])])   
                ax11.set_ylim(0,maxshortdata)
                fig.savefig(figurename+'.png', bbox_inches='tight', dpi = 150)
            f.write("\t".join(x for x in results)+ "\n")                  #we write in the file the result line
        f.close()
        visc_f.close()


#     ********* Execution *********
if __name__ == "__main__":
    main(sys.argv[1:])





