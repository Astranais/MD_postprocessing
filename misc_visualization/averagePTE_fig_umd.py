#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: average by Razvan, modified by Anaïs (figure plot)
Langage : Python3

                ****     Calculation of the arithmetic average and stdev    ****
          and plot of the variations of Pressure, energy  without entropy and temperature
       
"""


import numpy as np 
import sys
import getopt
import os
import subprocess
import matplotlib.pyplot as plt


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def grep_pattern(FileName, Pattern, SkipSteps): 
    data = []
    average = 0 
    stdev = 0 
    variance = 0
    anchor=Pattern.split()[0]
    patterns=subprocess.check_output(['grep',Pattern,FileName]).decode("utf-8")
    greps=patterns.split('\n')
    for isteps in range(SkipSteps+1,len(greps)):
        elems=greps[isteps].split()
        for ii in range(len(elems)):
            if elems[ii] == anchor:
                for jj in range(ii+1,len(elems)):
                    if is_number(elems[jj]):
                        data.append(float(elems[jj]))
                        break
    average = sum(data)/len(data)
    for ii in range(len(data)):
        variance = variance + (data[ii]-average)**2
    variance = variance/len(data)
    stdev = np.sqrt(variance)
    print('Averages over ',len(data),' steps for', Pattern, 'are: mean = ',average,' variance = ', variance, ' stdev = ', stdev)
    return data, average, stdev

def plot_function(name,ax,data,average,stdev):
    bleuspecial = '#00ccff'
    yellowcolors = ['#ffe680','#ffc900','#dbae00','#b28f00','#4e3f00'] #colors for impactor velocities above -  yellow shades
    x=np.arange(1,len(data)+1,1)
    if name == "Pressure":
        unit = "GPa"
        smallname = 'P'
    elif name == "Energy":
        unit = "eV"
        smallname = 'E'
    elif name == "Temperature":
        unit = "K"
        smallname = 'T'
    else:
        unit = "??"
    plt.minorticks_on()
    ax.plot(x,data,'.', color=yellowcolors[0])
    ax.plot((0,x[len(x)-1]),(average, average), 'r-', linewidth=4)
    txtaverage = smallname+' = ' + str(round(average)) + ' $\pm$ ' + str(round(stdev)).rstrip('0').rstrip('.') + ' ' + unit           #for temperature
    ax.text(0.02,0.01, txtaverage,transform=ax.transAxes, fontsize=12)
    ax.autoscale(enable=True,axis='both',tight=True)


def main(argv):
    FileName = ''
    SkipSteps=0
    try:
        opts, arg = getopt.getopt(argv,"hf:s:",["fFileName,sSkipSteps"])
    except getopt.GetoptError:
        print('averagePTE+fig.py -f <FileName> -s <SkipSteps>')
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print('averagePTE+fig.py program to extract and average numerical values, + plot its variations for the patterns "Pressure","Temperature" and "InternalEnergy"')
            print('averagePTE+fig.py -f <FileName> -s <SkipSteps>')
            sys.exit()
        elif opt in ("-f", "--fFileName"):
            FileName = str(arg)
        elif opt in ("-s", "--sSkipSteps"):
            SkipSteps = int(arg)
    print(FileName)
    if (os.path.isfile(FileName)):
        #****Calculations and creation of arrays
        Pressure, Average_P, stdev_P  = grep_pattern(FileName,"Pressure", SkipSteps)
        Temperature, Average_T, stdev_T  = grep_pattern(FileName,"Temperature", SkipSteps)
        Energy, Average_E, stdev_E  = grep_pattern(FileName,"InternalEnergy", SkipSteps)
        #****Figure Plot
        #Creation of the figure
        plt.close(1)
        fig, (ax1,ax2,ax3) = plt.subplots(3,figsize=(10,12), sharex=True)         #definition of our figure and subplots
        Titre='Fluctuations of P, T and E for ' + FileName
        #ax1.set_title(Titre, fontsize=15)
        ax3.set_xlabel('Step', fontweight = 'bold', fontsize=12)
        ax1.set_ylabel('Pressure (GPa)', fontweight = 'bold', fontsize=12)
        ax2.set_ylabel('Temperature (K)', fontweight = 'bold', fontsize=12)
        ax3.set_ylabel('Internal energy (eV)', fontweight = 'bold', fontsize=12)
        # Fine-tune figure; make subplots close to each other and hide x ticks for all but bottom plot.
        fig.subplots_adjust(hspace=0, left=0.1, bottom=0.07, right=0.97, top=0.94)
        major_ticks = np.arange(0, len(Pressure)+1, 500) 
        minor_ticks = np.arange(1, len(Pressure)+1, 100)
        for ax in [ax1,ax2,ax3]:
            ax.set_xticks(major_ticks)
            ax.set_xticks(minor_ticks, minor=True)                                           
            ax.xaxis.set_ticks_position('both')                                          
            ax.yaxis.set_ticks_position('both')
            #plt.setp(ax.get_yticklabels()[0], visible=False) 
            #plt.setp(ax.get_yticklabels()[-1], visible=False) 
            #plt.setp(ax.get_xticklabels()[-1], visible=False)
        #Plot data
        plot_function('Pressure',ax1,Pressure,Average_P,stdev_P)
        plot_function('Temperature',ax2,Temperature,Average_T,stdev_T)
        plot_function('Energy',ax3,Energy,Average_E,stdev_E)
        #Save Figure
        FigureName = FileName.split('/')[-1].split('.umd')[0]
        plt.savefig(FigureName+'.pdf', bbox_inches='tight', dpi=300)
        #******Histogram Plot
        fig2, ((ax1,ax2,ax3)) = plt.subplots(1,3, figsize=(18,6))
        ax1.hist(Pressure, bins=100, density=True)
        ax1.set_xlabel('Pressure (kbar)')
        ax2.hist(Temperature, bins=100, density=True)
        ax2.set_xlabel('Temperature (K)')
        ax3.hist(Energy, bins=100, density=True)
        ax3.set_xlabel('Energy (eV)')
        #plt.savefig(FigureName+'_hist'+'.png', bbox_inches='tight', dpi=90)
        #plt.show() # must be written AFTER savefig, otherwise nothing is saved
    else:
        print('No input file or file ',FileName,'does not exist')
        sys.exit()


if __name__ == "__main__":
   main(sys.argv[1:])




