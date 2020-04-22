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
    - 5 columns per atom type
        - Diffusivity (m^2/s)
        - stdev on diffusivity (m^2/s)
        - R^2  coefficient associated to the fit
        - slope of the fit
        - intercept (y(x=0))

and another .txt file with the filename of simulations corresponding to viscous fluids (msd of Si < 9 angstrom^2)
"""



"""     ********* Importation of the packages and modules used here *********     """
import sys
import getopt
import glob
import os
import numpy as np
from scipy import stats


def headerfile(firstfile, mineralfile):
    """creation of the newfile with correct header"""
    firstline = ['atom']  #beginning of the first line of the file gofr
    secondline = ['file']
    # creation of the header from the first line of the first file
    with open(firstfile, 'r') as f:
        line = f.readline()
    atoms = line.strip('time_(fs)').split() #list of atoms
    for atom in atoms:
        for i in range(0,5):
            firstline.append(atom)
        secondline.extend(['D(m2/s)','D_stdev','R_squared','slope','intercept'])
    newfilename = 'diffusivities.txt' 
    print('The file ',newfilename,' is created')
    f = open(newfilename,'w')
    if mineralfile != '':
        with open(mineralfile,'r') as hf:
            for line in hf:
                f.write(line)
    f.write("\t".join(x for x in firstline)+ "\n")
    f.write("\t".join(x for x in secondline)+ "\n")
    return f, atoms      #I return the newly created files f along with the list of element couples


def calculation_diffusivities(t,data,steps):
    """     ********* Calculation of the self diffusivities (using linear regression) *********     """
    i = 0
    while t[i] < steps:
        #print('time is',t[i])
        i = i+1
    NewTemps = t[i:]
    NewData = data[i:]
    #least squared linear regression giving the correlation coefficient
    slope, intercept, r_value, p_value, std_err = stats.linregress(NewTemps, NewData) #intercept = y(x=0), std_err = error on the slope
    R_squared = r_value**2
    diffusivity = slope * 10**(-5) / 6.
    diff_err = std_err * 10**(-5) / 6.
    return diffusivity , diff_err , R_squared , slope , intercept



def main(argv):
    """     ********* Main program *********     """
    SkipTemps = 0
    mineralfile = ''
    try:
        options,arg = getopt.getopt(argv,"hs:m:",["sSkipTemps","mineralfile"])
    except getopt.GetoptError:
        print('plot-msd.py -s <SkipTemps>(time in fs, default=0) -m <mineralfile (default none)>')
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('analyze_msd.py program to compute self diffusivities from msd.dat file and write them in a .txt file per subfolder')
            print('print also in viscous_simulations.txt all the files for which the diffusivity of Si is not large enough to consider the simulation as real liquid')
            print('analyze_msd.py -s <SkipTemps>(time in fs, default=0) -m <mineralfile (default none)>')
            print('')
            print('mineralfile is a .txt file containing the element and number (useful to compute the densities)')
            print('')
            print('This code produces a diffusivities.txt file with:')
            print('     - a column with filename')
            print('     - 5 columns per atom type')
            print('         - Diffusivity (m^2/s)')
            print('         - stdev on diffusivity (m^2/s)')
            print('         - R^2  coefficient associated to the fit')     
            print('         - slope of the fit') 
            print('         - intercept (y(x=0))') 
            sys.exit()
        elif opt in ("-s","--sSkipTemps"):
            SkipTemps = int(arg)
        elif opt in ('-m','--mineralfile'):
            mineralfile = str(arg)
    files = sorted(glob.glob('*.msd.dat')) #I list every msd files in alphabetic order
    if files != []:
        f, atoms = headerfile(files[0], mineralfile)                          #I create the first newfile for gofr and save the list of element couples 
        visc_f = open('viscous_simulations.txt','w')
        print("The viscous simulations are:")
        for file in files:
            #print("working on file",file)
            results = [file]
            for atom in atoms:
                data = np.loadtxt(file,skiprows=1,usecols=atoms.index(atom)+1,unpack=True)
                Temps = np.loadtxt(file,skiprows=1,usecols=0,unpack=True)
                try:
                    D, D_stdev, R, slope, intercept = calculation_diffusivities(Temps,data,SkipTemps) #we compute the diffusivities and stdev using linear fit
                    results.extend([str(D), str(D_stdev), str(R), str(slope), str(intercept)])
                    if atom == 'Si':
                        #print(data[-1])
                        if data[-1] <= 9.0: #if during the simu, the Si atoms do not move to the next Si site, then the simu is considered as viscous
                            visc_f.write(file.split('.outcar.msd.dat')[0]+'\t'+str(int(Temps[-1]))+'\t'+str(round(data[-1],1))+'\n')                        #we write in the visc_file the filename
                            print(file.split('.outcar.msd.dat')[0]+'\t'+str(int(Temps[-1]))+'\t'+str(round(data[-1],1)))
                except IndexError: #for too short simulation we just print the filename without any results
                    results.extend(['','','',''])
            f.write("\t".join(x for x in results)+ "\n")                  #we write in the file the result line
        f.close()
        visc_f.close()


#     ********* Execution *********
if __name__ == "__main__":
    main(sys.argv[1:])





