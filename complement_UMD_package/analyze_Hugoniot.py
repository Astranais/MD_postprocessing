#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3


        Produces the Hugoniot.txt file with the values of the Hugoniot
        

"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import numpy as np
from scipy.optimize import curve_fit
import scipy.odr as odr
import matplotlib.pyplot as plt
import crystallography as cr
import natsort



def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.outcar.umd.dat
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.strip('./')
    temperature = filename.split('_')[1]
    acell = filename.split('.outcar.umd.dat')[0].split('_')[3].strip('a')
    return temperature, acell

def Hugoniot(thermofile, ground_state, density_gs, impactor_velocities, impactor_density ):
    """creation of the newfile for hugoniot values"""
    #print(thermofile)
    newfilename = 'Hugoniot_'+thermofile.split('/')[-1].split('_')[1].split('.txt')[0]+'_'+'ground-state_'+str(density_gs)+'-'+str(ground_state[density_gs][2])+'.txt'
    nf = open(newfilename, 'w')
    string = 'Densities in kg/m3 for which the hugoniot equation equals 0, corresponding pressures in GPa from BM3 fit and Energies from Hugoniot equation\n'
    nf.write(string)
    string = "T(K)\trho(kg/m3)\tP(GPa)\trho0(kg/m3)\tK0(GPa)\tK'0\tChi2red\tE(J/kg)\tUp_crust(km/s)\tUp_Eimp1(km/s)\tUp_Eimp2(km/s)\tUp_Eimp3(km/s)\tUp_Mimp1(km/s)\tUp_Mimp2(km/s)\tUp_Mimp3(km/s)\n"
    nf.write(string)
    string = str(ground_state[density_gs][2])+'\t'+str(density_gs)+'\t'+str(ground_state[density_gs][1]*1E-9)+'\t-\t-\t-\t-\t'+str(ground_state[density_gs][0])
    #****** Calculation of particle velocities
    string = string + '\t0' 
    for Uimp in impactor_velocities:
        Up = -np.sqrt( abs((density_gs-impactor_density)) * ground_state[density_gs][1] / (density_gs*impactor_density)  )*1E-3 + Uimp
        string = string + '\t' + str(round(Up,2))
    nf.write(string+'\n')
    return nf

def convert_to_float(string):
    if string[0] == '>':
        value = float(string[1:])
    else:
        value = float(string[:])
    return value 


def poly3(x,a,b,c,d):
    """ 3rd order polynomial """
    xarr = np.asarray(x)
    yarr = a * xarr**3 + b * xarr**2 + c * xarr + d
    y = yarr.tolist()
    return y


def chi2red_3(rho,data,stdev_data,popt):
    """ function chi2 which has to be minimized"""
    y = poly3(rho,*popt)
    chi2 = 0.0
    for i in range(len(rho)):
        chi2 = chi2 + ( ( data[i] - y[i]  ) / stdev_data[i] )**2        
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

def chi2red(params,X,P,Err):
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
    chi2red = chi2 / (len(X) - len(params))
    return chi2red

def BM3_P_rho(parameters,rho):
    """definition of the function BM3 : 3rd order Birch-Murnaghan  P=f(V)"""
    rho0,K0,Kp0 = parameters
    BM3=3.*K0/2.*((rho/rho0)**(7./3.)-(rho/rho0)**(5./3.))*(1.+3./4.*(Kp0-4.)*((rho/rho0)**(2./3.)-1.))
    return BM3

def fonction_BM3(rho,rho0,K0,Kp0):
    """definition of the function BM3 : 3rd order Birch-Murnaghan  P=f(V)"""
    BM3=3.*K0/2.*((rho/rho0)**(7./3.)-(rho/rho0)**(5./3.))*(1.+3./4.*(Kp0-4.)*((rho/rho0)**(2./3.)-1.))
    return BM3




def creation_plot(filename):
    """     ********** Creation of the plot  **********    """
    plt.close()
    fig, ax1 = plt.subplots(figsize=(12,7))
    ax1.set_xlabel('Density ($kg/m^3$)',fontsize=12,fontweight='bold')
    for tl in ax1.get_xticklabels():
        tl.set_fontsize(11)
    for tl in ax1.get_yticklabels():
        tl.set_fontsize(11)
    ax1.set_ylabel("Pressure (GPa)",fontsize=12,fontweight='bold')
    ax1.grid(True, axis = 'x')
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.9, left = 0.1, hspace = 0, wspace = 0)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.title('P = f(rho) '+filename.split('.txt')[0] ,fontsize=13,fontweight='bold' )
    return fig,ax1

def main(argv):
    """     ********* Main program *********     """
    thermofile = '' #file for thermo values
    ground_state_file = '' #file for ground_state
    #extracted and created data dictionnaries
    rho = {}
    E = {}
    P = {}
    stdev_P = {}
    ground_state = {}
    hugoniot = {}
    #impact values
    impactor_velocities = [12.9,15.2,18.1,8.3,11.5,15.2] #in km/s, see bibliography excel file
    impactor_density = 3000 #kg/m3
    #figure dictionnary
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200','T5':'#ffcd01','T5.5':'#ff6e00','T6':'#ff0101','T6.5':'#ff00a2','T7':'#ff01de','T7.5':'#ffa6f4','T10':'#ffe86e','T15':'#ffbf90','T20':'#ff7788'}
    #other parameters
    Na=6.022*10**23
    eVtoJ = 1.6e-19        #1eV = 1.6e-19 J
    kb = 1.38064852e-23    #boltzmann constant
    try:
        options,arg = getopt.getopt(argv,"hf:g:t:",["filename","groundstate"])
    except getopt.GetoptError:
        print("analyze_Hugoniot.py -f <thermo_filename (from analyze fullaverage)> -t <type of the file: 'short' or 'all'> -g <ground-state_filename>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('analyze_Hugoniot.py program to find the Hugoniot point from a thermo file and the ground state file associated')
            print("analyze_Hugoniot.py -f <thermo_filename (from analyze fullaverage)> -t <type of the file: 'short' or 'all'> -g <ground-state_filename>")
            print('')
            print('requires the .txt file with ground states condition: two lines of header,  a first column containing the density in kg/m3, a second column with the specific energy in J/kg, a third column with the pressure in Pa and a last column with the temperature in K')
            print('')
            sys.exit()
        elif opt in ("-g","--groundstate"):
            ground_state_file = str(arg)
        elif opt in ('-f','--filename'):
            thermofile = str(arg)
        elif opt in ('-t','--type'):
            filetype = str(arg)
    #selection of the column depending on the type of the thermo file
    if filetype == 'all':
        column_number = {'rho':2,'P':5,'T':8,'E':14}
    else:
        column_number = {'rho':2,'P':5,'T':8,'E':11}
    #************ Extraction of elements info and ground state data
    #**creation of elements and number lists and initialization of T
    with open(thermofile,'r') as f:
        line = f.readline()
        print(line)
        entry=line.split()
        elements = entry[1:]
        line = f.readline()
        entry=line.split()
        number = entry[1:]
        f.readline()
        line = f.readline()
        entry=line.split()
        temperature0, acell0 = split_name(entry[0]) 
    #**calculation of M*N (in kg.part/mol) nedded for the calculation of densities and mass
    MN = 0
    ntot = 0
    for i in range(len(elements)):
        ntot = ntot + float(number[i])
        MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]*1E-3
    #extraction of ground-state data into a dictionnary
    with open(ground_state_file, 'r') as g:
        line = g.readline()
        line = g.readline()
        while True:
            line = g.readline()
            if not line: break
            else:
                entry = line.split()
                ground_state[int(entry[0])] = (float(entry[1]),float(entry[2]),round(float(entry[3])) ) #key = density value in kg/m3, value = (energy/mass of the cell,pressure, temperature)
                hugoniot[int(entry[0])] = {} #nested dictionnary of hugoniot equation for each T and each value of rho for ground state
    #************ Extraction of P-rho-E data
    with open(thermofile,'r') as f:
        [f.readline() for i in range(3)]
        while True:
            line = f.readline()
            if not line: break
            else:
                entry=line.split()
                temperature, acell = split_name(entry[0])
                if acell <= '15.0':
                    try:
                        rho[temperature].append(float(entry[column_number['rho']])*1000 )    #calculation density in kg/m3
                        E[temperature].append(float(entry[column_number['E']]) ) #in eV/atom
                        P[temperature].append(float(entry[column_number['P']]) ) #in GPa
                        stdev_P[temperature].append(convert_to_float(entry[column_number['P']+2]) ) #in GPa    
                    except KeyError:
                        rho[temperature] = [float(entry[column_number['rho']])*1000 ]    #calculation density in kg/m3
                        E[temperature] = [float(entry[column_number['E']]) ] #in eV/atom
                        P[temperature] = [float(entry[column_number['P']]) ] #in GPa
                        stdev_P[temperature] = [convert_to_float(entry[column_number['P']+2]) ] #in GPa
                             
    
    print("******")
    print("*********")
    print("************")
    print("****************")
    #*********** For each ground state we create the file with Hugoniot point and we fill it with data
    for density_gs in sorted(ground_state):
        print("******************************* For reference density",density_gs, 'kg/m3 and reference temperature',ground_state[density_gs][2], ' K')
        #creation of the file with hugoniot values and BM3 fit parameters
        nf = Hugoniot(thermofile, ground_state, density_gs,impactor_velocities,impactor_density  )
        #creation of the figure to see the  BM3 fit
        fig,ax = creation_plot(thermofile)
        #for each T,  calculation of Hugoniot values
        for temperature in natsort.natsorted(rho):
            print("************** For ",temperature)
            if temperature == 'T1.932' or temperature == 'T4.5'  or temperature == 'T5.5' or temperature == 'T6.5' or temperature == 'T2' or temperature == 'T7':
                print("we skip this T")
                continue
            else:
                #calculation of Hugoniot density
                string = str(int(float(temperature.strip('T'))*1000))
                hugoniot[temperature] = []
                previous_len = len(string)
                #we fill this array with data of the hugoniot equation
                for i in range(len(rho[temperature])):
                    Hg = E[temperature][i]*eVtoJ*Na*ntot/MN - ground_state[density_gs][0] + 0.5 * (P[temperature][i]*1E9+ground_state[density_gs][1]) * (1/rho[temperature][i] - 1/density_gs) 
                    hugoniot[temperature].append(Hg)
                    print(i, int(round(rho[temperature][i],0)), Hg , i-1, int(round(rho[temperature][i-1],0)), hugoniot[temperature][i-1] )
                    #and we compute the value of the rho for which Hg = 0
                    if np.sign(Hg) != np.sign(hugoniot[temperature][i-1]):
                        #print('change of sign --> we use the straight line between these 2 points to find the rho for which Hg = 0')
                        #print('the equation of such a line is Hg = (Hg[i]-Hg[i-1])/(rho[i]-rho[i-1]) * rho + Hg[i] - (Hg[i]-Hg[i-1])/(rho[i]-rho[i-1]) * rho[i]')
                        #print('Hg = 0 for rho = rho[i] - Hg[i]/( (Hg[i]-Hg[i-1])/(rho[i]-rho[i-1]) ) ')
                        rho0 = rho[temperature][i] - hugoniot[temperature][i]/( (hugoniot[temperature][i]-hugoniot[temperature][i-1])/(rho[temperature][i]-rho[temperature][i-1]) )
                        print('change of sign --> rho0 =' , rho0)
                        #******We fit a BM3 to this isotherm in order to find the correct pressure corresponding to the Hugoniot density
                        if len(rho[temperature]) >= 4:
                            #Fit of the equation and calculation of the corresponding chi2
                            m0 = {'T3':[2600,15,4],'T4':[2600,15,4],'T5':[2600,15,4],'T6':[2600,15,4],'T10':[2600,15,4],'T15':[2600,15,4],'T20':[2600,15,4]} #rho0 in kg/m3, K0 in GPa, Kp0
                            drho = 0.1
                            try: #odr fit
                                eos_model_BM3 = odr.Model(BM3_P_rho)  # model for fitting
                                stdev_rho = [i/1000 for i in rho[temperature]] #creation of stdev for rho
                                data_for_ODR = odr.RealData(rho[temperature], P[temperature], stdev_rho, stdev_P[temperature]) # RealData object using our data
                                odr_process = odr.ODR(data_for_ODR, eos_model_BM3, beta0=m0[temperature], maxit = 10000) # Set up orth-dist-reg with the model and data
                                odr_results = odr_process.run() # run the regression
                                print('***************',temperature, "Fit of BM3")
                                odr_results.pprint() # use the pprint method to display results
                                DensityX=np.arange(rho[temperature][-1]-drho,rho[temperature][0]+drho,drho) 
                                ax.plot(DensityX,BM3_P_rho(odr_results.beta,DensityX),'-', color=colors_T[temperature], label= temperature + ' BM3 fit')
                                ax.plot(rho[temperature],P[temperature],'o', color=colors_T[temperature], label= temperature + ' data')
                                #print(str(temperature) + '\t' + str('{:1.2e}'.format(odr_results.res_var)) + '\t' + str('{:1.2e}'.format(odr_results.beta[0])) + '\t' + str('{:1.2e}'.format(odr_results.beta[1])) + '\t' + str('{:1.2e}'.format(odr_results.beta[2])) + '\t' + str('{:1.2e}'.format(odr_results.sd_beta[0])) + '\t' + str('{:1.2e}'.format(odr_results.sd_beta[1])) + '\t' + str('{:1.2e}'.format(odr_results.sd_beta[2]))  + '\n')
                                chi2 = chi2red(odr_results.beta,rho[temperature],P[temperature],stdev_P[temperature])
                                #print('residual variance',odr_results.res_var)
                                P0 = fonction_BM3(rho0,odr_results.beta[0],odr_results.beta[1],odr_results.beta[2])
                                string = string + '\t' + str(round(rho0)) + '\t' + str(round(P0,2)) + '\t' + str(round(odr_results.beta[0])) + '\t' + str(round(odr_results.beta[1],2)) + '\t' + str(round(odr_results.beta[2],2)) + '\t' + str(round(chi2,2))
                            #This is the theoretical model (initial guesses), based on common values 
#                            m0=np.array([2600, 19, 6]) #rho0 in kg/m3, K0 in GPa, Kp0
#                            try:
#                                m, pcov = curve_fit(fonction_BM3, rho[temperature], P[temperature], p0=m0, sigma=stdev_P[temperature], absolute_sigma=False )
#                                chi2 = chi2red(m,rho[temperature],P[temperature],stdev_P[temperature])
#                                DensityX=np.arange(3200,6300,10) 
#                                ax.plot(DensityX,fonction_BM3(DensityX,m[0],m[1],m[2]),'-', color=colors_T[temperature], label= temperature + ' BM3 fit')
#                                ax.plot(rho[temperature],P[temperature],'o', color=colors_T[temperature], label= temperature + ' data')
#                                print("Least-square fit of BM3")
                                #print('T, Chi2red, rho0(kg/m3), K0(GPa), Kp0(no-units)')
                                #print(temperature,chi2,m[0],m[1],m[2])
#                                P0 = fonction_BM3(rho0,m[0],m[1],m[2])
#                                string = string + '\t' + str(round(rho0)) + '\t' + str(round(P0,2)) + '\t' + str(round(m[0])) + '\t' + str(round(m[1],2)) + '\t' + str(round(m[2],2)) + '\t' + str(round(chi2,2))
                            #if least-square minimization fails
                            except RuntimeError:
                                try:
                                    popt, pcov = curve_fit(poly3, rho[temperature], P[temperature], p0=None, sigma=stdev_P[temperature], absolute_sigma=False )
                                    chi2 = chi2red_3(rho[temperature],P[temperature],stdev_P[temperature],popt)
                                    DensityX=np.arange(3200,6300,10)
                                    ax.plot(DensityX,poly3(DensityX,*popt),'--', color=colors_T[temperature], label= temperature + ' poly3 fit')
                                    ax.plot(rho[temperature],P[temperature],'o', color=colors_T[temperature], label= temperature + ' data')
                                    print("Least-square fit of poly3")
                                    P0 = poly3(rho0,*popt)
                                    string = string + '\t' + str(round(rho0)) + '\t' + str(round(P0,2)) + '\t-\t-\t-\t' + str(round(chi2,2))
                                except RuntimeError:
                                    print("Fail of both BM3 a,nd poly3 fit")
                                    P0 = P[temperature][i] - hugoniot[temperature][i]/( (hugoniot[temperature][i]-hugoniot[temperature][i-1])/(P[temperature][i]-P[temperature][i-1]) )   
                                    string = string + '\t' + str(round(rho0)) + '\t' + str(round(P0,2)) + '\t-\t-\t-\t-'
                                    ax.plot(rho[temperature],P[temperature],'o', color=colors_T[temperature], label= temperature + ' data')
                        #if we don't have enough data, we find P using same method as for finding rho
                        else:
                            print("Not enough data for BM3 fit")
                            ax.plot(rho[temperature],P[temperature],'o', color=colors_T[temperature], label= temperature + ' data')
                            P0 = P[temperature][i] - hugoniot[temperature][i]/( (hugoniot[temperature][i]-hugoniot[temperature][i-1])/(P[temperature][i]-P[temperature][i-1]) )   
                            string = string + '\t' + str(round(rho0)) + '\t' + str(round(P0,2)) + '\t-\t-\t-\t-'
                        #Calculation of E0 in J/kg (E hugoniot) using Hugoniot equation and rho,P Hugoniot
                        E0 =  ground_state[density_gs][0] - 0.5 * (P0*1E9+ground_state[density_gs][1]) * (1/rho0 - 1/density_gs) 
                        string = string + '\t' + str(round(E0,2))
                        #****** Calculation of particle velocities
                        Up = np.sqrt( (P0*1E9 - ground_state[density_gs][1])*(rho0 - density_gs ) / (rho0*density_gs) )*1E-3
                        string = string + '\t' + str(round(Up,2))
                        for Uimp in impactor_velocities:
                            Up = -np.sqrt( (rho0-impactor_density) * P0*1E9 / (rho0*impactor_density)  )*1E-3 + Uimp
                            string = string + '\t' + str(round(Up,2))
                        break
                if len(string) == previous_len:
                    print('Hg = 0 not found')
                    string = string  + '\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-'
                        
                nf.write(string+'\n') 
        print('the file ', nf.name, ' is created')
        plt.legend(loc='best')
        figurename = 'BM3_fit_'+thermofile.split('/')[-1].split('_')[1].split('.txt')[0]+'_rho'+str(density_gs)+'T'+str(ground_state[density_gs][2])+'.png'
        fig.savefig(figurename, bbox_inches = 'tight', dpi = 150)
        print( figurename, ' is created')      
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



