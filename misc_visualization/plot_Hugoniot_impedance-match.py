#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3


        program to plot impedance match diagram
        requires the files Hugoniot-....txt  created by the script analyze_Hugoniot.py
         *************  ARTICLE VERSION  *************

"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D #useful to create a custom legend
import natsort
import glob
import re
from scipy.stats import linregress
from scipy.optimize import curve_fit




def creation_plot_Hug(plot_parameters,variable):
    """     ********** Creation of the plot  **********    """
    #print("I use creation plot")
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
    #Adjustment of ticks                
    if variable == 'T':                  
        major_xticks = np.arange(0, 21000, 5000) 
        minor_xticks = np.arange(0, 21000, 1000)
    else:
        major_xticks = np.arange(2500, 7000, 1000) 
        minor_xticks = np.arange(2500, 7000, 500)
        
    ax.yaxis.set_ticks_position('both')
                     
    ax.set_xticks(major_xticks)
    ax.set_xticks(minor_xticks, minor=True) 
    ax.xaxis.set_ticks_position('both') 

    plt.autoscale(enable=True,axis='x',tight=False)

    #ax.grid(True, which='both',axis = 'x', )
    ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
    ax.tick_params(axis='x', rotation=45)
    #labels
    ax.set_ylabel(r"Pressure (GPa)", fontweight = 'bold', fontsize = size_fonts)
    ax.set_ylim(-5,400)
    if variable == 'T':
        ax.set_xlabel(r"Temperature (K)", fontweight = 'bold', fontsize = size_fonts)
    else:
        ax.set_xlabel(r'Density (kg.m$^{-3}$)', fontweight = 'bold', 
                      fontsize = size_fonts, labelpad = shift_labelpad)
    return fig, ax



def creation_plot(plot_parameters):
    """     ********** Creation of the plot  **********    """
    #print("I use creation plot")
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
    #Adjustment of ticks
    major_xticks = np.arange(0, 20, 5) 
    minor_xticks = np.arange(0, 20, 1)

    ax.yaxis.set_ticks_position('both')
                     
    ax.set_xticks(major_xticks)
    ax.set_xticks(minor_xticks, minor=True) 
    ax.xaxis.set_ticks_position('both') 
    
    plt.autoscale(enable=True,axis='y',tight=False)
    plt.xlim(0,20)
    plt.ylim(0,600)

    #ax.grid(True, which='both',axis = 'x', )
    ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
    #labels
    ax.set_xlabel(r'Particle velocity (km/s)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad)
    ax.set_ylabel(r"Pressure (GPa)", fontweight = 'bold', fontsize = size_fonts)
    return fig, ax


def Shock_impact(file):
    """creation of the newfile for hugoniot values"""
    #print(file)
    newfilename = 'Shock-state'+file.split('.txt')[0].strip('Hugoniot')+'.txt'
    #print(newfilename)
    nf = open(newfilename, 'w')
    string = 'Shock state obtained from impedance match from impact velocities selected in analyze_Hugoniot script\n'
    nf.write(string)
    string = "rho(kg/m3)\tT(K)\tP(GPa)\tUp(km/s)\tU_impactor(km/s)\ta_imp\tb_imp\tc_imp\ta_crust\tb_crust\tc_crust\n"
    nf.write(string)
    return nf,newfilename

def impedance_match_linear(ax, ii, Up_imp, P_imp, colors, XUp, int_crust, s_crust, Up_shock, P_shock):
    #P-Up line for impactor and linear fit
    ax.plot(Up_imp,P_imp,'o',color = colors[ii])
    if len(Up_imp[4:]) > 3:
        s, intercept, r, p, stderr = linregress(Up_imp[4:],P_imp[4:])
    else: #when there are missing data we take all the remaining data
        s, intercept, r, p, stderr = linregress(Up_imp[1:],P_imp[1:])
    ax.plot(XUp, s*XUp+intercept, '--', color = colors[ii])
    #Up of shock = intersection of linear regression
    UpShock = (intercept-int_crust)/(s_crust-s)
    Up_shock.append(UpShock)
    PShock = s*UpShock+intercept
    P_shock.append(PShock)
    return Up_shock, PShock, P_shock

def impedance_match_poly(ax, ii, Up_imp, P_imp, colors, XUp, popt_crust, Up_shock, P_shock,  a_imp, b_imp, c_imp ):
    #P-Up cruve for impactor and poly2 fit
    ax.plot(Up_imp,P_imp,'o',color = colors[ii])
    popt, pcov = curve_fit(poly2, Up_imp, P_imp, p0=None, sigma=None, absolute_sigma=False)
    ax.plot(XUp, poly2(XUp, *popt), '-', color = colors[ii])
    #Up of shock = intersection of curves = intersection of polynomials = smallest solution (-b-sqrt(delta))/2a
    Delta = (popt[1]-popt_crust[1])**2 - 4*(popt[0]-popt_crust[0])*(popt[2]-popt_crust[2])
    UpShock = (-(popt[1]-popt_crust[1]) - np.sqrt(Delta))/(2*(popt[0]-popt_crust[0]))
    Up_shock.append(UpShock)
    PShock = popt[0]*UpShock**2+popt[1]*UpShock +popt[2]
    P_shock.append(PShock)
    a_imp.append(popt[0])
    b_imp.append(popt[1])
    c_imp.append(popt[2])
    return Up_shock, PShock, P_shock, a_imp, b_imp, c_imp



def poly2(x,a,b,c):
    """ 2nd order polynomial """
    xarr = np.asarray(x)
    yarr = a * xarr**2 + b * xarr + c 
    y = yarr.tolist()
    return y

    
def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),
                       "size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    #other dictionnaries and parameters for the figure
    Up_shock = [] #Up of shock = intersection of linear regression or poly2 fits
    P_shock = []  #P of shock
    rho_shock = []
    T_shock = []
    a_imp, b_imp, c_imp, a_crust, b_crust, c_crust = [],[],[],[],[],[] #coeff of poly2 fits for impactor and crust
    U_impactor = [12.9,15.2,18.1,8.3,11.5,15.2] #in km/s, impactor velocities, see bibliography excel file
    #colors = ['#5d9f00','#215f00','#003c00','#b8ff00','#90d500','#215f00'] #colors for impactor velocities above -  green shades
    #colors = ['#a13f9c','#691971','#440055','#ff80e5','#d260c2','#691971'] #colors for impactor velocities above - purple/pink shades
    colors = ['#dbae00','#b28f00','#4e3f00','#ffe680','#ffc900','#b28f00'] #colors for impactor velocities above -  yellow shades
    markertype = ['>','^','D','v','<'] #one marker per impactor velocity
    try:
        options,arg = getopt.getopt(argv,"hf:l:",["filename","letters"])
    except getopt.GetoptError:
        print("plot_Hugoniot_impedance-match.py -f <mineral name 1> -l <letters> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_Hugoniot_impedance-match.py program to plot impedance match diagram')
            print("plot_Hugoniot_impedance-match.py -f <mineral name 1> -l <letters (at least 2, ex: a,b)> ")
            print('')
            sys.exit()
        elif opt in ('-f','--filename'):
            mineral1 = str(arg)
        elif opt in ('-l','--letters'):
            letters = arg.split(',')
    print('**********************')
    print('***************')
    print('*********')
    print('*****')    
    #figure creation for impedance match -- ARTICLE VERSION
    fig4, ax4 = creation_plot(plot_parameters)
    figurename4 = 'impedance-match_'+mineral1
    #We take all the Hugoniot files created by the script analyze_hugoniot.py (one file per reference state)
    files = sorted(glob.glob('Hugoniot_'+mineral1+'_ground-state_*.txt'))
    #We remove the file with density = 3.0 which is used to simulate the impactor only
    selected_files = []
    for file in files:
        if re.search('ground-state_3000-0',file):
            impactor_Hugoniot_filename = file
        else:
            selected_files.append(file)
    #For each reference state we fit and plot 
    for file in selected_files:
        print('****** File',file)
        #initialization of lists
        exiting = 0
        Up_shock = []
        P_shock = []
        rho_shock = []
        T_shock = []
        #figure creation for impedance match
        fig, ax = creation_plot(plot_parameters)
        figurename = 'impedance-match'+file.split('.txt')[0].strip('Hugoniot')
        #figure creation for Hugoniot P-rho
        variable = 'rho'
        fig2, ax2 = creation_plot_Hug(plot_parameters, variable)    
        figurename2 = 'Hugoniot-P-rho'+file.split('.txt')[0].strip('Hugoniot')
        #figure creation for Hugoniot P-T
        variable = 'T'
        fig3, ax3 = creation_plot_Hug(plot_parameters, variable)                
        figurename3 = 'Hugoniot-P-T'+file.split('.txt')[0].strip('Hugoniot')
        #file creation for shock impact point
        nf,newfilename = Shock_impact(file)
        #creation of new X array to plot linear regression
        XUp = np.arange(0,18.2,0.1)
        Xrho = np.arange(2500,7000,500)
        XT = np.arange(0,20000,1000)
        #extraction, plot, fit of data and determination of intersection
        try:
            T, rho,P, Up_crust = np.loadtxt(file,delimiter = '\t', skiprows = 2,
                                            usecols = (0,1,2,8), unpack =True)
        except ValueError:
            print("Check your file",file,"there are missing data, skipping this file")
            continue        
        #we extract values for T >= 10 000 K
        smallP = P[-3:]
        smallrho = rho[-3:]
        smallT = T[-3:]

        #**** P-rho plot for the crust and linear fit using only T >= 10000 K
        ax2.plot(rho,P,'o',color = 'k')
        s_rho, int_rho, r_rho, p_rho,stderr_rho = linregress(smallrho,smallP)
        ax2.plot(Xrho,s_rho*Xrho+int_rho,'--',color = 'b')
        
        #**** P-T plot for the crust and linear fit using only T >= 10000 K
        ax3.plot(T,P,'o',color = 'k')
        s_T, int_T, r_T, p_T,stderr_T = linregress(smallT,smallP)
        ax3.plot(XT,s_T*XT+int_T,'--',color = 'b')
        
        #**** P-Up line for crust 
        ground_state = file.split('ground-state_')[-1].split('.txt')[0].split('-')
        if ground_state[0] == '2600':
            color = '#0080ed' #blue
        else: 
            color = '0.5'            
        ax.plot(Up_crust, P, 'o', color = color)
        ax.text(0.015,0.945, letters[0] , transform=ax.transAxes, horizontalalignment = 'left', 
                fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                bbox=dict(facecolor='w', edgecolor='k', pad=3.0))
        #linear fit
        #s_crust, int_crust, r_crust, p_crust, stderr_crust = linregress(Up_crust[4:],P[4:])
        #ax.plot(XUp, s_crust*XUp+int_crust, '--', color = '0.5')
        #2nd order poly fit
        popt_crust, pcov_crust = curve_fit(poly2, Up_crust, P, p0=None, sigma=None, absolute_sigma=False)
        ax.plot(XUp,poly2(XUp,*popt_crust),'-', color=color)
        a_crust = np.full(6, popt_crust[0])
        b_crust = np.full(6, popt_crust[1])
        c_crust = np.full(6, popt_crust[2])
        for ii in range(6):
            try: 
                P_imp, Up_imp =  np.loadtxt(impactor_Hugoniot_filename,delimiter = '\t',
                                            skiprows = 2, usecols = (2,ii+9), unpack = True)
                #print("Using the Hugoniot of the impactor for impedance match")
            except ValueError:
                print("Check your file",impactor_Hugoniot_filename ,"there are missing data, skipping this file")
                exiting = 1
                break
            except UnboundLocalError:
                print("Hugoniot of the impactor not found, we use the Hugoniot of the same material")
                P_imp, Up_imp =  np.loadtxt(file,delimiter = '\t', skiprows = 2, 
                                            usecols = (2,ii+9), unpack = True)
            #impedance match using linear regression of the last data
            #Up_shock, PShock, P_shock = impedance_match_linear(ax, ii, Up_imp, P_imp, colors, XUp, int_crust, s_crust, Up_shock, P_shock)
            #impedance match using 2nd order polynomial of all data
            XUp = np.arange(0, U_impactor[ii]+0.1,0.1)
            Up_shock, PShock, P_shock,  a_imp, b_imp, c_imp = impedance_match_poly(ax, ii, Up_imp, P_imp, 
                                                                                   colors, XUp, popt_crust, 
                                                                                   Up_shock, P_shock, 
                                                                                   a_imp, b_imp, c_imp)
            #***** We obtain corresponding rho and T shock using the Hugoniot values (we draw a line between the two points around the Pshock)
            #first we check of the shock pressure is among the hugoniot values
            index = -1
            for jj in range(len(P)):
                if PShock < P[jj]:
                    index = jj
                    break
            print('PShock is:',PShock)
            #if it is above the last hugoniot pressure, then we use a fit of the Hugoniot using the last 4 values to estimate the corresponding T and rho
            if index == -1:
                print('PShock is above last value',P[index])
                rho_shock.append((PShock - int_rho)/s_rho)
                T_shock.append((PShock - int_T)/s_T)
            #if PShock is surrounded by values of the Hugoniot, then we use them to bracket and estimate the corresponding T and rho
            else:
                print('index of P is:',index, 'with PShock between', P[index-1], 'and', P[index])
                slope_rho = ( P[index] - P[index-1] ) / (rho[index] - rho[index-1])
                inter_rho = P[index] - slope_rho * rho[index]
                rho_shock.append((PShock - inter_rho)/slope_rho)
                slope_T = ( P[index] - P[index-1] ) / (T[index] - T[index-1])
                inter_T = P[index] - slope_T * T[index]
                T_shock.append((PShock - inter_T)/slope_T)
            
        if exiting == 1:
            continue
        extension =  '.pdf'
        figurename = figurename + extension
        figurename2 = figurename2 + extension        
        figurename3 = figurename3 + extension
        for figname in [figurename,figurename2,figurename3]:
            print("figure saved with name ",figname)
        for ii in range(len(rho_shock)):
            nf.write(str(round(rho_shock[ii]))+'\t'+str(round(T_shock[ii]))+'\t'+str(round(P_shock[ii]))+'\t'+str(round(Up_shock[ii],2))+'\t'+str(round(U_impactor[ii],2))+'\t'+str(a_imp[ii])+'\t'+str(b_imp[ii])+'\t'+str(c_imp[ii])+'\t'+str(a_crust[ii])+'\t'+str(b_crust[ii])+'\t'+str(c_crust[ii])+'\n' )
        print("file saved with name ",newfilename)
        fig.savefig(figurename, bbox_inches = 'tight', dpi = 150)
        fig2.savefig(figurename2, bbox_inches = 'tight', dpi = 150)
        fig3.savefig(figurename3, bbox_inches = 'tight', dpi = 150)
        
        #*********** Plot the figure for the article
        #** extraction of data
        ground_state = file.split('ground-state_')[-1].split('.txt')[0].split('-')
        if ground_state[1] == '3000':
            #color_3000 = 'k'
            #color_3000= '#800000' #dark red
            color_3000= '#ff0000' #red
            colorfill_3000 = color_3000
            params_3000 = popt_crust
            P_shock_3000 = P_shock
            Up_shock_3000 = Up_shock
        elif ground_state[1] == '1932':
            #color_1932 = '#00bc00' #green
            #color_1932 = '#d260c2' #pink
            color_1932 = '#ff8a00' #orange 
            colorfill_1932 = color_1932
            params_1932 = popt_crust
            P_shock_1932= P_shock
            Up_shock_1932 = Up_shock
        elif ground_state[0] == '2500': 
            #color_2500 = '#0078dd' #blue
            color_2500 = '#b1dbff' #light blue
            colorfill_2500 =  color_2500
            params_2500 = popt_crust
            P_shock_2500 = P_shock
            Up_shock_2500 = Up_shock
        elif ground_state[0] == '2600':
            #color_2600 = '#ffc900' #yellow
            #color_2600 = '#cccccc' #grey
            color_2600 = '#0080ed' #blue
            colorfill_2600 =  color_2600
            params_2600 = popt_crust
            P_shock_2600 = P_shock
            Up_shock_2600 = Up_shock
        elif ground_state[0] == '2700':
            #color_2700 = '#cccccc' #grey
            #color_2700 = '#00bc00' #green
            color_2700 = '#004f92' #dark blue
            colorfill_2700 =  color_2700
            params_2700 = popt_crust
            P_shock_2700 = P_shock
            Up_shock_2700 = Up_shock
        if mineral1 == 'CaAl2Si2O8':
            line = '-'
        elif mineral1 == 'KAlSi3O8':
            line = '--'
            colorfill_3000 = colorfill_3000+'7f'
            colorfill_1932 = colorfill_1932+'7f'
            colorfill_2700 = colorfill_2700+'7f'
            colorfill_2600 = colorfill_2600+'7f'
            colorfill_2500 = colorfill_2500+'7f'
        else:
            colorfill_3000 = 'w'
            colorfill_1932 = 'w'
            colorfill_2700 = 'w'
            colorfill_2600 = 'w'
            colorfill_2500 = 'w'
            line = ':'
    #** plot the article figure
    #plot the crust
    ax4.plot(XUp, poly2(XUp, *params_2700), linestyle = line, 
             linewidth = plot_parameters['size_lines'], color = color_2700)
    ax4.plot(XUp, poly2(XUp, *params_2600), linestyle = line, 
             linewidth = plot_parameters['size_lines'], color = color_2600)
    ax4.plot(XUp, poly2(XUp, *params_2500), linestyle = line, 
             linewidth = plot_parameters['size_lines'], color = color_2500)
    ax4.plot(XUp, poly2(XUp, *params_1932), linestyle = line, 
             linewidth = plot_parameters['size_lines'], color = color_1932)
    ax4.plot(XUp, poly2(XUp, *params_3000), linestyle = line, 
             linewidth = plot_parameters['size_lines'], color = color_3000)
    for ii in range(len(markertype)):
        XUp = np.arange(0, U_impactor[ii]+0.1,0.1)
        #plot the impactors curves
        ax4.plot(XUp, poly2(XUp, a_imp[ii], b_imp[ii], c_imp[ii]), linestyle = line, color = colors[ii])
        #plot the peak shock
        ax4.plot(Up_shock_2700[ii], P_shock_2700[ii], ls = '',  marker = markertype[ii], 
                 markersize = plot_parameters["size_markers"]+2, markeredgecolor = color_2700, 
                 color = 'none',  markeredgewidth = 0.5,  markerfacecolor = colorfill_2700)
        ax4.plot(Up_shock_2600[ii], P_shock_2600[ii], ls = '',  marker = markertype[ii], 
                 markersize = plot_parameters["size_markers"]+2, markeredgecolor = color_2600, 
                 color = 'none',  markeredgewidth = 0.5,  markerfacecolor = colorfill_2600)
        ax4.plot(Up_shock_2500[ii], P_shock_2500[ii], ls = '',  marker = markertype[ii], 
                 markersize = plot_parameters["size_markers"]+2, markeredgecolor = color_2500, 
                 color = 'none',  markeredgewidth = 0.5,  markerfacecolor = colorfill_2500)
        ax4.plot(Up_shock_1932[ii], P_shock_1932[ii], ls = '',  marker = markertype[ii], 
                 markersize = plot_parameters["size_markers"]+2, markeredgecolor = color_1932, 
                 color = 'none',  markeredgewidth = 0.5,  markerfacecolor = colorfill_1932)
        ax4.plot(Up_shock_3000[ii], P_shock_3000[ii], ls = '',  marker = markertype[ii], 
                 markersize = plot_parameters["size_markers"]+2, markeredgecolor = color_3000, 
                 color = 'none',  markeredgewidth = 0.5,  markerfacecolor = colorfill_3000)
    #letter
    ax4.text(0.015,0.945, letters[1] , transform=ax4.transAxes, horizontalalignment = 'left', 
             fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
             bbox=dict(facecolor='w', edgecolor='k', pad=3.0))
    #save
    figurename4 = figurename4 + extension
    print("figure saved with name ",figurename4)
    fig4.savefig(figurename4, bbox_inches = 'tight', dpi = 300)
#    plt.show() # à mettre après savefig sinon on ne sauvegarde rien !!  

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



