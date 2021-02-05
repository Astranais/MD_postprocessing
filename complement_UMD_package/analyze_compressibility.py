#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3


                 Compute the thermal pressure coefficient,
        the isothermal compressibility and the isobaric expansivity
                 
         
         For use after the fullaverages.py script !!!!!!!!!!

"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import natsort
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import natsort
import glob
from matplotlib.colors import LinearSegmentedColormap


def create_colors():
    """ function to create the colors_T dictionnary from color map """
    colors_T = {}
    cm_data = np.loadtxt("/Users/akobsch/Dropbox/Recherche/PhD-MD-silicates/simulations/scripts/Analysis/ScientificColourMaps6/romaO/romaO.txt")
    cm_data = cm_data[::-1] #for reverse colors
    new_map = LinearSegmentedColormap.from_list('new', cm_data)
    temperatures = ['T2','T2.5','T3','T3.5','T4','T4.5','T5','T5.5','T6','T6.5','T7','T7.5','T7.7']
    color = iter(new_map(np.linspace(0,1,len(temperatures)))) #Creation of the color list    
    for T in temperatures:
        c = next(color)
        colors_T[T] = c    
    return colors_T

    
def split_name(filename): #ConstantType = 'T' for temperature 
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.outcar.umd.dat
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.strip('./')
    constant = filename.split('_')[1] #temperature
    #constant = filename.split('.outcar.umd.dat')[0].split('_')[3].strip('a') #acell
    return constant 
    
def convert_to_float(string):
    if string[0] == '>':
        value = float(string[1:])
    else:
        value = float(string[:])
    return value 

def creation_plot_rho_P(plot_parameters):
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
    fig2, ax2 = plt.subplots(1,1, figsize = size_figure, linewidth = size_lines /2)  
    #Adjustment of ticks
    major_xticks = np.arange(0, 350, 50) 
    minor_xticks = np.arange(0, 350, 10) 
    ax2.set_xticks(major_xticks)
    ax2.set_xticks(minor_xticks, minor=True) 
    ax2.xaxis.set_ticks_position('both') 
    major_yticks = np.arange(0,8,1)
    minor_yticks = np.arange(0,8,0.2)
    ax2.set_yticks(major_yticks)
    ax2.set_yticks(minor_yticks, minor=True) 
    ax2.yaxis.set_ticks_position('both')
    #ax2.set_ylim(yymin,yymax)
    #ax2.set_xlim(xxmin,xxmax)
    ax2.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
    plt.autoscale()
    #labels
    ax2.set_xlabel(r'Pressure (GPa)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    ax2.set_ylabel(r"Density (g.cm$^{-3}$)", fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    return fig2, ax2


def creation_plot(plot_parameters):
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
    fig, ax1 = plt.subplots(1,1, figsize = size_figure, linewidth = size_lines /2)  
    #Adjustment of ticks
    major_xticks = np.arange(0,10000,1000)
    minor_xticks = np.arange(0,10000,500)
    ax1.set_xticks(major_xticks)
    ax1.set_xticks(minor_xticks, minor=True) 
    ax1.xaxis.set_ticks_position('both')
    major_yticks = np.arange(0, 350, 50) 
    minor_yticks = np.arange(0, 350, 10)
    ax1.set_yticks(major_yticks)
    ax1.set_yticks(minor_yticks, minor=True) 
    ax1.yaxis.set_ticks_position('both') 
    #ax1.set_ylim(yymin,yymax)
    #ax1.set_xlim(xxmin,xxmax)
    ax1.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
    plt.autoscale()
    #labels
    ax1.set_xlabel(r'Temperature (K)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    ax1.set_ylabel(r"Pressure (GPa)", fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    return fig, ax1


def append_function(data, constant,entry, column_number, propY, propX):
    data[constant].append([  float(entry[column_number[propX]]) , float(entry[column_number[propY]]) , float(entry[column_number[propY]+1]) , convert_to_float(entry[column_number[propY]+2])   ])
def extractdata(file,data,column_number,ConstantType,propY, propX):
    #********* fill the dictionnaries with data
    with open(file,'r') as f:
        [f.readline() for i in range(3)]
        #extract data
        while True:
            line = f.readline()
            if not line: break
            else:    
                entry=line.split('\n')[0].split('\t')
                if ConstantType == 'T':
                    constant = split_name(entry[0]) #ConstantType = 'T' for temperature
                else:
                    constant = float(entry[column_number['rho']]) #density in g/cm3
                try:
                    append_function(data,constant,entry, column_number,propY, propX)
                except KeyError:                      
                    data[constant] = []
                    append_function(data,constant,entry, column_number,propY, propX)


def plot_fit_save(ax,rho, X, Y, plot_parameters,TPC,colors_rho):
    #transform to arrays
    X_data = np.asarray(X[rho])
    Y_data = np.asarray(Y[rho])
    #plot data
    try:
        ax.plot(X_data,Y_data, linestyle = '', marker='o', markersize = plot_parameters["size_markers"], 
                color = colors_rho[rho], label = rho)
        slope, intercept, rvalue, pvalue, stderr = linregress(X_data,Y_data)
        ax.plot(X_data,slope*X_data+intercept, linestyle = '--', color = colors_rho[rho])
        TPC[rho] = slope      
    except KeyError:
        print('pass plot')
        pass

    
def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),
                       "size_markers" : 6,"size_lines" : 1,"shift_labelpad" : 10}
    #other dictionnaries and parameters for the figure
    #colors_T = {'T3':'#297fff','T4':'#00ff00','T4.5':'#bae200','T5':'#ffcd01','T5.5':'#ff6e00','T6':'#ff0101','T6.5':'#ff00a2','T7':'#ff01de','T7.5':'#ffa6f4'}#,'T10':'#ffe86e','T15':'#ffbf90','T20':'#ff7788'}
    colors_T = create_colors()
    filetype = 'all'
    try:
        options,arg = getopt.getopt(argv,"hf:t:",["filename","type"])
    except getopt.GetoptError:
        print("analyze_compressibility.py -f <thermo_filename>  -t <type of the file: 'short' or 'all'> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('analyze_compressibility.py program to compute the thermal pressure coefficient, isothermal compressibility and isobaric expansivity.')
            print("analyze_compressibility.py  -f <thermo_filename>  -t <type of the file: 'short' or 'all'> ")
            print('')
            sys.exit()
        elif opt in ('-f','--filename'):
            thermofile = str(arg)
        elif opt in ('-t','--type'):
            filetype = str(arg)
    #************ initialization of the column number
    if filetype == 'all':
        column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,
                         'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,
                         'Cvm':25,'stdev_Cvm':26,'testCv':31,'stdev_testCv':32}
    else:
        column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,
                         'E':11,'stdev_E':12,'err_E':13,'Cvm_Nkb':14,'stdev_Cvm_Nkb':15,
                         'testCv':16,'stdev_testCv':17}
    #************ 1) compute thermal pressure coefficient: TPC = dP/dT at each rho
    print('*************** Thermal Pressure Coefficient (in GPa/K)')
    data = {}
    TPC = {}
    colors_rho = {}
    #extract all P vs T isochores (rho as keys of data dictionnary)
    extractdata(thermofile,data,column_number,'D','P','T')
    #refine data: remove the data we do not want to use
    X = {}
    Y = {}
    for rho in natsort.natsorted(data):
        if float(rho) < 1.5: #we do not take values in the liquid-gas part (density in g/cm3)
            continue
        else:
            #first we sort the data by temperature from smallest to biggest
            data[rho] = sorted(data[rho], key=lambda x: x[0], reverse=False) #sort by the first column: X = T
            #second we create X,Y dictionnaries
            for ii in range(len(data[rho])):
                if data[rho][ii][0] > 8000:#we do not take very high T
                    #if the temperature is too high then we stop the loop to not take the others
                    break
                else:
                    try:
                        X[rho].append(data[rho][ii][0] )
                        Y[rho].append(data[rho][ii][1] )
                    except KeyError:
                        X[rho] = []
                        Y[rho] = []
                        X[rho].append(data[rho][ii][0] )
                        Y[rho].append(data[rho][ii][1] )
    for rho in Y:
        if len(Y[rho]) > 1:
            colors_rho[rho] = ''            
    #create colors
    namecolor = 'lajolla'
    cm_data = np.loadtxt("/Users/akobsch/Dropbox/Recherche/PhD-MD-silicates/simulations/scripts/Analysis/ScientificColourMaps6/"+namecolor+"/"+namecolor+".txt")
    #cm_data = cm_data[::-1] #for reverse colors
    new_map = LinearSegmentedColormap.from_list('new', cm_data[20:]) # cm_data[20:]) #for devon
    color = iter(new_map(np.linspace(0,1,len(colors_rho)))) #Creation of the color list  
    for key in natsort.natsorted(colors_rho):
        c = next(color)
        colors_rho[key] = c
    #plot and fit P vs T for each isochor --> obtain TPC
    fig,ax = creation_plot(plot_parameters)
    for rho in natsort.natsorted(colors_rho,reverse=True):
        plot_fit_save(ax,rho, X, Y, plot_parameters,TPC, colors_rho)
    #print TPC
    print(TPC) #thermal pressure coefficient in GPa/K
    #legend
    ncols = 1
    rho_list = [label for label in natsort.natsorted(colors_rho, reverse=True)]
    custom_patch = [mpatches.Patch(color=colors_rho[key]) for key in rho_list]
    legend = ax.legend([col for col in custom_patch],[str(label) for label in rho_list],title = '$\\bf{Density}$ \n (g.cm$^{-3}$)', bbox_to_anchor=(1.02, 1), loc="upper left", fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol=ncols)
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])
 
    
    #************ 2) compute isothermal compressibility: beta = -1/rho drho/dP at each T
    print('*************** Isothermal Compressibility (in GPa-1)')
    data = {}
    beta = {}
    allT = {}
    #extract all P vs rho isotherms (T as keys of data dictionnary)
    extractdata(thermofile,data,column_number,'T','P','rho')
    #refine data: remove the data we do not want to use
    X = {} #rho
    Y = {} #P
    for temperature in data:
        if temperature == 'T10' or temperature == 'T15' or temperature == 'T20' or temperature == 'T1.932': #we do not take values at too high T
            continue
        else:
            allT[temperature] = '' #dictionnary for the legend
            #first we sort the data by density from smallest to biggest
            data[temperature] = sorted(data[temperature], key=lambda x: x[0], reverse=False) #sort by the first column: X = rho
            #second we create X,Y dictionnaries
            for ii in range(len(data[temperature])):
                try:
                    X[temperature].append(data[temperature][ii][0] )
                    Y[temperature].append(data[temperature][ii][1] )
                except KeyError:
                    X[temperature] = []
                    Y[temperature] = []
                    X[temperature].append(data[temperature][ii][0] )
                    Y[temperature].append(data[temperature][ii][1] )
    #plot rho vs P
    fig2,ax2 = creation_plot_rho_P(plot_parameters)
    #compute beta  
    for temperature in allT:
        beta[temperature] = {'rho':[],'P':[],'beta':[]}
        ax2.plot(Y[temperature], X[temperature],linestyle = '', marker='o', markersize = plot_parameters["size_markers"], color = colors_T[temperature], label = temperature )
        for ii in range(1,len(X[temperature])-1):
            #only for rho for which we have TPC
            if X[temperature][ii] in rho_list:
                beta[temperature]['rho'].append(X[temperature][ii])
                beta[temperature]['P'].append(Y[temperature][ii])
                #central finite difference
                slope = (X[temperature][ii+1]-X[temperature][ii-1])/(Y[temperature][ii+1]-Y[temperature][ii-1]) 
                beta[temperature]['beta'].append( -1/X[temperature][ii] * slope)
                ax2.plot(np.array([Y[temperature][ii-1],Y[temperature][ii+1]]),slope*np.array([Y[temperature][ii-1],Y[temperature][ii+1]])+(X[temperature][ii]-slope*Y[temperature][ii]), linestyle = '-', linewidth = 0.5, color = colors_T[temperature]  )
        if X[temperature][-1] in rho_list:
            beta[temperature]['rho'].append(X[temperature][-1])
            beta[temperature]['P'].append(Y[temperature][-1])
            #backward finite difference
            slope = (X[temperature][-1]-X[temperature][-2])/(Y[temperature][-1]-Y[temperature][-2]) 
            beta[temperature]['beta'].append( -1/X[temperature][-1] * slope)
            ax2.plot(np.array([Y[temperature][-2],Y[temperature][-1]]),slope*np.array([Y[temperature][-2],Y[temperature][-1]])+(X[temperature][-1]-slope*Y[temperature][-1]), linestyle = '-', linewidth = 0.5, color = colors_T[temperature]  )
        print('******',temperature, '\n',beta[temperature]['rho'],'\n', ['{:1.1e}'.format(beta[temperature]['beta'][ii]) for ii in range(len(beta[temperature]['beta']))]   ) #isothermal compressibility in GPa-1
    #legend
    ncols = 1
    T_list = [str(int((float(label.strip('T'))*1000))) for label in natsort.natsorted(allT)]
    custom_patch = [mpatches.Patch(color=colors_T[key]) for key in allT]
    legend = ax2.legend([col for col in custom_patch],[str(label) for label in T_list],title = '$\\bf{Temperature}$ \n         (K)', bbox_to_anchor=(1.02, 1), loc="upper left", fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol=ncols)
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])

    
    #************ 3) compute isobaric expansivity: alpha = 1/rho drho/dT at each T-rho point used in steps 1) and 2)
    print('*************** Isobaric Expansivity (in K-1)')
    alpha = {}
    for temperature in beta:
        alpha[temperature] = {'rho':[],'alpha':[]}
        for ii in range(len(beta[temperature]['rho'])):
            alpha[temperature]['rho'].append(beta[temperature]['rho'][ii])
            alpha[temperature]['alpha'].append(beta[temperature]['beta'][ii] * TPC[beta[temperature]['rho'][ii]] )
        print('******',temperature, '\n',alpha[temperature]['rho'],'\n',['{:1.1e}'.format(alpha[temperature]['alpha'][ii]) for ii in range(len(alpha[temperature]['alpha']))]) #isobaric expansivity in K-1
        
    #************ 4) save results
    newfilename = 'thermoelasticity_coefficients_'+thermofile.split('.txt')[0].split('_')[1]+'.txt'
    with open(newfilename,'w') as f:
        f.write('rho\ttemperature\tP\talpha\tbeta\tTPC\n')
        f.write('(g/cm3)\t(K)\t(GPa)\t(K-1)\t(GPa-1)\t(GPa/K)\n')
        for temperature in beta:
            for ii in range(len(beta[temperature]['rho'])):
                f.write(str(beta[temperature]['rho'][ii])+'\t'+str(int(float(temperature.strip('T'))*1000))+'\t'+str(beta[temperature]['P'][ii])+'\t'+'{:1.6e}'.format(alpha[temperature]['alpha'][ii])+'\t'+'{:1.6e}'.format(beta[temperature]['beta'][ii])+'\t'+'{:1.6e}'.format(TPC[beta[temperature]['rho'][ii]])+'\n' )
    print(newfilename,'created')
    
    figurename = 'P-T_'+thermofile.split('.txt')[0].split('_')[1]+'.png'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
#    plt.show() # à mettre après savefig sinon on ne sauvegarde rien !! 
    
    figurename2 = 'rho-P_'+thermofile.split('.txt')[0].split('_')[1]+'.png'
    print("figure saved with name ",figurename2)
    fig2.savefig(figurename2, bbox_inches = 'tight', dpi = 300)   
    
    
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



