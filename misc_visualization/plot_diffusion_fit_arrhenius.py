#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Arrhénius fit of self diffusion coefficient       ****
                and create .txt file with arhenius fit D=f(P) and D=f(T)
                
                ******* ARTICLE VERSION *****


"""


#     ********* Importation of the packages and modules used here *********
import sys
import getopt
import numpy as np
from scipy import stats
import natsort
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend import Legend
from matplotlib.lines import Line2D #useful to create a custom legend





    

#____                             ____________                             ____ 
#    '.                         .'            '.                         .' 
#      '.                     .'     FITS       '.                     .'
#        '.                 .'     & PLOTS        '.                 .'
#          '.             .'                        '.             .'
#            '._________.'                            '._________.'
def fit_n_plot(P,temperature,data, ax, nf, element, colors_T,data_P):
    print('************** ', temperature, element)
    P = np.asarray(P) * 1E-9 #converted to GPa
    data = np.asarray(data)
    #least squared linear regression giving the correlation coefficient
    slope, intercept, r_value, p_value, std_err = stats.linregress(P, np.log10(data)) #intercept = y(x=0), std_err = error on the slope
    R_squared = r_value**2
    #compute diffusivities at 25, 50 and 100 GPa from this fit
    diff25 = np.exp((slope * 25 + intercept)*np.log(10))
    diff50 = np.exp((slope * 50 + intercept)*np.log(10))
    diff100 = np.exp((slope * 100 + intercept)*np.log(10))
    nf.write(element+'\t'+temperature+'\t'+str(slope)+'\t'+str(intercept)+'\t'+str(round(R_squared,3))+'\t'+str(diff25)+'\t'+str(diff50)+'\t'+str(diff100)+'\n')
    #complete dictionnaries for future use in fit ln(D)=f(1/T)
    data_P[element][25].append(diff25)
    data_P[element][50].append(diff50)
    data_P[element][100].append(diff100)
    #plot
    ax.plot(P,data, 'o', color=colors_T[temperature])
    X = np.array([0.5,max(P)])
    ax.plot(X, np.exp((slope * X + intercept)*np.log(10)), '--', color=colors_T[temperature])
    ax.plot([25,50,100],[diff25,diff50,diff100], ls = '', marker = 'P', markersize=7,
            markeredgewidth = 1, markeredgecolor = 'k', markerfacecolor=colors_T[temperature] )
    
def fit_n_plotbis(T,P_val,data, ax, nf, element, color):
    T = np.asarray(T)
    data = np.asarray(data)
    #least squared linear regression giving the correlation coefficient
    slope, intercept, r_value, p_value, std_err = stats.linregress(1/T, np.log10(data)) #intercept = y(x=0), std_err = error on the slope
    R_squared = r_value**2
    nf.write(element+'\t'+str(P_val)+'\t'+str(slope)+'\t'+str(intercept)+'\t'+str(round(R_squared,3))+'\n')
    
    #plot
    ax.plot(1/T,data, 'o', color=color)#, label = element)
    X = np.array([1/7000,1/2000])
    ax.plot(X, np.exp((slope * X + intercept)*np.log(10)), '--', color=color)
    return slope, intercept

def plot_data_transp(P, temperature, data, ax, colors_T, plot_parameters):
    P = np.asarray(P) *1E-9 #converted to GPa
    data = np.asarray(data)
    ax.plot(P,data, linestyle = '', marker='o', markeredgewidth = 0, 
            markerfacecolor = colors_T[temperature]+'7f',
            markersize = 7)
    
#____                             ____________                             ____ 
#    '.                         .'            '.                         .' 
#      '.                     .'   CREATION     '.                     .'
#        '.                 .'     PLOTS          '.                 .'
#          '.             .'                        '.             .'
#            '._________.'                            '._________.'
def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.outcar.msd.dat
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.strip('./')
    temperature = str(int(float(filename.split('_')[1].strip('T'))*1000))
    acell = filename.split('.outcar.msd.dat')[0].split('_')[3].strip('a')
    return temperature, acell

def save_plots(fig, figurename):
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 150)


def creation_plot_4(plot_parameters,elements):
    """     ********** Creation of the plot D=f(P) **********    """
    print("I use creation plot 4")
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    #plot
    plt.close()
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, sharex=True, sharey=True, 
       figsize = (size_figure[0]*2,size_figure[1]*2))
    axis = [ax1,ax2,ax3,ax4]
    major_xticks = np.arange(0, 450, 50) 
    minor_xticks = np.arange(0, 450, 10) 
    for ax in [ax1,ax2,ax3,ax4]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_yscale('log')                                          
        ax.yaxis.set_ticks_position('both')
        plt.autoscale()
        ax.text(0.85,0.85, elements[axis.index(ax)], transform=ax.transAxes,
                                    horizontalalignment='left', fontweight = 'bold', fontsize = 12)

    #for ax in [ax1,ax2]:
    #    plt.setp(ax.get_xticklabels()[-1], visible=False) 
        #ax.grid(True, axis = 'x')

    ax1.set_xlim(0.5,210) 
    ax1.set_ylim(9e-12,1e-7)


    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax0 = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax0.set_xlabel(r'Pressure (GPa)', fontweight = 'bold', fontsize = size_fonts, 
                   labelpad = shift_labelpad*2)
    ax0.set_ylabel(r'Diffusion coefficient (m$^2$.s$^{-1}$)', fontweight = 'bold', 
                   fontsize = size_fonts, labelpad = shift_labelpad*3+shift_labelpad/2)
    return f, ax1, ax2, ax3, ax4, ax0



def creation_plot(plot_parameters):
    """     ********** Creation of the plot D=f(1/T)  **********    """
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
    fig, ax = plt.subplots(1,1, figsize = (size_figure[0]*2,size_figure[1]), linewidth = size_lines /2)
    Tlabels = np.array([1073,1473,1673,1873,2000,3000,4000,4500,5000,6000,7000])
    #creation of custom xticks labels
    labels = []
    for i in range(len(Tlabels)):
        labels.append(r'1/'+str(int(Tlabels[i])))
    plt.xticks(1/Tlabels, labels, rotation = 45)
    
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
    ax.set_yscale('log')                                                                                    
    ax.yaxis.set_ticks_position('both') 
    ax.set_xlim(1/9500,1/1050)
    ax.set_ylim(3e-13,3e-7)
    
    ax.set_xlabel(r'1/T (K$^{-1}$)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    ax.set_ylabel(r'Diffusion coefficient (m$^2$.s$^{-1}$)', fontweight = 'bold', 
                  fontsize = size_fonts, labelpad = shift_labelpad)
    return fig, ax

def creation_plot_4_bis(plot_parameters,elements):
    """     ********** Creation of the plot D=f(1/T) **********    """
    print("I use creation plot 4")
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    #plot
    plt.close()
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, sharex=True, sharey=True, 
       figsize = (size_figure[0]*2,size_figure[1]*2))
    axis = [ax1,ax2,ax3,ax4]
    #creation of custom xticks labels
    Tlabels = np.array([3000,4000,4500,5000,6000,7000])
    labels = []
    for i in range(len(Tlabels)):
        labels.append(r'1/'+str(int(Tlabels[i])))    
    for ax in axis:
        plt.xticks(1/Tlabels, labels)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.tick_params(axis = 'x',  rotation = 45)
        ax.set_yscale('log')                                                                                  
        ax.yaxis.set_ticks_position('both') 
        ax.text(0.85,0.85, elements[axis.index(ax)], transform=ax.transAxes, 
                                    horizontalalignment='left', fontweight = 'bold', fontsize = 12)


    ax.set_xlim(1/7500,1/2900)
    ax.set_ylim(9e-12,1e-7)
    
    
 
    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax0 = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax0.set_xlabel(r'1/T (K$^{-1}$)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*4)
    ax0.set_ylabel(r'Diffusion coefficient (m$^2$.s$^{-1}$)', fontweight = 'bold', 
                   fontsize = size_fonts, labelpad = shift_labelpad*3+shift_labelpad/2)
    return f, ax1, ax2, ax3, ax4, ax0

    
#____                             ____________                             ____ 
#    '.                         .'            '.                         .' 
#      '.                     .'    COMPUTE     '.                     .'
#        '.                 .'     FROM EXPE      '.                 .'
#          '.             .'                        '.             .'
#            '._________.'                            '._________.'
def Freda_n_Baker_1998(nfbis, element, slope1, intercept1):
    """ calculation of the diffusivity at low T from Freda and Baker 1998"""
    expeT = []
    try:
        with open('Expe_Freda1998.txt', 'r') as f:
            [f.readline() for i in range(2)]
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry = line.split('\n')[0].split('\t')
                    expeT.append(float(entry[2])+273)
        #[800,1200,1300,1400,1492,1600]) +  273
        nfbis.write('\nEstimation of diffusivities from fit at expe temperatures\n***** Expe Freda and Baker 1998\n')
        nfbis.write('elements\t'+'\t'.join(str(int(x)) for x in expeT)+'\n')
        expeT = np.asarray(expeT)
        Dex = np.zeros(len(expeT))
        for i in range(len(expeT)):
            Dex[i] = np.exp((slope1 * 1/expeT[i] + intercept1)*np.log(10))
        nfbis.write(element+'\t'+'\t'.join('{:.2e}'.format(d) for d in Dex)+'\n')
        return Dex
    except FileNotFoundError:
        print("File Freda and Baker 1998 not found")
        return 0
    

#____                             ____________                             ____ 
#    '.                         .'            '.                         .' 
#      '.                     .'     PLOT       '.                     .'
#        '.                 .'    EXPE DATA       '.                 .'
#          '.             .'                        '.             .'
#            '._________.'                            '._________.'
def plot_Freda_n_Baker_1998(ax,element, Dex, colors_elem, handles, legendlabels):
        columns_Freda1998 = {'Na':[9,10,11],'K':[18,19,20]}
        T, D, err = [], [] , []
        try:
            with open('Expe_Freda1998.txt', 'r') as f:
                [f.readline() for i in range(2)]
                while True:
                    line = f.readline()
                    if not line: break
                    else:
                        entry = line.split('\n')[0].split('\t')
                        T.append(float(entry[2])+273)
                        D.append(float(entry[columns_Freda1998[element][0]])*float(entry[columns_Freda1998[element][2]]))
                        err.append(float(entry[columns_Freda1998[element][1]]))
            #print(D)
            #print(err)
            Trev = 1/np.asarray(T)
            #Trev = Trev.tolist()
            D = np.asarray(D)
            ax.semilogy(Trev,D, 'D', markeredgecolor = 'k', markeredgewidth = 0.25, 
                        markerfacecolor=colors_elem[element])
            
            ax.semilogy(Trev,Dex, 'o', markeredgecolor=colors_elem[element],
                        markerfacecolor = 'w', markeredgewidth = 1 )
            
            #creation of custom xticks labels
            locations = ax.get_xticks()
            labels = [item.get_text() for item in ax.get_xticklabels()] #retrieve the previous labels
            newlabels = ['1/1873','','1/1673','','1/1473','1/1073']
            locations = np.append(locations,Trev)
            labels = labels+newlabels
            plt.xticks(locations, labels, rotation = 45)
            
            #add to legend
            handles.append(Line2D([0],[0], color = 'k', ls = '', marker = 'D',  markeredgewidth = 0.25))
            legendlabels.append('Freda and Baker, 1998')
            return handles, legendlabels
        except FileNotFoundError:
            print("File Freda and Baker 1998 not found")


def plot_Neilson2016(ax,elements, colors_elem, handles, legendlabels):
    columns_Neilson2016 = {'Na':8, 'Al':9, 'Si':10, 'O':11}
    Tlabel, T, D = [], [] , {'Na':[],'Al':[],'Si':[],'O':[]}
    try:
        with open('Simu_Neilson2016.txt', 'r') as f:
            [f.readline() for i in range(2)]
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry = line.split('\n')[0].split('\t')
                    P = round(float(entry[5]),0)
                    if P == 1.0:
                        Tlabel.append(float(entry[1]))
                        T.append(float(entry[3]))
                        for element in elements:
                            try:
                                D[element].append(float(entry[columns_Neilson2016[element]]))
                            except KeyError:
                                continue
        Trev = 1/np.asarray(T)
        for element in elements:
            try:
                D[element] = np.asarray(D[element])
                ax.semilogy(Trev,D[element], '*',  markeredgecolor = 'k', markeredgewidth = 0.25, 
                            markerfacecolor=colors_elem[element])
            except KeyError:
                continue
        
        #add to legend
        handles.append(Line2D([0],[0], color = 'k', ls = '', marker = '*',  markeredgewidth = 0.25))
        legendlabels.append('Neilson $et$ $al.$, 2016')
        return handles, legendlabels
    except FileNotFoundError:
        print("File Neilson2016 not found")   


def plot_Spera2009(ax,elements, colors_elem, handles, legendlabels):
    columns_Spera2009 = {'Ca':8, 'Al':9, 'Si':10, 'O':11}
    Tlabel, T, D = [], [] , {'Ca':[],'Al':[],'Si':[],'O':[]}
    try:
        with open('Simu_Spera2009.txt', 'r') as f:
            f.readline()
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry = line.split('\n')[0].split('\t')
                    P = round(float(entry[5]),0)
                    if P == 1.0:
                        Tlabel.append(float(entry[1]))
                        T.append(float(entry[3]))
                        for element in elements:
                            try:
                                D[element].append(float(entry[columns_Spera2009[element]]))
                            except KeyError:
                                continue
        Trev = 1/np.asarray(T)
        for element in elements:
            try:
                D[element] = np.asarray(D[element])
                ax.semilogy(Trev,D[element], 'P',  markeredgecolor = 'k', markeredgewidth = 0.25, 
                            markerfacecolor=colors_elem[element])
            except KeyError:
                continue
        
        #add to legend
        handles.append(Line2D([0],[0], color = 'k', ls = '', marker = 'P',  markeredgewidth = 0.25))
        legendlabels.append('Spera $et$ $al.$, 2009')
        return handles, legendlabels
    except FileNotFoundError:
        print("File Spera2009 not found")    
     
    
#____                             ____________                             ____ 
#    '.                         .'            '.                         .' 
#      '.                     .'     MAIN       '.                     .'
#        '.                 .'     PROGRAM        '.                 .'
#          '.             .'                        '.             .'
#            '._________.'                            '._________.'
def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),
                       "size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    colors_T = {'2000':'#800080','3000':'#297fff','4000':'#00ff00','4500':'#bae200',
                '5000':'#ffcd01','6000':'#ff0101','6500':'#ff00a2','7000':'#ff01de'}
    colors_elem = {'Al':'palevioletred','C':'0.25','Ca':'c','H':'w','K':'m','Na':'b','O':'r','Si':'y'} 
    colors_P = {1:'#2c1a4c',25:'#284073',50:'#3367a9',100:'#bdb7f1'}
    letter = ''
    Temperatures = []
    #variables for the data storage
    diffusivities = {} #dictionnary containing the diffusivities for each file
    thermo = {} #dictionnary containing the thermo data for each file
    try:
        options,arg = getopt.getopt(argv,"hf:g:t:l:",["filename","gthermofilename", "type", "letter"])
    except getopt.GetoptError:
        print("plot_fit_arrhenius.py -f <diffusion filename> -g <thermo filename>  -t <type of the file: 'short' or 'all'> -l <letter for article, default = ''>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_fit_arrhenius.py program to fit Arrhenius as a function of P for each T and as a function of T for 1 GPa (approximation)')
            print("plot_fit_arrhenius.py -f <diffusion filename> -g <thermo filename> -t <type of the file: 'short' or 'all'> -l <letter for article, default = ''>")
            print('')
            sys.exit()
        if opt in ('-f','--filename'):
            filename = str(arg)
        elif opt in ('-g','--gthermofilename'):
            thermofilename = str(arg)
        elif opt in ('-t','--type'):
            filetype = str(arg)
        elif opt in ('-l','--letter'):
            letter = str(arg)
    #************ initialization of the column number
    if filetype == 'all':
        column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,
                         'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,
                         'Cvm':25,'stdev_Cvm':26}
    else:
        column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,
                         'E':11,'stdev_E':12,'err_E':13,'Cvm_Nkb':14,'stdev_Cvm_Nkb':15}
    
    #**** Extraction of data
    with open (filename, 'r') as f:
        line = f.readline()
        elements = line.split('\n')[0].split('\t')[1:]
        [f.readline() for i in range(5)]
        while True:
            line = f.readline()
            if not line: break
            else:
                entry=line.split('\n')[0].split('\t')
                diffusivities[entry[0].split('/')[-1].split('.outcar')[0]] = [float(entry[1]),float(entry[7]),float(entry[13]),float(entry[19])] #diffusivities for cation, Al, Si, and O
    with open (thermofilename, 'r') as tf:
        [tf.readline() for i in range(3)]
        while True:
            line = tf.readline()
            if not line: break
            else:
                entry=line.split('\n')[0].split('\t')
                thermo[entry[0].split('/')[-1].split('.outcar')[0]] = [float(entry[column_number['rho']]),float(entry[column_number['P']])*1E9,float(entry[column_number['T']])] #thermo data: rho in g/cm3, P in Pa, T in K
    
    
    
    #****** Classic Arrhenius fit: ln(D) = f(P) for each T
    #** Creation of plots and .txt file
    fig1, ax1,ax2,ax3, ax4, ax0 = creation_plot_4(plot_parameters,elements)
    figurename1 = 'ArrheniusP_'+filename.split('/')[-1].split('_')[-1].split('.txt')[0]
    newfilename1  = 'ArrheniusP_'+filename.split('/')[-1].split('_')[-1]
    nf1 = open(newfilename1,'w')
    nf1.write('Arrhenius linear fit Y = aX+b with Y = log(D) in log(m2/s) and X = P in GPa and estimation of D for 3 P from these fits\nElement\tT(K)\ta(GPa-1)\tb(nounit)\tR2\tD_25GPa\tD_50GPa\tD_100GPa\n')
    #** We fit only between 5 and 150 GPa (in the approximately linear part)
    sortedfiles = []
    for file in natsort.natsorted(diffusivities):
        sortedfiles.append(file)
    temperature0, acell = split_name(sortedfiles[0]) 
    P,data1,data2,data3,data4 = [],[],[],[],[]
    T_P = []
    data_P = {elements[0]:{25:[],50:[],100:[]},
              elements[1]:{25:[],50:[],100:[]},
              elements[2]:{25:[],50:[],100:[]},
              elements[3]:{25:[],50:[],100:[]}}
    for file in sortedfiles:
        temperature, acell = split_name(file) 
        if temperature != temperature0 : 
            print('********** new T --> plot and fit the last one = ', temperature0)
            Temperatures.append(temperature0)
            try:
                #sort data with increassing P
                P,data1,data2,data3,data4 = zip(*sorted(zip(P,data1,data2,data3,data4), key=lambda x: x[0]))
                #for i in range(len(P)):
                #    print(P[i],'\t','{:.2e}'.format(data1[i]),'\t','{:.2e}'.format(data2[i]),'\t','{:.2e}'.format(data3[i]),'\t','{:.2e}'.format(data4[i]))      
                #** linear fit, write and plot for each data set
                fit_n_plot(P,temperature0,data1, ax1, nf1, elements[0], colors_T,data_P)
                fit_n_plot(P,temperature0,data2, ax2, nf1, elements[1], colors_T,data_P)
                fit_n_plot(P,temperature0,data3, ax3, nf1, elements[2], colors_T,data_P)
                fit_n_plot(P,temperature0,data4, ax4, nf1, elements[3], colors_T,data_P)
                T_P.append(float(temperature0))
            except ValueError:
                pass
            #re-initialization
            P,data1,data2,data3,data4 = [],[],[],[],[]
            temperature0= temperature
        else:
            #** extraction of data: diffusivities and pressure for separated T
            if thermo[file][1] < 25E9:  
                print(file)
                continue
            elif thermo[file][1] > 150E9:  
                print(file)
                continue
            else:
                P.append(thermo[file][1])
                data1.append(diffusivities[file][0])
                data2.append(diffusivities[file][1])
                data3.append(diffusivities[file][2])
                data4.append(diffusivities[file][3])
    print('plot and fit the last T = ',temperature)
    Temperatures.append(temperature)
    try:
        #sort data with increassing P
        P,data1,data2,data3,data4 = zip(*sorted(zip(P,data1,data2,data3,data4), key=lambda x: x[0]))
        #for i in range(len(P)):
        #    print(P[i],'\t','{:.2e}'.format(data1[i]),'\t','{:.2e}'.format(data2[i]),'\t','{:.2e}'.format(data3[i]),'\t','{:.2e}'.format(data4[i]))      
        #** linear fit, write and plot for each data set
        fit_n_plot(P,temperature,data1, ax1, nf1, elements[0], colors_T)
        fit_n_plot(P,temperature,data2, ax2, nf1, elements[1], colors_T)
        fit_n_plot(P,temperature,data3, ax3, nf1, elements[2], colors_T)
        fit_n_plot(P,temperature,data4, ax4, nf1, elements[3], colors_T)
        T_P.append(float(temperature0))
    except ValueError:
        pass
    #** We plot with transparency the other data
    temperature0, acell = split_name(sortedfiles[0]) 
    P,data1,data2,data3,data4 = [],[],[],[],[]
    for file in sortedfiles:
        temperature, acell = split_name(file) 
        if temperature != temperature0 : 
            #sort data with increassing P
            P,data1,data2,data3,data4 = zip(*sorted(zip(P,data1,data2,data3,data4), key=lambda x: x[0]))
            #** plot for each data set
            plot_data_transp(P, temperature0, data1, ax1, colors_T, plot_parameters)
            plot_data_transp(P, temperature0, data2, ax2, colors_T, plot_parameters)
            plot_data_transp(P, temperature0, data3, ax3, colors_T, plot_parameters)            
            plot_data_transp(P, temperature0, data4, ax4, colors_T, plot_parameters)
            #re-initialization
            P,data1,data2,data3,data4 = [],[],[],[],[]
            temperature0= temperature
        else:
            #** extraction of data: diffusivities and pressure for separated T
            P.append(thermo[file][1])
            data1.append(diffusivities[file][0])
            data2.append(diffusivities[file][1])
            data3.append(diffusivities[file][2])
            data4.append(diffusivities[file][3])
    #** plot last data set
    plot_data_transp(P, temperature, data1, ax1, colors_T, plot_parameters)
    plot_data_transp(P, temperature, data2, ax2, colors_T, plot_parameters)
    plot_data_transp(P, temperature, data3, ax3, colors_T, plot_parameters)            
    plot_data_transp(P, temperature, data4, ax4, colors_T, plot_parameters)
    
    
    
    #**** Approximate Arrhenius fit: ln(D)=f(1/T) for P ~ 1 GPa
    #** Creation of plots and .txt file
    figbis, axbis = creation_plot(plot_parameters)
    figurenamebis = 'ArrheniusT_'+filename.split('/')[-1].split('_')[-1].split('.txt')[0]+'_1GPa'   
    newfilenamebis  = 'ArrheniusT_'+filename.split('/')[-1].split('_')[-1]
    nfbis = open(newfilenamebis,'w')
    nfbis.write('Arrhenius linear fit Y = aX+b with Y = log(D) in log(m2/s) and X = 1/T in K-1\nElement\tP(GPa)\ta(K)\tb(nounit)\tR2\n')
    if letter != '':
        plt.text(0.006,0.945, letter , transform=axbis.transAxes, horizontalalignment = 'left',
                 fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                 bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
    #** extraction of data: diffusivities and temperatures for 1 GPa
    T,data1,data2,data3,data4 = [],[],[],[],[]
    for file in sortedfiles:
        rounded = round(thermo[file][1],-9)
        if rounded == 1.0E9:
            T.append(thermo[file][2])
            data1.append(diffusivities[file][0])
            data2.append(diffusivities[file][1])
            data3.append(diffusivities[file][2])
            data4.append(diffusivities[file][3])
    #sort data with increassing T
    T,data1,data2,data3,data4 = zip(*sorted(zip(T,data1,data2,data3,data4), key=lambda x: x[0]))
    #** linear fit, write and plot for each data set
    P_val = 1
    slope1, intercept1 = fit_n_plotbis(T,P_val,data1, axbis, nfbis, elements[0], colors_elem[elements[0]])
    slope2, intercept2 = fit_n_plotbis(T,P_val,data2, axbis, nfbis, elements[1], colors_elem[elements[1]])
    slope3, intercept3 = fit_n_plotbis(T,P_val,data3, axbis, nfbis, elements[2], colors_elem[elements[2]])
    slope4, intercept4 = fit_n_plotbis(T,P_val,data4, axbis, nfbis, elements[3], colors_elem[elements[3]])
        
    slopes = np.array([slope1,slope2,slope3,slope4])
    intercepts = np.array([intercept1,intercept2,intercept3,intercept4])
    
    
    
    #**** Arrhenius fit: ln(D)=f(1/T) for P ~ 25, 50 and 100 GPa obtained with the fit ln(D) = f(P)
    #** Creation of plots and .txt file
    fig, ax1b,ax2b,ax3b,ax4b,ax0b = creation_plot_4_bis(plot_parameters,elements)
    figurename = 'ArrheniusT_'+filename.split('/')[-1].split('_')[-1].split('.txt')[0]
    P_val = 1
    T = np.asarray(T)
    ax1b.plot(1/T, np.exp((slope1 * 1/T + intercept1)*np.log(10)), '--', color=colors_P[P_val])
    ax2b.plot(1/T, np.exp((slope2 * 1/T + intercept2)*np.log(10)), '--', color=colors_P[P_val])
    ax3b.plot(1/T, np.exp((slope3 * 1/T + intercept3)*np.log(10)), '--', color=colors_P[P_val])
    ax4b.plot(1/T, np.exp((slope4 * 1/T + intercept4)*np.log(10)), '--', color=colors_P[P_val])
    ax1b.plot(1/T, data1, 'o', color=colors_P[P_val])
    ax2b.plot(1/T, data2, 'o', color=colors_P[P_val])
    ax3b.plot(1/T, data3, 'o', color=colors_P[P_val])
    ax4b.plot(1/T, data4, 'o', color=colors_P[P_val])
    for P_val in [25,50,100]:
        #sort data with increassing T
        T,data1,data2,data3,data4 = [],[],[],[],[]
        T,data1,data2,data3,data4 = zip(*sorted(zip(T_P,data_P[elements[0]][P_val],
                                                    data_P[elements[1]][P_val],
                                                    data_P[elements[2]][P_val],
                                                    data_P[elements[3]][P_val]), key=lambda x: x[0]))
        #** linear fit, write and plot for each data set
        slope1, intercept1 = fit_n_plotbis(T,P_val,data1, ax1b, nfbis, elements[0], colors_P[P_val])
        slope2, intercept2 = fit_n_plotbis(T,P_val,data2, ax2b, nfbis, elements[1], colors_P[P_val])
        slope3, intercept3 = fit_n_plotbis(T,P_val,data3, ax3b, nfbis, elements[2], colors_P[P_val])
        slope4, intercept4 = fit_n_plotbis(T,P_val,data4, ax4b, nfbis, elements[3], colors_P[P_val])
    
    
    
    #**** Calculation of the diffusivity at low T (from expe)
    handles, legendlabels =  [] ,[]
    #from Freda and Baker 1998 (K and Na feldspars)
    Dex = Freda_n_Baker_1998(nfbis, elements[0], slopes[0], intercepts[0])
    #** Plot expe data    
    mineral = thermofilename.split('_')[1]
    if mineral == 'NaAlSi3O8':
        #from Freda and Baker 1998
        handles, legendlabels = plot_Freda_n_Baker_1998(axbis,elements[0],Dex, 
                                                        colors_elem, handles, legendlabels)      
        #from simu Neilson 2016
        handles, legendlabels = plot_Neilson2016(axbis,elements, colors_elem, 
                                                 handles, legendlabels)
    elif mineral  == 'KAlSi3O8':
        #from Freda and Baker 1998
        handles, legendlabels = plot_Freda_n_Baker_1998(axbis,elements[0],Dex,
                                                        colors_elem, handles, legendlabels)      
        #from simu Neilson 2016
        handles, legendlabels = plot_Neilson2016(axbis,elements, colors_elem, 
                                                 handles, legendlabels)        
    elif mineral  == 'CaAl2Si2O8':
        #from simu Spera2009
        handles, legendlabels = plot_Spera2009(axbis,elements, colors_elem, 
                                               handles, legendlabels)
            
    
    
    #********* Legend               
    #legend for plot D=f(P)
    Temperatures = list(set(Temperatures)) #get elements only once in the list
    custom_patch = [mpatches.Patch(color=colors_T[key]) for key in natsort.natsorted(Temperatures)]
    legend1 = ax0.legend([col for col in custom_patch],[label for label in natsort.natsorted(Temperatures)],
                         title = '$\\bf{Temperature~(K)}$', bbox_to_anchor=(0.5, 1.02), 
                         loc="lower center", fontsize = plot_parameters["size_font_ticks"], 
                         borderaxespad=0., ncol=len(Temperatures))
    plt.setp(legend1.get_title(),fontsize= plot_parameters["size_font_ticks"])
    
    #legend for plot D=f(1/T)
    custom_patch2 = [mpatches.Patch(color=colors_P[key]) for key in natsort.natsorted(colors_P)]
    legend2 = ax0b.legend([col for col in custom_patch2],[label for label in natsort.natsorted(colors_P)],
                          title = '$\\bf{Pressures~(GPa)}$', bbox_to_anchor=(0.5, 1.02), 
                          loc="lower center", fontsize = plot_parameters["size_font_ticks"], 
                          borderaxespad=0., ncol=len(Temperatures))
    plt.setp(legend2.get_title(),fontsize= plot_parameters["size_font_ticks"])
    
    #legend for plot D=f(1/T) at P=1GPa
    axbis.legend(handles, legendlabels, loc='best', fontsize = plot_parameters['size_fonts'])
    
    
    
    #**** Save plots
    save_plots(fig1, figurename1+'.pdf')
    save_plots(figbis, figurenamebis+'.pdf')
    save_plots(fig, figurename+'.pdf')
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












