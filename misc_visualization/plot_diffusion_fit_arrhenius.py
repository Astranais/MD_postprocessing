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
def fit_n_plot(P,temperature,data, ax, nf, element, colors_T):
    P = np.asarray(P) *1E-9 #converted to GPa
    data = np.asarray(data)
    #least squared linear regression giving the correlation coefficient
    slope, intercept, r_value, p_value, std_err = stats.linregress(P, np.log10(data)) #intercept = y(x=0), std_err = error on the slope
    R_squared = r_value**2
    nf.write(element+'\t'+temperature+'\t'+str(slope)+'\t'+str(intercept)+'\t'+str(round(R_squared,3))+'\n')
    #plot
    ax.semilogy(P,data, 'o', color=colors_T[temperature])
    ax.text(0.85,0.85, element, transform=ax.transAxes, horizontalalignment='left', fontweight = 'bold', fontsize = 12)
    X = np.array([0,max(P)])
    ax.semilogy(X, np.exp((slope * X + intercept)*np.log(10)), '--', color=colors_T[temperature])
    
def fit_n_plotbis(T,data, ax, nf, element, colors_elem):
    T = np.asarray(T)
    data = np.asarray(data)
    #least squared linear regression giving the correlation coefficient
    slope, intercept, r_value, p_value, std_err = stats.linregress(1/T, np.log10(data)) #intercept = y(x=0), std_err = error on the slope
    R_squared = r_value**2
    nf.write(element+'\t'+str(1E9)+'\t'+str(slope)+'\t'+str(intercept)+'\t'+str(round(R_squared,3))+'\n')
    
    #plot
    ax.semilogy(1/T,data, 'o', color=colors_elem[element])#, label = element)
    X = np.array([1/7000,1/2000])
    ax.semilogy(X, np.exp((slope * X + intercept)*np.log(10)), '--', color=colors_elem[element])
    return slope, intercept
    
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


def creation_plot_4(plot_parameters):
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
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, sharex=True, sharey=True, figsize = (size_figure[0]*2,size_figure[1]*2))
    major_xticks = np.arange(0, 450, 50) 
    minor_xticks = np.arange(0, 450, 10) 
    ax1.set_yscale('log')
    for ax in [ax1,ax2,ax3,ax4]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_yscale('log')                                          
    #        ax.set_yticks(major_yticks)
    #        ax.set_yticks(minor_yticks, minor=True)                                           
        ax.yaxis.set_ticks_position('both')
        plt.autoscale()
    #for ax in [ax1,ax2]:
    #    plt.setp(ax.get_xticklabels()[-1], visible=False) 
        #ax.grid(True, axis = 'x')
    ax1.set_xlim(0,120)
    ax1.set_ylim(4e-10,3e-7)


    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax0 = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    ax0.set_xlabel(r'Pressure (GPa)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    ax0.set_ylabel(r'Diffusion coefficient (m$^2$.s$^{-1}$)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*3+shift_labelpad/2)
    return f, ax1, ax2, ax3, ax4, ax0



def creation_plot(plot_parameters):
    """     ********** Creation of the plot D=f(T)  **********    """
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
#        ax.set_yticks(major_yticks)
#        ax.set_yticks(minor_yticks, minor=True)                                           
    ax.yaxis.set_ticks_position('both') 
    #ax.grid(True, axis = 'x')
    ax.set_xlim(1/9500,1/1050)
    ax.set_ylim(3e-13,3e-7)
    #plt.autoscale()
    
    ax.set_xlabel(r'1/T (K$^{-1}$)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    ax.set_ylabel(r'Diffusion coefficient (m$^2$.s$^{-1}$)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    return fig, ax

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
            ax.semilogy(Trev,D, 'D', markeredgecolor = 'k', markeredgewidth = 0.25, markerfacecolor=colors_elem[element])
            
            ax.semilogy(Trev,Dex, 'o', markeredgecolor=colors_elem[element], markerfacecolor = 'w', markeredgewidth = 1 )
            
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
                ax.semilogy(Trev,D[element], '*',  markeredgecolor = 'k', markeredgewidth = 0.25, markerfacecolor=colors_elem[element])
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
                ax.semilogy(Trev,D[element], 'P',  markeredgecolor = 'k', markeredgewidth = 0.25, markerfacecolor=colors_elem[element])
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
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    colors_T = {'2000':'#800080','3000':'#297fff','4000':'#00ff00','4500':'#bae200','5000':'#ffcd01','6000':'#ff0101','6500':'#ff00a2','7000':'#ff01de'}
    colors_elem = {'Al':'pink','C':'0.25','Ca':'c','H':'w','K':'m','Na':'b','O':'r','Si':'y'} 
    letter = ''
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
        column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,'Cvm':25,'stdev_Cvm':26}
    else:
        column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,'E':11,'stdev_E':12,'err_E':13,'Cvm_Nkb':14,'stdev_Cvm_Nkb':15}
    #**** Creation of plots and .txt files
    fig1, ax1,ax2,ax3, ax4, ax0 = creation_plot_4(plot_parameters)
    figurename1 = 'ArrheniusP_'+filename.split('/')[-1].split('_')[-1].split('.txt')[0]
    newfilename1  = 'ArrheniusP_'+filename.split('/')[-1].split('_')[-1]
    nf1 = open(newfilename1,'w')
    nf1.write('Arrhenius linear fit Y = aX+b with Y = log(D) in m2/s and X = P in Pa\nElement\tT(K)\ta\tb\tR2\n')
    figbis, axbis = creation_plot(plot_parameters)
    figurenamebis = 'ArrheniusT_'+filename.split('/')[-1].split('_')[-1].split('.txt')[0]    
    newfilenamebis  = 'ArrheniusT_'+filename.split('/')[-1].split('_')[-1]
    nfbis = open(newfilenamebis,'w')
    nfbis.write('Arrhenius linear fit Y = aX+b with Y = log(D) in m2/s and X = 1/T in K\nElement\tP(Pa)\ta\tb\tR2\n')
    #**** Extraction of data
    with open (filename, 'r') as f:
        line = f.readline()
        elements = line.split('\n')[0].split('\t')[1:]
        [f.readline() for i in range(3)]
        while True:
            line = f.readline()
            if not line: break
            else:
                entry=line.split('\n')[0].split('\t')
                diffusivities[entry[0].split('/')[-1].split('.outcar')[0]] = [float(entry[1]),float(entry[6]),float(entry[11]),float(entry[16])] #diffusivities for cation, Al, Si, and O
    with open (thermofilename, 'r') as tf:
        [tf.readline() for i in range(3)]
        while True:
            line = tf.readline()
            if not line: break
            else:
                entry=line.split('\n')[0].split('\t')
                thermo[entry[0].split('/')[-1].split('.outcar')[0]] = [float(entry[column_number['rho']]),float(entry[column_number['P']])*1E9,float(entry[column_number['T']])] #thermo data: rho in g/cm3, P in Pa, T in K
    #**** Classic Arrhenius fit: D = f(P) for each T
    sortedfiles = []
    for file in natsort.natsorted(diffusivities):
        sortedfiles.append(file)
    temperature0, acell = split_name(sortedfiles[0]) 
    P,data1,data2,data3,data4 = [],[],[],[],[]
    for file in sortedfiles:
        temperature, acell = split_name(file) 
        if temperature != temperature0 : 
            print('********** new T --> plot and fit the last one = ', temperature0)
            #sort data with increassing P
            P,data1,data2,data3,data4 = zip(*sorted(zip(P,data1,data2,data3,data4), key=lambda x: x[0]))
            #for i in range(len(P)):
            #    print(P[i],'\t','{:.2e}'.format(data1[i]),'\t','{:.2e}'.format(data2[i]),'\t','{:.2e}'.format(data3[i]),'\t','{:.2e}'.format(data4[i]))      
            #** linear fit, write and plot for each data set
            fit_n_plot(P,temperature0,data1, ax1, nf1, elements[0], colors_T)
            fit_n_plot(P,temperature0,data2, ax2, nf1, elements[1], colors_T)
            fit_n_plot(P,temperature0,data3, ax3, nf1, elements[2], colors_T)
            fit_n_plot(P,temperature0,data4, ax4, nf1, elements[3], colors_T)                
            #re-initialization
            P,data1,data2,data3,data4 = [],[],[],[],[]
            temperature0= temperature
        else:
            #** extraction of data: diffusivities and pressure for separated T
            if thermo[file][0] <= 1.4: 
                print(file)
                continue
            else:
                P.append(thermo[file][1])
                data1.append(diffusivities[file][0])
                data2.append(diffusivities[file][1])
                data3.append(diffusivities[file][2])
                data4.append(diffusivities[file][3])
    print('plot and fit the last T = ',temperature)
    #**** Approximate Arrhenius fit: D=f(T) for P ~ 1 GPa
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
    #sort data with increassing P
    T,data1,data2,data3,data4 = zip(*sorted(zip(T,data1,data2,data3,data4), key=lambda x: x[0]))
    #** linear fit, write and plot for each data set
    slope1, intercept1 = fit_n_plotbis(T,data1, axbis, nfbis, elements[0], colors_elem)
    slope2, intercept2 = fit_n_plotbis(T,data2, axbis, nfbis, elements[1], colors_elem)
    slope3, intercept3 = fit_n_plotbis(T,data3, axbis, nfbis, elements[2], colors_elem)
    slope4, intercept4 = fit_n_plotbis(T,data4, axbis, nfbis, elements[3], colors_elem)
        
    slopes = np.array([slope1,slope2,slope3,slope4])
    intercepts = np.array([intercept1,intercept2,intercept3,intercept4])
    
    handles, legendlabels =  [] ,[]
    #** calculation of the diffusivity at low T (from expe)
    #from Freda and Baker 1998 (K and Na feldspars)
    Dex = Freda_n_Baker_1998(nfbis, elements[0], slope1, intercept1)
    #** Plot expe data    
    mineral = thermofilename.split('_')[1]
    if mineral == 'NaAlSi3O8':
        #from Freda and Baker 1998
        handles, legendlabels = plot_Freda_n_Baker_1998(axbis,elements[0],Dex, colors_elem, handles, legendlabels)      
        #from simu Neilson 2016
        handles, legendlabels = plot_Neilson2016(axbis,elements, colors_elem, handles, legendlabels)
    elif mineral  == 'KAlSi3O8':
        #from Freda and Baker 1998
        handles, legendlabels = plot_Freda_n_Baker_1998(axbis,elements[0],Dex, colors_elem, handles, legendlabels)      
    elif mineral  == 'CaAl2Si2O8':
        #from simu Spera2009
        handles, legendlabels = plot_Spera2009(axbis,elements, colors_elem, handles, legendlabels)
    
    
    if letter != '':
        plt.text(0.006,0.945, letter , transform=axbis.transAxes, horizontalalignment = 'left', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
    #** legend
    plt.legend(handles, legendlabels, loc='best', fontsize = plot_parameters['size_fonts'])
    #**** Save plots
    #save_plots(fig1, figurename1+'.png')
    save_plots(figbis, figurenamebis+'.svg')
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












