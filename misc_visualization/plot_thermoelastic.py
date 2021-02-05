#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the thermoelastic coefficient      ****
                   and compare with data from Spera2009 and Neilson2016
                      *************  ARTICLE VERSION  *************

"""


#     ********* Importation of the packages and modules used here *********
import glob
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import crystallography as cr
import re
import natsort





def creation_plot_2(plot_parameters,xvariable,letters):
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
    f, (ax1, ax2) = plt.subplots(1,2, sharex=True, sharey=False, 
                               figsize = (size_figure[0]*2.3,size_figure[1]))
    
    if xvariable == 'rho':
        major_xticks = np.arange(0, 5.0, 0.5) 
        minor_xticks = np.arange(0, 5.1, 0.1) 
        for ax in [ax1,ax2]:
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True)
        ax1.set_xlim(2.0,4.5)
        label = r'Density (g.cm$^{-3}$)'
    else:
        ax1.set_xlim(1,275)
        ax1.set_xscale('log')
        label = r'Pressure (GPa)'

    for ax in [ax1,ax2]:
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)                                         
        ax.yaxis.set_ticks_position('both')
        #ax.autoscale(axis='y',tight=True)
        ax.set_xlabel(label, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad/2)
    
    ax1.set_ylim(0,10)
    ax2.set_ylim(0,0.2)

    if letters != '':
        ax1.text(0.985,0.95, letters[0] , transform=ax1.transAxes, 
                 horizontalalignment = 'right', fontweight = 'bold', 
                 fontsize = plot_parameters["size_fonts"], 
                 bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
        ax2.text(0.985,0.95, letters[1] , transform=ax2.transAxes, 
                 horizontalalignment = 'right', fontweight = 'bold', 
                 fontsize = plot_parameters["size_fonts"], 
                 bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
    

    ax1.set_ylabel(r'$\alpha$ (.10$^{-5}$ K$^{-1}$)', fontweight = 'bold',
                  fontsize = size_fonts, labelpad = shift_labelpad)
    ax2.set_ylabel(r'$\beta$ (GPa$^{-1}$)', fontweight = 'bold',
                  fontsize = size_fonts, labelpad = shift_labelpad)    

    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.3)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False,
                    labelleft=False, left=False, labelright = False, right=False)
    return f, ax1, ax2, ax





def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),
                       "size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    #other dictionnaries and parameters for the figure
    legend_labels = {}
    colors_T = {'2000':'#800080','3000':'#297fff','4000':'#00ff00','4500':'#bae200',
                '5000':'#ffcd01','6000':'#ff0101','6500':'#ff00a2','7000':'#ff01de'}
    filename = 'all'
    filename2 = ''
    filename3 = ''
    letters = ''
    xvariable = 'rho'
    #other dictionnaries and parameters for the figure
    Temperatures = []
    try:
        options,arg = getopt.getopt(argv,"hf:g:j:v:l:",["filename","gfilename","jfilename",'variable','letters'])
    except getopt.GetoptError:
        print("plot_thermoelastic.py -f <thermoelasticfilename> -g <comparisonfile1> -j <comparisonfile2 (option)>  -v <variable x axis (rho or P)> -l <2 letters for subplots, ex. 'a,b', default ''> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_thermoelastic.py program to plot thermoelastic coefficients as a function of density or pressure for each T and compare with other data')
            print("plot_thermoelastic.py -f <thermoelasticfilename> -g <comparisonfile1> -j <comparisonfile2 (option)>  -v <variable x axis (rho or P)> -l <2 letters for subplots, ex. 'a,b', default ''>")
            print('')
            sys.exit()
        if opt in ('-f','--filename'):
            filename = str(arg)
        elif opt in ('-g','--gfilename'):
            filename2 = str(arg)
        elif opt in ('-j','--jfilename'):
            filename3 = str(arg)
        elif opt in ('-v','--variable'):
            xvariable = str(arg)
        elif opt in ('-l','--letters'):
            letters = arg.split(',')
    #************ Creation of the figure to plot (with correct number of lines)
    mineral=filename.split('_')[2].split('.txt')[0]
    print(mineral)
    figurename = 'thermoelastic_'+ mineral
    files = [filename]
    for file in [filename2,filename3]:
        if file != '':
            files.append(file)
            figurename += '-'+file.split('.txt')[0].split('_')[-1]
    figurename +='_'+xvariable    
    typeline = {filename : '-', filename2: '*', filename3: '+'}
    fig,ax1,ax2, ax0 = creation_plot_2(plot_parameters,xvariable,letters)
    print("axis are ax1",ax1,"ax2",ax2,"ax0",ax0)
    #************** Extract our data and plot 
    #*** Initialization of data dictionnaries
    alpha = {}  #big dictionnary with inside alpha  all coresponding to different T 
    beta = {} #big dictionnary with inside beta  all coresponding to different T 
    X = {}  #idem for rho or P  
    if xvariable == 'P':
        Xcol = 2
    else:
        Xcol = 0
    #*** Extract data
    with open(filename,'r') as f:
        [f.readline() for i in range(2)]
        while True:
            line = f.readline()
            if not line:break
            else:
                entry = line.split('\n')[0].split('\t')
                try:
                    alpha[entry[1]].append(-float(entry[3])*1E5)
                    beta[entry[1]].append(-float(entry[4]))
                    X[entry[1]].append(float(entry[Xcol]))
                except KeyError:
                    alpha[entry[1]] = [-float(entry[3])*1E5]
                    beta[entry[1]] = [-float(entry[4])]
                    X[entry[1]] = [float(entry[Xcol])]
                    
    #*** Plot data
    for T in X:
        if T == '2000' or T == '4500' or T == '6500' or T == '7000':
            pass
        else:
            print('************ for ',T,' K')
            Temperatures.append(T)
            #Sort data by X values
            Xlist, Ylistalpha, Ylistbeta = zip(*sorted(zip(X[T],alpha[T],beta[T]), key=lambda x: x[0], reverse=False))    
            print(Xlist)
            print(Ylistalpha)
            print(Ylistbeta)  
            #Plot
            ax1.plot(Xlist,Ylistalpha, typeline[filename],  markersize = plot_parameters["size_markers"],
                     color = colors_T[T], linewidth = plot_parameters["size_lines"])
            ax2.plot(Xlist,Ylistbeta, typeline[filename],  markersize = plot_parameters["size_markers"], 
                     color = colors_T[T], linewidth = plot_parameters["size_lines"])
    #************** Extract the collected data and plot
    mineral_sources = {'thermoelastic_Neilson2016.txt':'NaAlSi3O8','thermoelastic_Spera2009.txt':'CaAl2Si2O8'}
    sources = {'thermoelastic_Neilson2016.txt':'Neilson $et~al.$ (2016)','thermoelastic_Spera2009.txt':'Spera $et~al.$ (2009)'}
    for file in files[1:]:
        print('For file ',file)
        if mineral_sources[file] == mineral:
            #***** Initialization of data dictionnaries
            alpha = {}  #big dictionnary with inside alpha  all coresponding to different T 
            beta = {} #big dictionnary with inside beta  all coresponding to different T 
            X = {}  #idem for rho or P  
            with open(file,'r') as f:
                f.readline()
                while True:
                    line = f.readline()
                    if not line:break
                    else:
                        entry = line.split('\n')[0].split('\t')
                        try:
                            alpha[entry[1]].append(float(entry[4])*1E5)
                            beta[entry[1]].append(float(entry[5]))
                            if xvariable == 'P':
                                X[entry[1]].append(float(entry[Xcol]))
                            else:
                                X[entry[1]].append(float(entry[Xcol])/1000)
                        except KeyError:
                            alpha[entry[1]] = [float(entry[4])*1E5]
                            beta[entry[1]] = [float(entry[5])]
                            if xvariable == 'P':
                                X[entry[1]] = [float(entry[Xcol])]
                            else:
                                X[entry[1]] = [float(entry[Xcol])/1000]
            for T in X:
                if T in Temperatures:
                    print('************ for ',T,' K')
                    #Sort data by X values
                    Xlist, Ylistalpha, Ylistbeta = zip(*sorted(zip(X[T],alpha[T],beta[T]), key=lambda x: x[0], reverse=False))
                    print(Xlist)
                    print(Ylistalpha)
                    print(Ylistbeta)  
                    #Plot
                    ax1.plot(Xlist,Ylistalpha, typeline[file],  markersize = plot_parameters["size_markers"], 
                             color = colors_T[T], linewidth = plot_parameters["size_lines"])
                    ax2.plot(Xlist,Ylistbeta, typeline[file],  markersize = plot_parameters["size_markers"], 
                             color = colors_T[T], linewidth = plot_parameters["size_lines"])
    #************** Create legend from custom artist/label lists
    custom_patch = [mpatches.Patch(color=colors_T[key]) for key in natsort.natsorted(Temperatures)]
    legend = ax0.legend([col for col in custom_patch],[label for label in natsort.natsorted(Temperatures)],
                        title = '$\\bf{Temperature~(K)}$', bbox_to_anchor=(0.5, 1.12), 
                        loc="lower center", fontsize = plot_parameters["size_font_ticks"],
                        borderaxespad=0., ncol=len(Temperatures))    
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])
    for file in files[1:]: 
        legend_labels[sources[file]] =  plt.Line2D((0,1),(0,0), markeredgecolor='k', 
                     markerfacecolor = 'k', linestyle = 'None', 
                     markersize = plot_parameters["size_markers"],
                     marker=typeline[file][0])
    legend_labels['this study'] =  plt.Line2D((0,1),(0,0), color='k', 
                 linestyle=typeline[filename][:], linewidth = plot_parameters["size_lines"])
    s = [(k, legend_labels[k]) for k in sorted(legend_labels.keys(),reverse = False)]      
    ax0.legend([v for k,v in s],[k for k,v in s], bbox_to_anchor=(0.5, 1.02), 
               loc='lower center', ncol=len(legend_labels), 
               fontsize = plot_parameters["size_font_ticks"], borderaxespad=0.) 
    ax0.add_artist(legend)
    figurename = figurename + '.pdf'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












