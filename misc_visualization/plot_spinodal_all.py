#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the spinodal diagrams      ****
               *************  ARTICLE VERSION  *************

"""


#     ********* Importation of the packages and modules used here *********
import glob
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import natsort
from matplotlib.lines import Line2D
from matplotlib.colors import to_hex
import re


def creation_plot_3_horizontal(plot_parameters, xvariable, yvariable):
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
    f, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, sharey=True, figsize = (size_figure[0]*2.8,size_figure[1]))
     #tick definition
    if xvariable =='rho':
        major_xticks = np.arange(0, 2.6, 0.5) 
        minor_xticks = np.arange(0, 2.1, 0.1)    
        (xxmin,xxmax) = (0.25,2)
        xlabel = r'Density (g.cm$^{-3}$)'
    elif  xvariable =='T':
        major_xticks = np.arange(1000, 8000, 1000) 
        minor_xticks = np.arange(1000, 8000, 500)
        (xxmin, xxmax) = (1500, 7500)
        xlabel = r'Temperature (K)'
    else:
        major_xticks = np.arange(0,1,1)
        minor_xticks = np.arange(0,1, 0.5)    
        (xxmin,xxmax) = (0,1)
        xlabel = 'Not defined yet'
        
    if yvariable == 'T':
        major_yticks = np.arange(1000, 8000, 1000) 
        minor_yticks = np.arange(1000, 8000, 500)
        (yymin, yymax) = (1500, 7700)
        ylabel = r'Temperature (K)'
    elif yvariable == 'P':
        major_yticks = np.arange(-3, 1.5, 0.5) 
        minor_yticks = np.arange(-3, 1.1, 0.1)
        (yymin, yymax) = (-1,1)
        ylabel = r'Pressure (GPa)'
    else:
        major_yticks = np.arange(0, 1, 1) 
        minor_yticks = np.arange(0, 1, 0.5) 
        (yymin, yymax) = (0, 1)
        ylabel = 'Not defined yet'
    
    for ax in [ax1,ax2,ax3]:
        #apply setup
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)                                            
        ax.set_yticks(major_yticks)
        ax.set_yticks(minor_yticks, minor=True)
        ax.yaxis.set_ticks_position('both') 
        ax.set_xlim(xxmin,xxmax)
        ax.set_ylim(yymin,yymax)
        ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad/2)

    #for ax in [ax1,ax2]:
    #    plt.setp(ax.get_xticklabels()[-1], visible=False)     

    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax0 = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
        
    
    ax0.set_ylabel(ylabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*4)
               
    return f, ax1, ax2, ax3, ax0


def creation_plot(plot_parameters, xvariable, yvariable):
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
    fig, ax = plt.subplots(1,1, figsize = size_figure, linewidth = size_lines /2)
    
    #tick definition
    if xvariable =='rho':
        major_xticks = np.arange(0, 2.6, 0.5) 
        minor_xticks = np.arange(0, 2.1, 0.1)    
        (xxmin,xxmax) = (0.25,2)#1.8)
        xlabel = r'Density (g.cm$^{-3}$)'
    elif  xvariable =='T':
        major_xticks = np.arange(1000, 8000, 1000) 
        minor_xticks = np.arange(1000, 8000, 500)
        (xxmin, xxmax) = (1500, 7500)
        xlabel = r'Temperature (K)'
    else:
        major_xticks = np.arange(0,1,1)
        minor_xticks = np.arange(0,1, 0.5)    
        (xxmin,xxmax) = (0,1)
        xlabel = 'Not defined yet'
        
    if yvariable == 'T':
        major_yticks = np.arange(1000, 8000, 1000) 
        minor_yticks = np.arange(1000, 8000, 500)
        (yymin, yymax) = (1500, 7500)
        ylabel = r'Temperature (K)'
    elif yvariable == 'P':
        major_yticks = np.arange(-3, 1, 0.5) 
        minor_yticks = np.arange(-3, 1, 0.1)
        (yymin, yymax) = (-2.5,0.6)
        ylabel = r'Pressure (GPa)'
    else:
        major_yticks = np.arange(0, 1, 1) 
        minor_yticks = np.arange(0, 1, 0.5) 
        (yymin, yymax) = (0, 1)
        ylabel = 'Not defined yet'
    
    
    #apply setup
    ax.set_xticks(major_xticks)
    ax.set_xticks(minor_xticks, minor=True)
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)                                            
    ax.set_yticks(major_yticks)
    ax.set_yticks(minor_yticks, minor=True)
    ax.yaxis.set_ticks_position('both') 
    #ax.grid(True, axis = 'x')
    ax.set_xlim(xxmin,xxmax)
    ax.set_ylim(yymin,yymax)
#    plt.autoscale()
    
    ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    ax.set_ylabel(ylabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad)
    ax0 = ax
    return fig, ax, ax0

def format1label(label):
    """formatage of compound label """
    i=0
    while True:
        if i <= len(label)-1:
            if re.match('[0-9]',label[i]):
                num=0
                try:
                    while re.match('[0-9]',label[i+1+num]):
                        num +=1
                except IndexError: #index error when we arrive at the end of the cluster name
                    pass
                    #print('end of the cluster')
                label = label[:i]+'$_{'+label[i:i+1+num]+'}$' + label[i+1+num:]
                i = i+5
            i = i+1
        else:break
    return label

def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    #plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),"size_markers" : 6,"size_lines" : 1,"shift_labelpad" : 10}
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),
                       "size_markers" : 8,"size_lines" : 2,"shift_labelpad" : 10}
    letters = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    letter = ''
    mineralname = 'all'
    critical = {'NaAlSi3O8':{'P':(0.03,0.14),'T':(5500,6000),'rho':(0.51,0.83)},
                'KAlSi3O8':{'P':(0.02,0.11),'T':(5000,5500),'rho':(0.56,0.87)},
                'CaAl2Si2O8':{'P':(0.11,0.21),'T':(7000,7500),'rho':(0.57,0.79)}  }
    # with constrained fit for Na
    critical_constrained = {'NaAlSi3O8':{'P':(0.1,0.2),'T':(6000,6500),'rho':(0.44,0.62)}, 
            'KAlSi3O8':{'P':(0.1,0.15),'T':(5500,6000),'rho':(0.36,0.64)} , 
            'CaAl2Si2O8':{'P':(0.1,0.17),'T':(7000,7500),'rho':(0.34,0.77)}  }
    try:
        options,arg = getopt.getopt(argv,"hf:t:l:",["filename1",'type', 'letter'])
    except getopt.GetoptError:
        print("plot_spinodal_all.py -f <mineralname or 'all'> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_spinodal_all.py program to plot the spinodal ')
            print("plot_spinodal_all.py -f <mineralname or 'all'> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> ")
            print('')
            sys.exit()
        if opt in ('-f','--filename1'):
            mineralname = str(arg)
        elif opt in ('-t','--type'):
            type_plot = str(arg)
            xvariable = type_plot.split('-')[1]
            yvariable = type_plot.split('-')[0]
        elif opt in ('-l','--letter'):
            letter = str(arg)
    if mineralname == 'all':
        minerals = ['KAlSi3O8','NaAlSi3O8','CaAl2Si2O8']
        figurename = 'spinodal_all'
        fig,ax1,ax2,ax3, ax0 = creation_plot_3_horizontal(plot_parameters,xvariable, yvariable)
    else:
        minerals = [mineralname]        
        figurename = 'spinodal_'+mineralname
        fig, ax1, ax0 = creation_plot(plot_parameters, xvariable, yvariable) 
    #************ creation of arrays and plot at the same time
    for ii in range(len(minerals)):
        mineralname = minerals[ii]
        print('*************** for mineralname', mineralname)
        files = sorted(glob.glob('spinodal_*_'+mineralname+'_*.txt')) #I list only the files we want
        print('files are', files)
        #** choice of subplot
        if ii == 0:
            ax = ax1
        elif ii == 1:
            ax = ax2
        elif ii == 2:    
            ax = ax3
        for file in files:
            print('***** for file', file)
            #**initialisation of data dicitonnaries for each column (marker type)  
            data = {'T' : [], 'rho':[],'P':[]}
            #********* fill the dictionnaries with data
            with open(file,'r') as f:
                f.readline()
                while True:
                    line = f.readline()
                    if not line: break
                    else:    
                        entry=line.split('\n')[0].split('\t')
                        data['T'].append(float(entry[0].strip('T'))*1000)
                        data['rho'].append(float(entry[1]))
                        data['P'].append(float(entry[2]))
            print(data['T'])
            print(data['rho'])
            #******* change of markers and colors style
            filename = file.split('_')
            if filename[1] == 'constrained':
                bordercolor =  '#ffc900' #'#ff0000'
                critical_area = critical_constrained
            else:
                bordercolor = '#000000'
                critical_area = critical
            if filename[-2] == 'soft':
                style = 'P-'
                fillcolor = bordercolor
            else:
                style = 'P:'
                fillcolor = '0.90'
            #******* plot the spinodal lines     
            ax.plot(data[xvariable],data[yvariable], style, markersize = plot_parameters["size_markers"]+3,
                    markeredgecolor = bordercolor, markeredgewidth = 1.5, markerfacecolor = fillcolor, 
                    linewidth = plot_parameters["size_lines"], color = bordercolor)
            #******* plot the critical area
            rect = mpatches.Rectangle((critical_area[mineralname][xvariable][0],
                               critical_area[mineralname][yvariable][0]), 
                width = critical_area[mineralname][xvariable][1]-critical_area[mineralname][xvariable][0], 
                height = critical_area[mineralname][yvariable][1]-critical_area[mineralname][yvariable][0], 
                linewidth = None, fill = True, color = bordercolor+'7f', joinstyle = 'round', zorder = -99 )
            ax.add_patch(rect)
               
        #text on plot
        if letter != '':
            ax.text(0.985,0.945,  letters[letters.index(letter)+ii] , transform=ax.transAxes, 
                                          horizontalalignment = 'right', verticalalignment = 'baseline',
                                          fontweight = 'bold', fontsize = plot_parameters["size_fonts"],
                                          bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
        else:
            ax.text(0.985,0.945, format1label(mineralname) , transform=ax.transAxes, 
                    horizontalalignment = 'right',verticalalignment = 'baseline',  
                    fontweight = 'bold', fontsize = plot_parameters["size_fonts"], 
                    bbox=dict(facecolor='none', edgecolor='k', pad=3.0))    
    
    
    #************ legend
    styles = ['-',':'] #low density dataset/ dataset for rho > 1
    fillcolors = ['k','0.90'] 
    labels1 = [r'low density dataset',r'$\rho$ > 1 g.cm$^{-1}$ dataset']
    custom_lines = [Line2D([0],[0], markeredgecolor = 'k', markeredgewidth = 1.5, 
                           markerfacecolor = fillcolors[ii], color = 'k', ls = styles[ii], 
                           marker = 'P', linewidth = plot_parameters["size_lines"], 
                           markersize = plot_parameters["size_markers"]+3) for ii in range(len(styles))]
    legend = ax0.legend([col for col in custom_lines],[label for label in labels1], 
                        bbox_to_anchor=(0.5, 1.09), loc="upper center", 
                        fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol=2)

    bordercolors = ['k','#ffc900']#'#ff0000']#not constrained/constrained fit
    labels2 = [r'ax$^3$+bx$^2$+cx+d fit', r'ax$^3$+bx$^2$+cx fit']
    custom_patch = [mpatches.Patch(color=key) for key in bordercolors]
    legend2 = ax0.legend([col for col in custom_patch],[label for label in labels2], 
                         bbox_to_anchor=(0.5, 1.09), loc="lower center", 
                         fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol=2)
    
    
    ax0.add_artist(legend)

    
    
    figurename = figurename+'_'+yvariable+'-'+xvariable + '.pdf'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












