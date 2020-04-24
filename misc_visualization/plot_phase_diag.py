#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the phase diagram in any projection     ****
               *************  ARTICLE VERSION  *************

"""


#     ********* Importation of the packages and modules used here *********
import glob
import sys
import getopt
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend import Legend
import re


def creation_plot_3(plot_parameters,  major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax, yvariable):
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
    
    
    #apply setup
    for ax in [ax1,ax2,ax3]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_yticks(major_yticks)                                        
        ax.set_yticks(minor_yticks, minor=True)
        ax.yaxis.set_ticks_position('both')
    
    #for ax in [ax2,ax3]:
    #    if yvariable != 'T':
    #        plt.setp(ax.get_yticklabels()[-1], visible=False)
    for ax in [ax1,ax2]:
        plt.setp(ax.get_xticklabels()[-1], visible=False)
    
    ax1.set_xlim(xxmin,xxmax)
    ax1.set_ylim(yymin,yymax)
#    plt.autoscale()
    
    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    ax.set_ylabel(ylabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*4)
    return f, ax1, ax2, ax3, ax


def creation_plot_2(plot_parameters,  major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax, yvariable):
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
    f, (ax1, ax2) = plt.subplots(2,1, sharex=True, sharey=True, figsize = (size_figure[0],size_figure[1]*2))  
    
    
    #apply setup
    for ax in [ax1,ax2]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_yticks(major_yticks)                                        
        ax.set_yticks(minor_yticks, minor=True)
        ax.yaxis.set_ticks_position('both')
    
    for ax in [ax2]:
        if yvariable != 'T':
            plt.setp(ax.get_yticklabels()[-1], visible=False)
    
    ax1.set_xlim(xxmin,xxmax)
    ax1.set_ylim(yymin,yymax)
#    plt.autoscale()
    
    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    ax.set_ylabel(ylabel, fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*4)
    return f, ax1, ax2, ax


def creation_plot(plot_parameters, major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax):
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
    return fig, ax


def append_value(data,entry,i):
    """      ******** append a float to the list, replacing empty str by NaN ****** """
    if entry[i] == '':
        data.append(float('nan'))
    elif entry[i] == '\n':
        data.append(float('nan'))
    else:
        data.append(int(entry[i]))
    return data


def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4)#(6,4)
    ,"size_markers" : 20,"size_lines" : 1,"shift_labelpad" : 10}
    #other dictionnaries and parameters for the figure
    filename = ''
    filename2 = ''
    init_letter = ''
    letters = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    #other dictionnaries and parameters for the figure
    compounds = []
    #location of the critical point
    critical_area = {'NaAlSi3O8':{'P':(0.1,0.2),'T':(6000,6500),'rho':(0.44,0.62)} , 'KAlSi3O8':{'P':(0,0.1),'T':(5000,5500),'rho':(0.56,0.87)} , 'CaAl2Si2O8':{'P':(0.1,0.2),'T':(7000,7500),'rho':(0.57,0.79)}  }
    try:
        options,arg = getopt.getopt(argv,"hf:g:t:l:v:",["filename","gfilename",'type', 'letter', "view"])
    except getopt.GetoptError:
        print("plot_phase_diag.py -f <filename> -g <filename 2 (if we want to plot 2 files)> -t <type of plot: y-x variables> -l <letter of first plot, default = ''> -v <view = zoom,classic,full>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_phase_diag.py program to plot the phase diagram ')
            print("plot_phase_diag.py -f <filename> -g <filename 2 (if we want to plot 2 files)> -t <type of plot: y-x variables, ex: T-rho or P-rho> -l <letter of first plot, default = ''> -v <view = zoom,classic,full>")
            print("plot_phase_diag.py requires to be lauched from the folder containing every phase_diag_compound.txt file")
            print('')
            sys.exit()
        if opt in ('-f','--filename'):
            filename = str(arg)
        elif opt in ('-g','--gfilename'):
            filename2 = str(arg)
        elif opt in ('-t','--type'):
            type_plot = str(arg)
            xvariable = type_plot.split('-')[1]
            yvariable = type_plot.split('-')[0]
        elif opt in ('-l','--letter'):
            init_letter = str(arg)
        elif opt in ('-v','--view'):
            view = str(arg)
    #************ tick definition
    tickdefinitions = {'P':{'zoom':(-1,1.6,0.5,0.1),'classic':(-3,6,1,0.5),'full':(0,550,50,25)} , 'rho':{'zoom':(0,5,0.5,0.1),'classic':(0,7,0.5,0.1),'full':(0,9,1,0.5)}, 'T':{'zoom':(0,9000,1000,500),'classic':(0,8000,1000,500),'full':(0,21000,5000,1000)} }   #begin/end/stepmajor/stepminor
    limitaxis =  {'P':{'zoom':(-1,1),'classic':(-3,5),'full':(0,500)} , 'rho':{'zoom':(0.3,2),'classic':(0.3,3),'full':(0.3,6.5)}, 'T':{'zoom':(1500,7700),'classic':(1500,7700),'full':(0,20500)} }  #begin/end
    axislabel =  {'P':r'Pressure (GPa)', 'rho':r'Density (g.cm$^{-3}$)', 'T':r'Temperature (K)' } 
    major_xticks = np.arange(tickdefinitions[xvariable][view][0], tickdefinitions[xvariable][view][1], tickdefinitions[xvariable][view][2]) 
    minor_xticks = np.arange(tickdefinitions[xvariable][view][0], tickdefinitions[xvariable][view][1], tickdefinitions[xvariable][view][3]) 
    (xxmin,xxmax) = (limitaxis[xvariable][view][0], limitaxis[xvariable][view][1])
    xlabel = axislabel[xvariable]
    major_yticks = np.arange(tickdefinitions[yvariable][view][0], tickdefinitions[yvariable][view][1], tickdefinitions[yvariable][view][2]) 
    minor_yticks = np.arange(tickdefinitions[yvariable][view][0], tickdefinitions[yvariable][view][1], tickdefinitions[yvariable][view][3]) 
    (yymin, yymax) = (limitaxis[yvariable][view][0], limitaxis[yvariable][view][1])
    ylabel = axislabel[yvariable]
    #************ initialization of the plot
    if filename == 'all':
        fig,ax1,ax2,ax3, ax0 = creation_plot_3(plot_parameters, major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax, yvariable)
        #print("axis are ax1",ax1,"ax2",ax2,"ax3",ax3,"ax0",ax0)
        files = sorted(glob.glob('phase_diag_*O8.txt'),reverse=True) #I list every phase diag files
        figurename = 'phase_diag_all'
    elif filename2 != '':
        fig,ax1,ax2, ax0 = creation_plot_2(plot_parameters, major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax, yvariable)
        print("axis are ax1",ax1,"ax2",ax2,"ax0",ax0)
        files = [filename, filename2] #we will plot 2 files
        figurename = filename.split('.txt')[0] +'_'+ filename2.split('.txt')[0].split('_')[2]
    else:
        fig, ax = creation_plot(plot_parameters, major_xticks,  minor_xticks, major_yticks,  minor_yticks, xlabel, ylabel, xxmin, xxmax, yymin, yymax)
        files = [filename] #I take only the file we want
        figurename = filename.split('.txt')[0]
    #************ creation of arrays and plot at the same time
    for file in files:
        print("****** For file",file)
        #********initialisations
        #**change of subplot
        if filename == 'all' or filename2 !='':
            if files.index(file) == 0:
                ax=ax1
                letter = init_letter                    
            if files.index(file) == 1:
                ax=ax2
                try:
                    #letter = letters[letters.index(init_letter)+3] #for all 3 projections
                    #letter = letters[letters.index(init_letter)+2] #for only 2 projections
                    letter = letters[letters.index(init_letter)+1] #for naming  from up to down
                except ValueError:
                    print("You haven't indicated any letter")
                    letter = ''
            if files.index(file) == 2:
                ax=ax3 
                try:
                    #letter = letters[letters.index(init_letter)+6] #for all 3 projections
                    #letter = letters[letters.index(init_letter)+4] #for only 2 projections
                    letter = letters[letters.index(init_letter)+2] #for naming  from up to down
                except ValueError:
                    print("You haven't indicated any letter")
            #print("I plot on axis",ax)
        #**extraction compound
        mineral = file.split('_')[2].split('.txt')[0]  
        compound = mineral
        #and formatage of label
        i=0
        while True:
            if i <= len(compound)-1:
                if re.match('[0-9]',compound[i]):
                    num=0
                    try:
                        while re.match('[0-9]',compound[i+num]):
                            num +=1
                    except IndexError: #index error when we arrive at the end of the cluster name
                        pass
                        #print('end of the cluster')
                    compound = compound[:i]+'$_{'+compound[i:i+num]+'}$' + compound[i+num:]
                    i = i+5
                i = i+1
            else:break
        compounds.append(compound)
        #**initialisation of data dicitonnaries for each column (marker type)
        not_started = {'rho':[],'T':[],'P':[]}
        ongoing = {'rho':[],'T':[],'P':[]}
        viscous_fluid = {'rho':[],'T':[],'P':[]}
        fluid = {'rho':[],'T':[],'P':[]}
        mix_liq_gaz = {'rho':[],'T':[],'P':[]}
        O2 = {'rho':[],'T':[],'P':[]}        
        liquid_spinodal = {'rho':[],'T':[],'P':[]}
        gas_spinodal = {'rho':[],'T':[],'P':[]}
        #********* fill the dictionnaries with data
        with open(file,'r') as f:
            line = f.readline()
            while True:
                line = f.readline()
                if not line: break
                else:    
                    entry=line.split('\t')
                    not_started['P'].append(float(entry[0]))
                    not_started['rho'].append(float(entry[1]))
                    not_started['T'] = append_value(not_started['T'],entry,2)
                    ongoing['P'].append(float(entry[0]))
                    ongoing['rho'].append(float(entry[1]))
                    ongoing['T'] = append_value(ongoing['T'],entry,3)
                    viscous_fluid['P'].append(float(entry[0]))
                    viscous_fluid['rho'].append(float(entry[1]))
                    viscous_fluid['T'] = append_value(viscous_fluid['T'],entry,4)
                    fluid['P'].append(float(entry[0]))
                    fluid['rho'].append(float(entry[1]))
                    fluid['T'] = append_value(fluid['T'],entry,5)
                    mix_liq_gaz['P'].append(float(entry[0]))
                    mix_liq_gaz['rho'].append(float(entry[1]))
                    mix_liq_gaz['T'] = append_value(mix_liq_gaz['T'],entry,6)
                    O2['P'].append(float(entry[0]))
                    O2['rho'].append(float(entry[1]))
                    O2['T'] = append_value(O2['T'],entry,7)
                    liquid_spinodal['P'].append(float(entry[0]))
                    liquid_spinodal['rho'].append(float(entry[1]))
                    liquid_spinodal['T'] = append_value(liquid_spinodal['T'],entry,8)
                    gas_spinodal['P'].append(float(entry[0]))
                    gas_spinodal['rho'].append(float(entry[1]))
                    gas_spinodal['T'] = append_value(gas_spinodal['T'],entry,9)
        #print('not_started',not_started)
        #print('fluid',fluid)
        #print('O2',O2)
        #print('liquid_spinodal',liquid_spinodal)
        
        
        #we remove data in case of P-rho plot
        if yvariable == 'P':
            dictionaries = [not_started,ongoing,viscous_fluid,fluid,mix_liq_gaz,O2,liquid_spinodal,gas_spinodal]
            for dicname in dictionaries:
                #print("*****************Change dicname")
                for ii in range(len(dicname['P'])):
                    #print(dicname['T'][ii])
                    if math.isnan(dicname['T'][ii]):
                        #print("NaN detected")
                        dicname['P'][ii] = float('nan')
        #create list for spinodal curve
        spinodal = {'P':[],'rho':[],'T':[]}
        for ii in range(len(gas_spinodal['T'])):
            if math.isnan(gas_spinodal['T'][ii]):
                continue
            else:
                spinodal[xvariable].append(gas_spinodal[xvariable][ii])
                spinodal[yvariable].append(gas_spinodal[yvariable][ii])
        for ii in range(len(liquid_spinodal['T'])):
            if math.isnan(liquid_spinodal['T'][ii]):
                continue
            else:
                spinodal[xvariable].append(liquid_spinodal[xvariable][ii])
                spinodal[yvariable].append(liquid_spinodal[yvariable][ii])
        #******* plot        
        rect = mpatches.Rectangle((critical_area[mineral][xvariable][0],critical_area[mineral][yvariable][0]), width = critical_area[mineral][xvariable][1]-critical_area[mineral][xvariable][0], height = critical_area[mineral][yvariable][1]-critical_area[mineral][yvariable][0], linewidth = None, fill = True, color = 'orange', joinstyle = 'round', zorder = -1 )
        ax.add_patch(rect)
        #ax.scatter(not_started[xvariable],not_started[yvariable],s=plot_parameters['size_markers'], facecolor = 'r', edgecolor='r', label = 'not started')
        #ax.scatter(ongoing[xvariable],ongoing[yvariable],s=plot_parameters['size_markers'], facecolor = 'orange', edgecolor='orange', label = 'ongoing')
        ax.scatter(fluid[xvariable],fluid[yvariable],s=plot_parameters['size_markers'], facecolor = 'b', edgecolor='b', label = 'fluid')
        ax.scatter(mix_liq_gaz[xvariable],mix_liq_gaz[yvariable],s=plot_parameters['size_markers'], facecolor = 'g', edgecolor='g', label = 'mix liquid + gaz')
        ax.scatter(viscous_fluid[xvariable],viscous_fluid[yvariable],marker = 'x', s=plot_parameters['size_markers']*2, facecolor = 'grey', edgecolor='grey', label = 'viscous fluid')
        #ax.scatter(O2[xvariable],O2[yvariable],s=plot_parameters['size_markers']*10, facecolor = 'none', edgecolor='r', label = 'O$_2$')
#        ax.scatter(liquid_spinodal[xvariable],liquid_spinodal[yvariable],s=plot_parameters['size_markers']*2, marker = 'P',  facecolor = 'k', edgecolor='k', label = 'liquid spinodal')
        ax.plot(spinodal[xvariable],spinodal[yvariable],linestyle = '--', linewidth =plot_parameters['size_lines']*2, color='k' )
        ax.plot(liquid_spinodal[xvariable],liquid_spinodal[yvariable],linestyle = '', linewidth =plot_parameters['size_lines']*2, color='k' , markersize=plot_parameters['size_markers']*0.5, marker = 'P',  markerfacecolor = 'k', markeredgecolor='k', label = 'liquid spinodal')
        ax.plot(gas_spinodal[xvariable],gas_spinodal[yvariable],linestyle = '', linewidth =plot_parameters['size_lines']*2, color='k' , markersize=plot_parameters['size_markers']*0.5, marker = 'P',  markerfacecolor = '0.90', markeredgecolor='k', label = 'liquid spinodal')
        

        
        
        #**text on plot        
        if yvariable == 'T' and xvariable == 'rho':
            positiontext = {'zoom':{'supercrit':(0.98,0.91),'liquid':(0.98,0.42),'liquidg':(0.9,3500),'+gas':(0.9,3200)} ,'classic':{'supercrit':(0.98,0.91),'liquid':(0.98,0.42),'liquidg':(0.7,3500),'+gas':(0.7,3200)}, 'full':{'supercrit':(0.4,0.91),'liquid':(0.5,0.21)}   }
            ax.text(positiontext[view]['supercrit'][0],positiontext[view]['supercrit'][1], "Supercritical" , transform=ax.transAxes, horizontalalignment = 'right', fontsize = plot_parameters["size_fonts"])
            ax.text(positiontext[view]['liquid'][0],positiontext[view]['liquid'][1], "Liquid" , transform=ax.transAxes, horizontalalignment = 'right', fontsize = plot_parameters["size_fonts"])
            if view != 'full':
                ax.text(positiontext[view]['liquidg'][0],positiontext[view]['liquidg'][1], "Liquid" ,  horizontalalignment = 'left', fontsize = plot_parameters["size_fonts"])
                ax.text(positiontext[view]['+gas'][0],positiontext[view]['+gas'][1], "+ Gas" , horizontalalignment = 'left', fontsize = plot_parameters["size_fonts"])                
        elif yvariable == 'P' and xvariable == 'rho':
            positiontext = {'zoom':{'supercrit':(0.4,0.85),'liquidg':(0.4,0.25),'+gas':(0.4,0.2)} ,'classic':{'supercrit':(0.4,0.85),'liquidg':(0.6,0.25),'+gas':(0.6,0.2)}, 'full':{'supercrit':(0.4,0.55)}   }
            ax.text(positiontext[view]['supercrit'][0],positiontext[view]['supercrit'][1], "Supercritical" , transform=ax.transAxes, horizontalalignment = 'right', fontsize = plot_parameters["size_fonts"])
            if view != 'full':
                ax.text(positiontext[view]['liquidg'][0],positiontext[view]['liquidg'][1], "Liquid" ,  transform=ax.transAxes, horizontalalignment = 'left', fontsize = plot_parameters["size_fonts"])
                ax.text(positiontext[view]['+gas'][0],positiontext[view]['+gas'][1], "+ Gas" , transform=ax.transAxes, horizontalalignment = 'left', fontsize = plot_parameters["size_fonts"])        
        elif yvariable == 'P' and xvariable == 'T':
            positiontext = {'zoom':{'supercrit':(7150,0.4),'liquid':(4000,0.6)} ,'classic':{'supercrit':(7150,1),'liquid':(4000,1)}, 'full':{'supercrit':(8000,10)}   }
            ax.text(positiontext[view]['supercrit'][0],positiontext[view]['supercrit'][1], "Supercritical" , horizontalalignment = 'left', rotation = 90, fontsize = plot_parameters["size_fonts"])
            if view != 'full':
                ax.text(positiontext[view]['liquid'][0],positiontext[view]['liquid'][1], "Liquid" , horizontalalignment = 'left', fontsize = plot_parameters["size_fonts"])
        if (filename == 'all' or filename2 != '') and letter != '':
            ax.text(0.02,0.95, letter , transform=ax.transAxes, horizontalalignment = 'left', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='w', edgecolor='k', pad=3.0))
        elif init_letter != '':
            ax.text(0.02,0.94, init_letter , transform=ax.transAxes, horizontalalignment = 'left', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='w', edgecolor='k', pad=3.0))
    #Create legend from custom artist/label lists  
#    legend_labels['not started'] = plt.Line2D((0,1),(0,0),  markersize = plot_parameters["size_markers"], marker='o', facecolor = 'r', edgecolor='r', linestyle='')
#    s = [(k, legend_labels[k]) for k in sorted(legend_labels.keys(),reverse = False)]
#    if filename == 'all':
#        ax=ax0
#    plt.legend([v for k,v in s],[k for k,v in s],loc='upper center', bbox_to_anchor=(0.46, 1.14),fancybox=True, fontsize = plot_parameters["size_fonts"], ncol=len(style_lines))
    figurename = figurename+'_'+yvariable+'-'+xvariable +'_'+ view + '.pdf'
    print("figure saved with name ",figurename)
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    #plt.show()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])












