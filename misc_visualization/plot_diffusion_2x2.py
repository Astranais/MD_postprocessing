#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the self diffusion coefficient of 4 different atoms       ****


"""


#     ********* Importation of the packages and modules used here *********
import glob
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
import crystallography as cr
import re



#____                             ____________                             ____ 
#    '.                         .'            '.                         .' 
#      '.                     .'     PLOT       '.                     .'
#        '.                 .'    EXPE DATA       '.                 .'
#          '.             .'                        '.             .'
#            '._________.'                            '._________.'
def plot_deKoker2010(ax1,ax2,ax3,ax4, handles, legendlabels):
    columns_deKoker2010 = {'Ca':8, 'Al':9, 'Si':10, 'O':11}
    rho, D = {'3000':[],'4000':[],'6000':[]} , {'Ca':{'3000':[],'4000':[],'6000':[]},'Al':{'3000':[],'4000':[],'6000':[]},'Si':{'3000':[],'4000':[],'6000':[]},'O':{'3000':[],'4000':[],'6000':[]}}
    colors_T_gray = {'3000':'0','4000':'0.2','5000':'0.5','6000':'0.85'} 
    try:
        with open('Simu_deKoker2010.txt', 'r') as f:
            f.readline()
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry = line.split('\n')[0].split('\t')
                    T=entry[1]
                    for elem in D:
                        D[elem][T].append(float(entry[columns_deKoker2010[elem]]))
                    rho[T].append(float(entry[2])/1000)
        for elem in D:
            #choice of axis
            if elem == 'O':
                ax = ax4
            elif elem == 'Si':
                ax = ax3
            elif elem == 'Al':
                ax = ax2
            else:
                ax = ax1
            for temp in D[elem]:
                ax.plot(rho[temp],D[elem][temp], '-', color = colors_T_gray[temp], linewidth = 2)      
        #add to legend
        handles.append(plt.Line2D([0],[0], color = 'k', ls = '-'))
        legendlabels.append('deKoker, 2010')
        return handles, legendlabels
    except FileNotFoundError:
        print("File deKoker2010 not found")            
    
def plot_Neilson2016(ax1,ax2,ax3,ax4, handles, legendlabels):
    columns_Neilson2016 = {'Na':8, 'Al':9, 'Si':10, 'O':11}
    rho, D = {'3000':[],'4000':[],'5000':[]} , {'Na':{'3000':[],'4000':[],'5000':[]},'Al':{'3000':[],'4000':[],'5000':[]},'Si':{'3000':[],'4000':[],'5000':[]},'O':{'3000':[],'4000':[],'5000':[]}}
    colors_T_gray = {'3000':'0','4000':'0.2','5000':'0.5','6000':'0.85'}
    try:
        with open('Simu_Neilson2016.txt', 'r') as f:
            [f.readline() for i in range(2)]
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry = line.split('\n')[0].split('\t')
                    T=entry[1]
                    if T in D['Na']:
                        for elem in D:
                            D[elem][T].append(float(entry[columns_Neilson2016[elem]]))
                        rho[T].append(float(entry[2])/1000)
        for elem in D:
            #choice of axis
            if elem == 'O':
                ax = ax4
            elif elem == 'Si':
                ax = ax3
            elif elem == 'Al':
                ax = ax2
            else:
                ax = ax1
            for temp in D[elem]:
                ax.plot(rho[temp],D[elem][temp], ':', color = colors_T_gray[temp], linewidth = 2)      
        #add to legend
        handles.append(plt.Line2D([0],[0], color = 'k', ls = ':'))
        legendlabels.append('Neilson $et$ $al.$, 2016')
        return handles, legendlabels
    except FileNotFoundError:
        print("File Neilson2016 not found")   

def plot_Spera2009(ax1,ax2,ax3,ax4, handles, legendlabels):
    columns_Spera2009 = {'Ca':8, 'Al':9, 'Si':10, 'O':11}
    rho, D = {'4000':[],'5000':[],'6000':[]} , {'Ca':{'4000':[],'5000':[],'6000':[]},'Al':{'4000':[],'5000':[],'6000':[]},'Si':{'4000':[],'5000':[],'6000':[]},'O':{'4000':[],'5000':[],'6000':[]}}
    colors_T_gray = {'3000':'0','4000':'0.2','5000':'0.5','6000':'0.85'}
    try:
        with open('Simu_Spera2009.txt', 'r') as f:
            f.readline()
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry = line.split('\n')[0].split('\t')
                    T=entry[1]
                    if T in D['Ca']:
                        for elem in D:
                            D[elem][T].append(float(entry[columns_Spera2009[elem]]))
                        rho[T].append(float(entry[2])/1000)
        for elem in D:
            #choice of axis
            if elem == 'O':
                ax = ax4
            elif elem == 'Si':
                ax = ax3
            elif elem == 'Al':
                ax = ax2
            else:
                ax = ax1
            for temp in D[elem]:
                ax.plot(rho[temp],D[elem][temp], '-', color = colors_T_gray[temp], linewidth = 2)
        #add to legend
        handles.append(plt.Line2D([0],[0], color = 'k', ls = '-'))
        legendlabels.append('Spera $et$ $al.$, 2009')
        return handles, legendlabels
    except FileNotFoundError:
        print("File Spera2009 not found")            

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
    temperature = filename.split('_')[1]
    acell = filename.split('.outcar.msd.dat')[0].split('_')[3].strip('a')
    mineral = filename.split('_')[0]
    return mineral, temperature, acell


def creation_plot_2x2(size_fonts,size_font_ticks,size_figure,shift_labelpad, size_lines):
    """     ********** Creation of the plot  **********    """
    plt.close()
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, sharex=True, sharey=True, figsize = size_figure)
    major_xticks = np.arange(0, 4.5, 0.5) 
    minor_xticks = np.arange(0, 4.1, 0.1) 
    for ax in [ax1,ax2,ax3,ax4]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = size_font_ticks, width = size_lines/2)
        ax.set_yscale('log')                                          
    #        ax.set_yticks(major_yticks)
    #        ax.set_yticks(minor_yticks, minor=True)                                           
        ax.yaxis.set_ticks_position('both')
    for ax in [ax3]:
        plt.setp(ax.get_xticklabels()[-1], visible=False) 
        #ax.grid(True, axis = 'x')
    ax1.set_xlim(1.0,4.1)
    ax1.set_ylim(4e-10,3e-7)

    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    ax.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*2)
    ax.set_ylabel(r'Diffusion coefficient (m$^2$.s$^{-1}$)', fontweight = 'bold', fontsize = size_fonts, labelpad = shift_labelpad*3+shift_labelpad/2)
    return f, ax1, ax2, ax3, ax4, ax


#____                             ____________                             ____ 
#    '.                         .'            '.                         .' 
#      '.                     .'     MAIN       '.                     .'
#        '.                 .'     PROGRAM        '.                 .'
#          '.             .'                        '.             .'
#            '._________.'                            '._________.'
def main(argv):
    """     ********* Main program *********     """
    #parameters for the figures depending on the output format (presentation or article)
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,8),"size_markers" : 5,"size_lines" : 1,"shift_labelpad" : 10}
    #other dictionnaries and parameters for the figure
    style_markers = {'CaAl2Si2O8':'o','KAlSi3O8':'s','NaAlSi3O8':'o'}
    style_lines = {'CaAl2Si2O8':'-','KAlSi3O8':'--','NaAlSi3O8':':'}
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200','T5':'#ffcd01','T6':'#ff0101','T6.5':'#ff00a2','T7':'#ff01de'}
    filename = 'all'
    ions = []
    compounds = []
    Temperatures = []
    #variables nedded fot the plot
    data = {} #dictionnary containing the data for each element, initialized for each T
    stdev = {} #same for stdev
    #other parameters
    Na=6.022*10**23
    try:
        options,arg = getopt.getopt(argv,"hf:t:",["filename","typefigure"])
    except getopt.GetoptError:
        print("plot_diffusion.py -f <filename>(default = all)")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_diffusion.py program to plot diffusion coefficient as a function of density for each T and the selected minerals containing 4 elements')
            print("plot_diffusion.py -f <filename>(default = all files of every compound)")
            print("plot_diffusion.py requires to be lauched from the folder containing every diffusivities file created by the script analyze_msd")
            print('')
            sys.exit()
        if opt in ('-f','--filename'):
            filename = str(arg)
    fig,ax1,ax2,ax3,ax4, ax0 = creation_plot_2x2(plot_parameters['size_fonts'],plot_parameters['size_font_ticks'],plot_parameters['size_figure'],plot_parameters['shift_labelpad'],plot_parameters['size_lines'])
    if filename == 'all':
        files = sorted(glob.glob('diffusivities*.txt'),reverse=True) #I list every diffusivities files
        figurename = 'diffusivities_2x2_all'
    else:
        files = [filename] #I take only the file we want
        figurename = filename.split('.txt')[0]+'_2x2'
    #** Plot expe data  
    handles, legendlabels =  [] ,[]
    handles, legendlabels = plot_deKoker2010(ax1,ax2,ax3,ax4, handles, legendlabels)
    handles, legendlabels = plot_Neilson2016(ax1,ax2,ax3,ax4, handles, legendlabels)
    #handles, legendlabels = plot_Spera2009(ax1,ax2,ax3,ax4, handles, legendlabels)
    #** Plot our data
    for file in files:
        #**extraction compound
        compound = file.split('_')[1].split('.txt')[0]  #we need the compound for the legend
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
        #**creation of elements and number lists and initialization of T
        with open(file,'r') as f:
            line = f.readline()
            entry=line.split()
            elements = entry[1:]
            ions.append(elements[0])
            line = f.readline()
            entry=line.split()
            number = entry[1:]
            [f.readline() for i in range(2)]
            line = f.readline()
            entry=line.split()
            mineral, temperature0, acell0 = split_name(entry[0])
        #**calculation of M*N nedded for the calculation of densities
        MN = 0
        for i in range(len(elements)):
            MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
        #**initialisation of data
        rho = [] #initialization of the array containing the densities
        for i in range(len(elements)):
            data[elements[i]] = []
            stdev[elements[i]] = []
        #**creation of arrays and plot at the same time
        with open(file,'r') as f:
            [f.readline() for i in range(4)]
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry=line.split()
                    mineral, temperature, acell = split_name(entry[0]) 
                    if temperature != temperature0:     #if we change T:
                        #print('**************** for ',temperature0)
                        #print('data is:',data)
                        #we add this T to the list used for the legend
                        Temperatures.append(temperature0)
                        #attribution of fill color
                        if mineral == 'NaAlSi3O8':
                            fillcolor = 'w'
                        elif mineral == 'CaAl2Si2O8':
                            fillcolor = colors_T[temperature0]+'ff'
                        elif mineral == 'KAlSi3O8':
                            fillcolor = colors_T[temperature0]+'7f'
                        #we plot the data
                        for elem in elements:
                            #choice of axis
                            if elem == 'O':
                                ax = ax4
                            elif elem == 'Si':
                                ax = ax3
                            elif elem == 'Al':
                                ax = ax2
                            else:
                                ax = ax1
                            #print("we plot for",temperature0, "and ",elem)
                            ax.plot(rho,data[elem], style_markers[mineral]+style_lines[mineral], markersize = plot_parameters["size_markers"], markeredgewidth = 0.5, color = colors_T[temperature0], markerfacecolor = fillcolor,  markeredgecolor = colors_T[temperature0], linewidth = plot_parameters["size_lines"])
                            #ax.errorbar(rho,data[elem], yerr=stdev[elem], fmt=style_markers[mineral], markersize = plot_parameters[size_markers], markeredgewidth = 0.5, edgecolor colors_T[temperature0], facecolor = fillcolor, linestyle = style_lines[mineral],  linewidth = plot_parameters["size_lines"] )
                        #we re-initialize the arrays
                        temperature0 = temperature
                        for i in range(len(elements)):
                            data[elements[i]] = []
                            stdev[elements[i]] = [] 
                        rho = []
                    if (temperature == 'T2') : continue
                    rho.append(MN/(Na*float(acell)**3*10**(-24)) )    #calculation density
                    for i in range(len(elements)):
                        data[elements[i]].append(float(entry[i*5+1]))
                        stdev[elements[i]].append(float(entry[i*5+2]))
            #we plot the last set of data (last T)
            #we add this T to the list used for the legend
            Temperatures.append(temperature0)
            #attribution of fill color
            if mineral == 'NaAlSi3O8':
                fillcolor = 'w'
            elif mineral == 'CaAl2Si2O8':
                fillcolor = colors_T[temperature0]+'ff'
            elif mineral == 'KAlSi3O8':
                fillcolor = colors_T[temperature0]+'7f'
            #we plot the data
            for elem in elements:
                #choice of axis
                if elem == 'O':
                    ax = ax4
                elif elem == 'Si':
                    ax = ax3
                elif elem == 'Al':
                    ax = ax2
                else:
                    ax = ax1
                ax.plot(rho,data[elem], style_markers[mineral]+style_lines[mineral], markersize = plot_parameters["size_markers"],  markeredgewidth = 0.5, color = colors_T[temperature0], markerfacecolor = fillcolor,  markeredgecolor = colors_T[temperature0], linewidth = plot_parameters["size_lines"])
                #ax.errorbar(rho,data[elem], yerr=stdev[elem], fmt=style_markers[mineral], markersize = plot_parameters[size_markers],  edgecolor colors_T[temperature0], facecolor = fillcolor, linestyle = style_lines[mineral],  linewidth = plot_parameters["size_lines"] ) 
    
    #********* Legend               
    #Create legend from custom artist/label lists
    #definition of line and marker style
    fillcolor = {}
    for mineral in style_markers:
        if mineral == 'CaAl2Si2O8':
            fillcolor[mineral] = '#000000ff'
        elif mineral == 'KAlSi3O8':
            fillcolor[mineral] = '#0000007f'
        else:
            fillcolor[mineral] = 'w'
    custom_lines = [plt.Line2D([0],[0], marker = style_markers[key], linestyle = style_lines[key], markeredgecolor = 'k', markerfacecolor = fillcolor[key],  color = 'k', markersize = plot_parameters["size_markers"], markeredgewidth = 0.5, linewidth = plot_parameters["size_lines"]) for key in sorted(style_markers)]
    ax0.legend([line for line in custom_lines],[label for label in sorted(compounds)], bbox_to_anchor=(0.5, 1.02), loc='lower center', fontsize = plot_parameters["size_fonts"],  borderaxespad=0., ncol =3)
    #Add elements on subplots
    if filename == 'all':
        string = ions[0]
        for i in range(1,len(ions)):
            string = string + ' ' + ions[i]
        ax1.text(0.95,0.9, string, transform=ax1.transAxes, horizontalalignment = 'right', fontweight = 'bold', fontsize = plot_parameters['size_fonts'])
    else:
        ax1.text(0.95,0.9, elements[0], transform=ax1.transAxes, horizontalalignment = 'right', fontweight = 'bold', fontsize = plot_parameters['size_fonts'])
    ax2.text(0.95,0.9, elements[1], transform=ax2.transAxes, horizontalalignment = 'right', fontweight = 'bold', fontsize = plot_parameters['size_fonts'])
    ax3.text(0.95,0.9, elements[2], transform=ax3.transAxes, horizontalalignment = 'right', fontweight = 'bold', fontsize = plot_parameters['size_fonts'])
    ax4.text(0.95,0.9, elements[3], transform=ax4.transAxes, horizontalalignment = 'right', fontweight = 'bold', fontsize = plot_parameters['size_fonts'])
    figurename = figurename + '.pdf'
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 150)
    print(figurename,' created')

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])
