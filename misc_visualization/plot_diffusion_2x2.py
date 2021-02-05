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
import natsort
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import crystallography as cr
import re
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



def format1label(label):
    """formatage of label with removing _1 """
    i=0
    while True:
        if i <= len(label)-1:
            if re.match('_',label[i]):
                num=0
                try:
                    while re.match('[0-9]',label[i+1+num]):
                        num +=1
                except IndexError: #index error when we arrive at the end of the cluster name
                    pass
                    #print('end of the cluster')
                if label[i+1:i+1+num] == '1':
                    label = label[:i] + label[i+1+num:]
                else:
                    label = label[:i]+'$_{'+label[i+1:i+1+num]+'}$' + label[i+1+num:]
                    i = i+5
            i = i+1
        else:break
    return label

def extract_TP(thermofile, column_number, TP, addtxt):
    """extract temperature and pressure from thermofile"""
    with open(thermofile, 'r') as f:
        [f.readline() for i in range(3)]
        #extract data
        while True:
            line = f.readline()
            if not line: break
            else:    
                entry=line.split('\n')[0].split('\t')
                TP[entry[0].split('outcar.umd.dat')[0].split('/')[-1]+addtxt] = (int(entry[column_number['T']]),float(entry[column_number['P']]))
    return TP



#____                             ____________                             ____ 
#    '.                         .'            '.                         .' 
#      '.                     .'     PLOT       '.                     .'
#        '.                 .'    EXPE DATA       '.                 .'
#          '.             .'                        '.             .'
#            '._________.'                            '._________.'
def plot_deKoker2010(ax1,ax2,ax3,ax4, handles, legendlabels,xvariable):
    columns_deKoker2010 = {'Ca':8, 'Al':9, 'Si':10, 'O':11}
    X, D = {'3000':[],'4000':[],'6000':[]} , {'Ca':{'3000':[],'4000':[],'6000':[]},
            'Al':{'3000':[],'4000':[],'6000':[]},'Si':{'3000':[],'4000':[],'6000':[]},
            'O':{'3000':[],'4000':[],'6000':[]}}
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
                    if xvariable == 'P':
                            X[T].append(float(entry[5]))
                    else:
                        X[T].append(float(entry[2])/1000)
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
                ax.plot(X[temp],D[elem][temp], ':', color = colors_T_gray[temp], linewidth = 2)      
        #add to legend
        handles.append(plt.Line2D([0],[0], color = 'k', ls = ':'))
        legendlabels.append('deKoker, 2010')
        return handles, legendlabels
    except FileNotFoundError:
        print("File deKoker2010 not found")            
    
def plot_Neilson2016(ax1,ax2,ax3,ax4, handles, legendlabels, xvariable):
    columns_Neilson2016 = {'Na':8, 'Al':9, 'Si':10, 'O':11}
    X, D = {'3000':[],'4000':[],'5000':[]} , {'Na':{'3000':[],'4000':[],'5000':[]},
            'Al':{'3000':[],'4000':[],'5000':[]},'Si':{'3000':[],'4000':[],'5000':[]},
            'O':{'3000':[],'4000':[],'5000':[]}}
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
                        if xvariable == 'P':
                            X[T].append(float(entry[5]))
                        else:
                            X[T].append(float(entry[2])/1000)
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
                ax.plot(X[temp],D[elem][temp], '-', color = colors_T_gray[temp], linewidth = 2)      
        #add to legend
        handles.append(plt.Line2D([0],[0], color = 'k', ls = '-'))
        legendlabels.append('Neilson $et$ $al.$, 2016')
        return handles, legendlabels
    except FileNotFoundError:
        print("File Neilson2016 not found")   

def plot_Spera2009(ax1,ax2,ax3,ax4, handles, legendlabels,xvariable):
    columns_Spera2009 = {'Ca':8, 'Al':9, 'Si':10, 'O':11}
    X, D = {'4000':[],'5000':[],'6000':[]} , {'Ca':{'4000':[],'5000':[],'6000':[]},
            'Al':{'4000':[],'5000':[],'6000':[]},'Si':{'4000':[],'5000':[],'6000':[]},
            'O':{'4000':[],'5000':[],'6000':[]}}
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
                        if xvariable == 'P':
                            X[T].append(float(entry[5]))
                        else:
                            X[T].append(float(entry[2])/1000)
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
                ax.plot(X[temp],D[elem][temp], ':', color = colors_T_gray[temp], linewidth = 2)
        #add to legend
        handles.append(plt.Line2D([0],[0], color = 'k', ls = ':'))
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


def creation_plot_2x2(plot_parameters, prop, xvariable):
    """     ********** Creation of the plot  **********    """
    plt.close()
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, sharex=True, sharey=True, figsize = plot_parameters['size_figure'])
    
    if xvariable == 'rho':
        major_xticks = np.arange(0, 4.5, 0.5) 
        minor_xticks = np.arange(0, 4.1, 0.1) 
        for ax in [ax1,ax2,ax3,ax4]:
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True)
        ax1.set_xlim(1.0,4.1)
        xlabel = r'Density (g.cm$^{-3}$)'
    else:
        ax1.set_xlim(1,275)
        ax1.set_xscale('log')
        xlabel = r'Pressure (GPa)'
        
    for ax in [ax1,ax2,ax3,ax4]:
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(which = 'both', labelsize = plot_parameters['size_font_ticks'], width = plot_parameters['size_lines']/2)                                            
    #        ax.set_yticks(major_yticks)
    #        ax.set_yticks(minor_yticks, minor=True)                                           
        ax.yaxis.set_ticks_position('both')
    for ax in [ax3]:
        plt.setp(ax.get_xticklabels()[-1], visible=False) 
        #ax.grid(True, axis = 'x')

    if prop == 'D':
        ax1.set_ylim(4e-10,3e-7)
        ylabel = r'Diffusion coefficient (m$^2$.s$^{-1}$)'
        ax1.set_yscale('log')
    else:
        ax1.set_ylim(0,500)
        ylabel = r'Transition time (fs)'
        

    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    f.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax = f.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, 
                    labelleft=False, left=False, labelright = False, right=False)
    ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = plot_parameters['size_fonts'],
                  labelpad = plot_parameters['shift_labelpad']*2)
    ax.set_ylabel(ylabel, fontweight = 'bold', fontsize = plot_parameters['size_fonts'], 
                  labelpad = plot_parameters['shift_labelpad']*3+plot_parameters['shift_labelpad']/2)
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
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,8),
                       "size_markers" : 5,"size_lines" : 1,"shift_labelpad" : 10}
    #other dictionnaries and parameters for the figure
    style_markers = {'CaAl2Si2O8':'d','KAlSi3O8':'s','NaAlSi3O8':'o'}
    style_lines = {'CaAl2Si2O8':':','KAlSi3O8':'--','NaAlSi3O8':'-'}
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200',
                'T5':'#ffcd01','T6':'#ff0101','T6.5':'#ff00a2','T7':'#ff01de'}
    #colors_T = create_colors()
    #variables nedded fot the plot
    filename = 'all'
    filename2 = ''
    ions = []
    compounds = []
    Temperatures = []
    data = {} #dictionnary containing the data for each element, initialized for each T
    stdev = {} #same for stdev
    #other parameters
    Na=6.022*10**23
    #dictionnary for fullaverages file in full version
    column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,
                     'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,
                     'Cvm':25,'stdev_Cvm':26,'testCv':31,'stdev_testCv':32}
    try:
        options,arg = getopt.getopt(argv,"hf:g:p:v:",["filename","gfilename2",'property','variable'])
    except getopt.GetoptError:
        print("plot_diffusion_2x2.py -f <filename>(default = all) -g <filename2(option)> -p <property to plot (D or t)> -v <variable x axis (rho or P)> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_diffusion_2x2.py program to plot diffusion coefficient (or time of ballistic to diffusive regime change) as a function of density or pressure for each T and the selected minerals containing 4 elements')
            print("plot_diffusion_2x2.py -f <filename>(default = all files of every compound) -g <filename2(option)>  -p <property to plot (D or t)> -v <variable x axis (rho or P)>")
            print("plot_diffusion_2x2.py requires to be lauched from the folder containing every diffusivities file created by the script analyze_msd")
            print('')
            print('For plots as function of P, make sure the files from fullaverages.py are in the current folder')
            sys.exit()
        if opt in ('-f','--filename'):
            filename = str(arg)
        elif opt in ('-g','--gfilename2'):
            filename2 = str(arg)
        elif opt in ('-p','--property'):
            prop = str(arg)
        elif opt in ('-v','--variable'):
            xvariable = str(arg)
    #** Figure creation 
    fig,ax1,ax2,ax3,ax4, ax0 = creation_plot_2x2(plot_parameters, prop, xvariable)
    if prop == 'D':
        name = 'diffusivities'
        index = 1
    else:
        name = 'regimechange'
        index = 6
    if filename == 'all':
        files = sorted(glob.glob('diffusivities_*.txt'),reverse=True) #I list every diffusivities files
        figurename = name + '_2x2_all_'+xvariable
    elif filename2 != '':
        files = [filename,filename2] #I list every diffusivities files
        figurename = name + '_2x2_'+filename.split('.txt')[0].split('_')[1]+'_'+filename2.split('.txt')[0].split('_')[1]+'_'+xvariable    
    else:
        files = [filename] #I take only the file we want
        figurename = name + '_2x2_'+filename.split('.txt')[0].split('_')[1]+'_'+xvariable
    #** Plot expe data  
    handles, legendlabels =  [] ,[]
    if prop == 'D':
        handles, legendlabels = plot_deKoker2010(ax1,ax2,ax3,ax4, handles, legendlabels, xvariable)
        handles, legendlabels = plot_Neilson2016(ax1,ax2,ax3,ax4, handles, legendlabels, xvariable)
        #handles, legendlabels = plot_Spera2009(ax1,ax2,ax3,ax4, handles, legendlabels, xvariable)
    #** Plot our data
    for file in files:
        #print("************************ for diffusivity file named",file)
        #**extraction compound        
        compounds.append(file.split('_')[1].split('.txt')[0]) #we need the compound for the legend
        #******* Extract P and T for each thermo file
        if xvariable == 'P':
            TP = {}
            mineralname = file.split('_')[1].split('.txt')[0]
            thermofiles = sorted(glob.glob('thermo_'+mineralname+'*.txt'))
            for thermofile in thermofiles:
                TP = extract_TP(thermofile, column_number, TP,'')
            if TP == {}:
                print('ERROR!!!! TP dictionnary empty')
        #**creation of elements and number lists and initialization of T
        with open(file,'r') as f:
            skiplines = 0
            while True:
                line = f.readline()
                skiplines += 1
                entry=line.split()
                if entry[0] == 'elements':
                    elements = entry[1:]
                    ions.append(elements[0])
                elif entry[0] == 'number':
                    number = entry[1:]
                elif entry[0] == 'file':
                    line = f.readline()
                    entry=line.split()
                    mineral, temperature0, acell0 = split_name(entry[0])
                    break
        #**calculation of M*N nedded for the calculation of densities
        MN = 0
        for i in range(len(elements)):
            MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
        #**initialisation of data
        data = {} #big dictionnary with inside self diffusion coefficient or time of regime change, all coresponding to different T and atom pairs
        stdev = {} #idem for stdev
        X = {}  #idem for rho or P
        for elem in elements:
            stdev[elem] = {}
            data[elem] = {}
            X[elem] = {}
        #****** Extraction of data
        with open(file,'r') as f:
            [f.readline() for i in range(skiplines)]
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    mineral, temperature, acell = split_name(entry[0]) 
                    #print(entry[0])
                    for i in range(len(elements)): 
                        elem = elements[i]
                        try:
                            data[elem][temperature].append(float(entry[i*6+index]))
                            stdev[elem][temperature].append(float(entry[i*6+2]))    
                            if xvariable == 'P':
                                X[elem][temperature].append(TP[entry[0].split('outcar.msd.dat')[0].split('/')[-1]][1])
                            else:
                                X[elem][temperature].append( MN/(Na*float(acell)**3*10**(-24)) )    #calculation density
                        except KeyError:
                            data[elem][temperature] = [ float(entry[i*6+index])]
                            stdev[elem][temperature] = [ float(entry[i*6+2]) ]
                            if xvariable == 'P':
                                X[elem][temperature] = [TP[entry[0].split('outcar.msd.dat')[0].split('/')[-1]][1]]
                            else:
                                X[elem][temperature] = [ MN/(Na*float(acell)**3*10**(-24)) ]    #calculation density
        #***** Plot each element on the correct subplot 
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
            for temperature in data[elem]:
                #print("******** ", temperature)
                #print(data[elem][temperature])
                Temperatures.append(temperature)
                #attribution of fill color
                if mineral == 'NaAlSi3O8':
                    fillcolor = colors_T[temperature]+'ff'#with custom dict    #np.array([colors_T[temperature][0],colors_T[temperature][1],colors_T[temperature][2],1])   #with dic from create_colors() function
                elif mineral == 'CaAl2Si2O8':
                    fillcolor = 'w'
                elif mineral == 'KAlSi3O8':
                    fillcolor =  colors_T[temperature]+'7f'#with custom dict    #np.array([colors_T[temperature][0],colors_T[temperature][1],colors_T[temperature][2],0.5])  #with dic from create_colors() function
                #remove data with nan 
                ThisData = np.array(data[elem][temperature])
                ThisX= np.array(X[elem][temperature])
                data_mask = np.isfinite(ThisData)                
                ax.plot(ThisX[data_mask],ThisData[data_mask], style_markers[mineral]+style_lines[mineral],
                        markersize = plot_parameters["size_markers"], markeredgewidth = 0.5, 
                        color = colors_T[temperature], markerfacecolor = fillcolor, 
                        markeredgecolor = colors_T[temperature], linewidth = plot_parameters["size_lines"])
                #if prop == 'D':
                    #ax.errorbar(X[elem][temperature],data[elem][temperature], yerr=stdev[elem][temperature], fmt=style_markers[mineral], markersize = plot_parameters[size_markers], markeredgewidth = 0.5, edgecolor colors_T[temperature], facecolor = fillcolor, linestyle = style_lines[mineral],  linewidth = plot_parameters["size_lines"] )
    #********* Legend               
    #Create legend from custom artist/label lists
    Temperatures = list(set(Temperatures)) #get elements only once in the list
    custom_patch = [mpatches.Patch(color=colors_T[key]) for key in natsort.natsorted(Temperatures)]
    legend = ax0.legend([col for col in custom_patch],[str(int(float(label.strip('T'))*1000)) for label in natsort.natsorted(Temperatures)],title = '$\\bf{Temperature~(K)}$', bbox_to_anchor=(0.5, 1.07), loc="lower center", fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol=len(Temperatures))
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])
    #definition of line and marker style
    fillcolor = {}
    for mineral in style_markers:
        if mineral == 'CaAl2Si2O8':
            fillcolor[mineral] = 'w'
        elif mineral == 'KAlSi3O8':
            fillcolor[mineral] = '#0000007f'
        else:
            fillcolor[mineral] = '#000000ff' 
    custom_lines = [plt.Line2D([0],[0], marker = style_markers[key], linestyle = style_lines[key],
                               markeredgecolor = 'k', markerfacecolor = fillcolor[key],  
                               color = 'k', markersize = plot_parameters["size_markers"], 
                               markeredgewidth = 0.5, linewidth = plot_parameters["size_lines"]) for key in sorted(compounds)]
    ax0.legend([line for line in custom_lines],[format1label(label) for label in sorted(compounds)],
               bbox_to_anchor=(0.5, 1.02), loc='lower center', fontsize = plot_parameters["size_fonts"], 
               borderaxespad=0., ncol =3)
    
    ax0.add_artist(legend)
    #Add elements on subplots
    if filename == 'all' or filename2 != '':
        string = ions[0]
        for i in range(1,len(ions)):
            string = string + ' ' + ions[i]
        ax1.text(0.95,0.9, string, transform=ax1.transAxes, horizontalalignment = 'right',
                 fontweight = 'bold', fontsize = plot_parameters['size_fonts'])
    else:
        ax1.text(0.95,0.9, elements[0], transform=ax1.transAxes, horizontalalignment = 'right', 
                 fontweight = 'bold', fontsize = plot_parameters['size_fonts'])
    ax2.text(0.95,0.9, elements[1], transform=ax2.transAxes, horizontalalignment = 'right', 
             fontweight = 'bold', fontsize = plot_parameters['size_fonts'])
    ax3.text(0.95,0.9, elements[2], transform=ax3.transAxes, horizontalalignment = 'right', 
             fontweight = 'bold', fontsize = plot_parameters['size_fonts'])
    ax4.text(0.95,0.9, elements[3], transform=ax4.transAxes, horizontalalignment = 'right', 
             fontweight = 'bold', fontsize = plot_parameters['size_fonts'])
    figurename = figurename + '.pdf'
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 150)
    print(figurename,' created')

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])
