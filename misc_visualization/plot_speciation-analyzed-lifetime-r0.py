#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: anais
Langage : Python3


        Plot the median or max lifetime of selected species as a function of T or rho for r0
        color scheme = coordination number from 1 to 11

"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D #useful to create a custom legend
from matplotlib.colors import LinearSegmentedColormap
import crystallography as cr
import natsort
import re
import glob





def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.outcar.umd.dat.r0.popul.dat_L*
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.split('/')[-1]
    temperature = str(int(float(filename.split('_')[1].strip('T'))*1000))
    acell = filename.split('.outcar')[0].split('_')[3].strip('a')[:]
    #print(acell)
    return temperature, acell


#____                             ____________                             ____ 
#    '.                         .'            '.                         .' 
#      '.                     .'     FITS       '.                     .'
#        '.                 .'     & PLOTS        '.                 .'
#          '.             .'                        '.             .'
#            '._________.'                            '._________.'
def color_change(color_scheme, color_list, color_dict):
    """apply colors to values of a dictionnary by folowing the color_list and color_scheme """
    for key in natsort.natsorted(color_list):
        c = next(color_scheme)
        color_dict[key] = c
    return color_dict


def creation_plot(max_den,xvariable,yvariable,plot_parameters,letter):
    """     ********** Creation of the plot  **********    """
    print("I use creation plot")
    #parameters for article format
    size_fonts =  plot_parameters["size_fonts"]
    size_font_ticks = plot_parameters["size_font_ticks"]
    size_figure = plot_parameters["size_figure"]
    size_markers = plot_parameters["size_markers"]
    size_lines = plot_parameters["size_lines"]
    shift_labelpad = plot_parameters["shift_labelpad"]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7,5))
    
    for ax in [ax1,ax2]:
        #Adjustment of ticks
        if xvariable == 'rho':
            if max_den <= 1.5:
                major_xticks = np.arange(0, max_den+0.1, 0.1) 
                minor_xticks = np.arange(0, max_den+0.05, 0.05)    
            else:
                major_xticks = np.arange(0, max_den+0.1, 0.5) 
                minor_xticks = np.arange(0, max_den+0.05, 0.1)    
            ax.set_xticks(major_xticks)
            ax.set_xticks(minor_xticks, minor=True)
            ax.xaxis.set_ticks_position('both')
            ax.tick_params(axis = 'x', which = 'both', labelbottom = False, width = size_lines/2)
            ax.tick_params(axis = 'y', which = 'both',labelsize = size_font_ticks, width = size_lines/2)
            ax.set_xlim(1,max_den)
        else:
            major_xticks = np.arange(0, 10000, 1000) 
            ax.set_xticks(major_xticks)
            ax.xaxis.set_ticks_position('both')
            ax.tick_params(axis = 'x', which = 'both',labelbottom = False, width = size_lines/2)
            ax.tick_params(axis = 'y', which = 'both',labelsize = size_font_ticks, width = size_lines/2)
            plt.autoscale(enable=True,axis='x',tight=False) 
        plt.autoscale(enable=True,axis='y',tight=False)
   
        ax.grid(True, which='major',axis = 'y' , linestyle=':', linewidth=size_lines/2 )
    
    ax2.tick_params(axis = 'x', which = 'both',labelbottom = True, labelsize = size_font_ticks, width = size_lines/2)
    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    fig.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax0 = fig.add_subplot(111, frameon=False)
    ax0.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    if xvariable == 'rho':
        ax0.set_xlabel(r"Density (g.cm$^{-3}$)", fontweight = 'bold', fontsize = size_fonts , labelpad = shift_labelpad)    
    if xvariable == 'T':
        ax0.set_xlabel(r"Temperature (K)", fontweight = 'bold', fontsize = size_fonts , labelpad = shift_labelpad)
    
    if yvariable == 'median':
        ax0.set_ylabel('Median cluster lifetime (fs)', fontweight = 'bold', fontsize = size_fonts , labelpad = shift_labelpad*2)
    else:
        ax0.set_ylabel('Maximum cluster lifetime (fs)', fontweight = 'bold', fontsize = size_fonts , labelpad = shift_labelpad*2)

    
    if letter != '':
        ax1.text(-0.12,0.98, letter , transform=ax1.transAxes, horizontalalignment = 'left', verticalalignment = 'top', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='none', edgecolor='k', pad=3.0))  
    return fig,ax1, ax2, ax0

def plot_and_write(xdata, xvariable, datadic,  cluster, colors_species, plot_parameters, nf, ax):
    #sort data using the first column = x
    x, y = zip(*sorted(zip( xdata[xvariable], datadic[cluster])))
    #obtain coordination # of polyhedra
    CN = int(cluster.split('_')[-1])
    #plot
    line, = ax.plot(x,y, '.-', color = colors_species[CN], markersize = plot_parameters["size_markers"], linewidth = plot_parameters["size_lines"]) #plot median
    #write
    nf.write(format_1label(cluster) + '\t' + '\t'.join(str(round(val,0)) for val in datadic[cluster]) + '\n')





#____                             ____________                             ____ 
#    '.                         .'            '.                         .' 
#      '.                     .'   FORMAT       '.                     .'
#        '.                 .'                    '.                 .'
#          '.             .'                        '.             .'
#            '._________.'                            '._________.'

def format_label(list_species):
    """ formatage of all labels of a list """
    for j in range(len(list_species)):
        list_species[j] = format_1label(list_species[j])
    return list_species

def format_1label(label):
    """formatage of label with removing _1 """
    if label != '':
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


def select_labels(datadic, xdata, colored_clusters, xvariable,ax):
    #we select the max for each cluster
    max_datadic = []
    max_datadic_melt = 0
    for cluster in datadic:
        if cluster != 'melt-like':
            max_datadic.append( (cluster, max(datadic[cluster]), xdata[xvariable][datadic[cluster].index(max(datadic[cluster]))]  )  )
        else:
            max_datadic_melt = (cluster, max(datadic[cluster]), xdata[xvariable][datadic[cluster].index(max(datadic[cluster]))]  ) 
    #we sort by the max
    max_datadic = sorted(max_datadic, key=lambda tup: tup[1], reverse=True)
    #We add the label for the 80% first cluster
    for i in range(int(len(colored_clusters)*0.8)):
        x_max = max_datadic[i][2]
        y_max = max_datadic[i][1]
        label = max_datadic[i][0]
        #formatage of label
        label = format_1label(label)
        ax.text(x_max, y_max, label, color='0.35')
    return max_datadic_melt




#____                             ____________                             ____ 
#    '.                         .'            '.                         .' 
#      '.                     .'     MAIN       '.                     .'
#        '.                 .'     PROGRAM        '.                 .'
#          '.             .'                        '.             .'
#            '._________.'                            '._________.'
def main(argv):
    """     ********* Main program *********     """
    atoms = 'all'
    max_den = 6
    size = {}
    lifetime = {}
    datadic = {}
    xdata = {'T':[],'rho':[]} #dictionnary containing the x data
    selected_clusters = []
    selected_files = [] #list containing the files we use for the plot (depends on the density limit)
    #other dictionnaries and parameters for the figure for article version
    list_colored_species = []
    list_colored_red = []
    list_colored_green = []
    list_colored_blue = []
    list_colored_purple = []
    list_minor_species = []
    colors_species = {"melt-like":'0'}
    letter = '' 
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (6,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 20}
    #other parameters
    Na=6.022*10**23
    sizelimit = 13
    try:
        options,arg = getopt.getopt(argv,"hf:v:m:l:a:t:d:",["filenameoutput","variable","mineralfile","letter","atoms",'type',"density_max"])
    except getopt.GetoptError:
        print("plot_speciation-analyzed-lifetime-r0.py  -f <output filename> -v <variable (rho,T)> -m <mineralfile with elements> -l <letter for article subplot, default = ''> -t <type (max or median)> -d <maximum density to plot in g/cm3> -a <list of first atoms determining the type of cluster to plot (e.g. 'Si_1O,Al_1O,melt,O')>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_speciation-analyzed-lifetime-r0.py program to  Plot the average and median lifetime of selected species as a function of T or rho')
            print("plot_speciation-analyzed-lifetime-r0.py -f <output filename> -v <variable (rho,T)> -m <mineralfile with elements> -l <letter for article subplot, default = ''> -t <type (max or median)>  -d <maximum density to plot in g/cm3> -a <list of first atoms determining the type of cluster to plot (e.g. 'Si_1O,Al_1O,melt,O')>")
            print('')
            print("all the species with more than 13 atoms are grouped together under the name 'melt-like'. To display the corresponding lifetime, indicate melt in the list of atoms")
            print('')
            print('requires to be launched from subfolder T or acell containing the requires the files popul.dat ')
            print('requires also the file containing elements and number (in order to compute the densities)')
            print('')
            sys.exit()
        elif opt in ('-f','--filename'):
            outputfilename = str(arg)
        elif opt in ("-v", "--variable"):
            xvariable = str(arg)
        elif opt in ("-m","--mineralfile"):
            mineralfile = str(arg)
        elif opt in ("-l","--letter"):
            letter = str(arg)
        elif opt in ("-a","--atoms"):
            atoms = arg.split(',')
        elif opt in ("-t","--type"):
            yvariable = str(arg)
        elif opt in ("-d","--density_max"):
            max_den = float(arg)
    #*****************************
    #*******************
    #*******
    #***** Calculation of the molecular mass
    with open(mineralfile,'r') as mf:
        entry = mf.readline()
        elements = entry.split()[1:]
        entry = mf.readline()
        number = entry.split()[1:]
    MN = 0
    for i in range(len(elements)):
        MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
    #*****************************
    #*******************
    #*******
    #***** Extraction of the lifetimes and creation of the xdata lists
    files = sorted(glob.glob('*r0.popul.dat')) #I list every popul.dat files in alphabetic order
    for file in files:
        #creation of the x data list
        temperature, acell = split_name(file)
        density = MN/(Na*float(acell)**3*10**(-24))
        if density <= max_den:    
            xdata['rho'].append(density )    #calculation density
            xdata['T'].append(int(temperature))   
            selected_files.append(file)
            lifetime[file] = {} #dictionnary with the lifetime of all the cluster in each file
            with open(file,'r') as f:     
                f.readline()
                f.readline()
                while True:
                    line = f.readline()
                    if not line: break
                    else:
                        entry=line.split('\t')
                        natoms = len(entry[4].split(','))
                        if natoms > sizelimit:
                            #print('for ',entry[0], 'natoms = ',natoms, '>', sizelimit)
                            try:
                                lifetime[file]['melt-like'].append(float(entry[3]))
                            except KeyError :
                                lifetime[file]['melt-like'] = [float(entry[3])]
                                size['melt-like'] = '>'+str(sizelimit)
                        else:
                            try:
                                lifetime[file][entry[0]].append(float(entry[3]))
                            except KeyError :
                                lifetime[file][entry[0]] = [float(entry[3])]
                                size[entry[0]] = natoms
                #print('lifetimes melt-like',lifetime[file]['melt-like'])
    print('for files', selected_files, len(selected_files))
    #print('we have xdata',xdata['rho'], len(xdata['rho']))
    #print('we have xdata',xdata['T'], len(xdata['T']))
    #**** Creation of the dictionnaries and lists of max or median lifetime for each cluster
    #initialisation of the dictionnary containing a list for each cluster selected
    if atoms == ['all']:
        print("we use all clusters")
        for file in selected_files:
            for cluster in lifetime[file]:
                datadic[cluster] = []
    else:
        print("we use the selected clusters")
        for file in selected_files:
            for cluster in lifetime[file]:
                for atom in atoms:
                    if cluster[:len(atom)] == atom:
                        datadic[cluster] = []                    
    #choice of the function to use
    if yvariable == 'median':
        function = np.median
    else:
        function = np.max
    for cluster in datadic:
        selected_clusters.append(cluster)
        for file in selected_files:
            try:
                datadic[cluster].append(function(lifetime[file][cluster]))
            except KeyError:
                datadic[cluster].append(0)
    #print('medianlife is',medianlife)
    #print('averagelife is',averagelife)
    print("The selected clusters are",selected_clusters)
    #*****************************
    #*******************
    #*******
    #***** Color attribution
    #**** Now the dictionary has every possible cluster, we attribute the colors to each cluster
    #check the lifetime not display the species with an average lifetime smaller than 50 fs and we add name of the species  
    for cluster in datadic : #loop on each cluster
        if max(datadic[cluster]) < 0:
            list_minor_species.append(cluster)
        elif cluster == 'melt-like':
            continue
        else:
            #this is for basic automatic color change
            list_colored_species.append(cluster)
    #**** we attribute the colors to each cluster 
    #this is for basic automatic color change
    color = iter(plt.cm.jet(np.linspace(0,1,11))) #Creation of the color list for 11 coordination #
    for i in range(1,12,1):
        c = next(color)
        colors_species[i] = c
    #*****************************
    #*******************
    #*******
    #*** Creation of the plots
    #one subplot per type without interstitial --> 2 subplots
    fig, ax1, ax2, ax0 = creation_plot(max_den,xvariable,yvariable,plot_parameters,letter)
    #plot for each cluster and write file with median and average lifetime
    figurename = outputfilename+'_'+xvariable+'_r0_'+yvariable
    newfilename = figurename + '.txt'
    nf = open(newfilename,'w')
    if xvariable == 'rho':
        nf.write('Density(g/cm3)\t' + '\t'.join(str(round(density,2)) for density in xdata[xvariable]) + '\n')
    else:
        nf.write('Temperature(K)\t' + '\t'.join(str(temp) for temp in xdata[xvariable]) + '\n')
    #plot lines
    for cluster in list_colored_species:
        print("for cluster",cluster)
        if cluster[:5] == 'Al_1O' :
            #print('plot on ax1 for cluster', cluster)
            plot_and_write(xdata, xvariable,  datadic, cluster, colors_species, plot_parameters, nf, ax1)
        elif cluster[:5] == 'Si_1O' :
            #print('plot on ax2 for cluster', cluster)
            plot_and_write(xdata, xvariable,  datadic, cluster, colors_species, plot_parameters, nf, ax2)
    #*****************************
    #*******************
    #*******
    #legend           
    list_minor_species = format_label(natsort.natsorted(list_minor_species))       
    custom_lines = [Line2D([0],[0],color = colors_species[CN], ls = '--', marker = '.', markersize = plot_parameters["size_markers"], linewidth = plot_parameters["size_lines"]) for CN in range(1,12,1)]
    list_colored_species = [str(CN) for CN in range(1,12,1)]
    print(list_colored_species)
    legend = ax0.legend([line for line in custom_lines],[label for label in list_colored_species],title = '$\\bf{x}$', bbox_to_anchor=(0.5, 1.04), loc='lower center', fontsize = plot_parameters["size_fonts"],  borderaxespad=0., ncol = 11)
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_fonts"])
    ax1.text(0.99,0.99, r'AlO$_{x}$' , transform=ax1.transAxes, horizontalalignment = 'right', verticalalignment = 'top', fontweight = 'bold', fontsize = plot_parameters["size_fonts"])
    ax2.text(0.99,0.99, r'SiO$_{x}$' , transform=ax2.transAxes, horizontalalignment = 'right', verticalalignment = 'top', fontweight = 'bold', fontsize = plot_parameters["size_fonts"])
    #*****************************
    #*******************
    #*******
    #save plots
    extension = '.svg'
    figurename = figurename+ extension
    fig.savefig(figurename, bbox_inches='tight', dpi=300)
    print(figurename, 'is created')    
    nf.write('\n')
    nf.write('Species_with_medianlife<50fs\t'+ '\t'.join(str(species) for species in list_minor_species) + '\n')
    nf.close()
    print(newfilename, 'is created')
#    plt.show() # à mettre après savefig sinon on ne sauvegarde rien !!  

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



