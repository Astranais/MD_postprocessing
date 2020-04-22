#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on  Oct   2019

@author: anais
Langage : Python3


        Horizontal Bar plot of clusters life as a function of simulation time
        indicating where are the selected atom index (in the melt or in which gas species)

"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re
import natsort



def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.split('/')[-1]
    temperature = str(int(float(filename.split('_')[1].strip('T'))*1000))
    acell = filename.split('_')[3].strip('a')
    return temperature, acell


def creation_plot():
    """     ********** Creation of the plot  **********    """
    plt.close()
    fig, ax = plt.subplots(figsize=(12,12))
    ax.set_ylabel(r'Species formula', fontweight = 'bold', fontsize = 12)
    ax.set_xlabel(r"Step in the simulation", fontweight = 'bold', fontsize = 12, labelpad = 10)
    #Adjustment of ticks
    ax.tick_params(which = 'both', labelsize = 10, width = 0.5)
    plt.tick_params(left = False, right = False, labelbottom = True)
    plt.autoscale(enable=True,axis='both',tight=True)
    plt.margins(y=0)
    return fig,ax


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

def main(argv):
    """     ********* Main program *********     """
    data = {} #dictionnary containing list of begin and end lifetime for each cluster of each species
    data_atoms = {} #dictionnary used to plot data
    compospecies = {} #dictionnary with all atom indices in the simu
    sizes = {} #dictionnary containing size of species
    colors_indices = {} #dict for colors
    label_clusters = [] #list for label
    begintime = 0
    endtime = 0
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    try:
        options,arg = getopt.getopt(argv,"hf:b:e:i:",["file","begintime","endtime",'indices'])
    except getopt.GetoptError:
        print("plot_speciation-lifetime-species_loc-index.py -f <.popul.dat filename>  -b <(option) begin of time window> -e <(option) end of time window>  -i <atoms indices, either list of selected indices or range min-max>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_speciation-lifetime-species_loc-index.py horizontal Bar plot of clusters life as a function of simulation time for one simulation')
            print('plot_speciation-lifetime-species_loc-index.py show whether the atom indices are in the melt or in which gas phase')
            print("plot_speciation-lifetime-species_loc-index.py -f <.popul.dat filename>  -b <(option) begin of time window> -e <(option) end of time window>  -i <atoms indices, either list of selected indices or range min-max>")
            print('')
            sys.exit()
        elif opt in ("-f", "--file"):
            populfile = str(arg)
        elif opt in ("-b","--begintime"):
            begintime = int(arg)-1
        elif opt in ("-e","--endtime"):
            endtime = int(arg) +1   
        elif opt in ("-i","--indices"):
            entry = str(arg)
            if '-' in entry:
                minindex, maxindex = entry.split('-')
                indices = [str(i) for i in range(int(minindex), int(maxindex)+1)]
            else:
                indices = natsort.natsorted(entry.split(','))
    print("working with indices", indices)
    #***** initialization of data dictionnary
    with open(populfile,'r') as f:
        [f.readline() for i in range(2)]
        #we read and store all the data into lists
        while True: 
            line = f.readline()
            if not line: break
            else:
                entry=line.split('\t')
                data[entry[0]] = []
                compospecies[entry[0]] = []
    data['melt'] = []
    compospecies['melt'] = []
    #***** Extraction of maxsimu
    maxsimu = max(np.loadtxt(populfile, usecols = 2, skiprows = 2))
    if endtime == 0:        
        endtime = maxsimu+1
    print("max of simu is",endtime-1)
    #***** Extraction of data in the selected time window
    with open(populfile,'r') as f:
        [f.readline() for i in range(2)]
        #we read and store all the data into lists
        while True: 
            line = f.readline()
            if not line: break
            else:
                entry=line.split('\n')[0].split('\t')
                if int(entry[1]) > begintime and int(entry[2]) < endtime:
                    composition = re.findall('\d+',entry[4])#list of atom indices  for this cluster
                    sizes[entry[0]] = len(composition)
                    if sizes[entry[0]] > 13:
                        data['melt'].append( (float(entry[1]),float(entry[2]), composition)  ) #begin lifetime, end lifetime, composition
                        compospecies['melt'].extend(composition)
                    else:
                        data[entry[0]].append( (float(entry[1]),float(entry[2]), composition)  ) #begin lifetime, end lifetime, composition
                        compospecies[entry[0]].extend(composition)
    sizes['melt'] = 208
    #***** Selection of species to plot and initialization of the dictionnary used to plot data
    selected_species0 = []
    for species in data:
        for atomindex in indices:
            if atomindex in compospecies[species]:
                selected_species0.append(species)
                break
    selected_species0 = natsort.natsorted(selected_species0)
    #ordering of species based now on sizes
    selected_species = []
    for i in range(0,13):
        for species in selected_species0:
            if sizes[species] == i:
                selected_species.append(species)
    selected_species.append('melt')
    print("selected species are", selected_species)
    #***** Creation of color scheme for indices, y axis and plot
    color = iter(plt.cm.jet(np.linspace(0,1,len(indices)))) #Creation of the color list
    for key in natsort.natsorted(indices):
        c = next(color)
        colors_indices[key] = c
    #***** Extraction of data - Creation of the dictionnary with the begin and end lifetime of each cluster for each atom index wanted
    #initialization
    for species in data:
        data_atoms[species] = {}
        for ii in range(208):
            data_atoms[species][str(ii)] = []
    #fill with all the data for every atom index of the selected species
    for species in selected_species:
        for ii in range(len(data[species])):
            for atomindex in data[species][ii][2]:
                data_atoms[species][atomindex].append( (data[species][ii][0],data[species][ii][1]) )
    #***** Creation of the plot
    # Compute the number of lines (one line per atomindex per species)
    Nvertical = [0 for i in range(len(selected_species))]
    for ii in range(len(selected_species)):
        for index in indices:
            if data_atoms[selected_species[ii]][index] != []:
                Nvertical[ii] +=1
    print("each having # different atoms",Nvertical)
    fig,ax = creation_plot()
    y = -0.5
    for ii in range(len(selected_species)):
        label_clusters.append('')
        #label_clusters.append(format1label(selected_species[ii]))
        #count the number of lines in order to place the label of the species in the middle
        count = 0
        for index in indices:
            if data_atoms[selected_species[ii]][index] != []:
                count+=1
        #plot the bars and create the label list
        count2 = 0
        for index in indices:
            if data_atoms[selected_species[ii]][index] != []:
                y+=1
                count2 +=1
                for jj in range(len( data_atoms[selected_species[ii]][index]  )):
                    ax.barh(y,data_atoms[selected_species[ii]][index][jj][1]-data_atoms[selected_species[ii]][index][jj][0]+1, left = data_atoms[selected_species[ii]][index][jj][0] ,height = 1, color =  colors_indices[index]   )
                if count2 == count//2+1:
                    label_clusters.append(format1label(selected_species[ii]))
                else:
                    label_clusters.append('')
        y+=1
        try:
            if sizes[selected_species[ii+1]] == sizes[selected_species[ii]]:
                ax.axhline(y,xmin=0,xmax = 1, color = '0', ls='-', linewidth = 0.75) #space between different species only
            else:
                ax.axhline(y,xmin=0,xmax = 1, color = '0', ls='-', linewidth = 0.75) #space between different sizes 
        except IndexError:#when we arrive at the end of selected_species list
            pass
    
    #**** Legend
    plt.yticks(np.arange(0,y+1,1), label_clusters, rotation = 25, fontsize = 12)   
    plt.ylim(0,y-0.5)
    custom_patch = [mpatches.Patch(color=colors_indices[key]) for key in natsort.natsorted(colors_indices, reverse = True)]
    legend = plt.legend([col for col in custom_patch], [label for label in natsort.natsorted(colors_indices, reverse = True)],title = '$\\bf{atomic}$\n$\\bf{index}$', bbox_to_anchor=(1.01, 0), loc='lower left', fontsize = plot_parameters["size_font_ticks"],  borderaxespad=0., ncol = 1)
    plt.setp(legend.get_title(),fontsize= plot_parameters["size_font_ticks"])



    figurename = populfile.split('/')[-1]+'_species_loc_index.pdf'
    plt.savefig(figurename, bbox_inches='tight', dpi=300)
    print(figurename, 'is created')
    #plt.show()


#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



