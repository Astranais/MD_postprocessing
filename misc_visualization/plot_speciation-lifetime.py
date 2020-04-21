#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on  Jun  13 2018

@author: anais
Langage : Python3


        Bar plot of clusters lifetime in the gas phase 

"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator,AutoMinorLocator
import crystallography as cr
import natsort



def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.strip('./')
    temperature = str(int(float(filename.split('_')[1].strip('T'))*1000))
    acell = filename.split('_')[3].strip('a')
    return temperature, acell


def creation_plot():
    """     ********** Creation of the plot  **********    """
    plt.close()
    fig, ax = plt.subplots(figsize=(12,7))
    ax.set_ylabel(r'Cluster absolute lifetime (steps)', fontweight = 'bold', fontsize = 10)
    ax.set_xlabel(r"Cluster formula", fontweight = 'bold', fontsize = 10, labelpad = 10)
    #Adjustment of ticks
    ymajorLocator = AutoLocator()
    yminorLocator = AutoMinorLocator()
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    plt.yticks(fontsize = 8)
    plt.tick_params(bottom = False, top = False, labelbottom = True)
    plt.autoscale(enable=True,axis='both',tight=False)
    ax.grid(True, which='major',axis = 'y', linestyle=':', linewidth=1/2 )
    return fig,ax


def main(argv):
    """     ********* Main program *********     """
    clusters = [] #list containing all clusters
    lifetime = [] #list containing all clusters lifetime
    size = [] #list containing all clusters size
    shortpop = [] #list contianing all data below length
    label_clusters = [] #list for label
    total_lifetime = [] #list for all final lifetime sorted
    total_size = []
    all_data = {} #dictionnary containing all lifetime for each type of cluster
    try:
        options,arg = getopt.getopt(argv,"hf:l:",["file","length"])
    except getopt.GetoptError:
        print("plot_speciation-lifetime.py -f <.popul.dat filename>  -l <maximum length of gas cluster>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_speciation-lifetime.py program to  Plot cluster lifetime for one simulation and create txt files with only lifetimes of small clusters inside')
            print("plot_speciation-lifetime.py -f <.popul.dat filename> -l <maximum length of gas cluster>")
            print('')
            sys.exit()
        elif opt in ("-f", "--file"):
            populfile = str(arg)
        elif opt in ("-l","--length"):
            length = int(arg)
    #print(timesteps)
    #***** Creation of the lists with cluster name, lifetime and cluster length
    with open(populfile,'r') as f:
        [f.readline() for i in range(2)]
        #we read and store all the data into lists
        while True: 
            line = f.readline()
            if not line: break
            else:
                entry=line.split('\t')
                clusters.append(entry[0])
                lifetime.append(entry[3])
                composition = entry[4].split(',')
                size.append(len(composition))
    #we sort the population with the cluster size
    population = zip(size, clusters, lifetime) 
    population = natsort.natsorted(population)
    #we remove all the clusters larger than the length indicated
    for i in range(len(population)): 
        if population[i][0] <= length: 
            shortpop.append(population[i])
    try:
        data = [shortpop[0]] #we initialize the data for the plot
    except IndexError:
        print('No clusters smaller than ',length,' atoms in file',populfile)
        sys.exit()
    #we create  final big lists (label_clusters, total_lifetime, total_size) used to plot data that are sorted by size AND lifetime, with correct labels 
    maxlen = 1 #maximum length of lifetime list among all clusters, needed to fill the gaps in order to print out the .txt file
    counter = 1  #initialization of the counter (count the number of time we see this cluster)
    print("(cluster size, cluster name, lifetime)")
    for i in range(1,len(shortpop)):
        print(shortpop[i-1])
        #if we do not change of cluster then we add the current data set to the data list which contain all lifetime for the current cluster type 
        if shortpop[i][1] == shortpop[i-1][1]:
            data.append(shortpop[i]) 
            counter += 1 
        #if we do change of cluster type, then we sort data by lifetime and ad them to the final list
        else:
            #print(shortpop[i][1], '!=', shortpop[i-1][1])
            #print('the species', shortpop[i-1][1], 'appeared', counter, 'times')
            #we should not write "if counter > 0:" to save only clusters we see more than one time because at low T we don't see the big polymer anymore
            size, clusters, lifetime = zip(*data) #we extract the lifetime of the selected data (i.e. only for the cluster type we worked until now)
            lifetime = [float(x) for x in lifetime] #convert str to int
            lifetime = [lifetime[i] for i in range(len(lifetime))] #multiply by the timestep to obtain lifetime in fs
            lifetime = sorted(lifetime, reverse=True) # we sort data by lifetime
            clusters = list(clusters) #convert to list
            size = list(size) #convert to list
            #print(size)
            #all_data[str(clusters[0])]=lifetime
            all_data[str(clusters[0])]= {'lifetime':lifetime,'size':size[0]}
            if counter > maxlen:
                maxlen = counter
            for j in range(1,len(clusters)): #we remove all labels but the first (for the graph)
                clusters[j] = ''
            [label_clusters.append(clusters[j]) for j in range(len(clusters))] #we make the final label list
            [total_lifetime.append(lifetime[j]) for j in range(len(lifetime))] #we make the final lifetime list
            [total_size.append(int(size[j])) for j in range(len(size))]
            data = [shortpop[i]] #we re-initialize data (other set of clusters)
            counter = 1 #initialization of the counter
    #same with last set of data
    size, clusters, lifetime = zip(*data) #we extract the lifetime of the selected data (i.e. only for the cluster type we worked until now)
    lifetime = [float(x) for x in lifetime] #convert str to int
    lifetime = [lifetime[i] for i in range(len(lifetime))] #multiply by the timestep to obtain lifetime in fs
    lifetime = sorted(lifetime, reverse=True) # we sort data by lifetime
    clusters = list(clusters) #convert to list
    size = list(size) #convert to list
    #all_data[str(clusters[0])]=lifetime
    all_data[str(clusters[0])]= {'lifetime':lifetime,'size':size[0]}
    if counter > maxlen:
        maxlen = counter
    for j in range(1,len(clusters)): #we remove all labels but the first (for the graph)
        clusters[j] = ''
    [label_clusters.append(clusters[j]) for j in range(len(clusters))] #we make the final label list
    [total_lifetime.append(lifetime[j]) for j in range(len(lifetime))] #we make the final lifetime list
    [total_size.append(int(size[j])) for j in range(len(size))]
    
    #***** Creation of the txt file for plotting elsewhere or for use in other scripts as plot_speciation-lifetime-comp-1fig3T
    try:
        print('max size of clusters is ', max(total_size))
        #print(label_clusters)
        #print(total_lifetime)
        
        filename = populfile.split('/')[-1]+'_L'+str(length)+'.txt'
        for key in all_data:
            while len(all_data[key]['lifetime']) < maxlen:
                all_data[key]['lifetime'].append('')
        with open(filename,'w') as nf:
            nf.write("\t".join(str(all_data[x]['size']) for x in sorted(all_data))+ "\n")    
            nf.write("\t".join(x for x in sorted(all_data))+ "\n")    
            for i in range(maxlen):
                nf.write("\t".join(str(all_data[cluster]['lifetime'][i]) for cluster in sorted(all_data))+ "\n")
        print(filename, 'is created')
    except ValueError:
        print("No cluster smaller than ",length, "atoms in file",populfile)
        sys.exit()
    
    #***** Creation of the plot
    #fig, ax = creation_plot()
    print("tot len",len(total_lifetime))
    #x = np.arange(len(total_lifetime))
    #plt.bar(x,total_lifetime, width = 1 )
    #print(label_clusters)
    #plt.xticks(x, label_clusters, rotation = 45, fontsize = 8)    
    #figurename = populfile.split('/')[-1]+'_L'+str(length)+'.png'
    #plt.savefig(figurename, bbox_inches='tight', dpi=150)
    #print(figurename, 'is created')
    #plt.show()

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



