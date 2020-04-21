#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 12:32:18 2019

@author: akobsch



        program to  Plot lifetime barchart for all clusters and T,  one subfigure per cluster, one figure per speciation_r1.popul.dat file
        
        need to be launched from a folder containing all the T files we want to plot
"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import re
import crystallography as cr



def max_axis(number):
    number = int(number)
    if number < 10:
        maxnumber = 10
    else:
        if number < 100:
            ndigits = 1
        else:
            ndigits = len(str(number))-2
        maxnumber = np.ceil( number /  10**ndigits ) * 10**ndigits
    return maxnumber


def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.strip('./')
    temperature = str(int(float(filename.split('_')[1].strip('T'))*1000))
    acell = filename.split('_')[3].split('.outcar')[0].strip('a')
    return temperature, acell

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

def format_label(list_species):
    """ formatage of all labels of a list """
    for j in range(len(list_species)):
        if list_species[j] != '':
            list_species[j] = format_1label(list_species[j])
    return list_species


def main(argv):
    """     ********* Main program *********     """
    files = []
    all_lifetimes = {} #dictionnary containing lifetimes for each cluster 
    all_length = {}#dictionnary containing length of cluster for each cluster 
    all_clusters_sizes = [] #list with all the cluster sizes accross all files
    all_species = [] #list with all the species accross all files
    nsubfig = np.zeros(13,dtype=int) #list of number of subfig per natoms cluster
    order_species={}  #dictionnary to order all the species in each column
    for ii in range(1,14):
        order_species[ii]=[]
    #parameters for the figures depending on the output format (presentation or article)
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 20}
    colors_T = {'2000':'#800080','3000':'#297fff','4000':'#00ff00','4500':'#bae200','5000':'#ffcd01','6000':'#ff0101','6500':'#ff00a2','7000':'#ff01de'}
    letters = ['a','b','c','d','e','f','g','h','i']
    #other parameters
    Na=6.022*10**23
    try:
        options,arg = getopt.getopt(argv,"hm:f:",["mineralfile","filename"])
    except getopt.GetoptError:
        print("plot_speciation-lifetime-r1.py  -m <mineralfile with elements> -f <one filename of the same format name than all we want to plot>  ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_speciation-lifetime-r1.py program to  Plot lifetime barchart for all clusters smaller then 13 atoms,  one subfigure per cluster,  one figure per speciation_r1.popul.dat file created by the script speciation_lifetime.py')
            print("plot_speciation-lifetime-r1.py  -m <mineralfile with elements> -f <one filename of the same format name than all we want to plot>  ")
            print("")
            print('requires the file containing elements and number (in order to compute the densities)')
            print('')
            sys.exit()
        elif opt in ("-m","--mineralfile"):
            mineralfile = str(arg)
        elif opt in ("-f","--filename"):
            firstfilename = str(arg)
    #***** Calculation of the molecular mass
    with open(mineralfile,'r') as mf:
        entry = mf.readline()
        elements = entry.split()[1:]
        entry = mf.readline()
        number = entry.split()[1:]
    MN = 0
    for i in range(len(elements)):
        MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
    #***** 1st step : List of all the files concerned by the  T and acell
    cell = firstfilename.split('.outcar.umd.dat')[0].split('_')[3].strip('a')
    speciation = firstfilename.split('.outcar.umd.dat.r')[-1].split('.popul.dat')[0]
    length = firstfilename.split('.popul.dat_L')[-1].split('.txt')[0]
    filename = '*_a'+cell+'*outcar.umd.dat.r'+speciation+'.popul.dat_L'+length+'*.txt'
    print('I search for files of type', filename)
    allfiles = sorted(glob.glob(filename))                             #I list every population file created by plot_speciation-lifetime.py in alphabetic order
    #print('all files are:',allfiles)
    #we remove the empty files (the ones created by hand) from the file list
    for file in allfiles:
        if os.stat(file).st_size != 0:
            files.append(file)
    #print('selected files are:',files)
    #***** 2nd step : Extraction of all cluster type from all files  in order to have max number of rows per columns
    for file in files:
        print('************* for file',file)
        with open(file,'r') as f:
            #**** Extract species and lengths
            line = f.readline() #we read the first line with cluster sizes
            clusters_sizes = line.split('\n')[0].split('\t')
            line = f.readline() #we read the second line with cluster names
            clusters=line.split('\n')[0].split('\t')
            #print('clusters in file:',clusters)
            for ii in range(len(clusters)):
                all_length[clusters[ii]] = int(clusters_sizes[ii])
        #***** Count the number of subfigure per natoms cluster and Create dictionnary to order all the species in each column
        #update of the list containing each cluster size
        for species in clusters:
            if species not in all_species:
                all_clusters_sizes.append(clusters_sizes[clusters.index(species)])
        #update of the newcount based on all the cluster that have been seen up to now
        for ii in range(1,14):
            newcount = all_clusters_sizes.count(str(ii))
            if nsubfig[ii-1] < newcount:
                nsubfig[ii-1] = newcount
        for species in clusters:
            if species in order_species[all_length[species]]: continue
            else:
                order_species[all_length[species]].append(species)
            all_species.append(species)
    #creation of the list containing each species ONCE
    all_species = []
    for ii in range(1,14):
        for species in order_species[ii]:
            all_species.append(species)
    #print('total number of subfig:',nsubfig)
    print('order of species',order_species)
    #print('all the species are',all_species)
    #**** 3rd step : Core of the script = extraction of data and plot     
    for file in files:
        print('************* for file',file)
        letter = letters[allfiles.index(file)]
        #**** 3.1) Extraction of T and acell
        temperature, acell = split_name(file)
        #**** 3.2) Extract all data
        with open(file,'r') as f:
            #**** Extract species and lengths
            line = f.readline() #we read the first line with cluster sizes
            clusters_sizes = line.split('\n')[0].split('\t')
            line = f.readline() #we read the second line with cluster names
            clusters=line.split('\n')[0].split('\t')
            for ii in range(len(clusters)):
                all_length[clusters[ii]] = int(clusters_sizes[ii])
                all_lifetimes[clusters[ii]] = []
            #**** Extract lifetimes
            while True:
                line = f.readline() #we read all the other lines with lifetime for each apparition
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    for ii in range(0,len(clusters)): #we store the correct lifetime in the dictionnary for the corresponding key
                        if entry[ii] != '':
                            all_lifetimes[clusters[ii]].append(float(entry[ii]))
                        else: continue
        #**** 3.3) plot
        #fig 1 has the clusters up to 7 atoms
        print("creation of plot 1 with",max(nsubfig[0:7]),'subplots in rows')
        fig1 = plt.figure(figsize=(7*2,max(nsubfig[0:7])*1.9), constrained_layout=False)
        gs1 = fig1.add_gridspec(nrows=max(nsubfig[0:7]), ncols=7, width_ratios =np.ones(7), height_ratios =np.ones(max(nsubfig[0:7])), wspace = 0.5, hspace = 0.5, figure=fig1)
        #fig 2 has the clusters from 8 to 13 atoms
        print("creation of plot 2 with",max(nsubfig[7:]),'subplots in rows')
        fig2 = plt.figure(figsize=(7*2,max(nsubfig[7:])*1.9), constrained_layout=False)
        gs2 = fig2.add_gridspec(nrows=max(nsubfig[7:]), ncols=7,  width_ratios =np.ones(7), height_ratios =np.ones(max(nsubfig[7:])), wspace = 0.5, hspace = 0.5,  figure=fig2)
        for species in all_species:
            if all_length[species] < 8:                
                #print(species, order_species[all_length[species]].index(species), all_length[species]-1 )
                #plot only if there is something to draw
                if species in clusters:
                    #creation subplot
                    ax = fig1.add_subplot(gs1[order_species[all_length[species]].index(species), all_length[species]-1])
                    #text
                    ax.text(0.85,0.85, format_1label(species) , transform=ax.transAxes, horizontalalignment = 'right', fontsize = plot_parameters["size_fonts"])
                    if order_species[all_length[species]].index(species) ==  0:
                        ax.text(0.5,1.1, str(all_length[species]) , transform=ax.transAxes, horizontalalignment = 'center', fontsize = plot_parameters["size_fonts"], fontweight = 'bold')            
                    #plot
                    ax.bar(np.arange(len(all_lifetimes[species])),all_lifetimes[species],  width = 1, align = 'edge', color = colors_T[temperature])
                    #axis limits (pretty, rounded up)
                    maxnumberx = max_axis(len(all_lifetimes[species]))
                    ax.set_xticks([0, maxnumberx])
                    ax.set_xlim([0, maxnumberx])
                    maxnumbery = max_axis(max(all_lifetimes[species]))
                    ax.set_yticks([0,maxnumbery]) 
                    ax.set_ylim([0,maxnumbery])          
                    ax.tick_params(which = 'both', labelsize = plot_parameters["size_font_ticks"], width = plot_parameters["size_lines"]/2)               
                #empty subplot otherwise
                else:
                    ax = fig1.add_subplot(gs1[order_species[all_length[species]].index(species), all_length[species]-1], frameon=False)                    
                    ax.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)     
            else:
                #plot only if there is something to draw
                if species in clusters:
                    #creation subplot
                    ax = fig2.add_subplot(gs2[order_species[all_length[species]].index(species), all_length[species]-8])
                    #text
                    ax.text(0.85,0.85, format_1label(species) , transform=ax.transAxes, horizontalalignment = 'right', fontsize = plot_parameters["size_fonts"])
                    if order_species[all_length[species]].index(species) ==  0:
                        ax.text(0.5,1.1, str(all_length[species]) , transform=ax.transAxes, horizontalalignment = 'center', fontsize = plot_parameters["size_fonts"], fontweight = 'bold')
                    #plot
                    ax.bar(np.arange(len(all_lifetimes[species])),all_lifetimes[species],  width = 1, align = 'edge', color = colors_T[temperature])
                    #axis limits (pretty, rounded up)
                    maxnumberx = max_axis(len(all_lifetimes[species]))
                    ax.set_xticks([0, maxnumberx])
                    ax.set_xlim([0, maxnumberx])
                    maxnumbery = max_axis(max(all_lifetimes[species]))
                    ax.set_yticks([0,maxnumbery]) 
                    ax.set_ylim([0,maxnumbery])          
                    ax.tick_params(which = 'both', labelsize = plot_parameters["size_font_ticks"], width = plot_parameters["size_lines"]/2) 
                else:
                    ax = fig2.add_subplot(gs2[order_species[all_length[species]].index(species), all_length[species]-8], frameon=False)
                    ax.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False) 
        #**** 4th step : Add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
        ax0 = fig1.add_subplot(gs1[:,:], frameon=False)
        ax0bis = fig2.add_subplot(gs2[:,:], frameon=False)
        for ax in [ax0,ax0bis]:
            ax.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
            ax.set_xlabel(r'Number of species occurence', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], labelpad = plot_parameters["shift_labelpad"]*1.2)
            ax.set_ylabel(r'Lifetime (fs)', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], labelpad = plot_parameters["shift_labelpad"]*2.5)
            ax.text(-0.1,1, letter , transform=ax.transAxes, horizontalalignment = 'left', fontweight = 'bold', fontsize = 12, bbox=dict(facecolor='none', edgecolor='k', pad=3.0))  
        #**** 5th step: Save plot
        figurename1 = file.split('/')[-1].split('.txt')[0]+'_matrix1.pdf'
        fig1.savefig(figurename1, dpi=300, bbox_inches='tight') #remove 'tight' to take into account options of subplots_adjust
        figurename2 = file.split('/')[-1].split('.txt')[0]+'_matrix2.pdf'
        fig2.savefig(figurename2, dpi=300, bbox_inches='tight') #remove 'tight' to take into account options of subplots_adjust
        print(figurename1, 'is created')
        print(figurename2, 'is created')
        #plt.show()  
    #**** 6th step: For files that does not exist and for which we created a corresponding empty file, we create an empty figure and remove them of the list of files
    for file in allfiles:
        if os.stat(file).st_size == 0:
            temperature, acell = split_name(file)
            letter = letters[allfiles.index(file)]
            #fig 1 has the clusters up to 7 atoms
            fig1 = plt.figure(figsize=(7*2,max(nsubfig[0:7])*1.9), constrained_layout=False)
            gs1 = fig1.add_gridspec(nrows=max(nsubfig[0:7]), ncols=7, width_ratios =np.ones(7), height_ratios =np.ones(max(nsubfig[0:7])), wspace = 0.5, hspace = 0.5, figure=fig1)
            #fig 2 has the clusters from 8 to 13 atoms
            fig2 = plt.figure(figsize=(7*2,max(nsubfig[7:])*1.9), constrained_layout=False)
            gs2 = fig2.add_gridspec(nrows=max(nsubfig[7:]), ncols=7,  width_ratios =np.ones(7), height_ratios =np.ones(max(nsubfig[7:])), wspace = 0.5, hspace = 0.5,  figure=fig2)
            #**** 3rd step : Add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
            ax0 = fig1.add_subplot(gs1[:,:], frameon=False)
            ax0bis = fig2.add_subplot(gs2[:,:], frameon=False)
            for ax in [ax0,ax0bis]:
                ax.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
                ax.set_xlabel(r'Number of species occurence', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], labelpad = plot_parameters["shift_labelpad"]*1.2)
                ax.set_ylabel(r'Lifetime (fs)', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], labelpad = plot_parameters["shift_labelpad"]*2.5)
                ax.text(-0.1,1, letter , transform=ax.transAxes, horizontalalignment = 'left', fontweight = 'bold', fontsize = 12, bbox=dict(facecolor='none', edgecolor='k', pad=3.0))  
            
            
    
    
    
    
    
    
    
    
    
    
    
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:]) 
    