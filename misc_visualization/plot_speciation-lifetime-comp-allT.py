#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on  Jun  13 2018

@author: anais
Langage : Python3


       program to  Plot lifetime barchart for all temperatures at a selected acell,  one figure per file
       need to be launched from a folder containing all the T files we want to plot

"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator,AutoMinorLocator
import matplotlib.patches as mpatches
import crystallography as cr
import re
import os
import glob
import natsort




def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./CaAl2Si2O8_T3_nvt_a12.5.
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.strip('./')
    temperature = str(int(float(filename.split('_')[1].strip('T'))*1000))
    acell = filename.split('_')[3].split('.outcar')[0].strip('a')
    return temperature, acell

def empty_plot(fig, ax, file, ymax, letter, colors_T, speciation, length, atomslist):
    legend_labels = {}
    temp, acell = split_name(file)
    label= str(temp)+' K'
    ax.set_ylim(0,ymax)
    ax.yaxis.set_ticks_position('both')
    plt.text(0.01,0.95, letter , transform=ax.transAxes, horizontalalignment = 'left', fontweight = 'bold', fontsize = 12, bbox=dict(facecolor='none', edgecolor='k', pad=3.0))      
    legend_labels[label] =  mpatches.Patch(color=colors_T[temp])
    s = [(k, legend_labels[k]) for k in sorted(legend_labels.keys(),reverse = False)]
    plt.legend([v for k,v in s],[k for k,v in s],loc='upper right', fontsize = 12)
    plt.subplots_adjust(top = 0.97, bottom = 0.26, right = 0.94, left = 0.12, hspace = 0, wspace = 0)
    plt.tick_params(#
        axis='both',
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        left = False, 
        labelbottom=False, # labels along the bottom edge are off
        labelleft = False) 
    figurename = 'lifetime_allT_'+file.split('.outcar.umd.dat')[0]+'_r'+speciation+'_L'+length+'_'+atomslist+'.png'
    plt.savefig(figurename, dpi=150)#, bbox_inches='tight') #remove 'tight' to take into account options of subplots_adjust
    print(figurename, 'is created')


def creation_plot(speciation):
    """     ********** Creation of the plot  **********    """
    plt.close()
    fig, ax = plt.subplots(figsize=(8,5))
    if speciation == '1':
        xlabel = 'Chemical species'
    else:
        xlabel = 'Coordinating polyhedra'
    ax.set_ylabel(r'Cluster absolute lifetime (fs)', fontweight = 'bold', fontsize = 12)
    ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = 12)
    ax.xaxis.set_label_coords(0.5, -0.3)
    ax.yaxis.set_label_coords(-0.1, 0.5)
    #Adjustment of ticks
    ymajorLocator = AutoLocator()
    yminorLocator = AutoMinorLocator()
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which = 'both', labelsize = 10, width = 0.5)
    #plt.autoscale(enable=True,axis='y',tight=True)
    #ax.set_ylim(0,770) #for NaAlSi3O8 a19.0
    #ax.set_ylim(0,60) #for NaAlSi3O8 a15.0
    plt.tick_params(bottom = False, top = False, labelbottom = True)
    ax.grid(True, which='major',axis = 'y', linestyle=':', linewidth=0.5  )
    return fig,ax



def main(argv):
    """     ********* Main program *********     """
    files = []
    data = {} #dictionnary containing lifetimes for each cluster 
    cluster_lengths = {}#dictionnary containing length of cluster for each cluster 
    sizes = {}
    labels = {}
    note = {} #dictionnary to add a note on the figure in case of max lifetime = simulation time
    #other dictionnaries and parameters for the figure
    colors_T = {'2000':'#800080','3000':'#297fff','4000':'#00ff00','4500':'#bae200','5000':'#ffcd01','6000':'#ff0101','6500':'#ff00a2','7000':'#ff01de'}
    atoms = 'all'
    atomslist = 'all'
    letters = ['a','b','c','d','e','f','g','h','i']
    #other parameters
    Na=6.022*10**23
    try:
        options,arg = getopt.getopt(argv,"hm:a:f:l:",["mineralfile","atom","filename","lifetime"])
#        options,arg = getopt.getopt(argv,"hc:m:a:r:l:",["cell","mineralfile",'atoms','rspeciation',"length"])
    except getopt.GetoptError:
        print("plot_speciation-lifetime-comp-allT.py  -m <mineralfile with elements> -a <list of first atom determining the type of cluster to plot  (default = 'all')>  -f <one filename of the same format name than all we want to plot> -l <list of max lifetime of all the simu>  ")
#        print("plot_speciation-lifetime-comp-allT.py  -m <mineralfile with elements> -a <list of first atom determining the type of cluster to plot  (default = 'all')>  -c <cell we want to plot>  -r <speciation type (0 or 1)> -l <max length of clusters>  ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_speciation-lifetime-comp-allT.py program to  Plot lifetime barchart for all temperatures at a selected acell for some cluster,  one figure per file')
            print("plot_speciation-lifetime-comp-allT.py  -m <mineralfile with elements> -a <list of first atom determining the type of cluster to plot  (default = 'all')>  -f <one filename of the same format name than all we want to plot> -l <list of max lifetime of all the simu>  ")
#            print("plot_speciation-lifetime-comp-allT.py  -m <mineralfile with elements> -a <list of first atom determining the type of cluster to plot (default = 'all')>  -c <cell we want to plot>  -r <speciation type (0 or 1)> -l <max length of clusters> ")
            print("")
            print('requires the file containing elements and number (in order to compute the densities)')
            print('requires to launch the script from the directory containing all the files .txt created with the plot_speciation_lifetime.py script')
            print('')
            sys.exit()
#        elif opt in ("-c", "--cell"):
#            cell = str(arg)
        elif opt in ("-m","--mineralfile"):
            mineralfile = str(arg)
        elif opt in ("-a","--atoms"):
            atomslist = str(arg)
            atoms = arg.split(',')                      #list of atom couples we want to analyze here
        elif opt in ("-f","--filename"):
            firstfilename = str(arg)
        elif opt in ("-l","--lifetime"):
            max_lifetimes = list(map(int,arg.split(',')))
#        elif opt in ("-r","--rspeciation"):
#            speciation = str(arg)
#        elif opt in ("-l","--length"):
#            length = str(arg)
    #***** Calculation of the molecular mass
    with open(mineralfile,'r') as mf:
        entry = mf.readline()
        elements = entry.split()[1:]
        entry = mf.readline()
        number = entry.split()[1:]
    MN = 0
    for i in range(len(elements)):
        MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
    #***** List of all the files concerned by the  T and acell
    cell = firstfilename.split('.outcar.umd.dat')[0].split('_')[3].strip('a')
    speciation = firstfilename.split('.outcar.umd.dat.r')[-1].split('.popul.dat')[0]
    length = firstfilename.split('.popul.dat_L')[-1].split('.txt')[0]
    filename = '*_a'+cell+'*outcar.umd.dat.r'+speciation+'.popul.dat_L'+length+'*.txt'
    print('I search for files of type', filename)
    allfiles = sorted(glob.glob(filename))                             #I list every population file created by plot_speciation-lifetime.py in alphabetic order
    #print(files)
    print('I use the corresponding limit lifetimes', max_lifetimes)
    #we remove the empty files (the ones created by hand) from the file list
    for file in allfiles:
        if os.stat(file).st_size != 0:
            files.append(file)
    #***** Extraction of all cluster type from all files along with max size of data per cluster
    print("*********************** 1st step: calculation of total number of bars in order to have the same x axis")
    for file in files:
        print('********** for file',file)
        with open(file,'r') as f:
            line = f.readline() #we read the first line with cluster sizes
            all_length = line.split('\n')[0].split('\t')
            print(all_length)
            line = f.readline() #we read the second line with cluster names
            clusters=line.split('\n')[0].split('\t')
            print('clusters in file:',clusters)
            #we take all the clusters
            if atoms == 'all':
                print("I use all atomic clusters")
                for i in range(len(clusters)):
                    data[clusters[i]] = [] #initialization of data 
                    cluster_lengths[clusters[i]] = all_length[i]
            else: #or we select only the cluster names we want 
                print('I use only clusters starting by ', atoms)
                for atom in atoms:
                    for i in range(len(clusters)):
                        if clusters[i][:len(atom)] == atom:
                            print(clusters[i])
                            data[clusters[i]] = [] #initialization of data 
                            cluster_lengths[clusters[i]] = all_length[i]
            #extraction of data, all in the same data dictionnary in order to know the total (cumulative) number of bars in the bar plot
            while True:
                line = f.readline() #we read all the other lines with lifetime for each apparition
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    for i in range(0,len(clusters)): #we store the correct lifetime in the dictionnary for the corresponding key
                        if entry[i] != '':
                            try:
                                data[clusters[i]].append(float(entry[i]))
                            except KeyError:
                                continue
                        else: continue
        #print('data',data) 
        #test to find the maximum value for data sizes (for each cluster)
        if atoms == 'all':
            for cluster in clusters:
                try:
                    if len(data[cluster]) > sizes[cluster]:
                        sizes[cluster] = len(data[cluster])
                except KeyError:
                    sizes[cluster] = len(data[cluster])
        else:
            for atom in atoms:
                for cluster in clusters:
                    if cluster[:len(atom)] == atom:
                        try:
                            if len(data[cluster]) > sizes[cluster]:
                                sizes[cluster] = len(data[cluster])
                        except KeyError:
                            sizes[cluster] = len(data[cluster])
        #print("sizes",sizes)
    #************ Extraction of the data and determination of the max of all figures (to keep the y axis identical)
    print("*********************** 2nd step: Extraction of data per file and deterination of ymax for all fig")
    ymax = 0
    for file in files: 
        print('********** for file',file)
        letter = letters[allfiles.index(file)]
        print('letter is', letter)
        data = {} #initialization of data
        with open(file,'r') as f:
            line = f.readline() #we read the first line with cluster sizes
            line = f.readline() #we read the second line with cluster names
            clusters=line.split('\n')[0].split('\t')
            #print('clusters in file:',clusters)
            #we take all the clusters
            if atoms == 'all':
                for i in range(len(clusters)):
                    data[clusters[i]] = [] #initialization of data 
            else: #or we select only the cluster names we want    
                for atom in atoms:
                    for i in range(len(clusters)):
                        if clusters[i][:len(atom)] == atom:
                            data[clusters[i]] = [] #initialization of data 
            #extraction of data
            while True:
                line = f.readline() #we read all the other lines with lifetime for each apparition
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    for i in range(0,len(clusters)): #we store the correct lifetime in the dictionnary for the corresponding key
                        if entry[i] != '':
                            try:
                                data[clusters[i]].append(float(entry[i]))
                            except KeyError:
                                continue
                        else: continue
        #we complete the data lists with nan in  order to have the same x axis
        for cluster in sizes:
            try:
                if len(data[cluster]) < sizes[cluster]:
                    for i in range(sizes[cluster]-len(data[cluster])):
                        data[cluster].append(float('nan'))
            except KeyError:
                data[cluster] = []
                for i in range(sizes[cluster]):
                    data[cluster].append(float('nan'))
        #********* Creation of the dictionnary containing the lists of labels (with the cluster name and then voids)   
        for cluster in sizes: 
            print('cluster',cluster)
            print(data[cluster])
            if str(data[cluster][0]) != 'nan': #we replace labels by voids when the specie is not there
                labels[cluster]=[cluster]
                for i in range(sizes[cluster]-1):
                    labels[cluster].append('') 
            else:
                labels[cluster]=['']
                for i in range(sizes[cluster]-1):
                    labels[cluster].append('') 
        #********* Creation of the big list of data (concatenation) and labels sorted by length of clusters
        totdata = [] 
        totlabels = []
        sorted_cluster_lengths = natsort.natsorted(cluster_lengths.items(), key=lambda kv: kv[1])
        print(sorted_cluster_lengths)
        for i in range(len(sorted_cluster_lengths)):
            cluster = sorted_cluster_lengths[i][0]
            totdata.extend(data[cluster])
            totlabels.extend(labels[cluster])
        print('maximum of data',np.nanmax(totdata))
        if np.nanmax(totdata) < max_lifetimes[files.index(file)] :
            if np.nanmax(totdata) > ymax:
                ymax = np.nanmax(totdata)
            note[file] = ''
        else:
            note[file] = 'maximum lifetime \n = simulation time' 
    print("*********************************** ymax =",ymax)   
    #************ Same as before with plot, each on a different figure (one fig per T)
    print("*********************** 3rd step: Extraction of data per file and plot each graph on a different figure")    
    for file in files: 
        print('********** for file',file)
        letter = letters[allfiles.index(file)]
        print('letter is', letter)
        data = {} #initialization of data
        with open(file,'r') as f:
            line = f.readline() #we read the first line with cluster sizes
            line = f.readline() #we read the second line with cluster names
            clusters=line.split('\n')[0].split('\t')
            #print('clusters in file:',clusters)
            #we take all the clusters
            if atoms == 'all':
                for i in range(len(clusters)):
                    data[clusters[i]] = [] #initialization of data 
            else: #or we select only the cluster names we want    
                for atom in atoms:
                    for i in range(len(clusters)):
                        if clusters[i][:len(atom)] == atom:
                            data[clusters[i]] = [] #initialization of data 
            #extraction of data
            while True:
                line = f.readline() #we read all the other lines with lifetime for each apparition
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    for i in range(0,len(clusters)): #we store the correct lifetime in the dictionnary for the corresponding key
                        if entry[i] != '':
                            try:
                                data[clusters[i]].append(float(entry[i]))
                            except KeyError:
                                continue
                        else: continue
        #we complete the data lists with nan in  order to have the same x axis
        for cluster in sizes:
            try:
                if len(data[cluster]) < sizes[cluster]:
                    for i in range(sizes[cluster]-len(data[cluster])):
                        data[cluster].append(float('nan'))
            except KeyError:
                data[cluster] = []
                for i in range(sizes[cluster]):
                    data[cluster].append(float('nan'))
        #********* Creation of the dictionnary containing the lists of labels (with the cluster name and then voids)   
        for cluster in sizes: 
            if str(data[cluster][0]) != 'nan': #we replace labels by voids when the specie is not there
                labels[cluster]=[cluster]
                for i in range(sizes[cluster]-1):
                    labels[cluster].append('') 
            else:
                labels[cluster]=['']
                for i in range(sizes[cluster]-1):
                    labels[cluster].append('') 
        #********* Creation of the big list of data (concatenation) and labels sorted by length of clusters
        totdata = [] 
        totlabels = []
        sorted_cluster_lengths = natsort.natsorted(cluster_lengths.items(), key=lambda kv: kv[1])
        print(sorted_cluster_lengths)
        for i in range(len(sorted_cluster_lengths)):
            cluster = sorted_cluster_lengths[i][0]
            totdata.extend(data[cluster])
            totlabels.extend(labels[cluster])     
        #formatage of labels
        for j in range(len(totlabels)):
            if totlabels[j] != '':
                i=0
                while True:
                    #print('i=',i, 'len totlabels -1 = ',len(totlabels[j])-1 )
                    if i <= len(totlabels[j])-1:
                        #print('totlabels j i = ', totlabels[j][i])
                        if re.match('_',totlabels[j][i]):
                            num=0
                            try:
                                while re.match('[0-9]',totlabels[j][i+1+num]):
                                    num +=1
                            except IndexError: #index error when we arrive at the end of the cluster name
                                pass
                                #print('end of the cluster')
                            totlabels[j] = totlabels[j][:i]+'$_{'+totlabels[j][i+1:i+1+num]+'}$' + totlabels[j][i+1+num:]
                            #print('new totlabels j =',totlabels[j])
                            i = i+5
                        i = i+1
                    else:break
        #print(totlabels)
        #*********** Plot
        #Creation of the plot
        fig, ax = creation_plot(speciation)
        if (len(totdata) == 0) :
            print("There is nothing to plot for clusters of type ",atoms, 'in file ', file )
            #save empty plot
            empty_plot(fig, ax, file, ymax, letter, colors_T, speciation, length, atomslist)
            continue
        #creation of x vector
        totsize = 0
        for key in sizes:
            totsize = totsize + sizes[key]
        x = np.arange(totsize)
        print('len of x= ', len(x))
        #creation of custom labels
        temp, acell = split_name(file)
        label= str(temp)+' K'
        fixedvar = str(round(MN/(Na*float(acell)**3*10**(-24)),2))
        title = file.split('_')[0]
        #formatage of the mineral name
        i=0
        while True:
            if i <= len(title)-1:
                if re.match('[0-9]',title[i]):
                    num=0
                    try:
                        while re.match('[0-9]',title[i+num]):
                            num +=1
                    except IndexError: #index error when we arrive at the end of the cluster name
                        pass
                        #print('end of the cluster')
                    title = title[:i]+'$_{'+title[i:i+num]+'}$' + title[i+num:]
                    i = i+5
                i = i+1
            else:break
        title = title+' at '+fixedvar+' g.cm$^{-3}$'
        if speciation == '0' or len(x)>1000: #because of resolution problems (too many bars so they are too thin to appear) we should use a line plot filled for r0
            plt.fill_between(x,totdata, y2 =0, color=colors_T[temp], alpha = 1, label=label)  
        else: #but for r1 we keep the standard bar plot (less bars to plot)
            plt.bar(x,totdata, width = 1, color=colors_T[temp], alpha = 1, label=label)  
        plt.xticks(x, totlabels, rotation = 45, fontsize = 10)
        width=20
        ax.set_xlim(-width,len(x)+width)
        ax.set_ylim(0,ymax)
        ax.yaxis.set_ticks_position('both')
        #plt.title(title, fontsize = 12, fontweight = 'bold')
        plt.text(0.01,0.95, letter , transform=ax.transAxes, horizontalalignment = 'left', fontweight = 'bold', fontsize = 12, bbox=dict(facecolor='none', edgecolor='k', pad=3.0))  
        plt.legend(loc='upper right', fontsize = 12)
        plt.subplots_adjust(top = 0.97, bottom = 0.26, right = 0.94, left = 0.12, hspace = 0, wspace = 0)
        #addition of the note if it exists
        if note[file] != '':
            plt.text(0.7,0.7, note[file] , transform=ax.transAxes, horizontalalignment = 'left', fontsize = 12)  
        #save plot
        figurename = 'lifetime_allT_'+file.split('.outcar.umd.dat')[0]+'_r'+speciation+'_L'+length+'_'+atomslist+'.png'
        plt.savefig(figurename, dpi=150)#, bbox_inches='tight') #remove 'tight' to take into account options of subplots_adjust
        print(figurename, 'is created')
        #plt.show()    

    #For files that does not exist and for which we created a corresponding empty file, we create an empty figure and remove them of the list of files
    for file in allfiles:
        if os.stat(file).st_size == 0:
            letter = letters[allfiles.index(file)]
            fig, ax = creation_plot(speciation)
            empty_plot(fig, ax, file, ymax, letter, colors_T, speciation, length, atomslist)
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])





