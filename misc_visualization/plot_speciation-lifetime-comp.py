#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on  Jun  13 2018

@author: anais
Langage : Python3


        Bar plot of clusters lifetime for selected type of cluster !
        All the files data on the same figure
         *************  ARTICLE VERSION  *************

"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator,AutoMinorLocator
import crystallography as cr
import re
import natsort



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



def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format CaAl2Si2O8_T3_nvt_a12.5.
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.split('/')[-1]
    temperature = str(int(float(filename.split('_')[1].strip('T'))*1000))
    acell = filename.split('_')[3].split('.outcar')[0].strip('a')
    return temperature, acell



def creation_plot(xlabel):
    """     ********** Creation of the plot  **********    """
    plt.close()
    fig, ax = plt.subplots(figsize=(10,4))
    ax.set_ylabel(r'Cluster absolute lifetime (fs)', fontweight = 'bold', fontsize = 12)
    ax.set_xlabel(xlabel, fontweight = 'bold', fontsize = 12, labelpad = 25)
    #Adjustment of ticks
    ymajorLocator = AutoLocator()
    yminorLocator = AutoMinorLocator()
    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)
    ax.tick_params(which = 'both', labelsize = 10, width = 0.5)
    plt.tick_params(bottom = False, top = False, labelbottom = True)
    plt.autoscale(enable=True,axis='y',tight=True)
    ax.grid(True, which='major',axis = 'y', linestyle=':', linewidth=0.5 )
    return fig,ax



def main(argv):
    """     ********* Main program *********     """
    data1 = {} #dictionnary containing lifetimes for each cluster (file 1)
    data2 = {} #dictionnary containing lifetimes for each cluster (file 2)
    data3 = {} #dictionnary containing lifetimes for each cluster (file 3)
    cluster_lengths = {}#dictionnary containing length of cluster for each cluster of all files
    filename3 = ''
    atoms = 'all'
    letter = ''
    #other dictionnaries and parameters for the figure
    colors_T = {'2000':'#800080','3000':'#297fff','4000':'#00ff00','5000':'#ffcd01','6000':'#ff0101','7000':'#ff01de'}
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 20}
    #other parameters
    Na=6.022*10**23
    try:
        options,arg = getopt.getopt(argv,"hf:g:j:a:v:m:l:",["file1","gfile2","jfile3","atom","variable","mineralfile",'letter'])
    except getopt.GetoptError:
        print("plot_speciation-lifetime-comp.py -m <mineralfile with elements>  -a <list of first atoms determining the type of cluster to plot (default = 'all')> -v <variable (rho,T)> -l <letter, default = ''> -f <small_popul.txt filename> -g <small_popul.txt filename n°2> -j <small_popul.txt filename n°3 (option)>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_speciation-lifetime-comp.py program to  Plot lifetime barchart for 3 files (.txt files obtained from plot_speciation-lifetime.py) for selected clusters')
            print("plot_speciation-lifetime-comp.py -m <mineralfile with elements>  -a <list of first atoms determining the type of cluster to plot (default = 'all')> -v <variable (rho,T)> -l <letter, default = ''> -f <small_popul.txt filename> -g <small_popul.txt filename n°2> -j <small_popul.txt filename n°3 (option)> ")
            print('')
            print('requires the file containing elements and number (in order to compute the densities)')
            print('requires the files .txt created with the plot_speciation_lifetime.py script')
            print('')
            sys.exit()
        elif opt in ("-f", "--file1"):
            filename1 = str(arg)
        elif opt in ("-g", "--gfile2"):
            filename2 = str(arg)
        elif opt in ("-j", "--jfile3"):
            filename3 = str(arg)
        elif opt in ("-a","--atom"):
            atoms = arg.split(',')
        elif opt in ("-v", "--variable"):
            variable = str(arg)
        elif opt in ("-m","--mineralfile"):
            mineralfile = str(arg)
        elif opt in ('-l','letter'):
            letter = str(arg)
    #***** Calculation of the molecular mass
    with open(mineralfile,'r') as mf:
        entry = mf.readline()
        elements = entry.split()[1:]
        entry = mf.readline()
        number = entry.split()[1:]
    MN = 0
    for i in range(len(elements)):
        MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
    #***** Extraction of all cluster type from all files along with max size of data per cluster
    print("*********************** 1st step: extraction of data separately")
    with open(filename1,'r') as f:
        print('********** for file',filename1)
        line = f.readline() #we read the first line with cluster sizes
        all_length = line.split('\n')[0].split('\t')
        line = f.readline() #we read the second line with cluster names
        clusters=line.split('\n')[0].split('\t')
        #print('clusters in file:',clusters)
        #we take all the clusters
        if atoms == 'all':
            print("I use all atomic clusters")
            for i in range(len(clusters)):
                data1[clusters[i]] = [] #initialization of data 
                cluster_lengths[clusters[i]] = int(all_length[i])
        else: #or we select only the cluster names we want 
            print('I use only clusters starting by ', atoms)
            for atom in atoms:
                for i in range(len(clusters)):
                    if clusters[i][:len(atom)] == atom:
                        data1[clusters[i]] = [] #initialization of data 
                        cluster_lengths[clusters[i]] = int(all_length[i])
        #print(data1)
        while True:
            line = f.readline() #we read all the other lines with lifetime for each apparition
            if not line: break
            else:
                entry=line.split('\n')[0].split('\t')
                for i in range(0,len(clusters)): #we store the correct lifetime in the dictionnary for the corresponding key
                    if entry[i] != '':
                        try:
                            data1[clusters[i]].append(float(entry[i]))
                        except KeyError:
                            continue
                    else: continue
    with open(filename2,'r') as f:
        print('********** for file',filename2)
        line = f.readline() #we read the first line with cluster sizes
        all_length = line.split('\n')[0].split('\t')
        line = f.readline() #we read the second line with cluster names
        clusters=line.split('\n')[0].split('\t')
        #print('clusters in file:',clusters)
        #we take all the clusters
        if atoms == 'all':
            print("I use all atomic clusters")
            for i in range(len(clusters)):
                data2[clusters[i]] = [] #initialization of data 
                cluster_lengths[clusters[i]] = int(all_length[i])
        else: #or we select only the cluster names we want 
            print('I use only clusters starting by ', atoms)
            for atom in atoms:
                for i in range(len(clusters)):
                    if clusters[i][:len(atom)] == atom:
                        data2[clusters[i]] = [] #initialization of data 
                        cluster_lengths[clusters[i]] = int(all_length[i])
        #print(data1)
        while True:
            line = f.readline() #we read all the other lines with lifetime for each apparition
            if not line: break
            else:
                entry=line.split('\n')[0].split('\t')
                for i in range(0,len(clusters)): #we store the correct lifetime in the dictionnary for the corresponding key
                    if entry[i] != '':
                        try:
                            data2[clusters[i]].append(float(entry[i]))
                        except KeyError:
                            continue
                    else: continue
    if filename3 != '':
        with open(filename3,'r') as f:
            print('********** for file',filename3)
            line = f.readline() #we read the first line with cluster sizes
            all_length = line.split('\n')[0].split('\t')
            line = f.readline() #we read the second line with cluster names
            clusters=line.split('\n')[0].split('\t')
            #print('clusters in file:',clusters)
            #we take all the clusters
            if atoms == 'all':
                print("I use all atomic clusters")
                for i in range(len(clusters)):
                    data3[clusters[i]] = [] #initialization of data 
                    cluster_lengths[clusters[i]] = int(all_length[i])
            else: #or we select only the cluster names we want 
                print('I use only clusters starting by ', atoms)
                for atom in atoms:
                    for i in range(len(clusters)):
                        if clusters[i][:len(atom)] == atom:
                            data3[clusters[i]] = [] #initialization of data 
                            cluster_lengths[clusters[i]] = int(all_length[i])
            #print(data1)
            while True:
                line = f.readline() #we read all the other lines with lifetime for each apparition
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    for i in range(0,len(clusters)): #we store the correct lifetime in the dictionnary for the corresponding key
                        if entry[i] != '':
                            try:
                                data3[clusters[i]].append(float(entry[i]))
                            except KeyError:
                                continue
                        else: continue
    #creation of the list containing all the cluster once and only once
    clusters = []
    [clusters.append(key) for key in data1]
    for key in data2:
        if key not in clusters:
            clusters.append(key)
    if filename3 != '':
        for key in data3:
            if key not in clusters:
                clusters.append(key)
    clusters.sort()
    print(clusters)
    print("*********************** 2nd step: Addition of '0' data in order to have the same x axis ")
    #*********creation of the dictionnaries containing:
    #-sizes of the data
    #-lists of labels (with the cluster name and then voids)
    #-data completed by voids
    labels = {}
    sizes = {}
    size3 = 0
    for cluster in clusters:
        #***test to obtain the size of data list for current cluster
        try:
            size1 = len(data1[cluster])
            #print("#1:",cluster,size1)
        except KeyError:
            size1 = 0
            #print("#1:",cluster,size1)
        try:
            size2 = len(data2[cluster])
            #print("#2:",cluster,size2)
        except KeyError:
            size2 = 0
            #print("#2:",cluster,size2)
        if filename3 != '':
            try:
                size3 = len(data3[cluster])
                #print("#2:",cluster,size2)
            except KeyError:
                size3 = 0
        #the number we keep is the max of the three previous values
        sizes[cluster]=max(size1,size2,size3)
        #***now we can create the labels
        #either with simple method of 'label'+ voids
        labels[cluster]=[cluster]
        for i in range(sizes[cluster]-1):
            labels[cluster].append('') 
        #or by putting the label in the middle of the arrays
        #labels[cluster]=[]
        #if sizes[cluster] > 1:
        #    if sizes[cluster] > 2:
        #        for i in range(round(sizes[cluster]/2)-1):
        #            labels[cluster].append('')
        #        labels[cluster].append(cluster)
        #        for i in range(round(sizes[cluster]/2)-1):
        #            labels[cluster].append('')
        #    else:
        #        labels[cluster].extend([cluster,''])
        #else:
        #    labels[cluster]=[cluster]
        #***Now we add 'nan' to data array in order to have the same x axis
        try:
            if size1 < sizes[cluster]:
                for i in range(sizes[cluster]-size1-1):
                    data1[cluster].append(float('nan'))
                data1[cluster].append(0)
        except KeyError:
            data1[cluster] = []
            for i in range(sizes[cluster]-1):
                data1[cluster].append(float('nan'))
            data1[cluster].append(0)
        try:
            if size2 < sizes[cluster]:
                for i in range(sizes[cluster]-size2-1):
                    data2[cluster].append(float('nan'))
                data2[cluster].append(0)
        except KeyError:
            data2[cluster] = []
            for i in range(sizes[cluster]-1):
                data2[cluster].append(float('nan'))
            data2[cluster].append(0)
        if filename3 != '':
            try:
                if size3 < sizes[cluster]:
                    for i in range(sizes[cluster]-size3-1):
                        data3[cluster].append(float('nan'))
                    data3[cluster].append(0)
            except KeyError:
                data3[cluster] = []
                for i in range(sizes[cluster]-1):
                    data3[cluster].append(float('nan'))
                data3[cluster].append(0)
    #print(labels)
    #print("sizes of lists",sizes)
    print("*********************** 3rd step: Creation of labels and arrays to plot")
    #creation of the big list of data (concatenation) and labels
    totdata1 = [] 
    totdata2 = []
    totdata3 = []
    totlabels = []    
    #if we want to sort clusters by size and name
    sorted_cluster_lengths = [(v[0],v[1]) for v in natsort.natsorted(cluster_lengths.items(), key=lambda  kv: (kv[1], kv[0]))] #[v[0] for v in sorted(d.items(), key=lambda kv: (-kv[1], kv[0]))] #  natsort.natsorted(cluster_lengths.items(), key=lambda kv: kv[1])
    print(sorted_cluster_lengths)
    for i in range(len(sorted_cluster_lengths)):
        cluster = sorted_cluster_lengths[i][0]
        totdata1.extend(data1[cluster])
        totdata2.extend(data2[cluster])
        if filename3 != '':
            totdata3.extend(data3[cluster])
        totlabels.extend(labels[cluster])
        #print(cluster, sizes[cluster])
    #formatage of labels
    for j in range(len(totlabels)):
        if totlabels[j] != '':
            totlabels[j] = format1label(totlabels[j])
    #print(totlabels)
    if (len(totdata1) == 0) and (len(totdata2) == 0) and (len(totdata3) == 0):
        print("There is nothing to plot for clusters of type ",atoms )
        sys.exit()
    #creation of x vector
    totsize = 0
    for key in sizes:
        totsize = totsize + sizes[key]
    x = np.arange(totsize)
    #creation of custom labels
    temp1, acell1 = split_name(filename1)
    temp2, acell2 = split_name(filename2)
    if filename3 != '':
        temp3, acell3 = split_name(filename3)
    if variable == 'rho':
        var1='a'+acell1
        label1 = str(round(MN/(Na*float(acell1)**3*10**(-24)),2))+' g.cm$^{-3}$'
        print(label1)
        var2='a'+acell2
        label2 = str(round(MN/(Na*float(acell2)**3*10**(-24)),2))+' g.cm$^{-3}$'
        print(label2)
        if filename3 != '':
            var3='a'+acell3
            label3 = str(round(MN/(Na*float(acell3)**3*10**(-24)),2))+' g.cm$^{-3}$'
            print(label3)
        fixedvar = temp1
        title = filename1.split('_')[0]
        #formatage of the mineral name
        i=0
        while True:
            if i <= len(title)-1:
                if re.match('[0-9]',title[i]):
                    title = title[:i]+'$_{'+title[i]+'}$' + title[i+1:]
                    i = i+5
                i = i+1
            else:break
        title = title +' at '+fixedvar+' K'
    else:
        var1=temp1+'K'
        label1= str(temp1)+' K'
        var2=temp2+'K'
        label2= str(temp2)+' K'
        if filename3 != '':
            var3=temp3+'K'
            label3= str(temp3)+' K'
        fixedvar = str(round(MN/(Na*float(acell2)**3*10**(-24)),2))
        title = filename1.split('_')[0]
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
    print("*********************** 4th step: Plot everything on the same figure")
    #*********** Plot
    #Creation of the plot
    speciation_type = filename1.split('umd.dat.')[1].split('.popul.dat')[0]
    if speciation_type == 'r1':
        xlabel = 'Chemical species'
    else:
        xlabel = 'Coordinating polyhedra'
    fig, ax = creation_plot(xlabel)
    #if speciation_type == 'r0': #because of resolution problems (too many bars so they are too thin to appear) we should use a line plot filled for r0
    #    print("speciation 0")
    #    plt.fill_between(x,totdata1, y2 =0, color=colors_T[temp1], alpha = 1, label=label1)
    #    plt.fill_between(x,totdata2, y2 =0, color=colors_T[temp2], alpha = 0.5, label=label2)
    #    if filename3 != '':
    #        plt.fill_between(x,totdata3, y2=0, color=colors_T[temp3], alpha = 0.35, label=label3)   
    #else: #but for r1 we keep the standard bar plot (less bars to plot)    
    #    print("speciation 1")
    #    plt.bar(x,totdata1, width = 1, color=colors_T[temp1], alpha = 1, label=label1)
    #    plt.bar(x,totdata2, width = 1, color=colors_T[temp2], alpha = 0.5, label=label2)    
    #    if filename3 != '':
    #        plt.bar(x,totdata3, width = 1, color=colors_T[temp3], alpha = 0.35, label=label3)    
    plt.plot(x, totdata1, '-', color=colors_T[temp1],  label=label1)
    plt.plot(x, totdata2, '-', color=colors_T[temp2],  label=label2)
    if filename3 != '':
        plt.plot(x, totdata3, '-', color=colors_T[temp3],  label=label3)
    plt.xticks(x, totlabels, rotation = 45, fontsize = 10)
    ax.yaxis.set_ticks_position('both')
    if letter != '':
        ax.text(0.007,0.94, letter , transform=ax.transAxes, horizontalalignment = 'left', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], bbox=dict(facecolor='none', edgecolor='k', pad=3.0))
    #plt.title(title, fontsize = 12, fontweight = 'bold')
    #plt.legend(loc='upper right', fontsize = 12)
    #save plot
    if filename3 != '':
        figurename = 'comparison_lifetime_'+filename1.split('/')[-1].split('_')[0]+'_'+speciation_type+'_'+fixedvar+'_'+var1+'-'+var2+'-'+var3
    else:
        figurename = 'comparison_lifetime_'+filename1.split('/')[-1].split('_')[0]+'_'+speciation_type+'_'+fixedvar+'_'+var1+'-'+var2
    figurename = figurename + '.svg'
    plt.savefig(figurename, bbox_inches='tight', dpi=300)
    print(figurename, 'is created')
    #plt.show()    
  
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])





