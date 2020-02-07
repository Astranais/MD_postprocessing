#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 13:58:46 2020

@author: akobsch
"""
#    ********* Importation of the packages and modules used here *********    
import sys
import getopt
import glob
import os
import numpy as np


def headerfile(allpairs_ordered, dirpath):
    """creation of the newfile with correct header"""
    firstline = ['pair']  #beginning of the first line of the file gofr
    secondline = ['file']
    for pair in allpairs_ordered:
        for i in range(0,4):
            firstline.append(pair)
        secondline.extend(['xmax','average','median','bondr3'])
    newfilename = dirpath+'_bonds.txt' 
    print('The file ',newfilename,' is created')
    f = open(newfilename,'w')
    f.write("\t".join(x for x in firstline)+ "\n")
    f.write("\t".join(x for x in secondline)+ "\n")
    return f   #I return the newly created files f

def average_bond(radius,gofr,xmin_gofr,ymax_gofr):
    """compute weighted average, same as [int(r gofr(r)]/[int(gofr(r))] up to the 1st xmin"""
    ii = 0
    average = 0
    sumweight = 0
    while radius[ii] < xmin_gofr:
        average += gofr[ii]/ymax_gofr * radius[ii]
        sumweight += gofr[ii]/ymax_gofr
        ii+=1
    return average/sumweight

def compute_bondr1(radius,gofr,xmin_gofr):
    """compute [int(r gofr(r)]/[int(gofr(r))] up to the 1st xmin"""
    #this is the weighted average
    ii = 0
#    dr = radius[2]-radius[1]
    rgofr = 0.0
    intgofr=0.0
    while radius[ii] < xmin_gofr:
        rgofr +=  radius[ii] * gofr[ii] # * dr
        intgofr += gofr[ii]   # * dr
        ii += 1
    return rgofr/intgofr

def compute_bondr2(radius,gofr,xmin_gofr):
    """compute [int(r2 gofr(r)]/[int(r gofr(r))] up to the 1st xmin"""
    ii = 0
#    dr = radius[2]-radius[1]
    r2gofr = 0.0
    rgofr=0.0
    while radius[ii] < xmin_gofr:
        r2gofr +=  radius[ii]**2 * gofr[ii] # * dr
        rgofr +=  radius[ii] * gofr[ii] # * dr
        ii += 1
    return r2gofr/rgofr

def compute_bondr3(radius,gofr,xmin_gofr):
    """compute [int(r3 gofr(r)]/[int(r2 gofr(r))] up to the 1st xmin"""
    ii = 0
#    dr = radius[2]-radius[1]
    r3gofr = 0.0
    r2gofr=0.0
    while radius[ii] < xmin_gofr:
        r3gofr +=  radius[ii]**3 * gofr[ii] # * dr
        r2gofr +=  radius[ii]**2 * gofr[ii] # * dr
        ii += 1
    return r3gofr/r2gofr

def median_bond(radius,gofr,xmin_gofr):
    """compute the median bond from the gofr up to the 1st xmin"""
    ii = 0
    dr = radius[2]-radius[1]
    intgofr=0.0
    cumulint=[]
    while radius[ii] < xmin_gofr:
        intgofr += gofr[ii] * dr
        cumulint.append(intgofr)
        ii += 1
    cumulint = np.asarray(cumulint)
    return radius[np.abs(cumulint - intgofr/2).argmin()]



def main(argv):
    """     ********* Main program *********     """
    atoms = []
    allpairs = {}
    allpairs_ordered = []
    data = {}
    directory = os.curdir
    try:
        options,arg = getopt.getopt(argv,"hg:a:d:",["gofrsfile","atoms","directory"])
    except getopt.GetoptError:
        print("bond_analysis.py -g <gofrs.txt filename> -a <pair of atoms>(ex: 'Ca-O,Ca-Ca') -d <directory where the gofr.dat files are located (option)>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('bond_analysis.py program to compute several estimation of the average bond length and write them in a bonds.txt file')
            print("bond_analysis.py -g <gofrs.txt filename> -a <pair of atoms>(ex: 'Ca-O,Ca-Ca') -d <directory where the gofr.dat files are located (option but required if you have different gofr.dat file with same name in different subfolders)>")
            print('')
            sys.exit()
        elif opt in ("-g", "--gofrsfile"):
            gofrfile = str(arg)
        elif opt in ("-a", "--atoms"):
            atoms = arg.split(',')                      #list of atom pairs we want to analyze here)
        elif opt in ('-d','--directory'):
            directory = str(arg)
    #********** 1st step: extraction of pairs and data from the gofrs.txt file
    #read the header and extract the pairs
    skip_head = 0
    with open(gofrfile, 'r') as f:
        while True:
            line = f.readline()
            if not line: break
            else:
                entry = line.split('\n')[0].split('\t')
                skip_head +=1
                if entry[0] == 'pair':
                    allpairs_ordered = [entry[1]]        #we initialize the pairs list
                    for ii in range(2,len(entry)):       #we fill the pairs list with atom pairs occuring only once
                        if entry[ii] != entry[ii-1]:
                            allpairs_ordered.append(entry[ii])
                    for ii in range(len(allpairs_ordered)):
                        allpairs[allpairs_ordered[ii]] = ii #and we add the index for each pair in the dictionnay allpairs with atom pairs as key
                elif entry[0] == 'file':
                    break
    #fill the data dictionnary with all the data from gofrs
    with open(gofrfile, 'r') as f:
        [f.readline() for ii in range(skip_head)]
        while True:
            line = f.readline()
            if not line: break
            else:
                entry = line.split('\n')[0].split('\t')
                data[entry[0].split('/')[-1]] = entry[1:]
    print('all pairs available are ', allpairs)
    print("we analyze", atoms)
    #********** 2nd step: for each gofr.dat file in the directories below the directory indicated we test if they are in the gofr.txt file
    # if yes, then we compute the different bonds
    # if no we move to the next one
    for dirpath, dirnames, filenames in os.walk(directory):        
        files = sorted(glob.glob(dirpath+'/*.gofr.dat')) #I list every gofr files in alphabetic order
        print('in folder ', dirpath, 'the matching files are', files)
        if files != []:
            f = headerfile(allpairs_ordered, dirpath)     
            for file in files:
                if file.split('/')[-1] in data:
                    print(file)
                    allbonds = {}
                    #extraction of data
                    for pair in allpairs:
                        for atom_pair in atoms:
                            if (pair == atom_pair):
                                #we load the data of the gofr.dat files
                                if pair == 'O2':
                                    distance, gofr, int_gofr = np.loadtxt(file, usecols=(0, allpairs['O-O']*2 + 1, allpairs['O-O'] *2 + 2), skiprows=1, unpack=True) 
                                    print('analysis', pair)
                                else:
                                    distance, gofr, int_gofr = np.loadtxt(file, usecols=(0, allpairs[pair]*2 + 1, allpairs[pair] *2 + 2), skiprows=1, unpack=True) 
                                    #we compute also the reverse pair
                                    reversepair = pair.split('-')[1]+'-'+pair.split('-')[0]
                                    print('analysis', pair, 'and', reversepair)
                                #we remove the last 0 value
                                distance = np.delete(distance,[len(distance)-1]) 
                                gofr = np.delete(gofr,[len(gofr)-1])
                                int_gofr = np.delete(int_gofr,[len(int_gofr)-1])  
                                #we load the computed data from the gofrs.txt file
                                try:
                                    xmax = float(data[file.split('/')[-1]][allpairs[pair]*5])
                                    xmin = float(data[file.split('/')[-1]][allpairs[pair]*5 + 2])
                                    bondr3 = float(data[file.split('/')[-1]][allpairs[pair]*5 + 4])
                                    median = median_bond(distance,gofr,xmin)
                                    try: 
                                        average = compute_bondr1(distance,gofr,xmin)
                                    except ZeroDivisionError: #happens if xmin = 0
                                        average = 0
                                except ValueError:
                                    print(pair, '\t --> nan encountered')
                                    xmax= 0
                                    bondr3 = 0
                                    average = 0
                                    median = 0
                                allbonds[pair]=[str(xmax),str(average),str(median),str(bondr3)]
                                if pair != 'O2':
                                    allbonds[reversepair]=[str(xmax),str(average),str(median),str(bondr3)]
                    #we write in the file the results
                    #first we complete by X the values not computed 
                    for pair in allpairs_ordered:
                        #print(pair)
                        if pair not in allbonds:
                            print(pair, '\t --> not computed')
                            allbonds[pair] = ['0','0','0','0']
                    #then we write all the results        
                    f.write(file.split('/')[-1])
                    for pair in allpairs_ordered:
                        f.write("\t")
                        f.write("\t".join(x for x in allbonds[pair]))
                    f.write('\n')
            f.close()
            print(f,'written')
 
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])
 