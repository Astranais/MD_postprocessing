#!/usr/bin/env python
# -*- coding: utf-8 -*-

###
#AUTHORS: ANAIS KOBSCH, NATALIA SOLOMATOVA
###

#*********** Importation of the packages and modules used here ************
import glob
import os
import sys
import getopt
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button




def headerfile(firstfile, dirpath):
    """creation of the newfile with correct header"""
    allpairs = {}
    allpairs_ordered = []
    firstline = ['pair']  #beginning of the first line of the file gofr
    secondline = ['file']
    # creation of the header from the first line of the first file
    with open(firstfile, 'r') as f:
        line = f.readline()
    line = line.strip('dist')
    line = re.sub('(Int\([A-Za-z-]*\))', ' ',line).split()        #I substitute all Int(...)  by a blankspace and split using spaces
    for ii in range(len(line)):
        allpairs[line[ii]] = ii*2 + 1
        allpairs_ordered.append(line[ii])
    for pair in allpairs_ordered:
        for i in range(0,5):
            firstline.append(pair)
        secondline.extend(['xmax','ymax','xmin','coord','bond'])
    newfilename = dirpath+'_gofrs.txt' 
    print('The file ',newfilename,' is created')
    f = open(newfilename,'w')
    f.write("\t".join(x for x in firstline)+ "\n")
    f.write("\t".join(x for x in secondline)+ "\n")
    return f, allpairs, allpairs_ordered     #I return the newly created files f along with the list of element allpairs

def average_bond(radius,gofr,xmin_gofr):
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


def interactive_fit(data,file,allpairs, guess_xmax, guess_xmin, distance, bonds, pair):
    """fit max and min of data using interactive plot for one pair of atoms"""
    #print(pair, allpairs[pair], allpairs[pair] + 1)
    #********* First extraction of data
    gofr, intgofr = np.loadtxt(file,
                               usecols=(allpairs[pair], allpairs[pair] +1),
                               skiprows=1, unpack=True)
    #we compute also the Int(g(r)) of the reverse pair
    reversepair = pair.split('-')[1]+'-'+pair.split('-')[0]
    intgofr_reverse = np.loadtxt(file,
                   usecols=(allpairs[reversepair] + 1),
                   skiprows=1, unpack=True)
                    
    #********** Second we cut the data and work separately on the max and the min
    orig_size = np.size(distance)  # size of our data
    nonzeros = np.count_nonzero(gofr)
    max_index = np.argmax(gofr)  # find the index of max gofr
    last_index = -1
    cutdistance = distance[max_index:last_index]  # cut off beginning for finding minimum
    cutgofr = gofr[max_index:last_index]  # cut off beginning for finding minimum
    cut_size = np.size(cutdistance)  # size of data
    cut_min_index = np.argmin(cutgofr)  # find the min of gofr in the cut segment
    min_index = orig_size - cut_size + cut_min_index - 1  # min index of the *uncut* gofr
    end_imin = int(min_index + min_index * 0.08)  # 8% to the right
    end_imax = int(max_index + max_index * 0.12)  # 12% to the left

    if nonzeros > 10  and end_imin < (orig_size - 1): #and end_imax < (orig_size - 1):
        plt.close()
        def creation_fig_xmax():
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.subplots_adjust(bottom=0.2)
            ax.set_title(file.split('.gofr.dat')[0] + '  ' + pair)
            plt.xlabel('Distance ($\AA$)')
            plt.ylabel('g(r)')
            plt.grid()
            ax.plot(distance, gofr, 'o', markersize=2)  # plot the fitted max
            return fig
        
        fig = creation_fig_xmax()
        
        def onclick(event):
            global xmax_click
            xmax_click = event.xdata  # record x value of click
            plt.close()
        
        fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        
        try:
            max_index = (np.abs(distance - xmax_click)).argmin()  # index of the closest value to clicked value
        except TypeError:
            print('Please do not click outside the plot area, try again to click on the max value')
            fig = creation_fig_xmax()
            fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()
            try:
                max_index = (np.abs(distance - xmax_click)).argmin()
            except TypeError:
                print('I told you, now I use 0 as clicked value')
                max_index = (np.abs(distance - 0)).argmin() 

        #****** determine start and end cols for fitting the polynomial for the max region:
        end_imax = int(max_index + max_index * 0.12)  # 12% to the left (arbitrary)
        start_imax = int(max_index - max_index * 0.08)  # 8% to the right (arbitrary)
        if end_imax > (orig_size - 1):  # for rare cases where the max distance is out of bounds
            end_imax = orig_size - 1  # subtract one because cols start at 0

        # define the regions that need to be fitted for the max area:
        fit_max_dist = distance[int(start_imax):int(end_imax + 1)]  # ranges do not include last element, so add one
        fit_max_gofr = gofr[int(start_imax):int(end_imax + 1)]

        order = 3  # order of the polynomials
        max_fit = np.poly1d(np.polyfit(fit_max_dist, fit_max_gofr, order))
        xaxis_max = np.linspace(distance[int(start_imax)], distance[int(end_imax)],
                                1000)  # create new x values for the fitted max region
        yaxis_max = max_fit(xaxis_max)  # create new y values for fitted max region
        poly_max_i = np.argmax(yaxis_max)  # determine index of max y value

        
        def creation_fig_xmin():
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.subplots_adjust(bottom=0.2)
            ax.set_title(file.split('.gofr.dat')[0] + '  ' + pair)
            ax.plot(distance, gofr, 'o', xaxis_max, max_fit(xaxis_max), '-', markersize=2)  # plot the fitted max
            plt.xlabel('Distance ($\AA$)')
            plt.ylabel('g(r)')
            plt.grid()
            return fig

        fig = creation_fig_xmin()

        def onclick(event):
            global xmin_click
            xmin_click = event.xdata  # record x value of click
            plt.close()

        fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()

        try:
            min_index = (np.abs(distance - xmin_click)).argmin()  # index of the closest value to clicked value
        except TypeError:
            print('Please do not click outside the plot area, try again to click on the min value')
            fig = creation_fig_xmin()
            fig.canvas.mpl_connect('button_press_event', onclick)
            plt.show()
            try:
                min_index = (np.abs(distance - xmin_click)).argmin()  # index of the closest value to clicked value
            except TypeError:
                print('I told you, now I use 0 as clicked value')
                min_index = (np.abs(distance - 0)).argmin() 


        #****** determine start and end cols for fitting the polynomial for the min region:
        start_imin = int(min_index - min_index * 0.08)  # 8% to the left (arbitrary)
        end_imin = int(min_index + min_index * 0.08)  # 8% to the right (arbitrary)
        if end_imin > (orig_size - 1):  # for rare cases where the min distance is out of bounds
            end_imin = orig_size - 1  # subtract one because cols start at 0
        # define the regions that need to be fitted for min area using the clicked min:
        fit_min_dist = distance[int(start_imin):int(end_imin + 1)]  # ranges do not include last element, so add one
        fit_min_gofr = gofr[int(start_imin):int(end_imin + 1)]
        intgofr = intgofr[int(start_imin):int(end_imin + 1)]
        intgofr_reverse = intgofr_reverse[int(start_imin):int(end_imin + 1)]

        min_fit = np.poly1d(np.polyfit(fit_min_dist, fit_min_gofr, order))
        intgofr_fit = np.poly1d(np.polyfit(fit_min_dist, intgofr, order))
        intgofr_fit_reverse = np.poly1d(np.polyfit(fit_min_dist, intgofr_reverse, order))
        xaxis_min = np.linspace(distance[int(start_imin)], distance[int(end_imin)],
                                1000)  # create new x values for the fitted min region
        yaxis_min = min_fit(xaxis_min)  # create new y values for fitted min region
        poly_min_i = np.argmin(yaxis_min)  # determine index of min y value
        yaxis_intgofr = intgofr_fit(xaxis_min)
        yaxis_intgofr_reverse = intgofr_fit_reverse(xaxis_min)

        xmax_gofr = xaxis_max[int(poly_max_i)]  # x value of the max (x,y) of the max region of gofr
        ymax_gofr = yaxis_max[int(poly_max_i)]  # max y value of the max region of gofr
        xmin_gofr = xaxis_min[int(poly_min_i)]  # x value of the min (x,y) of the min region of gofr
        y_intgofr = yaxis_intgofr[int(poly_min_i)]  # y value of the integral of gofr that corresponds to the min x value of gofr
        y_intgofr_reverse = yaxis_intgofr_reverse[int(poly_min_i)]  # y value of the integral of gofr that corresponds to the min x value of gofr


        class Decision(object):
            def good(self, event):
                global value_array
                plt.close(fig)  # close and move to next plot
                print('xmax=', round(xmax_gofr, 2), ' ymax=', round(ymax_gofr, 2), ' xmin=', round(xmin_gofr, 2), ' coord=', round(y_intgofr, 2))
                value_array = (xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, y_intgofr_reverse)
                return value_array

            def bad(self, event):
                global value_array
                plt.close(fig)  # close and move to next plot
                print('bad')
                xmax_gofr = 0
                ymax_gofr = 0
                xmin_gofr = 0
                y_intgofr = 0
                y_intgofr_reverse = 0
                value_array = (xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, y_intgofr_reverse)
                return value_array

        fig = plt.figure()  # show figure with the newly fitted minimum curve
        ax = fig.add_subplot(111)
        plt.subplots_adjust(bottom=0.2)
        ax.set_title(file.split('.gofr.dat')[0] + '  ' + pair)
        ax.plot(distance, gofr, 'o', xaxis_max, max_fit(xaxis_max), '-', xaxis_min, min_fit(xaxis_min), 'r-',
                markersize=2)
        plt.xlabel('Distance ($\AA$)')
        plt.ylabel('g(r)')

        callback = Decision()
        axbad = plt.axes([0.12, 0.05, 0.1, 0.075])
        axgood = plt.axes([0.8, 0.05, 0.1, 0.075])
        bgood = Button(axgood, 'Good')
        bgood.on_clicked(callback.good)
        bbad = Button(axbad, 'Bad')
        bbad.on_clicked(callback.bad)
        plt.show()

        xmax_gofr = value_array[0]
        ymax_gofr = value_array[1]
        xmin_gofr = value_array[2]
        y_intgofr = value_array[3]
        y_intgofr_reverse = value_array[4]

    else: # simply output 0's if the data is no good
        print('skip')
        xmax_gofr = 0
        ymax_gofr = 0
        xmin_gofr = 0
        y_intgofr = 0
        y_intgofr_reverse = 0
    
    try: 
        bond = average_bond(distance,gofr,xmin_gofr) #we compute the real bond length using the selected xmin value
    except ZeroDivisionError: #happens if xmin = 0
        bond = 0
        
    #*********** Last step, we save in the dictionnary the resultsfor the current pair of atoms
    data[pair] = [str(xmax_gofr),str(ymax_gofr),str(xmin_gofr),str(y_intgofr), str(bond)]
    data[reversepair] = [str(xmax_gofr),str(ymax_gofr),str(xmin_gofr),str(y_intgofr_reverse), str(bond)]
    bonds[pair] = str(xmin_gofr)
    bonds[reversepair] = str(xmin_gofr)
    guess_xmax[pair] = xmax_gofr                  #update the initial guesses for the next acell
    guess_xmin[pair] = xmin_gofr
    
    cutdistance = 0  # re initialization of variables for the next loop
    fit_min_dist = 0
    fit_max_dist = 0
    fit_min_gofr = 0
    fit_max_gofr = 0
    gofr = 0
    intgofr = 0
    intgofr_reverse
    return (data, bonds, guess_xmin, guess_xmax)



xmin_click = []
xmax_click = []
def analyze_gofrs_interactive(data,file,allpairs, guess_xmax, guess_xmin, atoms, bonds):
    """ extraction of the min, max and coordination number using fits and interactive plot """
    distance = np.loadtxt(file, usecols=(0,), skiprows = 1, unpack=True)
    for pair in allpairs:
        for atom_pair in atoms:
            if (pair == atom_pair):
                #interactive fit only for the selected allpairs of atoms
                print(pair, '\t --> interactive fit')
                data, bonds, guess_xmin, guess_xmax = interactive_fit(data,file,allpairs, guess_xmax, guess_xmin, distance, bonds, pair)
    return(data,bonds, guess_xmax, guess_xmin)






def analyze_gofrs_automatic(data,file,allpairs, guess_xmax, guess_xmin, atoms, bonds):
    """ extraction of the min, max and coordination number using fits and previous values as initial guesses """
    distance = np.loadtxt(file, usecols=(0,), skiprows = 1, unpack=True)
    for pair in allpairs:
        for atom_pair in atoms:
            if (pair == atom_pair):
                if (guess_xmax[pair] == 0 or guess_xmin[pair] == 0):  
                    #if we haven't succeeded to have initial guesses, then we try again with interactive plot
                    print(pair, '\t --> interactive fit')
                    data, bonds, guess_xmin, guess_xmax = interactive_fit(data,file,allpairs, guess_xmax, guess_xmin, distance, bonds, pair)
                else: #if we have correct initial guesses for this pair of atoms, then we use an automatic fitting process
                    print(pair, '\t --> automatic fit')
                    #print(pair,allpairs[pair],allpairs[pair]+1)
                    #********* First extraction of data
                    gofr, intgofr = np.loadtxt(file,
                                               usecols=(allpairs[pair], allpairs[pair] +1),
                                               skiprows=1, unpack=True)
                    #we compute also the Int(g(r)) of the reverse pair
                    reversepair = pair.split('-')[1]+'-'+pair.split('-')[0]
                    data[reversepair] = []
                    intgofr_reverse = np.loadtxt(file,
                                   usecols=(allpairs[reversepair] + 1),
                                   skiprows=1, unpack=True)
                    
    
                    #********** Second we cut the data and work separately on the max and the min
                    orig_size = np.size(distance)
                    max_index = (np.abs(distance - guess_xmax[pair])).argmin()  # index of the closest value to guessed value
                    min_index = (np.abs(distance - guess_xmin[pair])).argmin()  # index of the closest value to guessed value
        
                    # determine start and end cols for fitting the polynomial:
                    start_imin = int(min_index - min_index * 0.08) # 8% to the left (arbitrary)
                    end_imin = int(min_index + min_index * 0.08)  # 8% to the right (arbitrary)
                    if end_imin > (orig_size - 1):  # for rare cases where the max distance is out of bounds
                        end_imin = orig_size - 1  # subtract one because cols start at 0
                    end_imax = int(max_index + max_index * 0.12) # 12% to the left (arbitrary)
                    start_imax = int(max_index - max_index * 0.08) # 8% to the right (arbitrary)
                    if end_imax > (orig_size - 1):  # for rare cases where the max distance is out of bounds
                        end_imax = orig_size - 1 # subtract one because cols start at 0
        
                    # define the regions that need to be fitted:
                    # ranges do not include last element, so add one
                    fit_max_dist = distance[int(start_imax):int(end_imax + 1)]
                    fit_min_dist = distance[int(start_imin):int(end_imin + 1)]
                    fit_max_gofr = gofr[int(start_imax):int(end_imax + 1)]
                    fit_min_gofr = gofr[int(start_imin):int(end_imin + 1)]
                    intgofr = intgofr[int(start_imin):int(end_imin + 1)]
                    intgofr_reverse = intgofr_reverse[int(start_imin):int(end_imin + 1)]
        
                    nonzeros = np.count_nonzero(gofr)
                    # fit nth order polynomials to max and min regions
                    # only include gofrs with sufficient data and that doesn't go out of bounds
                    if nonzeros > 10 and end_imax < (orig_size - 1) and end_imin < (orig_size - 1):
                        order = 4
                        max_fit = np.poly1d(np.polyfit(fit_max_dist, fit_max_gofr, order))
                        min_fit = np.poly1d(np.polyfit(fit_min_dist, fit_min_gofr, order))
                        intgofr_fit = np.poly1d(np.polyfit(fit_min_dist, intgofr, order))
                        intgofr_fit_reverse = np.poly1d(np.polyfit(fit_min_dist, intgofr_reverse, order))
                        xaxis_max = np.linspace(distance[int(start_imax)], distance[int(end_imax)], 1000)  # create new x values for the fitted max region
                        xaxis_min = np.linspace(distance[int(start_imin)], distance[int(end_imin)],1000)  # create new x values for the fitted min region
                        yaxis_max = max_fit(xaxis_max)  # create new y values for fitted max region
                        yaxis_min = min_fit(xaxis_min)  # create new y values for fitted min region
                        yaxis_intgofr = intgofr_fit(xaxis_min)
                        yaxis_intgofr_reverse = intgofr_fit_reverse(xaxis_min)
                        poly_max_i = np.argmax(yaxis_max)  # determine index of max y value
                        poly_min_i = np.argmin(yaxis_min)  # determine index of min y value
                        
                        xmax_gofr = xaxis_max[int(poly_max_i)]  # x value of the max (x,y) of the max region of gofr
                        ymax_gofr = yaxis_max[int(poly_max_i)]  # max y value of the max region of gofr
                        xmin_gofr = xaxis_min[int(poly_min_i)]  # x value of the min (x,y) of the min region of gofr
                        y_intgofr = yaxis_intgofr[int(poly_min_i)]  # y value of the integral of gofr that corresponds to the min x value of gofr
                        y_intgofr_reverse = yaxis_intgofr_reverse[int(poly_min_i)]  # y value of the integral of gofr that corresponds to the min x value of gofr
                        
                        plt.plot(distance, gofr, 'o', xaxis_max, max_fit(xaxis_max), '-', xaxis_min, min_fit(xaxis_min), 'r-', markersize=2)
                               
                    else: # simply output 0's if the data is no good
                        xmax_gofr = 0
                        ymax_gofr = 0
                        xmin_gofr = 0
                        y_intgofr = 0
                        y_intgofr_reverse = 0
                    
                    bond = average_bond(distance,gofr,xmin_gofr) #we compute the real bond length using the selected xmin value
                    
                    #*********** Last step, we save in the dictionnary the results for the current pair of atoms
                    data[pair] = [str(xmax_gofr),str(ymax_gofr),str(xmin_gofr),str(y_intgofr),str(bond)]
                    data[reversepair] = [str(xmax_gofr),str(ymax_gofr),str(xmin_gofr),str(y_intgofr_reverse),str(bond)]
                    bonds[pair] = str(xmin_gofr)
                    bonds[reversepair] = str(xmin_gofr)
                    guess_xmax[pair] = xmax_gofr                  #update the initial guesses for the next acell
                    guess_xmin[pair] = xmin_gofr
        
                    #re initialization of variables for the next loop
                    fit_min_dist = 0
                    fit_max_dist = 0
                    fit_min_gofr = 0
                    fit_max_gofr = 0
                    gofr = 0
                    intgofr = 0 
                    intgofr_reverse = 0   
    return(data,bonds,guess_xmax, guess_xmin)



def main(argv):
    """     ********* Main program *********     """
    atoms = []
    elements = []
    bondfile = 1
    try:
        options,arg = getopt.getopt(argv,"ha:b:",["atoms","bondfile"])
    except getopt.GetoptError:
        print("analyze_gofr.py -a <pairs of atoms (option)>(ex: 'Ca-O,Ca-Ca', default includes all) -b < 0 or 1 (1= write bond file, 0 do not write bond file)>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('*******')
            print('analyze_gofr.py program to extract all relevant data from the gofr.dat files created by the script gofrs.py and write them into temperature_gofrs.txt files')
            print("analyze_gofr.py -a <pairs of atoms (option)>(ex: 'Ca-O,Ca-Ca', default includes all) -b < 0 or 1 (1= write bond file, 0 do not write bond file)>")
            print('WARNING! for efficiency, your data needs to be located in different folders for each T and you should launch this script from the folder containing every subfolder temperature')
            print('WARNING! you have to click FIRST on the maximum, and THEN on the minimum')
            print(' ')
            print('This code produces a gofrs.txt file with:')
            print('    - a column with filename')
            print('    - 5 column per A-B pair of atoms')
            print('        - x value of max(gofr)')
            print('        - y value associated')
            print('        - x value of the first min(gofr)')
            print('        - coordination # = y value of int(gofr) associated to previous x value')
            print('        - average bond length = [int(r gofr(r)]/[int(gofr(r))] up to the 1st xmin')
            print('Default -b = 1,  the code write a bond file per gofr.dat file with the x value of the first minimum for each pair of atoms in a format that can be used by the speciation scripts ')
            print(' ')
            sys.exit()
        elif opt in ("-a", "--atoms"):
            atoms = arg.split(',')                      #list of atom allpairs we want to analyze here
        elif opt in ('-b', '--bondfile'):
            bondfile = int(arg)
    for dirpath, dirnames, filenames in os.walk(os.curdir):
        files = sorted(glob.glob(dirpath+'/*.gofr.dat')) 
        print("**** files in the current directory are:", files)
        if files != []:
            f, allpairs, allpairs_ordered = headerfile(files[0], dirpath)                          #I create the first newfile for gofr and save the dictionnary of element allpairs with their column n°
            print("all pairs of atoms are",allpairs_ordered)
            #if we want to analyze all the atoms, then we create a list with only half the pair (no duplicates since g(r_A-B) = g(r_B-A) )
            if atoms == []:
                #first we extract the individual elements
                for pair in allpairs:
                    elem = pair.split('-')[0]
                    if elem not in elements:
                        elements.append(elem) 
                #second we create the A-B allpairs
                for i in range(len(elements)):
                    for j in range(i,len(elements)):
                        atoms.append(elements[i]+'-'+elements[j])
            print("selected pairs of atoms are",atoms)
            interactive = 0
            guess_xmax = {}                                                     #dictionnaries for initial guesses (key = pair of atoms, value = xmin or max value)
            guess_xmin = {}
            for file in files:
                print("analyze file",file)
                data = {}
                bonds = {}
                newbondsfile = file.split('.gofr.dat')[0]+'.bonds.inp' 
                if bondfile == 1:
                    b = open(newbondsfile,'w')  #I create the newfile for bonds
                if file == files[0]:
                    interactive = 1
                if interactive == 1:
                    #we compute the min,max etc. from the gofr and int using fit and interactive plot
                    data, bonds, guess_xmax, guess_xmin = analyze_gofrs_interactive(data,file,allpairs, guess_xmax, guess_xmin, atoms, bonds)
                    interactive = 0
                else:
                    #we compute the min,max etc. from the gofr and int using automatic fit based on previous values for initial guesses
                    data, bonds, guess_xmax, guess_xmin = analyze_gofrs_automatic(data,file,allpairs, guess_xmax, guess_xmin, atoms, bonds)       
                #we write in the file the results
                #first we complete by X the values not computed 
                for pair in allpairs:
                    #print(pair)
                    if pair not in data:
                        print(pair, '\t --> not computed')
                        data[pair] = ['0','0','0','0','0']
                        bonds[pair] = '0'               
                #then we write all the results        
                f.write(file.split('/')[-1])
                for pair in allpairs_ordered:
                    f.write("\t")
                    f.write("\t".join(x for x in data[pair]))
                f.write('\n')
                #we write the bond file
                if bondfile == 1:
                    for pair in allpairs_ordered:
                        b.write(pair.split('-')[0]+'\t'+pair.split('-')[1]+'\t'+bonds[pair]+'\n')
                    b.close()
                    print("bond file ",newbondsfile, ' created') 
            f.close()
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



