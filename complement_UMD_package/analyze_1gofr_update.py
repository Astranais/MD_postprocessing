#!/usr/bin/env python
# -*- coding: utf-8 -*-

###
#AUTHORS: ANAIS KOBSCH, NATALIA SOLOMATOVA
###

#*********** Importation of the packages and modules used here ************
import sys
import os
import getopt
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button


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


def interactive_fit(data,gofrfile,allpairs,distance, bonds, pair, method):
    """fit max and min of data using interactive plot for one pair of atoms"""
    skipdata = 0
    decision = 'c'
    #********* First extraction of data
    #we also consider the case when we added the columns O2 at the end of  the full gofrs data file 
    if pair == 'O2':
        #print('O2 = O-O', allpairs['O-O']*2 + 1, allpairs['O-O']*2 +  2)
        gofr, intgofr = np.loadtxt(gofrfile,
                               usecols=(allpairs['O-O']*2 + 1, allpairs['O-O']*2 + 2),
                               skiprows=1, unpack=True)

    else:
        #print(pair, allpairs[pair]*2 + 1, allpairs[pair]*2 +  2)
        gofr, intgofr = np.loadtxt(gofrfile,
                               usecols=(allpairs[pair]*2 + 1, allpairs[pair]*2 + 2),
                               skiprows=1, unpack=True)
        #we compute also the Int(g(r)) of the reverse pair
        reversepair = pair.split('-')[1]+'-'+pair.split('-')[0]
        intgofr_reverse = np.loadtxt(gofrfile,
                   usecols=(allpairs[reversepair]*2 + 2),
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

    if nonzeros > 10 : #and end_imin < (orig_size - 1): #and end_imax < (orig_size - 1):
        plt.close()
        def creation_fig_xmax():
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.subplots_adjust(bottom=0.2)
            ax.set_title(gofrfile.split('.gofr.dat')[0] + '  ' + pair)
            plt.xlabel('Distance ($\AA$)')
            plt.ylabel('g(r)')
            plt.grid()
            ax.plot(distance, gofr, 'o', markersize=2)  # plot the fitted max
            if pair == 'O2':
                major_xticks_zoom = np.arange(0, 3, 0.5) 
                minor_xticks_zoom = np.arange(0, 3, 0.1)
                major_yticks_zoom = np.arange(0, 0.6, 0.1) 
                minor_yticks_zoom = np.arange(0, 0.6, 0.05)
                
                ax.set_xticks(major_xticks_zoom)
                ax.set_xticks(minor_xticks_zoom, minor=True)
                ax.set_xticklabels(major_xticks_zoom)      #to choose to display only major tick labels
                ax.set_yticks(major_yticks_zoom)
                ax.set_yticks(minor_yticks_zoom, minor=True)                                           
            
                ax.set_xlim([0.5,2.5])
                ax.set_ylim([0,0.3])
            return fig
        
        fig = creation_fig_xmax()
        
        def onclick(event):
            global xmax_click, ymax_click
            xmax_click = event.xdata  # record x value of click
            ymax_click = event.ydata 
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
            ax.set_title(gofrfile.split('.gofr.dat')[0] + '  ' + pair)
            ax.plot(distance, gofr, 'o', xaxis_max, max_fit(xaxis_max), '-', markersize=2)  # plot the fitted max
            plt.xlabel('Distance ($\AA$)')
            plt.ylabel('g(r)')
            plt.grid()
            if pair == 'O2':
                major_xticks_zoom = np.arange(0, 3, 0.5) 
                minor_xticks_zoom = np.arange(0, 3, 0.1)
                major_yticks_zoom = np.arange(0, 0.6, 0.1) 
                minor_yticks_zoom = np.arange(0, 0.6, 0.05)
                
                ax.set_xticks(major_xticks_zoom)
                ax.set_xticks(minor_xticks_zoom, minor=True)
                ax.set_xticklabels(major_xticks_zoom)      #to choose to display only major tick labels
                ax.set_yticks(major_yticks_zoom)
                ax.set_yticks(minor_yticks_zoom, minor=True)                                           
            
                ax.set_xlim([0.5,2.5])
                ax.set_ylim([0,0.3])    
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

        y_intgofr2 = intgofr[min_index]     #the coordination number corresponding to the clicked position            


        #****** determine start and end cols for fitting the polynomial for the min region:
        start_imin = int(min_index - min_index * 0.08)  # 8% to the left (arbitrary)
        end_imin = int(min_index + min_index * 0.08)  # 8% to the right (arbitrary)
        if end_imin > (orig_size - 1):  # for rare cases where the min distance is out of bounds
            end_imin = orig_size - 1  # subtract one because cols start at 0
        # define the regions that need to be fitted for min area using the clicked min:
        fit_min_dist = distance[int(start_imin):int(end_imin + 1)]  # ranges do not include last element, so add one
        fit_min_gofr = gofr[int(start_imin):int(end_imin + 1)]
        intgofr = intgofr[int(start_imin):int(end_imin + 1)]
        
        min_fit = np.poly1d(np.polyfit(fit_min_dist, fit_min_gofr, order))
        intgofr_fit = np.poly1d(np.polyfit(fit_min_dist, intgofr, order))
            
        xaxis_min = np.linspace(distance[int(start_imin)], distance[int(end_imin)],1000)  # create new x values for the fitted min region
        yaxis_min = min_fit(xaxis_min)  # create new y values for fitted min region
        poly_min_i = np.argmin(yaxis_min)  # determine index of min y value
        yaxis_intgofr = intgofr_fit(xaxis_min)

        xmax_gofr = xaxis_max[int(poly_max_i)]  # x value of the max (x,y) of the max region of gofr
        ymax_gofr = yaxis_max[int(poly_max_i)]  # max y value of the max region of gofr
        xmin_gofr = xaxis_min[int(poly_min_i)]  # x value of the min (x,y) of the min region of gofr
        y_intgofr = yaxis_intgofr[int(poly_min_i)]  # y value of the integral of gofr that corresponds to the min x value of gofr

        
        #same for reverse pair if we are not analyzing O2 or other special additional column
        if pair != 'O2':
            y_intgofr_reverse2 = intgofr_reverse[min_index]     #the coordination number corresponding to the clicked position
            intgofr_reverse = intgofr_reverse[int(start_imin):int(end_imin + 1)]
            intgofr_fit_reverse = np.poly1d(np.polyfit(fit_min_dist, intgofr_reverse, order))
            yaxis_intgofr_reverse = intgofr_fit_reverse(xaxis_min)
            y_intgofr_reverse = yaxis_intgofr_reverse[int(poly_min_i)]  # y value of the integral of gofr that corresponds to the min x value of gofr
        else:
            y_intgofr_reverse = 0
            y_intgofr_reverse2 = 0
    
        #******* Third we save the fitted of clicked data depending on the choice made (good or bad)
        class Decision(object):
            def fit(self, event):
                global value_array, decision
                plt.close(fig)  # close and move to next plot
                print('fitted values:')
                print('xmax=', round(xmax_gofr, 2), ' ymax=', round(ymax_gofr, 2), ' xmin=', round(xmin_gofr, 2), ' coord=', round(y_intgofr, 2))
                value_array = (xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, y_intgofr_reverse)
                decision = 'f'
                return value_array

            def bad(self, event):
                global value_array, decision
                plt.close(fig)  # close and move to next plot
                print('bad values:')
                print('xmax=', 0, ' ymax=', 0, ' xmin=', 0, ' coord=', 0)
                value_array = (0, 0, 0, 0, 0)
                decision = 'b'
                return value_array
            
            def click(self, event):
                global value_array, decision
                plt.close(fig)  # close and move to next plot
                print('clicked values:')
                print('xmax=', round(xmax_click, 2), ' ymax=', round(ymax_click, 2), ' xmin=', round(xmin_click, 2), ' coord=', round(y_intgofr2, 2))
                value_array = (xmax_click, ymax_click, xmin_click, y_intgofr2, y_intgofr_reverse2)
                decision = 'c'
                return value_array
        
        fig = plt.figure()  # show figure with the newly fitted minimum curve
        ax = fig.add_subplot(111)
        plt.subplots_adjust(bottom=0.2)
        ax.set_title(gofrfile.split('.gofr.dat')[0] + '  ' + pair)
        ax.plot(distance, gofr, 'o', xaxis_max, max_fit(xaxis_max), '-', xaxis_min, min_fit(xaxis_min), 'r-', markersize=2)
        plt.xlabel('Distance ($\AA$)')
        plt.ylabel('g(r)')
        if pair == 'O2':
            major_xticks_zoom = np.arange(0, 3, 0.5) 
            minor_xticks_zoom = np.arange(0, 3, 0.1)
            major_yticks_zoom = np.arange(0, 0.6, 0.1) 
            minor_yticks_zoom = np.arange(0, 0.6, 0.05)
            
            ax.set_xticks(major_xticks_zoom)
            ax.set_xticks(minor_xticks_zoom, minor=True)
            ax.set_xticklabels(major_xticks_zoom)      #to choose to display only major tick labels
            ax.set_yticks(major_yticks_zoom)
            ax.set_yticks(minor_yticks_zoom, minor=True)                                           
        
            ax.set_xlim([0.5,2.5])
            ax.set_ylim([0,0.4])

        callback = Decision()
        axclick = plt.axes([0.12, 0.05, 0.1, 0.075])
        axbad = plt.axes([0.45, 0.05, 0.1, 0.075])
        axfit = plt.axes([0.8, 0.05, 0.1, 0.075])
        bfit = Button(axfit, 'Fit')
        bfit.on_clicked(callback.fit)
        bbad = Button(axbad, 'Bad')
        bbad.on_clicked(callback.bad)
        bclick = Button(axclick, 'Click')
        bclick.on_clicked(callback.click)
        plt.show()

        xmax_gofr = value_array[0]
        ymax_gofr = value_array[1]
        xmin_gofr = value_array[2]
        y_intgofr = value_array[3]
        if pair != 'O2':
            y_intgofr_reverse = value_array[4]

    else: # simply output 0's if the data is no good
        print('skip')
        xmax_gofr, ymax_gofr, xmin_gofr, y_intgofr, y_intgofr_reverse = 0,0,0,0,0
        skipdata = 1
        decision = 'b'
    try: 
        bond = average_bond(distance,gofr,xmin_gofr) #we compute the real bond length using the selected xmin value
    except ZeroDivisionError: #happens if xmin = 0
        bond = 0
    

    #*********** Last step, we save in the dictionnary the resultsfor the current pair of atoms
    #**** Update method
    if method[0] == 'u':
        print('data before change for', pair, data[gofrfile.split('/')[-1]][allpairs[pair]*5:allpairs[pair]*5+5])
        data[gofrfile.split('/')[-1]][allpairs[pair]*5:allpairs[pair]*5+5] = [str(xmax_gofr),str(ymax_gofr),str(xmin_gofr),str(y_intgofr), str(bond)]
        print('data after change for',pair, data[gofrfile.split('/')[-1]][allpairs[pair]*5:allpairs[pair]*5+5])
        if pair != 'O2':
            print('change also reverse pair')
            data[gofrfile.split('/')[-1]][allpairs[reversepair]*5:allpairs[reversepair]*5+5] = [str(xmax_gofr),str(ymax_gofr),str(xmin_gofr),str(y_intgofr_reverse), str(bond)]
    #**** Writing method
    else:
        #we compute the real bond length using the clicked xmin value
        try: 
            bond_click = average_bond(distance,gofr,xmin_click) 
        except ZeroDivisionError: #happens if xmin = 0
            bond_click = 0  
        #we write the click value
        #if skipdata == 1 or decision == 'b':
        if decision == 'b':
            data['clicked'][allpairs[pair]*5:allpairs[pair]*5+5] = ['0','0','0','0','0']
            if pair != 'O2':
                print('change also reverse pair')
                data['clicked'][allpairs[reversepair]*5:allpairs[reversepair]*5+5] = ['0','0','0','0','0']
        else:
            data['clicked'][allpairs[pair]*5:allpairs[pair]*5+5] = [str(xmax_click),str(ymax_click),str(xmin_click),str(y_intgofr2), str(bond_click)]
            if pair != 'O2':
                print('change also reverse pair')
                data['clicked'][allpairs[reversepair]*5:allpairs[reversepair]*5+5] = [str(xmax_click),str(ymax_click),str(xmin_click),str(y_intgofr_reverse2), str(bond_click)]
        #if the fitted value is good, we write it
        if decision == 'f':
            data['fitted'][allpairs[pair]*5:allpairs[pair]*5+5] = [str(xmax_gofr),str(ymax_gofr),str(xmin_gofr),str(y_intgofr), str(bond)]
            if pair != 'O2':
                print('change also reverse pair')
                data['fitted'][allpairs[reversepair]*5:allpairs[reversepair]*5+5] = [str(xmax_gofr),str(ymax_gofr),str(xmin_gofr),str(y_intgofr_reverse), str(bond)]
        #if not we write 0
        else:
            data['fitted'][allpairs[pair]*5:allpairs[pair]*5+5] = ['0','0','0','0','0']
            if pair != 'O2':
                print('change also reverse pair')
                data['fitted'][allpairs[reversepair]*5:allpairs[reversepair]*5+5]  = ['0','0','0','0','0']    
    #**** Bondfile
    if bonds != {}:
        bonds[pair] = str(xmin_gofr)
        if pair != 'O2':
            bonds[reversepair] = str(xmin_gofr)
    
    cutdistance = 0  # re initialization of variables for the next loop
    fit_min_dist = 0
    fit_max_dist = 0
    fit_min_gofr = 0
    fit_max_gofr = 0
    gofr = 0
    intgofr = 0
    intgofr_reverse = 0
    if bonds != {}:
        return (data, bonds)
    else:
        return(data)


decision = []
xmin_click = []
xmax_click = []
def analyze_gofrs_interactive(data,gofrfile,allpairs, atoms, bonds, method):
    """ extraction of the min, max and coordination number using fits and interactive plot """
    distance = np.loadtxt(gofrfile, usecols=(0,), skiprows = 1, unpack=True)
    for pair in allpairs:
        for atom_pair in atoms:
            if (pair == atom_pair):
                #interactive fit only for the selected allpairs of atoms
                if bonds != {}:
                    data, bonds  = interactive_fit(data,gofrfile,allpairs, distance, bonds, pair, method)
                else:
                    data  = interactive_fit(data,gofrfile,allpairs, distance, bonds, pair, method)
    if bonds != {}:
        return(data,bonds)
    else:
        return(data)



def main(argv):
    """     ********* Main program *********     """
    bondfile = ''
    datafile = ''
    gofrfile = ''
    bonds = {}
    add_X_O2 = 0
    try:
        options,arg = getopt.getopt(argv,"hf:g:a:b:",["fgofrfile","gofrsfile","atoms","bondfile"])
    except getopt.GetoptError:
        print("analyze_1gofr_update.py -f <gofr_filename.dat> -g <subfolder_gofrs.txt> -a <couples of atoms>(ex: 'Ca-O,Ca-Ca') -b <bond filename to update> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('*******')
            print('analyze_1gofr_update.py program to extract all relevant data from 1 gofr.dat file created by the script gofrs.py and write the corrected value (fitted or clicked) into the full gofrs.txt file and into the bonds.inp if indicated ')
            print("analyze_1gofr_update.py -f <gofr_filename.dat> -g <subfolder_gofrs.txt> -a <couples of atoms>(ex: 'Ca-O,Ca-Ca') -b <bond filename to update (option)>")
            print('WARNING! you have to click FIRST on the maximum, and THEN on the minimum')
            print('WARNING!!!! If you want precise clicked value, you have to click on the correct position of maximum and minimum')
            print(' ')
            print("If subfolder_gofrs.txt indicated, then this code replace the old data of the gofrs.txt file created by analyze_gofr_semi_automatique.py by the fitted value if you choose 'fit', by the clicked value if you choose 'click' or by 0 if you choose 'bad'")
            print('If no file is provided, then this code creates a gofr_....txt file with the clicked value and the fitted value if fit was selected.')
            print('*******')
            print(' ')
            sys.exit()
        elif opt in ("-f", "--fgofrfile"):
            gofrfile = str(arg)
        elif opt in ("-a", "--atoms"):
            atoms = arg.split(',')                      #list of atom couples we want to analyze here)
        elif opt in ("-b","--bondfile"):
            bondfile = str(arg)
        elif opt in ("-g","--gofrsfile"):
            datafile = str(arg)
    #******* -1st step: check for the presence of files
    if os.path.isfile(gofrfile):
        if bondfile != '':
            if os.path.isfile(bondfile):
                pass
            else:
                print("The bondfile '", bondfile, "' does not exist, continue without modifying bondfile")
                bondfile = ''
        if os.path.isfile(datafile):
            print("We take the file '", datafile, "' as datafile which will be updated with the new selected values")
            method='update'
        else:
            print("the datafile '", datafile, "' does not exist, this code will create a new .txt file with the values fitted and clicked")
            method='write'
    else:
        print("The gofrfile '", gofrfile, "' does not exist")
        sys.exit()
    #******* 1st step: extraction of pairs of atoms from gofr.dat file 
    allpairs = {}  #to store the index of each pair of atoms (in order to find the correct column number in gofr and data files)
    allpairs_ordered = []  #to store all the pair of atoms in the same order as in gofr file (--> always keep the same order)
    with open(gofrfile, 'r') as f:
        line = f.readline()
    line = line.strip('dist')
    line = re.sub('(Int\([A-Za-z-]*\))', ' ',line).split()        #I substitute all Int(...)  by a blankspace and split using spaces
    for ii in range(len(line)):
        allpairs[line[ii]] = ii
        allpairs_ordered.append(line[ii]) 
    if 'O2' in atoms:
        allpairs['O2'] = allpairs['O-O'] + 1
        allpairs_ordered.append('O2')
    print("all pairs of atoms in header", allpairs)
    
      
    #******* 2nd step: extraction of data from the full gofrs.txt file or initialization of data depending on the method used (update or writing)
    #**** Update method
    if method[0] == 'u':
        #First we check if the header contains additionnal "pairs", like O2, in order to add them if needed
        with open(datafile,'r') as d:
            while True:
                line = d.readline()
                entry = line.split('\n')[0].split('\t')
                if entry[0] == 'pair':
                    break
        #if O2 is selected in the atoms list but there is not yet header for it, then we add it to the header
        if 'O2' not in line and 'O2' in atoms:
            add_X_O2 = 1
            add_header_pairs = '\tO2\tO2\tO2\tO2\tO2\n'
            add_header_file = '\txmax\tymax\txmin\tcoord\tbond\n'
        else:
            add_header_pairs = '\n'
            add_header_file = '\n'    
        #Second we extract the header and complete it with new pairs if needed        
        data = {}
        allfiles_ordered = []
        with open(datafile,'r') as d:
            #save the header
            header = ''
            while True:
                line = d.readline()
                entry = line.split('\n')[0].split('\t')
                if entry[0] == 'pair':
                    header = header + line.split('\n')[0] + add_header_pairs
                elif entry[0] == 'file':
                    header = header + line.split('\n')[0] + add_header_file
                    break
                else:
                    header = header + line
            while True:
                line = d.readline()            
                if not line : break
                else:
                    entry = line.split('\n')[0].split('\t')
                    #data is the dictionnary with filenames as keys
                    data[entry[0].split('/')[-1]] = entry[1:]
                    allfiles_ordered.append(entry[0].split('/')[-1])
                    #add X,X,X,X,X for the O2 columns if we create them
                    if add_X_O2 == 1:
                        data[entry[0].split('/')[-1]].extend(['0','0','0','0','0'])
    #**** Writing method
    else:
        #save the header
        header = 'method'
        for pair in allpairs_ordered:
            for i in range(0,5):
                header = header + '\t' + pair
        header = header + '\nmethod'
        for pair in allpairs_ordered:
            header = header + '\txmax\tymax\txmin\tcoord\tbond'
        header = header + '\n'
        #initialize the data dictionnary with method as keys
        data = {}
        data['clicked'] = []
        data['fitted'] = []
        for pair in allpairs_ordered:
            data['clicked'].extend(['0','0','0','0','0'])
            data['fitted'].extend(['0','0','0','0','0'])
    #******* 3rd step: extraction of bond file informations if bondfile indicated
    if bondfile != '':
        with open(bondfile,'r') as b:
            while True:
                line = b.readline()
                if not line:break
                else:
                    entry = line.split('\n')[0].split('\t')
                    bonds[entry[0]+'-'+entry[1]] = entry[2]
    #******* 4th step: Analysis
    if bondfile != '':
        data, bonds = analyze_gofrs_interactive(data,gofrfile,allpairs, atoms, bonds, method)
    else:
        data  = analyze_gofrs_interactive(data,gofrfile,allpairs, atoms, bonds, method)
    #******* 5th step: overwrite the files
    #**** Update method
    if method[0] == 'u':
        with open(datafile,'w') as f:
            #print the header previously saved
            f.write(header)
            #then we write all the results        
            for file in allfiles_ordered:
                f.write(file+'\t')
                f.write("\t".join(x for x in data[file]))
                f.write('\n')
        print('The file ',datafile,' is updated')
    #**** Writing method
    else:
        newfilename = 'gofr_'+gofrfile.split('/')[-1].split('.dat')[0]+'.txt'
        with open(newfilename,'w') as f:
            #print the header previously saved
            f.write(header)
            #then we write all the results    
            f.write('fitted\t'+ "\t".join(x for x in data['fitted']) + '\n')
            f.write('clicked\t'+ "\t".join(x for x in data['clicked']) + '\n')
        print('The file ',newfilename,' is created')
    #**** Bondfile
    if bondfile != '':
        with open(bondfile, 'w') as b:
            for pair in allpairs_ordered:
                #we remove O2 from the loop because we only need elements pairs for bond files
                if pair == 'O2': continue
                #but if we need the xmin for O2, then when we compute it, it replaces the one for O-O
                elif pair == 'O-O' and 'O2' in atoms:
                    print("replace O-O xmin by O2 xmin in bond file")
                    b.write('O\tO\t'+bonds['O2']+'\n')
                else:
                    b.write(pair.split('-')[0]+'\t'+pair.split('-')[1]+'\t'+bonds[pair]+'\n')
            b.close()
            print("bond file ",bondfile, ' updated') 
            
        
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])



