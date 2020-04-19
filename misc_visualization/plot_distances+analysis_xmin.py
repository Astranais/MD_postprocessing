#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot one of the distance in .gofrs.txt       ****
                  and update/create the .bonds.inp files with data smoothing
                          one column figure version
                    
"""



#     ********* Importation of the packages and modules used here *********
import sys
import getopt
import glob
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import AutoLocator,AutoMinorLocator
import numpy as np
import crystallography as cr
from scipy.optimize import curve_fit



def poly3(x,a,b,c,d):
    """ 3rd order polynomial """
    xarr = np.asarray(x)
    yarr = a * xarr**3 + b * xarr**2 + c * xarr + d
    y = yarr.tolist()
    return y

def linear(x,a,b):
    """ linear equation """
    xarr = np.asarray(x)
    yarr = a * xarr + b 
    y = yarr.tolist()
    return y

def chi2red(popt,rho,data,function):
    """ function chi2 which has to be minimized"""
    y = function(rho,*popt)
    # compute chi-square using stdev data = 1
    chi2 = 0.0
    for ii in range(len(rho)):
        chi2 = chi2 + ( (data[ii] - y[ii]) )**2
    chi2red = chi2 / (len(rho) - len(popt))
    return chi2red



def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./XX/CaAl2Si2O8_T3_nvt_a12.5.outcar.gofr.dat
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.split('/')[-1]
    temperature = filename.split('_')[1]
    acell = filename.split('.outcar.gofr.dat')[0].split('_')[3].strip('a')
    return temperature, acell




def smoothing_xmin(data,selected_pairs):
    """Smooth xmin data using fitting"""
    #**** Fit of each curve for each atom pairs and T
    chi2_reduced = {}
    params = {}
    new_xmin = {}
    for temp in sorted(data):
        #print(temp)
        chi2_reduced[temp] = {}
        params[temp] = {}             
        for i in range(len(selected_pairs)):
            atom_pair = selected_pairs[i][0]
            #remove nan from arrays
            idx = np.isfinite(data[temp]['rho']) & np.isfinite(data[temp]['dist'][atom_pair])
            newrho = []
            newdata =[]
            for i in range(len(idx)):
                if idx[i]:
                    newrho.append(data[temp]['rho'][i])
                    newdata.append(data[temp]['dist'][atom_pair][i])
            #fit only if the size of the data is enough
            #print(atom_pair,"len of data is", len(newdata))
            if len(newdata) > 1:
                if atom_pair == 'O2' or len(newdata) < 5:
                    function=linear
                else:
                    function=poly3
                params[temp][atom_pair], pcov = curve_fit(function, newrho, newdata, p0=None, sigma=None, absolute_sigma=False )
                chi2_reduced[temp][atom_pair]=chi2red(params[temp][atom_pair],data[temp]['rho'],data[temp]['dist'][atom_pair],function)
            else:
                print("Not enough data to fit any curve for atom pair ", atom_pair, " and temperature ", temp)
    #**** Compute the new xmin for each file and atom pair
    maxrho = 0
    minrho = 9999
    for temp in sorted(data):
        #print(temp)
        for i in range(len(data[temp]['file'])):
            new_xmin[data[temp]['file'][i]] = {}
            for j in range(len(selected_pairs)):
                atom_pair = selected_pairs[j][0]
                rho=data[temp]['rho'][i]
                if rho > maxrho: maxrho =rho
                if rho < minrho: minrho =rho
                if np.isnan(data[temp]['dist'][atom_pair][i]):
                    new_xmin[data[temp]['file'][i]][atom_pair]= np.nan
                else:
                    try:
                        if len(params[temp][atom_pair]) == 2:
                            function=linear
                        else: 
                            function=poly3
                    except KeyError:  
                        params[temp][atom_pair] = [0,0]
                        function=linear
                    new_xmin[data[temp]['file'][i]][atom_pair]=function(rho,*params[temp][atom_pair])  
                #print("for", data[temp]['file'][i], "and atom pair", atom_pair, "new xmin is", round(new_xmin[data[temp]['file'][i]][atom_pair],3), "compared to old", round(data[temp]['dist'][atom_pair][i],3))
    return new_xmin, minrho, maxrho, params, chi2_reduced




def update_bonds(bondfile,selected_pairs,atoms,new_xmin,unique_pairs,file):
    print("bond file ",bondfile, ' updated')
    # Extraction of bond file informations 
    bonds = {}
    with open(bondfile,'r') as b:
        while True:
            line = b.readline()
            if not line:break
            else:
                entry = line.split('\n')[0].split('\t')
                bonds[entry[0]+'-'+entry[1]] = entry[2]
    # Update the bond dictionnary for this bondfile                
    for i in range(len(selected_pairs)):
        pair = selected_pairs[i][0]
        if np.isnan(new_xmin[file][pair]):pass
        else:
            bonds[pair] = str(new_xmin[file][pair])
        if pair != 'O2':
            reversepair = pair.split('-')[1]+'-'+pair.split('-')[0]
            bonds[reversepair] = str(new_xmin[file][pair])
    # Update the bondfile
    with open(bondfile, 'w') as b:
        for pair in unique_pairs:
            #we remove O2 from the loop because we only need elements pairs for bond files
            if pair == 'O2': continue
            #but if we need the xmin for O2, then when we compute it, it replaces the one for O-O
            elif pair == 'O-O' and 'O2' in atoms:
                try:
                    b.write('O\tO\t'+bonds['O2']+'\n')
                    print("replace O-O xmin by O2 xmin in bond file")
                except KeyError:
                    print("no O2 bond, nothing to change in O-O xmin")
                    b.write('O\tO\t'+bonds['O-O']+'\n')
            else:
                b.write(pair.split('-')[0]+'\t'+pair.split('-')[1]+'\t'+bonds[pair]+'\n')
        b.close()
    
def create_bonds(bondfile,selected_pairs,atoms,new_xmin,unique_pairs,file):
    print("bond file ",bondfile, ' created')
    # Create the bond dictionnary for this bondfile                
    bonds = {}
    selected_pairs_name= [] #we create this list instead of using atoms list directly in case we indicated a wrong atom pair in atoms list
    for i in range(len(selected_pairs)):
        selected_pairs_name.append(selected_pairs[i][0])
    for pair in unique_pairs:
        if pair in selected_pairs_name:
            #we do analyze this atom_pair so we indicate the computed value
            if np.isnan(new_xmin[file][pair]):
                bonds[pair] = str(0)
            else:
                bonds[pair] = str(new_xmin[file][pair])
            if pair != 'O2':
                reversepair = pair.split('-')[1]+'-'+pair.split('-')[0]
                bonds[reversepair] = str(new_xmin[file][pair])
        else:
            #we do not analyze this atom_pair so we fill with 0 except if it is a reverse pair
            if pair != 'O2':
                reversepair = pair.split('-')[1]+'-'+pair.split('-')[0]
            if reversepair in selected_pairs_name: pass
            else:
                bonds[pair] = str(0)
    # Write the bondfile
    with open(bondfile, 'w') as b:
        for pair in unique_pairs:
            #we remove O2 from the loop because we only need elements pairs for bond files
            if pair == 'O2': continue
            #but if we need the xmin for O2, then when we compute it, it replaces the one for O-O
            elif pair == 'O-O' and 'O2' in atoms:
                try:
                    b.write('O\tO\t'+bonds['O2']+'\n')
                    print("replace O-O xmin by O2 xmin in bond file")
                except KeyError:
                    print("no O2 bond, nothing to change in O-O xmin")
                    b.write('O\tO\t'+bonds['O-O']+'\n')
            else:
                b.write(pair.split('-')[0]+'\t'+pair.split('-')[1]+'\t'+bonds[pair]+'\n')
        b.close()   


def main(argv):
    """     ********* Main program *********     """
    #other dictionnaries and parameters for the figure
    markers = ['o','^','*','P','X','<','>','v','8','s','p','h','+','D'] #************************************Be sure you have enough markers for all the pairs you want to plot
    colors_T = {'T2':'#800080','T3':'#297fff','T4':'#00ff00','T4.5':'#bae200','T5':'#ffcd01','T5.5':'#ff6e00','T6':'#ff0101','T6.5':'#ff00a2','T7':'#ff01de','T7.5':'#ffa6f4','T10':'#ffe86e','T15':'#ffbf90','T20':'#ff7788'}
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),"size_markers" : 8,"size_lines" : 1,"shift_labelpad" : 20}
    label = {'xmax':r"1st $g(r)$ peak location ($\AA$)",'xmin':r"Coordination sphere radius ($\AA$)",'bond':r"Bond length ($\AA$)"}
    #initialization of variables
    bondanalysis = 0
    elements = ''
    number = ''
    atoms = []
    legend_labels = {} 
    data = {} #big dictionnary with inside dist dictionnary and rho, all coresponding to different T
    distance_columns = {'xmax':1,'xmin':3,'bond':5}
    #other parameters
    Na=6.022*10**23
    try:
        options,arg = getopt.getopt(argv,"hg:a:d:b:",["gofrsfilename","atom","distance",'bondanalysis'])
    except getopt.GetoptError:
        print("plot_distances+analysis_xmin.py -g <_gofrs.txt> -a <pairs of atoms>(ex: 'Ca-O,Ca-Ca,O2')  -d <distance type to print ('xmax' or 'xmin' or 'bond')> -b <=1 if analysis and update of bondfiles, default =0>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_distances+analysis_xmin.py program to plot xmin,xmax or bond length as a function of density or acell for each T and all the selected pairs of atoms')
            print("plot_distances+analysis_xmin.py -g <_gofrs.txt> -a <pairs of atoms>(ex: 'Ca-O,Ca-Ca,O2')  -d <distance type to print ('xmax' or 'xmin' or 'bond')> -b <=1 if analysis and update of bondfiles, default =0>")
            print("plot_distances+analysis_xmin.py requires _gofrs.txt file (use analyze_gofrs_semi_automatic.py)")
            print('')
            print('WARNING: this script use the filenames inside the _gofrs.txt file to extract the temperature and cell size for the plot. If you do have both information in your filenames at the same place, then update the function split_name to extract them correctly.')
            sys.exit()
        if opt in ('-g','--gofrsfilename'):
            filename = str(arg)
        elif opt in ('-a','--atoms'):
            atoms = arg.split(',')                      #list of atom pairs we want to analyze here
        elif opt in ('-d','--distance'):
            distance_type = str(arg)
        elif opt in ('-b','--bondanalysis'):
            bondanalysis = int(arg)
    #******* 1st step: read the header and extract all relevant informations
    #creation of elements and number lists and initialization of T
    skip_head = 0
    with open(filename,'r') as f:
        while True:
            line = f.readline()
            if not line: break
            else:
                entry = line.split('\n')[0].split('\t')
                skip_head +=1
                if entry[0] == 'elements':
                    elements = list(filter(None, entry[1:]))
                    print(elements)
                if entry[0] == 'number':
                    number = list(filter(None, entry[1:]))
                    print(number)
                if entry[0] == 'pair':
                    allpairs_ordered = entry[1:]
                if entry[0] == 'file':
                    line = f.readline()
                    entry=line.split()
                    temperature0, acell0 = split_name(entry[0]) 
                    break
    if elements != '' and number != '':
        #calculation of M*N nedded for the calculation of densities
        MN = 0
        for i in range(len(elements)):
            MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
    else:
        MN =0
    #creation of the list containing all the pairs only once and in the same order
    unique_pairs = [allpairs_ordered[0]]
    for i in range(1,len(allpairs_ordered)):
        if allpairs_ordered[i] != allpairs_ordered[i-1]:
            unique_pairs.append(allpairs_ordered[i])
    print('All the pairs availables in the file are:',unique_pairs)
    print('The pairs selected for the plot are:', atoms)
    print('The index for each pair is:')
    #initialization of the dictionnaries containing the data
    #we create sub dictionnaries of the data dictionnary
    data[temperature0] = {'file':[],'rho':[],'dist':{}}
    selected_pairs = [] #pairs to analyze along with their column number
    for i in range(0,len(unique_pairs)):
        for atom_pair in atoms:
            if (unique_pairs[i] == atom_pair): 
                data[temperature0]['dist'][atom_pair] = []
                selected_pairs.append((atom_pair,i*5+distance_columns[distance_type]))
                print(atom_pair, i*5+distance_columns[distance_type])
    #******* 2nd step: extraction of data from the fullgofrs.txt file 
    with open(filename,'r') as f:
        [f.readline() for i in range(skip_head)]
        while True:
            line = f.readline()
            if not line: break
            else:
                entry=line.split('\t')
                temperature, acell = split_name(entry[0]) 
                if temperature != temperature0:     #if we change T:
                    #we create sub dictionnaries of the data dictionnary  and we re-initialize the arrays
                    temperature0 = temperature
                    data[temperature0] = {'file':[],'rho':[],'dist':{}}
                    for i in range(len(selected_pairs)):
                        data[temperature0]['dist'][selected_pairs[i][0]] = []
                if MN ==0:
                    data[temperature0]['rho'].append(float(acell))     #rho is replaced by acell
                else:
                    data[temperature0]['rho'].append(MN/(Na*float(acell)**3*10**(-24)))     #calculation density
                data[temperature0]['file'].append(entry[0].split('/')[-1].split('.gofr.dat')[0])              #extract filename
                #we replace the 0 (no data) by NaN
                for i in range(len(selected_pairs)):
                    if float(entry[selected_pairs[i][1]]) == 0.0:
                        data[temperature0]['dist'][selected_pairs[i][0]].append(float('nan'))
                    else:
                        data[temperature0]['dist'][selected_pairs[i][0]].append(float(entry[selected_pairs[i][1]]))
    #******* 3rd step: bond analysis 
    #Smooth xmin data using fitting and get minrho, maxrho even if we don't want to update/create .bonds.inp
    new_xmin, minrho, maxrho, params, chi2_reduced = smoothing_xmin(data,selected_pairs)
    if bondanalysis == 1:    
        #**** Extract all gofrfiles available
        allgofrfiles = []
        for dirpath, dirnames, filenames in os.walk(os.curdir):
            allgofrfiles.extend(sorted(glob.glob(dirpath+'/*.gofr.dat')))
        print("**** gofrfiles in the current directory are:", allgofrfiles)
        gofrfileslocation = {}
        for gofrfile in allgofrfiles:
            gofrfileslocation[gofrfile.split('/')[-1].split('.gofr.dat')[0]] = gofrfile
        #**** For each line (gofrfile) in the .gofrs.txt file
        for temp in sorted(data):
            for file in  data[temp]['file']:
                #Check if the bondfile already exist: if yes we update it, if not we create it
                bondfile = gofrfileslocation[file].split('.gofr.dat')[0]+'.bonds.inp'
                if (os.path.isfile(bondfile)):
                    update_bonds(bondfile,selected_pairs,atoms,new_xmin,unique_pairs,file)
                else:
                    create_bonds(bondfile,selected_pairs,atoms,new_xmin,unique_pairs,file)           
    else:
        print("Bond analysis not done. If you want it, use option -b 1")
    #******* 4th step: plot of the data and write the deviation from mean in file
    plt.close(1)
    h = 2 * len(atoms) #height of the figure depends on the number of pairs we display
    fig = plt.figure(1,figsize = (4,h))
    plt.subplots_adjust(top = 0.97, bottom = 0.07, right = 0.89, left = 0.07, hspace = 0, wspace = 0)
    #Creation of ticks
    if MN != 0:
        major_ticks = np.arange(0.5, 7, 0.5) 
        minor_ticks = np.arange(0.5, 7, 0.1)
    else:
        major_ticks = AutoLocator()
        minor_ticks = AutoMinorLocator()                                                                 
    #plot & write percentage variation of distance
    ylim_allT = {}
    string = ''
    for atom_pair in atoms:
        string = string+'_'+atom_pair
    newfilename = 'distance_variation_'+distance_type+'_'+filename[:-4]+string+'.txt'
    f=open(newfilename,'w')
    for i in range(len(selected_pairs)):
        atom_pair = selected_pairs[i][0]
        ylim_allT[atom_pair] = [9999,0]
        #print('***************** for atom pair:',atom_pair)
        f.write(atom_pair+'\t'+"\t".join(temp for temp in sorted(data) )+ "\n")
        all_means='mean'
        all_range='range'
        all_percent='%var_from_mean'
        all_percent_rho='%var_from_mean_per_rho'
        all_params_a='fitted_parameter_a'
        all_params_b='fitted_parameter_b'
        all_params_c='fitted_parameter_c'
        all_params_d='fitted_parameter_d'
        all_chi2='chi2_reduced'
        plt.subplot(len(atoms), 1, atoms.index(atom_pair)+1)
        ax = plt.gca()
        #plot and prepare the strings to write in the new log file
        for temp in sorted(data): #loop over temperatures = key in data dictionnary
            #classic analysis (mean, range, %variation)
            #print('******* for temperature:',temp)
            mean_dist=np.nanmean(data[temp]['dist'][atom_pair])
            all_means = all_means+'\t'+str(np.round(mean_dist,2))
            #print('mean dist for',atom_pair,'=',np.round(mean_dist,2))
            range_var=(np.nanmax(data[temp]['dist'][atom_pair])- np.nanmin(data[temp]['dist'][atom_pair]) )
            all_range = all_range+'\t'+str(np.round(range_var,2))
            #print('dist range for', atom_pair, '=', np.round(range_var,2))
            percent_var_mean = (range_var /  mean_dist) * 100
            all_percent=all_percent+'\t'+str(np.round(percent_var_mean,1))
            #print('percent var mean for', atom_pair, '=', np.round(percent_var_mean,1))
            percent_var_mean_per_rho = percent_var_mean / ( np.nanmax(data[temperature0]['rho']) - np.nanmin(data[temperature0]['rho']) )
            all_percent_rho=all_percent_rho+'\t'+str(np.round(percent_var_mean_per_rho,1))
            #print('percent var mean  per rho for', atom_pair, '=', np.round(percent_var_mean_per_rho,1))
            #results from fit
            if bondanalysis == 1:
                newrhoX = np.arange(minrho, maxrho+0.1, 0.1)
                try:
                    all_params_a= all_params_a+'\t'+str(params[temp][atom_pair][0])
                    all_params_b= all_params_b+'\t'+str(params[temp][atom_pair][1])
                    if len(params[temp][atom_pair])==4:
                        all_params_c= all_params_c+'\t'+str(params[temp][atom_pair][2])
                        all_params_d= all_params_d+'\t'+str(params[temp][atom_pair][3])
                        function = poly3
                    else:
                        all_params_c= all_params_c+'\t-'
                        all_params_d= all_params_d+'\t-'
                        function = linear
                    all_chi2 = all_chi2+'\t'+str(chi2_reduced[temp][atom_pair])
                    #print('chi2 for', atom_pair, '=', chi2_reduced[temp][atom_pair])
                    ax.plot(newrhoX,function(newrhoX,*params[temp][atom_pair]), '-',color=colors_T[temp], linewidth = plot_parameters["size_lines"] )
                except KeyError:
                    all_chi2 = all_chi2+'\t-'
                    all_params_a= all_params_a+'\t-'
                    all_params_b= all_params_b+'\t-'
                    all_params_c= all_params_c+'\t-'
                    all_params_d= all_params_d+'\t-'
            #plot
            ax.plot(data[temp]['rho'],data[temp]['dist'][atom_pair], '--',  color=colors_T[temp], linewidth = plot_parameters["size_lines"], marker = markers[atoms.index(atom_pair)], markersize = plot_parameters["size_markers"])
#            ax.text(1.025,0.5, atom_pair ,transform=ax.transAxes, fontsize=plot_parameters["size_fonts"], fontweight='bold')
            ax.text(0.85,0.9, atom_pair, transform=ax.transAxes, fontsize=plot_parameters["size_fonts"], fontweight='bold')
            #define bounds for axis
            if ylim_allT[atom_pair][1] < np.nanmax(data[temp]['dist'][atom_pair]):
                ylim_allT[atom_pair][1] = np.nanmax(data[temp]['dist'][atom_pair])
            if ylim_allT[atom_pair][0] > np.nanmin(data[temp]['dist'][atom_pair]):
                ylim_allT[atom_pair][0] = np.nanmin(data[temp]['dist'][atom_pair])
        f.write(all_means+'\n')
        f.write(all_range+'\n')
        f.write(all_percent+'\n')
        f.write(all_percent_rho+'\n')
        f.write(all_chi2+'\n')
        f.write(all_params_a+'\n')
        f.write(all_params_b+'\n')
        f.write(all_params_c+'\n')
        f.write(all_params_d+'\n')
        f.write('\n')                           
        #Adjustment of ticks and make the graph prettier
        if MN != 0:
            ax.set_xticks(major_ticks)
            ax.set_xticks(minor_ticks, minor=True)
        else:
            ax.xaxis.set_major_locator(major_ticks)
            ax.xaxis.set_minor_locator(minor_ticks)    
        ax.xaxis.set_ticks_position('both')
        ax.yaxis.set_ticks_position('both')
        majorLocator = AutoLocator()
        minorLocator = AutoMinorLocator()
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_minor_locator(minorLocator)                
        ax.set_xlim(minrho,maxrho)
        ax.set_ylim(ylim_allT[atom_pair][0],ylim_allT[atom_pair][1])
        if  atoms.index(atom_pair) != len(atoms)-1:
            plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(which = 'both', labelsize = plot_parameters["size_font_ticks"], width = plot_parameters["size_lines"]/2)   
    f.close()
    # Fine-tune figure
    ax0 = fig.add_subplot(111, frameon=False) 
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright=False, right=False)
    if MN != 0:
        ax0.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], labelpad = plot_parameters["shift_labelpad"])
    else:
        ax0.set_xlabel(r'Cell size ($\AA$)', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], labelpad = plot_parameters["shift_labelpad"])    
    ax0.set_ylabel(label[distance_type],fontsize=plot_parameters["size_fonts"],fontweight='bold', labelpad = plot_parameters["shift_labelpad"]+plot_parameters["shift_labelpad"]/1.3)
    #Legend 
    for temp in data: 
        legend_labels[str(int(float(temp.strip('T'))*1000))] =  mpatches.Patch(color=colors_T[temp])
    s = [(k, legend_labels[k]) for k in sorted(legend_labels.keys(),reverse = False)]      
    #legend = ax0.legend([v for k,v in s],[k for k,v in s], title = ' $\\bf{Temperature}$ \n           (K)',loc='upper left',bbox_to_anchor=(1.2, 1), fontsize = plot_parameters["size_fonts"], borderaxespad=0.,ncol=1)
    #plt.setp(legend.get_title(),fontsize= plot_parameters["size_fonts"])
    #plt.title(filename.split('.txt')[0], fontsize = plot_parameters["size_fonts"], fontweight='bold')
    #save the figure
    string = ''
    for atom_pair in atoms:
        string = string+'_'+atom_pair
    figurename = 'distance_'+distance_type+'_'+filename[:-4]+string+'.png'
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 150)
    print(figurename, '   created')
    #plt.show()
   
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])

