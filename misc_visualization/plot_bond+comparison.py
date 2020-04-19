#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the different bond lenght measurements       ****
                    plot the enveloppe of the bond lenght variations
                    One subplot per pairs of atoms
                    comparison with data from litterature

            
                    
"""



#     ********* Importation of the packages and modules used here *********
import sys
import getopt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend import Legend
from matplotlib.ticker import AutoLocator,AutoMinorLocator
import numpy as np
import crystallography as cr








def split_name(filename):
    """ Function to extract compound name from filename """
    # *********** My filename are in the format ./XX/CaAl2Si2O8_T3_nvt_a12.5.outcar.gofr.dat
    # ******* so I can split their name with _ and take the compound and T from their name
    filename = filename.split('/')[-1]
    temperature = str(int(float(filename.split('_')[1].strip('T'))*1000))
    acell = filename.split('.outcar.gofr.dat')[0].split('_')[3].strip('a')
    return temperature, acell


def main(argv):
    """     ********* Main program *********     """
    #parameters for the figures depending on the output format (presentation or article)
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (8,4),"size_markers" : 8,"size_lines" : 1,"shift_labelpad" : 20}
    colors_T = {'0':'k','2000':'#800080','3000':'#297fff','4000':'#00ff00','4500':'#bae200','5000':'#ffcd01','5500':'#ff6e00','6000':'#ff0101','6500':'#ff00a2','7000':'#ff01de'}
    plot_area = 1 #0 to only plot lines of each colors for our data, 1 to plot instead area (enveloppe) of our data
    #initialization of parameters
    filename4 = ''
    filename3 = ''
    #other parameters
    Na=6.022*10**23
    try:
        options,arg = getopt.getopt(argv,"hf:g:j:i:",["filename2","gfilename","jfilename","ifilename"])
    except getopt.GetoptError:
        print("plot_coordinence_separated.py -f <bonds_compound.txt> -g <gofr_comparison1.txt> -j <gofr_comparison2.txt> -i <gofr_comparison3.txt> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('plot_coordinence_separated.py program to plot coordinence as a function of density for each T and all cations with O')
            print("plot_coordinence_separated.py-f <bonds_compound.txt> -g <gofr_comparison1.txt> -j <gofr_comparison2.txt> -i <gofr_comparison3.txt> ")
            print("plot_coordinence_separated.py requires bonds.txt file containing every bond distance of one compound for every acell (use bond_analysis.py script)")
            print('')
            sys.exit()
        if opt in ('-f','--filename'):
            filename = str(arg) #plot bond length --> column 5 and xmax --> column 1
        elif opt in ('-g','--gfilename'):
            filename2 = str(arg) # plot from deKoker 2010 
        elif opt in ('-j','--jfilename'):
            filename3 = str(arg) 
        elif opt in ('-i','--ifilename'):
            filename4 = str(arg)
    #************ creation of the figure to plot (with correct number of lines)
    if filename3 == '':    
        if filename4 == '':
            files = [filename,filename2]
            string = filename.split('_')[1]+'-'+filename2.split('.txt')[0].split('_')[-1]
        else:
            files = [filename,filename2, filename3]
            string = filename.split('_')[1]+'-'+filename2.split('.txt')[0].split('_')[-1]+'-'+filename3.split('.txt')[0].split('_')[-1]
    else:
        if filename4 == '':
            files = [filename,filename2, filename3]
            string = filename.split('_')[1]+'-'+filename2.split('.txt')[0].split('_')[-1]+'-'+filename3.split('.txt')[0].split('_')[-1]
        else:
            files = [filename,filename2,filename3,filename4]
            string = filename.split('_')[1]+'-'+filename2.split('.txt')[0].split('_')[-1]+'-'+filename3.split('.txt')[0].split('_')[-1]+'-'+filename4.split('.txt')[0].split('_')[-1]
    typeline = {filename : '-', filename2: '+', filename3: '*',filename4: 'x'}
    atoms = ['Na-O','Ca-O','Al-O','Si-O']
    nlines = 3
    size_figure = (6,2.3*nlines) #height of the figure depends on the number of pair we display
    offset = 0 #offset for printing text on Na-O K-O subplot
    plt.close(1)
    fig, (ax1,ax2,ax3) = plt.subplots(nrows = nlines, ncols = 1, sharex = True, sharey = False, figsize=size_figure)
    major_xticks = np.arange(0, 8, 0.5) 
    minor_xticks = np.arange(0, 8, 0.1)
    for ax in [ax1,ax2,ax3]:
        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)                                           
        ax.xaxis.set_ticks_position('both')
        majorLocator = AutoLocator()
        minorLocator = AutoMinorLocator()                                                        
        ax.yaxis.set_ticks_position('both')
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_minor_locator(minorLocator)
        ax.set_xlim(0.9,4.4)
        plt.autoscale(axis='y')
        ax.tick_params(which = 'both', labelsize = plot_parameters["size_font_ticks"], width = plot_parameters["size_lines"]/2)   
    # Fine-tune figure : 1) make subplots close to each other (+ less whitespace around), 2) hide x ticks for all but bottom plot, 3) add a big invisible subplot in order to center x and y labels (since the ticklabels are turned off we have to move the x and y labels with labelpad)
    fig.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1, hspace = 0.03, wspace = 0.02)
    ax0 = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    ax0.set_xlabel(r'Density (g.cm$^{-3}$)', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], labelpad = plot_parameters["shift_labelpad"])
    ax0.set_ylabel(r'Bond length ($\AA$)', fontweight = 'bold', fontsize = plot_parameters["size_fonts"], labelpad = plot_parameters["shift_labelpad"]*1.5)
    figurename = 'bondlength_'+ string
    #************** Extract our data and plot 
    #*** Initialization of data dictionnaries
    pair_columns = {}
    data = {} #big dictionnary with inside coord  all coresponding to different T and atom pairs
    rho = {}  #idem for rho
    enveloppe = {}
    rhoenveloppe = {}
    yvariable = 1 #first file --> xmax
    linestyle = '-'
    colorline = '#368400ff'
    colorfill = '#3684007f'
    for iteration in range(3):
        #*** Calculation of the molecular mass
        with open(filename, 'r') as f:
            line = f.readline()
            entry=line.split()
            elements = entry[1:]
            line = f.readline()
            entry=line.split()
            number = entry[1:]
            line = f.readline()
            entry=line.split()
            pair = entry[1:]
        MN = 0
        for i in range(len(elements)):
            MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]
        #******** extraction of data        
        pair_columns[pair[0]] =yvariable
        for i in range(1,len(pair)):
            if pair[i] != pair[i-1]:
                pair_columns[pair[i]] = i+yvariable
        print('All the pairs available in the file with their index are:', pair_columns)
        for atom_pair in pair_columns:
            for atom in atoms:
                if atom_pair == atom:
                    data[atom] = {}
                    rho[atom] = {}
                    enveloppe[atom] = {}
                    rhoenveloppe[atom]= {}
        with open(filename,'r') as f:
            [f.readline() for i in range(4)]
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    temperature, acell = split_name(entry[0]) 
                    for atom_pair in pair_columns:
                        for atom in atoms:
                            if atom_pair == atom:
                                try:
                                    data[atom][temperature].append( float(entry[pair_columns[atom]]) )
                                    rho[atom][temperature].append( MN/(Na*float(acell)**3*10**(-24)) )    #calculation density
                                except KeyError:
                                    data[atom][temperature] = [ float(entry[pair_columns[atom]]) ]
                                    rho[atom][temperature] = [ MN/(Na*float(acell)**3*10**(-24)) ]    #calculation density
        #***** Compute the enveloppe of our data
        if plot_area == 1:
            for atom in data:
                #take all the data together
                alldata = []
                allrho = []
                for temp in data[atom]:
                    alldata.extend(data[atom][temp])
                    allrho.extend(rho[atom][temp])
                allrho, alldata = zip(*sorted(zip(allrho, alldata)))
                #extract enveloppe
                enveloppe[atom]['up'] = [alldata[0]]
                enveloppe[atom]['down'] = [alldata[0]]
                rhoenveloppe[atom] = [allrho[0]]
                j = 0
                for i in range(1,len(allrho)):
                    if allrho[i] == allrho[i-1]:
                        if alldata[i] > enveloppe[atom]['up'][j]:
                            enveloppe[atom]['up'][j] = alldata[i]
                        if alldata[i] < enveloppe[atom]['down'][j]:
                            enveloppe[atom]['down'][j] = alldata[i]
                    else:
                        j+=1
                        rhoenveloppe[atom].append(allrho[i])
                        enveloppe[atom]['up'].append(alldata[i])
                        enveloppe[atom]['down'].append(alldata[i])
                #print('len enveloppe',len(enveloppe[atom]['up']),len(enveloppe[atom]['down']), len(rhoenveloppe[atom]))
        #***** Plot each atom pair on the correct subplot --> ax1 for Na/Ca-O, ax2 for Al-O, ax3 for Si-O
        for temperature in data['Si-O']:
            for atom in data:
                if atom == 'Ca-O' or atom == 'Na-O':
                    if plot_area == 1:
                        ax1.fill_between(rhoenveloppe[atom],enveloppe[atom]['up'],enveloppe[atom]['down'], facecolor = colorfill, linewidth=plot_parameters['size_lines'], alpha = 0.1)
                        ax1.plot(rhoenveloppe[atom],enveloppe[atom]['down'], linestyle = linestyle, color=colorline, linewidth=plot_parameters["size_lines"] )
                        ax1.plot(rhoenveloppe[atom],enveloppe[atom]['up'], linestyle = linestyle, color=colorline, linewidth=plot_parameters["size_lines"] )
                    else:
                        ax1.plot(rho[atom][temperature],data[atom][temperature], linestyle = linestyle, color=colors_T[temperature], markersize = plot_parameters["size_markers"], linewidth=plot_parameters["size_lines"] )
                elif atom == 'Si-O':
                    if plot_area == 1:
                        ax3.fill_between(rhoenveloppe[atom],enveloppe[atom]['up'],enveloppe[atom]['down'], facecolor = colorfill, linewidth=plot_parameters['size_lines'], alpha = 0.1)
                        ax3.plot(rhoenveloppe[atom],enveloppe[atom]['down'], linestyle = linestyle, color=colorline, linewidth=plot_parameters["size_lines"] )
                        ax3.plot(rhoenveloppe[atom],enveloppe[atom]['up'], linestyle = linestyle, color=colorline, linewidth=plot_parameters["size_lines"] )    
                    else:
                        ax3.plot(rho[atom][temperature],data[atom][temperature], linestyle = linestyle, color=colors_T[temperature], markersize = plot_parameters["size_markers"], linewidth=plot_parameters["size_lines"] )
                    #ax3.text( 0.02,0.91, atom ,transform=ax3.transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  ) #for coord
                    ax3.text( 0.89,0.91, atom ,transform=ax3.transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  ) #for bond
                elif atom == 'Al-O':
                    if plot_area == 1:                    
                        ax2.fill_between(rhoenveloppe[atom],enveloppe[atom]['up'],enveloppe[atom]['down'], facecolor = colorfill, linewidth=plot_parameters['size_lines'], alpha = 0.1)
                        ax2.plot(rhoenveloppe[atom],enveloppe[atom]['down'], linestyle = linestyle, color=colorline, linewidth=plot_parameters["size_lines"] )
                        ax2.plot(rhoenveloppe[atom],enveloppe[atom]['up'], linestyle = linestyle, color=colorline, linewidth=plot_parameters["size_lines"] )    
                    else:
                        ax2.plot(rho[atom][temperature],data[atom][temperature],linestyle = linestyle, color=colors_T[temperature], markersize = plot_parameters["size_markers"], linewidth=plot_parameters["size_lines"] )
                    #ax2.text(0.02,0.91, atom ,transform=ax2.transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  ) #for coord
                    ax2.text(0.89,0.91, atom ,transform=ax2.transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  ) #for bond
        if iteration == 0:
            linestyle = '--'
            colorline = '#368400ff'
            colorfill = '#3684007f'
            yvariable = 2 #second run --> average
        elif iteration == 1:
            linestyle = '-.'
            colorline = '#6fc933ff'
            colorfill = '#6fc9337f'
            yvariable = 3 #third run --> median
        elif iteration == 2:
            linestyle = ':'
            colorline = 'k'
            colorfill = 'gray'
            yvariable = 4 #fourth run --> bondr3
    #************** Extract the collected data and plot
    sources = {'gofrs-bond_deKoker2010.txt':'de Koker (2010)','gofrs-bond_Angel1990.txt':'Angel $et~al.$ (1990)','gofrs-bond_Taylor1979.txt':'Taylor and Brown (1979)'}
    print('files from other sources',files[1:])
    for file in files[1:]:
        print('For file ',file)
        #***** Initialization of data dictionnaries
        pair_columns={'Na-O':2,'Ca-O':2,'Al-O':3,'Si-O':4}
        data = {} #big dictionnary with inside coord  all coresponding to different T and atom pairs
        rho = {}  #idem for rho
        with open(file, 'r') as f:
            line = f.readline()
            entry=line.split()
            pairs = entry[2:9]
        print('All the pairs available in the file are:', pairs)
        #******** extraction of data        
        for atom_pair in pairs:
            for atom in atoms:
                if atom_pair == atom:
                    data[atom] = {}
                    rho[atom] = {}
        with open(file,'r') as f:
            f.readline()
            while True:
                line = f.readline()
                if not line: break
                else:
                    entry=line.split('\n')[0].split('\t')
                    temperature = str(int(float( entry[9])))
                    for atom_pair in pairs:
                        for atom in atoms:
                            if atom_pair == atom:
                                datapoint = float(entry[pair_columns[atom]])
                                #replace 0 by nan
                                if datapoint == 0.0:
                                    datapoint = float(np.nan)
                                #add datapoint to data dic
                                try:
                                    data[atom][temperature].append( datapoint )
                                    rho[atom][temperature].append( float(entry[0])/1000) 
                                except KeyError:
                                    data[atom][temperature] = [ datapoint ]
                                    rho[atom][temperature] = [ float(entry[0])/1000]   
        #***** Plot each atom pair on the correct subplot --> ax1 for Na/Ca-O, ax2 for Al-O, ax3 for Si-O
        for temperature in data['Si-O']:
            for atom in data:
                if atom == 'Ca-O' or atom == 'Na-O':
                    #print('plot M-O for ', temperature)
                    ax1.plot(rho[atom][temperature],data[atom][temperature], marker = typeline[file][0], linestyle = 'None', color=colors_T[temperature], markersize = plot_parameters["size_markers"] )
                elif atom == 'Si-O':
                    #print('plot Si-O for ', temperature)
                    ax3.plot(rho[atom][temperature],data[atom][temperature], marker = typeline[file][0], linestyle = 'None', color=colors_T[temperature], markersize = plot_parameters["size_markers"])
                elif atom == 'Al-O':
                    #print('plot Al-O for ', temperature)
                    ax2.plot(rho[atom][temperature],data[atom][temperature], marker = typeline[file][0], linestyle = 'None',  color=colors_T[temperature], markersize = plot_parameters["size_markers"])
    for atom in ['Ca-O']:#,'Na-O']:
        #ax1.text( 0.02,0.91-offset, atom,transform=ax1.transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  ) #for coord
        ax1.text( 0.89,0.91-offset, atom,transform=ax1.transAxes, fontsize=plot_parameters['size_fonts'], fontweight='bold'  ) #for bond
        offset = offset + 0.10        #we update the offset for printing text on Na-O K-O subplot after finishing the first file              
    #Legend 
    legend_labels = {}
    legend_labels2 = {}
    for file in files[1:]: 
        legend_labels[sources[file]] =  plt.Line2D((0,1),(0,0), color='k', linestyle = 'None', markersize = plot_parameters["size_markers"]-2, marker=typeline[file][0])
    if plot_area == 1:
        legend_labels2['xmax'] =  mpatches.Patch(facecolor='#3684007f',  linestyle = '-',edgecolor='#368400ff')
        legend_labels2['average'] = mpatches.Patch(facecolor='#3684007f', linestyle = '--',edgecolor='#368400ff')
        legend_labels2['median'] = mpatches.Patch(facecolor='#6fc9337f', linestyle = '-.',edgecolor='#6fc933ff')
    else:        
        legend_labels2['xmax'] =  plt.Line2D((0,1),(0,0), color='k', linestyle='-', linewidth = plot_parameters["size_lines"])
        legend_labels2['average'] =  plt.Line2D((0,1),(0,0), color='k', linestyle='--', linewidth = plot_parameters["size_lines"])
        legend_labels2['median'] =  plt.Line2D((0,1),(0,0), color='k', linestyle='-.', linewidth = plot_parameters["size_lines"])
    s = [(k, legend_labels[k]) for k in sorted(legend_labels.keys(),reverse = False)]      
    s2 = [(k, legend_labels2[k]) for k in sorted(legend_labels2.keys(),reverse = False)]      
    #legend1 = ax1.legend([v for k,v in s],[k for k,v in s], loc='lower left', bbox_to_anchor=(0.3, 1.02), ncol=1, fontsize = plot_parameters["size_font_ticks"], borderaxespad=0.)
    #legend2 = ax1.legend([v for k,v in s2],[k for k,v in s2], loc='lower left', bbox_to_anchor=(0.01, 1.02), ncol=1, fontsize = plot_parameters["size_font_ticks"], borderaxespad=0., title='This study')
    #plt.setp(legend2.get_title(),fontsize= plot_parameters["size_font_ticks"], fontweight = 'bold')
    legend1 = ax1.legend([v for k,v in s],[k for k,v in s], loc='lower left', bbox_to_anchor=(0.3, 0.02), ncol=1, fontsize = plot_parameters["size_font_ticks"]-2, borderaxespad=0.)
    legend2 = ax1.legend([v for k,v in s2],[k for k,v in s2], loc='lower left', bbox_to_anchor=(0.01, 0.02), ncol=1, fontsize = plot_parameters["size_font_ticks"]-2, borderaxespad=0.)
    
    ax1.add_artist(legend1)
    
    #save the figure
    figurename = figurename + '.pdf'
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    
    print(figurename, '   created')
    #plt.show()
   
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])
