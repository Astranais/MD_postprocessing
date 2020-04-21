#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot size of all clusters along with T and/or rho       ****
                         *************  ARTICLE VERSION  *************

"""



#     ********* Importation of the packages and modules used here *********     """
import sys
import getopt
import glob
import numpy as np
import matplotlib.pyplot as plt
import crystallography as cr
from matplotlib.colors import LinearSegmentedColormap




def split_name(filename):
    """ Function to extract T and acell from filename """
    # *********** My filename are in the format CaAl2Si2O8_T3_nvt_a12.0.outcar.umd.dat.r1.stat.dat
    # ******* so I can split their name with _ and tale the acell and T from their name
    # **** Then change these two lines in harmony with your filename
    filename = filename.strip('./')
    temperature = str(int(float(filename.split('_')[1].strip('T'))*1000))
    acell = filename.split('.outcar.umd.dat.r1.stat.dat')[0].split('_')[3].strip('a')
    return temperature, acell





    
def main(argv):
    """     ********* Main program *********     """
    xdata = {'T':[],'rho':[]} #dictionnary containing the values of T and rho for the subplots
    #other parameters
    Na=6.022*10**23
    mineralfile = ''
    letter = ''
    try:
        options,arg = getopt.getopt(argv,"hm:v:l:s:",["mineralfile","variable","lifetime","subfig"])
    except getopt.GetoptError:
        print("plot_speciation-r1-size.py  -m <mineralfile with elements> -v <variable (rho,T)> -l <max_lifetime for colorbar (in fs)> -s <subfigure letter (option, default = ' ')>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print("plot_speciation-r1-size.py  Plot size of all clusters along with T and/or rho    ")
            print("plot_speciation-r1-size.py  -m <mineralfile with elements> -v <variable (rho,T)> -l <max_lifetime for colorbar (in fs)> -s <subfigure letter (option, default = ' ')>")
            print('requires the file containing elements and number (in order to compute the densities)')
            print('requires to be in the folder T or acell containing all the *outcar.umd.dat.r1.stat.dat files')
            sys.exit()
        elif opt in ("-m","--mineralfile"):
            mineralfile = str(arg)
        elif opt in ("-v", "--variable"):
            variable = str(arg)
        elif opt in ("-l", "--lifetime"):
            max_lifetime = int(arg)
        elif opt in ("-s", "--subfig"):
            letter = str(arg)
    #****** Calulcation of the molecular mass
    with open(mineralfile,'r') as mf:
        entry = mf.readline()
        elements = entry.split()[1:]
        entry = mf.readline()
        number = entry.split()[1:]
    MN = 0
    for i in range(len(elements)):
        MN = MN + float(number[i]) * cr.Elements2rest(elements[i])[3]  #I compute the molecular mass    
    #****** Extraction of data
    files = sorted(glob.glob('*outcar.umd.dat.r1.stat.dat'))                             #I list every speciation files in alphabetic order
    #creation of the list containing all T and rho only once
    temperature0, acell0 = split_name(files[0])
    xdata['rho'] = [acell0]
    xdata['T'] = [temperature0]
    for file in files:
        temperature, acell = split_name(file)
        if temperature != temperature0:
            xdata['T'].append(temperature)      #xdata['T'] = ['3000','4000','5000' ...]
            temperature0 = temperature
        counter = 0 #counter to check the presence of the same acell in the xdata list
        for a in xdata['rho']:
            if a == acell: #if the current acell is already in the list, then the counter is changed to 1
                counter = 1
        if counter == 0:
            xdata['rho'].append(acell)   #xdata['rho'] = ['12.0','13.0','16.5', ...]
    xdata['rho'] = sorted(xdata['rho'], reverse=False)
    xdata['T'] = sorted(xdata['T'], reverse=True)
    #****** Creation_plot     
    #Adjustment of ticks, labels etc.
    if variable == 'rho':
        number = xdata['T'][:]     #number = ['6000','5000', ...]
        major_ticks = []
        #major_ticks = np.arange(0, 7, 0.5) 
        #minor_ticks = np.arange(0, 7, 0.1)
        xlimits = (0.9,2.3)
        xlabel = "Density (g.cm$^{-3}$)"
        for i in range(len(xdata['rho'])):      
            xdata['rho'][i] = MN/(Na* float(xdata['rho'][i]) **3*10**(-24))  #calculation density, xdata['rho'] = [2.19, 1.98, 1.65, ...]
            major_ticks.append(round(xdata['rho'][i],2))
        note = xdata['T'][:]
        for i in range(len(note)):
            note[i] = str(note[i]) + ' K'  #note = ['6000 K','5000 K', ...]
    else:
        number = xdata['rho'][:]                #number = ['12.0','13.0','16.5', ...]
        major_ticks = np.arange(0, 10000, 1000) 
        minor_ticks = np.arange(0, 10000, 500)
        xlimits = (1500,7500)   
        xlabel = "Temperature (K)"    
        for i in range(len(xdata['rho'])):      
            xdata['rho'][i] = MN/(Na* float(xdata['rho'][i]) **3*10**(-24))  #calculation density, xdata['rho'] = [2.19, 1.98, 1.65, ...]
        note = xdata['rho'][:]
        for i in range(len(note)):
            note[i] = str(round(note[i],2)) + r' g.cm$^{-3}$'       #note = ['2.19 g/cm3','1.98 g/cm3', ...]     
    #Colormap definition
    #cm = LinearSegmentedColormap.from_list("", ["blue","green","red"])
    norm=plt.Normalize(0,max_lifetime)
    #plot
    plt.close(1)  
    h = 4 * len(number) #height of the figure depends on the number of T or rho we display
    fig = plt.figure(1,figsize = (7,h))
    plt.subplots_adjust(top = 0.97, bottom = 0.07, right = 0.89, left = 0.07, hspace = 0, wspace = 0)
    for subplot in range(1,len(number)+1):
        plt.subplot(len(number), 1, subplot)
        ax = plt.gca()
        for file in sorted(files, reverse=False):
            #print(file)
            temperature, acell = split_name(file)
            if variable == 'rho':
                fixed = 'T'+temperature0
                if temperature == number[subplot-1]: 
                    lifetime, length = np.loadtxt(file, usecols = (1,3), skiprows = 2, unpack = True ) 
                    if np.size(length) != 1:
                        result=sorted(zip(lifetime,length),key=lambda x: x[0]) #we sort the cluster lifetime by lifetime
                        lifetime, length = zip(*result)
                    #xarray = np.linspace(-np.size(length)/2,+np.size(length)/2,num=np.size(length))*0.001 +  np.ones(np.size(length)) *  MN/(Na* float(acell) **3*10**(-24)) 
                    xarray = (np.random.rand(np.size(length))*2-1) * (0.02)   +  (MN/(Na* float(acell) **3*10**(-24)))  #random vector around 0, * spread (here 0.02g/cm3), + center of the x data (acell) 
                    #xarray = np.ones(np.size(length)) *  MN/(Na* float(acell) **3*10**(-24)) 
                    try:
                        s = ax.scatter(xarray, length, s=20,  c=lifetime, facecolor = 'none', lw = 1, norm=norm, cmap = plt.cm.get_cmap('jet'))#cmap = cm)
                    except IndexError: #problem if we have only one data on lifetime = only have the melt
                        s = ax.scatter(xarray, length, s=20,  edgecolor=(0.5,0,0), facecolor = 'none', lw = 1)
                    s.set_facecolor('none')
            else:
                fixed = 'a'+acell0
                if acell == number[subplot-1]:
                    lifetime, length = np.loadtxt(file, usecols = (1,3,), skiprows = 2, unpack = True ) 
                    if np.size(length) != 1:
                        result=sorted(zip(lifetime,length),key=lambda x: x[0]) #we sort the cluster lifetime by lifetime
                        lifetime, length = zip(*result)
                    #xarray = np.linspace(-np.size(length)/2,+np.size(length)/2,num=np.size(length)) + np.ones(np.size(length)) *  float(temperature)
                    xarray = (np.random.rand(np.size(length))*2-1) * (100)   +   float(temperature)  #random vector around 0, * spread (here 100K), + center of the x data (T) 
                    #xarray = np.ones(np.size(length)) *  float(temperature)
                    try:
                        s = ax.scatter(xarray, length, s=20,  c=lifetime, facecolor = 'none', lw = 1, norm=norm, cmap = plt.cm.get_cmap('jet'))#cmap = cm)
                    except IndexError: #problem if we have only one data on lifetime = only have the melt
                        s = ax.scatter(xarray, length, s=20,  edgecolor=(0.5,0,0), facecolor = 'none', lw = 1)
                    s.set_facecolor('none')
        cb = fig.colorbar(s)
        cb.set_label('$\\bf{Total}$ $\\bf{lifetime}$ $\\bf{(fs)}$')
        ax.set_xticks(major_ticks)
        #ax.set_xticks(minor_ticks, minor=True)  
        ax.xaxis.set_tick_params(which = 'both', direction='inout')
        major_yticks = np.arange(0,210,50)
        minor_yticks = np.arange(0,210,10)
        ax.set_yticks(major_yticks)
        ax.set_yticks(minor_yticks, minor=True)  
        ax.yaxis.set_tick_params(which = 'both', direction='inout')
        plt.xticks(rotation = 45, fontsize = 10)
        
        ax.set_xlim(xlimits)
        ax.set_ylim(0,208)
        ax.grid(True, which='major',axis = 'y', linestyle=':', linewidth=0.5)
        
        #we make the graph prettier
        if  subplot != len(number):
            plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(which = 'both', labelsize = 10, width = 0.5)
        
        #plt.title(note[subplot-1], fontsize=15, fontweight='bold')
        if letter != '':
            ax.text(0.01,0.95, letter , transform=ax.transAxes, horizontalalignment = 'left', fontweight = 'bold', fontsize = 12, bbox=dict(facecolor='none', edgecolor='k', pad=3.0)) 
            
        ax.text(1.195,0.985, '+' , transform=ax.transAxes, horizontalalignment = 'right', fontweight = 'bold', fontsize = 10) 
        
    
    
    # Fine-tune figure
    ax0 = fig.add_subplot(111, frameon=False) 
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright = False, right=False)
    ax0.set_xlabel(xlabel, fontweight = 'bold', fontsize = 12, labelpad = 35)    
    ax0.set_ylabel('Cluster size', fontweight = 'bold', fontsize = 12 ,labelpad = 30)

      
    figurename = 'speciation-size_' + files[0].split('_')[0]  + '_' + fixed + '-v' + variable + '-l' + str(max_lifetime) + 'fs' +'.png'
    print(figurename, 'is created')
    fig.savefig(figurename, bbox_inches = 'tight', dpi = 300)
    #plt.show()





#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])























