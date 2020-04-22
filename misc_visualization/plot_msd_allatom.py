#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anais
Langage : Python3

                ****     Plot the mean square displacement of each atoms at a selected T and rho       ****
                            each element on separated subplot
                            

"""



#     ********* Importation of the packages and modules used here *********     """
import sys
import getopt
import glob
import numpy as np
import matplotlib.pyplot as plt
import crystallography as cr
import re
import labellines as ll
#if you have ModuleNotFoundError: No module named 'labellines'
#then do in your terminal python3 -m pip install matplotlib-label-lines





def label_line(ax, line, label, halign, color='0.35', fs=12):
    """Add an annotation to the given line with appropriate placement and rotation.
    Based on code from:
        [How to rotate matplotlib annotation to match a line?]
        (http://stackoverflow.com/a/18800233/230468)
        User: [Adam](http://stackoverflow.com/users/321772/adam)
    Arguments
    ---------
    ax : `matplotlib.axes.Axes` object
        Axes on which the label should be added.
    line : `matplotlib.lines.Line2D` object
        Line which is being labeled.
    label : str
        Text which should be drawn as the label.
    ...
    Returns
    -------
    text : `matplotlib.text.Text` object
    """
    xdata, ydata = line.get_data()
#    x1 = xdata[0]
#    x2 = xdata[-1]
#    y1 = ydata[0]
#    y2 = ydata[-1]
    
    #formatage of label
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
                label = label[:i]+'$_{'+label[i+1:i+1+num]+'}$' + label[i+1+num:]
                i = i+5
            i = i+1
        else:break


    if halign.startswith('l'):
        x1 = xdata[0]
        x2 = xdata[1]
        y1 = ydata[0]
        y2 = ydata[1]
        xx = x1
        halign = 'left'
    elif halign.startswith('r'):
        x1 = xdata[-2]
        x2 = xdata[-1]
        y1 = ydata[-1]
        y2 = ydata[-1]
        xx = x2
        halign = 'right'
    elif halign.startswith('c'):        
        x1 = xdata[int(len(xdata)/2)-1]
        x2 = xdata[int(len(xdata)/2)]
        y1 = ydata[int(len(xdata)/2)-1]
        y2 = ydata[int(len(xdata)/2)]
        if ax.get_xscale() == 'log':
            xx = 10**(0.5*(np.log10(x1) + np.log10(x2)))
        else:
            xx = 0.5*(x1 + x2)
        halign = 'center'
    else:
        raise ValueError("Unrecogznied `halign` = '{}'.".format(halign))

    if ax.get_xscale() == 'log' and ax.get_yscale() == 'log':
        yy = 10**(np.interp(np.log10(xx), np.log10(xdata), np.log10(ydata)))
    elif ax.get_xscale() == 'log' and ax.get_yscale() != 'log':
        yy = np.interp(np.log10(xx), np.log10(xdata), ydata)
    else:
        yy = np.interp(xx, xdata, ydata)

    ylim = ax.get_ylim()
    # xytext = (10, 10)
    xytext = (0, 0)
    text = ax.annotate(label, xy=(xx, yy), xytext=xytext, textcoords='offset points',
                       size=fs, color=color, zorder=1,
                       horizontalalignment='left',#halign,
                       verticalalignment='center_baseline')

    sp1 = ax.transData.transform_point((x1, y1))
    sp2 = ax.transData.transform_point((x2, y2))

    rise = (sp2[1] - sp1[1])
    run = (sp2[0] - sp1[0])

    slope_degrees = np.degrees(np.arctan2(rise, run))
    text.set_rotation_mode('anchor')
    text.set_rotation(slope_degrees)
    ax.set_ylim(ylim)
    return text


def split_name(filename):
    """ Function to extract T and acell from filename """
    # *********** My filename are in the format CaAl2Si2O8_T3_nvt_a12.0.outcar.gofr.dat
    # ******* so I can split their name with _ and tale the acell and T from their name
    # **** Then change these two lines in harmony with your filename
    temperature = filename.split('_')[1]                                     #I extract the corresponding temperature
    acell = filename.split('.outcar.msd')[0].split('_')[3].strip('a')   #I extract the corresponding acell
    return temperature, acell



def main(argv):
    """     ********* Main program *********     """
    #parameters for the figure for article
    plot_parameters = {"size_fonts" : 12,"size_font_ticks":10,"size_figure" : (4,4),"size_markers" : 4,"size_lines" : 1,"shift_labelpad" : 10}
    #other dictionnaries and parameters for the figure
    lines = {}   #dictionnary for the lines for legend (the keys are acell or densities, and values are the color of lines)
    limits = {} #dictionnary for y limits  
    index = []
    Na=6.022*10**23  #avogadro constant
    try:
        options,arg = getopt.getopt(argv,"hm:f:a:i:",["mineralfile","file","atoms",'index'])
    except getopt.GetoptError:
        print("plot_msd_allatom.py  -m <mineralfile> -f <file> -a <atoms> -i <index atoms>")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print("plot_msd_allatom.py script to plot all the msd of all atoms for one file")
            print("plot_msd_allatom.py  -m <mineralfile> -f <file> -a <atoms> -i <index atoms>")
            sys.exit()
        elif opt in ("-m","--mineralfile"):
            mineralfile = str(arg)
        elif opt in ("-f","--file"):
            file = str(arg)
        elif opt in ("-a","--atoms"):
            elements = arg.split(',') 
        elif opt in ('-i','--index'):
            index = arg.split(',')
    #extraction of elements 
    with open(mineralfile,'r') as mf:
        entry = mf.readline()
        allelements = entry.split()[1:]
        entry = mf.readline()
        number = entry.split()[1:]
    MN = 0
    for i in range(len(allelements)):
        MN = MN + float(number[i]) * cr.Elements2rest(allelements[i])[3]  #I compute the molecular mass
    #calculation of densities in g/cm3
    temperature, acell = split_name(file)  
    Volume = float(acell)**3
    Density = MN/(Na*Volume*10**(-24))
    with open(file, 'r') as f:
        atoms = f.readline()
        atoms = atoms.strip('time_(fs)').split()                            #I extract the atoms
    #creation_plot
    plt.close(1)
    if len(elements) == 1: #width of the figure depends on the number of element we display
        w = 4
    elif len(elements) == 2:
        w = 7 
    else:
        w = 10
    fig = plt.figure(1,figsize = (w,4))
    for element in elements:
        #change of subplot
        plt.subplot(1,len(elements),elements.index(element)+1)
        ax = plt.gca() 
        plt.subplots_adjust(top = 0.91, bottom = 0.1, right = 0.98, left = 0.1, hspace = 0, wspace = 0.2)
        for atom in atoms:
            if atom == element: #if the current atom name is exactly the current element then draw the msd with a big line
                print(atom,'MSD elem')
                Temps,DataAtom = np.loadtxt(file,skiprows=1,usecols=(0,atoms.index(atom)+1),unpack=True)
                ax.plot(Temps,DataAtom,'-', color = 'r', linewidth = plot_parameters['size_lines']*2, label = 'mean msd')
                ax.set_xlabel(elements[elements.index(element)], fontweight = 'bold', fontsize = plot_parameters['size_fonts'])
                ax.xaxis.set_label_position('top')
            elif atom[:len(element)] == element: #if the first caracters of the current atom name matches the current element then draw msd
                print(atom,'MSD atom')
                Temps,DataAtom = np.loadtxt(file,skiprows=1,usecols=(0,atoms.index(atom)+1),unpack=True)
                if index == []:
                    line, = ax.plot(Temps,DataAtom,'-', color = '0.5',linewidth = plot_parameters['size_lines'], label = atom)
                    label_line(ax, line, atom, halign='right')  
                else:
                    indexatom = re.search('\d+',atom).group(0)
                    if indexatom in index:
                        print(atom,'MSD index')
                        line, = ax.plot(Temps,DataAtom,'-', color = '#068cff',linewidth = plot_parameters['size_lines'], label = atom)
                        label_line(ax, line, atom, halign='right', color = '#068cff')
                    else:
                        line, = ax.plot(Temps,DataAtom,'-', color = '0.5',linewidth = plot_parameters['size_lines'], label = atom)
        
        #xvals = []
        #[xvals.append(14000) for i in range(len(atoms))]
        #ll.labelLines(plt.gca().get_lines(),align=False, xvals=xvals,color='k')    
        
        #limitation of data along y
        if elements.index(element) == 0:
            ax.autoscale(axis='y', tight=True)
            limits[temperature] = ax.get_ylim()
        ax.set_ylim(limits[temperature])
        #we make the graph prettier
        if  elements.index(element) != 0:
            plt.setp(ax.get_yticklabels(), visible=False)
        #if elements.index(element) == len(elements)-1:
        #    ax.set_ylabel(temperature.strip('T')+'000 K' +'\n' + r'  $\rho$ = '+ str(round(Density,2)) + r' (g.cm$^{-3}$)' , fontweight = 'bold', fontsize = 12, rotation  = 0, labelpad = 20)
        #    ax.yaxis.set_label_position('right')    
    
   
    # Fine-tune figure
    ax0 = fig.add_subplot(111, frameon=False)
    plt.tick_params(labeltop=False, top=False, labelbottom=False, bottom=False, labelleft=False, left=False, labelright =False, right=False)
    ax0.set_xlabel(r'Time (fs)', fontweight = 'bold', fontsize = plot_parameters['size_fonts'], labelpad = 20)
    ax0.set_ylabel(r"Mean square displacement ($\mathregular{\AA^2}$)",fontsize=plot_parameters['size_fonts'],fontweight='bold', labelpad = 35)
    title = file.split('.outcar.msd.dat')[0]+r'  $\rho$ = '+ str(round(Density,2)) + r' (g.cm$^{-3}$)'   
    #ax0.set_title(title ,fontsize=plot_parameters['size_fonts'],fontweight='bold' , y = 1.1)
    
    figurename = 'msd_allatoms_'+file+'.png'
    plt.savefig(figurename, bbox_inches='tight', dpi = 150)
    print("the figure ",figurename, "is created")
    #plt.show()
            
 

#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])























