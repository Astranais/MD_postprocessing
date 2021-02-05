#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 14:11:32 2020

@author: akobsch
"""
import sys, glob, getopt



def extract_TP(thermofile, column_number, TP, addtxt):
    """extract temperature and pressure from thermofile"""
    with open(thermofile, 'r') as f:
        [f.readline() for i in range(3)]
        #extract data
        while True:
            line = f.readline()
            if not line: break
            else:    
                entry=line.split('\n')[0].split('\t')
                TP[entry[0].split('outcar.umd.dat')[0].split('/')[-1]+addtxt] = (int(entry[column_number['T']]),float(entry[column_number['P']]))
    return TP



def main(argv):
    """     ********* Main program *********     """
    #dictionnaries 
    TP = {}
    #dictionnary for fullaverages file in full version
    column_number = {'rho':2,'P':5,'stdev_P':6,'err_P':7,'T':8,'stdev_T':9,'err_T':10,'E':14,'stdev_E':15,'err_E':16,'Cvm_Nkb':23,'stdev_Cvm_Nkb':24,'Cvm':25,'stdev_Cvm':26,'testCv':31,'stdev_testCv':32}
    try:
        options,arg = getopt.getopt(argv,"hf:m:",["filename","mineralfile"])
    except getopt.GetoptError:
        print("vibr2diffusion.py -f <thermofilename> -m <mineralfile> ")
        sys.exit()
    for opt,arg in options:
        if opt == '-h':
            print('')
            print('vibr2diffusion.py program to create diffusion file from all the values obtained from vibr_sectrum_umd.py ')
            print("vibr2diffusion.py -f <thermofilename> -m <mineralfile> ")
            print("vibr2diffusion.py requires to be lauched from the folder containing every diff.dat file created by the script vibr_sectrum_umd.py")
            print('')
            sys.exit()
        if opt in ('-f','--filename'):
            thermofilename = str(arg)
        elif opt in ('-m','--mineralfile'):
            mineralfile = str(arg)
    #******* Extract TP values from thermofilename obtained from fullaverages.py
    TP = extract_TP(thermofilename, column_number, TP, '')
    #******* Initialisation of the newfile
    newfilename = 'diffusivities-vibr_'+thermofilename.split('/')[-1].split('_')[1]+'.txt'
    nf = open(newfilename,'w') 
    with open(mineralfile,'r') as mf:
        while True:
            line = mf.readline()
            if not line: break
            else:
                nf.write(line)
                entry = line.split()
                if entry[0] == 'elements':
                    elements = entry[1:]
    string = ''
    for elem in elements:
        string += '\tD_'+ elem+'(m^2/s)'
    nf.write("file\tP(GPa)"+string+'\n')
    #******* Extract all the data from diff.dat files and write in new file
    files = sorted(glob.glob('*.diff.dat'),reverse=False)
    for file in files:
        filename = file.split('/')[-1].split('outcar')[0]
        with open(file, 'r') as f:
            line = f.readline()
        nf.write(filename+'outcar\t'+str(TP[filename][1])+'\t'+line)
    print("diff file ", newfilename, "written")
    
#    ********* Execution *********  
if __name__ == "__main__":
    main(sys.argv[1:])
  

        
        
        
        
        
        
        
        
        
        
        
        
        
        