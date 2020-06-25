# Python 2
# Find Filament, Specie1, and Specie2 .vtk files
# Get the clean points
# Save the points as .csv files
# Regards, Sirmaxford

import os, glob, csv, sys
from paraview.simple import *
import numpy as np

paraview.simple._DisableFirstRenderCameraReset()

current_path = os.getcwd()
flmt = 'Filament_'
mtr1 = 'Specie1_'
mtr2 = 'Specie2_'
r = 0.8 
seed = '113'
Fcount = 0
M1count = 0
M2count = 0

dirlist = ['R0.8B13.0ATP2000.0MD3000.0']
#dirlist = glob.glob(current_path+'/*/') # comment out this if only using the above one folder
#dirlist = sorted(dirlist, key=lambda x:x[-28:]) # comment out this one also if only using the above one folder

for i in dirlist:
    os.chdir(i); print(i)
    filament_list = glob.glob('Filament_A001T**.vtk') # os.listdir()
    filaments = sorted(filament_list, key=lambda x:x[-11:])
    fnos = len(filaments)
    specie1_list = glob.glob('IntSpecie1_A001T**.vtk') 
    specie1 = sorted(specie1_list, key=lambda x:x[-11:])
    s1nos = len(specie1)
    specie2_list = glob.glob('IntSpecie2_A001T**.vtk') 
    specie2 = sorted(specie2_list, key=lambda x:x[-11:])
    s2nos = len(specie2)
    if fnos == s1nos == s2nos:
        print('Number of files = %s'%fnos)
        for i in filaments:
            f = LegacyVTKReader(FileNames=[i]) # create a new 'Legacy VTK Reader'
            try:
                os.remove(flmt+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(Fcount,1))+'.csv')
                print('Deleted: '+flmt+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(Fcount,1))+'.csv')
            except(OSError, RuntimeError, TypeError, NameError):
                pass
            SaveData(flmt+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(Fcount,1))+'.csv', proxy=f, Precision=6)
            print('Saved: '+flmt+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(Fcount,1))+'.csv')
            Fcount = Fcount+1
        for j in specie1:
            m1 = LegacyVTKReader(FileNames=[j])
            try:
                os.remove(mtr1+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(M1count,1))+'.csv')
                print('Deleted: '+mtr1+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(M1count,1))+'.csv')
            except(OSError, RuntimeError, TypeError, NameError):
                pass
            SaveData(mtr1+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(M1count,1))+'.csv', proxy=m1, Precision=6)
            print('Saved: '+mtr1+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(M1count,1))+'.csv')
            M1count = M1count+1
        for k in specie2:
            m2 = LegacyVTKReader(FileNames=[k])
            try:
                os.remove(mtr2+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(M2count,1))+'.csv')
                print('Deleted: '+mtr2+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(M2count,1))+'.csv')
            except(OSError, RuntimeError, TypeError, NameError):
                pass
            SaveData(mtr2+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(M2count,1))+'.csv', proxy=m2, Precision=6)
            print('Saved: '+mtr2+seed+'R'+str(np.round(r,1))+'Ts'+str(np.round(M2count,1))+'.csv')
            M2count = M2count+1
    else:
        sys.exit('Error! The number of files are different!')

    r=r+0.1
    Fcount = 1
    M1count = 1
    M2count = 1


