"""
Go inside various complete simulation data.
Search for bm*.csv
Copy and paste to some place with a new name.
Regards, Sam Sirmaxford.
"""

import os, glob, shutil
import numpy as np

current_path = os.getcwd()
file1 = 'bmAc.csv'; file2 = 'bminAc.csv'
change = 0.1 # change in the variable
name = '273s5'

dirlist = glob.glob(current_path+'/*/')
dirlist = sorted(dirlist, key=lambda x:x[-18:])

for i in dirlist:
    os.chdir(i); print(i)
    shutil.copy2(file1,current_path+'/'+'acBm'+str(np.round(change,1))+'_'+str(name)+'.csv')
    shutil.copy2(file2,current_path+'/'+'defBm'+str(np.round(change,1))+'_'+str(name)+'.csv')
    change=change+0.1
