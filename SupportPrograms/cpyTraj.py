"""
Go inside various complete simulation data.
Search for TipXY***.txt
Copy and paste to some place with a new name.
Regards, Sam Sirmaxford.
"""

import os, glob, shutil
import numpy as np

current_path = os.getcwd()
filename = 'TipXY_A001.txt'
md = 0.1 # smallest motor density
seed = '273ATP50'

dirlist = glob.glob(current_path+'/*/')
dirlist = sorted(dirlist, key=lambda x:x[-18:])

for i in dirlist:
    os.chdir(i)
    shutil.copy2(filename,current_path+'/'+'TipXY'+str(np.round(md,1))+'seed'+str(seed)+'.txt')
    md=md+0.1
