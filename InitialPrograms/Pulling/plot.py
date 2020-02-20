# prepared by Sam. feel free to consult (sirmaxford@gmail.com).
import fileinput, sys, shutil, os, time, socket, subprocess
import matplotlib; matplotlib.use('Agg') #set matplotlib to not use the Xwindows backend
import matplotlib.pyplot as plt; import pandas as pd; import numpy as np
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }

results_dir = 'PLOTS'
plot_name = 'Graph'

try:  
    os.mkdir(results_dir)
except OSError:  
    print ("=> Creation of the directory: %s failed" % results_dir)
else:  
    print ("=> Successfully created %s directory." % results_dir)

columns = ['time', 'x_tip', 'y_tip']; DtFile=0.01
df_load = pd.read_csv('TipXY_A001.txt', names=columns, delim_whitespace=True)
df_nice = df_load.drop(['time'], axis=1)
df = df_nice.iloc[0::10, :] #pick every 9th row starting from the first
df.to_csv(results_dir+'/SkippedTipXY_A001.csv', encoding='utf-8', index=False)
x_max = df['x_tip'].max(); x_min = df['x_tip'].min()
y_max = df['y_tip'].max(); y_min = df['y_tip'].min()
x1 = (0.01*(x_max-x_min))+x_min
y1 = (0.90*(y_max-y_min))+y_min
x2 = (0.15*(x_max-x_min))+x_min 
y2 = (0.85*(y_max-y_min))+y_min 

Dx_tip = np.diff(df['x_tip']); Dy_tip = np.diff(df['y_tip'])
DD=np.sqrt((Dx_tip**2)+(Dy_tip**2))
np.savetxt(results_dir+'/DD_A001.csv', DD, delimiter=',')
v=DD/(10*DtFile); Av_vel = np.mean(v)
np.savetxt(results_dir+'/v.csv', v, delimiter=',')
vSD=np.sum(((v-Av_vel)**2)/(np.size(v)-1)); vSD=np.sqrt(vSD)

plt.figure(figsize=(10,8))
plt.text(x1, y1, 'Av_Vel = ', fontdict=font); plt.text(x2, y1, '%.5f'%Av_vel, fontdict=font)
plt.text(x1, y2, 'Vel_SD = ', fontdict=font); plt.text(x2, y2, '%.5f'%vSD, fontdict=font)
plt.plot(df['x_tip'],df['y_tip'], label='Leading tip', color='green', marker='o', linestyle='dashed', linewidth=2, markersize=7)
plt.xlabel('X', fontdict=font); plt.ylabel('Y', fontdict=font)
plt.title('Actin Filament Movement '+plot_name); plt.legend(loc='upper left'); plt.grid()
plt.savefig(results_dir+'/'+plot_name+'.svg', bbox_inches='tight', format='svg', dpi=500)
plt.savefig(results_dir+'/'+plot_name+'.png', bbox_inches='tight', format='png', dpi=500)

print ("=> %s successfully saved. Done!\n" %results_dir)
