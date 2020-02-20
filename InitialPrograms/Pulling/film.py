# prepared by Sam. feel free to consult (sirmaxford@gmail.com)
import fileinput, sys, shutil, os, time, socket, subprocess, glob
from paraview.simple import *
import pandas as pd

columns = ['time', 'x_tip', 'y_tip']
df = pd.read_csv('TipXY_A001.txt', names=columns, delim_whitespace=True)
x_max = df['x_tip'].max(); x_min = df['x_tip'].min()
y_max = df['y_tip'].max(); y_min = df['y_tip'].min()

results_dir = 'FILM'
film_name = 'filament_film'

try:  
    os.mkdir(results_dir)
except OSError:  
    print ("=> Creation of the directory: %s failed" % results_dir)
else:  
    print ("=> Successfully created %s directory." % results_dir)

Filament_file_list = glob.glob('Filament**0.vtk') # all files starting with 'Filament' and ending with '.vtk'
MotorSpecie1_file_list = glob.glob('MotorSpecie1**0.vtk')
MotorSpecie2_file_list = glob.glob('MotorSpecie2**0.vtk')

filament_files = sorted(Filament_file_list, key=lambda x:x[-16:]) # sort by the last 16 characters of the file name
MotorSpecie1_files = sorted(MotorSpecie1_file_list, key=lambda x:x[-16:])
MotorSpecie2_files = sorted(MotorSpecie2_file_list, key=lambda x:x[-16:])

# create a new 'Legacy VTK Reader'
Filaments = LegacyVTKReader(FileNames=filament_files)
MotorSpecie1 = LegacyVTKReader(FileNames=MotorSpecie1_files)
MotorSpecie2 = LegacyVTKReader(FileNames=MotorSpecie2_files)

animationScene1 = GetAnimationScene() # get animation scene

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

renderView1 = GetActiveViewOrCreate('RenderView') # get active view
renderView1.ViewSize = [1546, 860] # set a specific view size

renderView1.AxesGrid.Visibility = 1
renderView1.UseTexturedBackground = 0

# set active source
SetActiveSource(Filaments)
FilamentsDisplay = Show(Filaments, renderView1) # show data in view
FilamentsDisplay.ColorArrayName = [None, ''] # trace defaults for the display properties.
FilamentsDisplay.GlyphType = 'Arrow'
FilamentsDisplay.ScalarOpacityUnitDistance = 3.0
FilamentsDisplay.LineWidth = 3.0
FilamentsDisplay.DiffuseColor = [0.24, 0.0, 0.0]

SetActiveSource(MotorSpecie1)
MotorSpecie1Display = Show(MotorSpecie1, renderView1)
MotorSpecie1Display.ColorArrayName = [None, '']
MotorSpecie1Display.GlyphType = 'Arrow'
MotorSpecie1Display.ScalarOpacityUnitDistance = 0.16
MotorSpecie1Display.Opacity = 0.3
MotorSpecie1Display.DiffuseColor = [0.0, 0.0, 1.0]

SetActiveSource(MotorSpecie2)
MotorSpecie2Display = Show(MotorSpecie2, renderView1)
MotorSpecie2Display.ColorArrayName = [None, '']
MotorSpecie2Display.GlyphType = 'Arrow'
MotorSpecie2Display.ScalarOpacityUnitDistance = 0.33
MotorSpecie2Display.Opacity = 0.8
MotorSpecie2Display.DiffuseColor = [0.69, 0.69, 0.0]

# create a new 'Annotate Time'
annotateTime1 = AnnotateTime()
SetActiveSource(annotateTime1)
annotateTime1Display = Show(annotateTime1, renderView1)
annotateTime1.Format = 'Time: 60 sec. [%.2f]'
annotateTime1Display = Show(annotateTime1, renderView1)

# current camera placement for renderView1
#renderView1.CameraPosition = [5.566566450844984, -3.271004214234651, 18.733021468562196]
#renderView1.CameraFocalPoint = [5.566566450844984, -3.271004214234651, 0.012500000186264515]
# reset view to fit data
renderView1.ResetCamera(0.0,x_max,y_min,0.0,x_min,y_max)
renderView1.CameraParallelScale = 1.5

# save animation images/movie
WriteAnimation(results_dir+'/'+film_name+'.avi', Magnification=1, FrameRate=7.0, Quality=70, Compression=True)
