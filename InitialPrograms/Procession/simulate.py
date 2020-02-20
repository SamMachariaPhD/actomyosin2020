# prepared by Sam. feel free to consult (sirmaxford@gmail.com).
import fileinput, sys, shutil, os, time, socket, subprocess

max_comp_simulations = 11

files = ['MotilityAssayActin2MotorsParameters_v5.f90', 'MotilityAssayConfinements_v1.f90', \
'MotilityAssayForceForceFunctions_v3.f90', 'MotilityAssaySubstrateDeformation_v2.f90', 'mt.f90', \
'MotilityAssayActin2MotorsMain_v9.f90', 'autocomp_py.sh', 'autorun_py.sh', 'analysis.py', 'film1.py', 'film2.py', 'tar.py', 'untar.py']

param_file = "MotilityAssayActin2MotorsParameters_v5.f90"
main_file = "MotilityAssayActin2MotorsMain_v9.f90"

pc_hostname = socket.gethostname()
date_today = time.strftime("%d-%m-%Y")
dir_name = date_today+pc_hostname #dynamic folder based on date and PC used

current_path = os.getcwd()  

run_range = range(1,max_comp_simulations+1)
simulations = int(input("\n=> Enter the total number (int) of simulations to be performed in this Computer today: "))
if simulations in run_range:
    print("\n=> %s simulations will be performed in this computer." %simulations)
else:
    sys.exit("\n=> Please enter a number (int) between 1 and %s" %max_comp_simulations)

print ("=> The current working directory is: %s" % current_path) 

try:  
    os.mkdir(dir_name)
except OSError:  
    print ("=> Creation of the directory: %s failed" % dir_name)
else:  
    print ("=> Successfully created today's simulation directory: %s\n" % dir_name)

simulation_runs = simulations
simulations_counter = 0

while simulation_runs > 0:
    params = open('param_set.txt')
    read_params = params.readlines()
    species_ratio, beads_number, ATP_value, MD_value = map(float,read_params[simulations_counter].split(',')) #dynamic parameters based on .txt file
    params.close()
    #species_ratio, beads_number, ATP_value, MD_value = map(float,input().split(',')) #Bad. I have to enter each param.
    new_dir = 'R'+str(species_ratio)+'B'+str(beads_number)+'ATP'+str(ATP_value)+'MD'+str(MD_value) #dynamic folder name
    os.chdir(dir_name)
    os.mkdir(new_dir)
    #print ("=> New folder successfully created: %s " % new_dir) #print outs for debugging
    os.chdir(current_path)
    for f in files:
        shutil.copy(f, dir_name+'/'+new_dir)
    #print ("=> Files successfully copied to the new folder:\n%s" % files)
    os.chdir(dir_name+'/'+new_dir)
    # replace all occurrences of 'ATP = 2000.0' with 'ATP = ATP_value'
    for i, line in enumerate(fileinput.input(param_file, inplace=1)):
        sys.stdout.write(line.replace('Species1Ratio = 0.60', 'Species1Ratio = '+str(species_ratio)))
    for i, line in enumerate(fileinput.input(param_file, inplace=1)):
        sys.stdout.write(line.replace('NumBeads = 13', 'NumBeads = '+str(beads_number)))
    for i, line in enumerate(fileinput.input(param_file, inplace=1)):
        sys.stdout.write(line.replace('ATP = 2000.0', 'ATP = '+str(ATP_value)))
    for i, line in enumerate(fileinput.input(param_file, inplace=1)):
        sys.stdout.write(line.replace('Motor_Density = 3000.0', 'Motor_Density = '+str(MD_value)))
    #print ("\n=> Param. file in %s successfully updated." % new_dir)
    #print ("\n=> Simulation program started!: %s" % new_dir)
    subprocess.call("ifort mt.f90 MotilityAssayActin2MotorsParameters_v5.f90 MotilityAssayConfinements_v1.f90 MotilityAssaySubstrateDeformation_v2.f90 MotilityAssayForceForceFunctions_v3.f90 MotilityAssayActin2MotorsMain_v9.f90", shell=True)
    #print ("=> Programs successfully compiled:\n%s " % files)
    tic=time.time()
    subprocess.call("ulimit -s unlimited;./a.out", shell=True)
    toc=time.time()
    tym=toc-tic
    print("\nTotal time is sec. = ", file=open('pyout.txt','a'))
    print(tym, file=open('pyout.txt','a'))
    print ("\n=> Programs in %s have successfully run complete!\n" % new_dir)
    simulation_runs = simulation_runs-1
    simulations_counter = simulations_counter+1 #simulation_counter for making dirs and prog. progress
    os.chdir(current_path)
    time.sleep(1)


print("\n=> All the %s simulations are successfully completed.\nDone!\n" %simulations)

"""
#-------------------------------PLOTTING THE OUTPUT DATA--------------------------------------------
plot_runs = simulations
plot_counter = 0
plot_file = ['plot.py']
ploting_file = "plot.py"

try:
    while plot_runs > 0:
        params = open('param_set.txt')
        read_params = params.readlines()
        species_ratio, beads_number, ATP_value, MD_value = map(float,read_params[plot_counter].split(','))
        params.close()
        new_dir = 'R'+str(species_ratio)+'B'+str(beads_number)+'ATP'+str(ATP_value)+'MD'+str(MD_value)
        print("\n=> Plotting results for %s " %new_dir)
        for f in plot_file:
            shutil.copy(f, dir_name+'/'+new_dir)
        os.chdir(dir_name+'/'+new_dir)
        for i, line in enumerate(fileinput.input(ploting_file, inplace=1)):
            sys.stdout.write(line.replace("'Graph'", "'"+new_dir+"'"))
        os.system('python3 plot.py')
        plot_runs = plot_runs-1
        plot_counter = plot_counter+1
        os.chdir(current_path)
        time.sleep(1)
except (OSError, RuntimeError, TypeError, NameError):
    print("Plotting script has an Error. ")
    print("Plotting not completed. ")
    pass

#------------------------------------------------------------------------------------------------------

#-------------------------------MAKE FILM FROM OUTPUT DATA--------------------------------------------
film_runs = simulations
film_counter = 0
film_file = ['film.py']
filming_file = "film.py"

try:
    while film_runs > 0:
        params = open('param_set.txt')
        read_params = params.readlines()
        species_ratio, beads_number, ATP_value, MD_value = map(float,read_params[film_counter].split(','))
        params.close()
        new_dir = 'R'+str(species_ratio)+'B'+str(beads_number)+'ATP'+str(ATP_value)+'MD'+str(MD_value)
        print("\n=> Making film for %s " %new_dir)
        for f in film_file:
            shutil.copy(f, dir_name+'/'+new_dir)
        os.chdir(dir_name+'/'+new_dir)
        for i, line in enumerate(fileinput.input(filming_file, inplace=1)):
            sys.stdout.write(line.replace("'filament_film'", "'"+new_dir+"'"))
        os.system('pvpython film.py')
        film_runs = film_runs-1
        film_counter = film_counter+1
        os.chdir(current_path)
        time.sleep(1)
except (OSError, RuntimeError, TypeError, NameError):
    print("Film script has an Error. ")
    print("Film not completed. ")
    pass

#------------------------------------------------------------------------------------------------------
"""

os.chdir(dir_name+'/'+new_dir)
#==============================================
try:
    os.system('python3 analysis.py')
except (Exception, e):
    print("Sorry, 'analysis.py' has an error.")
    pass
#==============================================
try:
    os.system('pvpython film1.py')
except (Exception, e):
    print("Sorry, 'film1.py' has an error.")
    pass
#==============================================
try:
    os.system('pvpython film2.py')
except (Exception, e):
    print("Sorry, 'film2.py' has an error.")
    pass
#==============================================
try:
    os.system('python3 tar.py')
except (Exception, e):
    print("Sorry, 'tar.py' has an error.")
    pass
#==============================================

