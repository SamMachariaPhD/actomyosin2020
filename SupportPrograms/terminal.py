import os,sys
#==============================================
try:
    os.system("gnome-terminal -e 'python3 analysis.py' ")
    #os.system('python3 analysis.py')
except (Exception, e):
    print("Sorry, 'analysis.py' has an error.")
    pass
#==============================================
try:
    os.system("gnome-terminal -e 'pvpython film1.py' ")
except (Exception, e):
    print("Sorry, 'film1.py' has an error.")
    pass
#==============================================
try:
    os.system("gnome-terminal -e 'pvpython film2.py' ")
except (Exception, e):
    print("Sorry, 'film2.py' has an error.")
    pass
#==============================================
try:
    os.system("gnome-terminal -e 'python3 tar.py' ")
except (Exception, e):
    print("Sorry, 'tar.py' has an error.")
    pass
#==============================================
