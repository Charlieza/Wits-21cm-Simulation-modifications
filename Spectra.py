import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os
import csv
import glob

#This code takes a spectrum file and making a signal out of it at various
#redshifts and adds the signal to a delta_T output box from 21cmFAST

#This version picks out only redshifts in the spectrum files

#This has the foreground included#

###### Lets start by bringing in the spectrum #######



wd = '/Volumes/CHARLES/Boxes/theprog/hunds'
os.chdir(wd)

#Lets work on the Slice and other signal
for i in os.listdir(os.getcwd()):
        if i.startswith('delta_T_v3_no_halos'):

                filename = i
                print ' ****************** THE BEGINING OF WORKING ON NEW FILE***********************'
                print ('The current file: '+filename)
                print '*************creating new directory***************'
                location = '/Users/charles/Desktop/pics/'

                if not os.path.exists(location):
                        os.makedirs(location)

                z = float(filename[21:-86])
                Tb=float(filename[21:-86])
                print ('redshift is '+ str('%06.2f'%z))
                print ('Tb is '+ str('%06.2f'%z))
                file = open(location +'log' +str(x) +'.txt', 'w')
                file.write(str(z))
                file.write(' '+str(Tb))
                file.close()
