"""Written Charles Takalani:

24/11/2016

Written to analyse data of 21cm boxes working with the SZE-21cm"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os
import csv
import glob
from random import randint

#This code takes a comptonised spectrum in flux units
#then finds a corresponding redshift, inserting the spectrum signal as a masked region i.e. "the cluster signal"
#The end result of this code gives us the results in Colafrancesco, Emritte and Paolo, July 2016 paper
#The sze-21cm

#This version picks out only redshifts in the spectrum file

k = 1.32*(10**(-23)) #Boltzman constant
h = 6.63*(10**(-34)) #Planck constant
c = 3.00*(10**(8)) #speed of light
To = 2.725 #CMB Temperature in Kelvin

#Bring in first the comptonised spectrum and put this into three arrays


f = open('2paolocmbcompt_21cmdt6_10kev_1em3.csv', 'rb')
reader = csv.reader(f)
freq = []
redshift = []
mod = []
for row in reader:
    freq.append(float(row[0]))
    redshift.append(float(row[1]))
    mod.append(float(row[2]))
f.close()

print freq
print redshift
print mod


#The directory that contains Boxes to work with

wd = '/home/charles/Documents/simulation_test_folder/data/hunds'
os.chdir(wd)

for i in os.listdir(os.getcwd()):
        if i.startswith('delta_T_v3_no_halos'):

                filename = i
                
                #Now we start working with the boxes
                
                print ' ****************** THE BEGINING OF WORKING ON NEW FILE***********************'
                print ('The current file: '+filename)
                
                #Create the directory where images and and output data will be stored
                
                print '*************creating new directory***************'
                location = '/home/charles/Documents/simulation_test_folder/pictures/'

                if not os.path.exists(location):
                        os.makedirs(location)

                z = float(filename[21:-86])
                
                #We get the redshift from the filename string
                
                print ('Printing redshift just to make sure' + str('%06.2f'%z))
                
                
                #Now this is the part we bring in the other patch of the sky that does not contain the simulated cluster

                s = '/home/charles/Documents/simulation_test_folder/data2/hunds'
                os.chdir(s)
                for j in os.listdir(os.getcwd()):
                    if j.startswith('delta_T_v3_no_halos_z'+str('%06.2f'%(z))):
                        myfile = j

                        print ('The other patch of the sky is from : '+myfile)

                        print ('The redshift from filename is: ' + str(z))
                        fi = float(1420/(1+z))
                        v = float("{0:.2f}".format(fi))
                        f = (v*(10**6))
                        print ('frequency in MHz is:' + str(v))
                        
                        #Lets pick only the boxes at the redshifts defined in the spectrum file

                        for i in range (len(redshift)):
                            #Working through redshifts
                            if redshift[i] == z:
                                x = mod[i]
                                I_mod=x
                                #p defined below is x in the paper mentioned
                                
                                p = ((h*f)/(k*To))
                                print ( 'x, the modification is:' +str(x))


                                #This is a box with the cluster signal, the unperturbed CMB, and a normal box


                                def box_with_cluster(ndim):
                                    shape = (ndim, ndim, ndim)

                                    print 'Now working with box with cluster...........'

                                    print 'creating unperturbed CMB box at this frequency..........'

                                    #The unperturbed CMB
                                    Io_st = (2*(((k*To)**3)/((h*c)**2))*((p**3)/((np.exp(p))-1))*(10**26)*(1.18*(10**-7)))*(np.ones(shape))

                                    print 'Done !!!!!!!'

                                    print 'Now bringing in the box.....'

                                    #The Box we want to work with

                                    fd = open('/home/charles/Documents/simulation_test_folder/data/hunds/'+filename, 'rb')

                                    dI = (((2*(k*(f**2))/(c**2)))*(10**(-3))*(10**26)*(1.18*(10**-7)))*(np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape))

                                    bc = dI + Io_st
                                    return bc
                                    
                                    print 'masking box with cluster'
                                
                                box = box_with_cluster(256)


                                print 'Calculating cluster temperature'
                                a =np.mean(box)
                                print 'Outputting values to file'

                                file = open(location +'SZE-21cmSpec' +str(x) +'.txt', 'w')
                                file.write(str(z))
                                file.write(' '+str(v))
                                file.write(' '+str(x))
                                file.write(' '+str(a))
                                file.close()

                                break
                            else:
                                x=0
                                print 'Now I hope this is gonna work and we are done'




                                
