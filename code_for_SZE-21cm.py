"""Written Charles Takalani:

24/11/2016

Written to analyse data of 21cm boxes working with the SZE-21cm"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os
import csv
import glob

#This code takes a comptonised spectrum in flux units
#then finds a corresponding redshift, inserting the spectrum signal as a masked region i.e. "the cluster signal"
#The end result of this code gives us the results in Colafrancesco, Emritte and Paolo, July 2016 paper
#The sze-21cm

#This version picks out only redshifts in the spectrum file

k = 1.32*(10**(-23)) #Boltzman constant
h = 6.63*(10**(-34)) #Planck constant
c = 3.00*(10**(8)) #speed of light
To = 2.725*(10**3) #CMB Temperature in Kelvin

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

wd = '/mnt/storage1/simulations/charles/data/hunds'
os.chdir(wd)

for i in os.listdir(os.getcwd()):
        if i.startswith('delta_T_v3_no_halos'):

                filename = i
                
                #Now we start working with the boxes
                
                print ' ****************** THE BEGINING OF WORKING ON NEW FILE***********************'
                print ('The current file: '+filename)
                
                #Create the directory where images and and output data will be stored
                
                print '*************creating new directory***************'
                location = '/home/charlest/Desktop/work/pics/'

                if not os.path.exists(location):
                        os.makedirs(location)

                z = float(filename[21:-86])
                
                #We get the redshift from the filename string
                
                print ('Printing redshift just to make sure' + str('%06.2f'%z))
                
                
                #Now this is the part we bring in the other patch of the sky that does not contain the simulated cluster

                s = '/mnt/storage1/simulations/charles/data/hunds'
                os.chdir(s)
                for j in os.listdir(os.getcwd()):
                    if j.startswith('delta_T_v3_no_halos_z'+str('%06.2f'%(z+0.07))):
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

                                    fd = open(filename, 'rb')

                                    dI = (((2*(k*(f**2))/(c**2)))*(10**26)*(1.18*(10**-7)))*(np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape))

                                    bc = dI + Io_st
                                    print 'masking box with cluster'
                                    
                                    bc[118:138,118:138,118:138] = I_mod
                                    fd.close()
                                    return bc
                                box_cluster = box_with_cluster(256)

                                print 'done with box containg cluster'

                                print 'now working on box containg the background'

                                print 'working on Unperturbed CMB'

                                def Unperturbed_CMB(ndim):
                                    #create a onesarray
                                    shape = (ndim, ndim, ndim)
                                    Io_st = (2*(((k*To)**3)/((h*c)**2))*((p**3)/((np.exp(p))-1))*(10**26)*(1.18*(10**-7)))*(np.ones(shape))
                                    array = (np.ones(shape))
                                    Io_st_array = array * Io_st
                                    return Io_st_array
                                fa = Unperturbed_CMB(256)

                                print 'working on the background box'


                                #Take slice of the box
                                def readslice(ndim):
                                    shape = (ndim,ndim,ndim)
                                    fd = open(filename, 'rb')
                                    print filename
                                    dI = (((2*(k*(f**2))/(c**2)))*(10**26)*(1.18*(10**-7))*(np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape)))
                                    fd.close()
                                    return dI
                                fb =  readslice(256)

                                print ' now adding background with unperturbed CMB'

                                background = fa + fb

                                print 'Now performing differential procedure'

                                SZEI =  box_cluster - background

                                print 'Now getting you a good temperature box'
                                SZET = (SZEI*(10**(-26)))*(c**2)*(1/((2*k*(f**2))))*(1.18*(10**7))*(10**3)

                                print 'Calculating cluster temperature'
                                a =np.mean(SZET[118:138,118:138,118:138])

                                print 'Calculating foreground temperature'

                                b = np.mean(SZET[139:256,139:256,139:256])

                                print 'Outputting values to file'

                                file = open(location +'SZE-21cmSpec' +str(x) +'.txt', 'w')
                                file.write(str(z))
                                file.write(' '+str(v))
                                file.write(' '+str(x))
                                file.write(' '+str(a))
                                file.write(' '+str(b))
                                #file.write(' '+str(c))
                                #file.write(' '+str(d))
                                file.close()

                                #The rest is just imaging and naming
                                v1 = ['f = ', str(v) , ' MHz ', ', z = ', str(z), ' with signal_T = ', str(x), ' mK']
                                v2= ['f_=_', str(v) , '_MHz_', '_and_z_=_', str(z),'_signal_x=_', str(x),'.pdf']

                                print 'Outputting Simulation Box'

                                plt.subplot(2,2,1)
                                #plt.title('21cmFAST BOX slice at' +'/n' + ''.join(v) )
                                plt.title('21cmFAST BOX slice at ' + ''.join(v1) )
                                plt.xlabel('Mpc')
                                plt.ylabel('Mpc')
                                plt.imshow(box_cluster[118,:,:], cmap='jet')
                                plt.colorbar(format='%.e')


                                plt.subplot(2,2,2)
                                #plt.title('21cmFAST BOX slice with signal at' +''.join(v) )
                                plt.xlabel('Mpc')
                                plt.ylabel('Mpc')
                                plt.imshow(background[118,:,:], cmap='jet')
                                plt.colorbar(format='%.e')


                                plt.subplot(2,2,3)
                                #plt.title('21cmFAST BOX slice foreground and signal at' + ''.join(v) )
                                plt.xlabel('Mpc')
                                plt.ylabel('Mpc')
                                plt.imshow(SZEI[118,:,:], cmap='jet')
                                plt.colorbar(format='%.e')

                                plt.subplot(2,2,4)
                                #plt.title('Masked signal without foreground at'+ ''.join(v) )
                                plt.xlabel('Mpc')
                                plt.ylabel('Mpc')
                                plt.imshow(SZET[118,:,:], cmap='jet')

                                plt.colorbar()
                                plt.savefig(location + 'All_plots_at_' +''.join(v2), )
                                #plt.show()
                                plt.close()

                                break
                            else:
                                x=0
                                print 'Now I hope this is gonna work and we are done'




                                
