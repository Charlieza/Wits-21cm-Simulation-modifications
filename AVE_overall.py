import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os
import csv

#This code takes a spectrum file and making a signal out of it at various
#redshifts and adds the signal to a delta_T output box from 21cmFAST

#This version picks out only redshifts in the spectrum files

#This has the foreground included#

###### Lets start by bringing in the spectrum #######


f = open('cmbcompt_21cmdt5_10kev_1em3.csv', 'rb')
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

wd = '/mnt/storage1/simulations/charles/21cmFAST/256simulation/Boxes'
os.chdir(wd)

#Lets work on the Slice and other signal
for i in os.listdir(os.getcwd()):
        if i.startswith('delta_T_v3_no_halos'):

                filename = i
                print ' ****************** THE BEGINING OF WORKING ON NEW FILE***********************'
                print ('The current file: '+filename)
                print '*************creating new directory***************'
                location = '/home/charlest/Desktop/multipictures/'

                if not os.path.exists(location):
                    os.makedirs(location)


                    
                    z = float(filename[21:-85])
                    print (str('%06.2f'%z))


                    print ('redshift from filename is: ' + str('%06.2f'%z))
                    fi = float(1420/(1+z))
                    f = float("{0:.2f}".format(fi))
                    print ('frequency is:' + str('%06.2f'%f))
                    #Lets pick only the boxes at the redshifts defined in the spectrum file


                    for i in range (len(redshift)):
                        #Working through redshifts
                        if redshift[i] == z:
                            x = mod[i]
                            print ( 'x, the modification is:' +str('%06.2f'%x))



                            #Lets build the signal to be added
                            def another(ndim):
                                #create a zeroarray
                                shape = (ndim, ndim, ndim)
                                array = np.zeros(shape)
                                #choose region to mask, and the value to mask it by
                                array[0:256, 0:256, 0:256] = x
                                return array
                            fa = another(256)
                            #Take slice of the box

                            def readslice(ndim):
                                shape = (ndim,ndim,ndim)
                                fd = open(filename, 'rb')
                                print filename
                                data = np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape)
                                fd.close()
                                return data
                            fb =  readslice(256)

                            fc = fa -fb 


                            a = np.mean(ff)
                            b = np.std(ff)
                            c = np.mean(fc)
                            d = np.std(fc)

                            print ('ave is    '+str('%06.2f'%c))
                            print ('std is  ' +str('%06.2f'%d))

                            file = open(location + 'temps' +str('%06.2f'%x)+',txt', 'w')
                            file.write(str('%06.2f'%z))
                            file.write(' '+str('%06.2f'%f))
                            file.write(' '+str('%06.2f'%x))
                            file.write(' '+str('%06.2f'%a))
                            file.write(' '+str('%06.2f'%b))
                            file.write(' '+str('%06.2f'%c))
                            file.write(' '+str('%06.2f'%d))
                            file.close()


                            #The rest is just imaging and naming
                            v = ['f = ', str(f) , ' MHz ', ', z = ', str(z), ' with signal_T = ', str(x), ' mK']
                            v2= ['f_=_', str(f) , '_MHz_', '_and_z_=_', str(z),'_signal_x=_', str(x),'.png']


                            print 'Outputting Simulation Box'

                            plt.subplot(2,2,1)
                            #plt.title('21cmFAST BOX slice at' +'/n' + ''.join(v) )
                            plt.title('21cmFAST BOX slice at ' + ''.join(v) )
                            plt.xlabel('Mpc')
                            plt.ylabel('Mpc')
                            plt.imshow(fa[118,:,:], cmap='jet')
                            plt.colorbar()


                            plt.subplot(2,2,2)
                            #plt.title('21cmFAST BOX slice with signal at' +''.join(v) )
                            plt.xlabel('Mpc')
                            plt.ylabel('Mpc')
                            plt.imshow(fb[118,:,:], cmap='jet')
                            plt.colorbar()

                            plt.subplot(2,2,3)
                            #plt.title('21cmFAST BOX slice foreground and signal at' + ''.join(v) )
                            plt.xlabel('Mpc')
                            plt.ylabel('Mpc')
                            plt.imshow(fc[118,:,:], cmap='jet')
                            plt.colorbar()

                            break
                        else:
                            x=0

                        print ( 'x, the modification is:' +str(x))
                        print ' *****************Done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ******************'
