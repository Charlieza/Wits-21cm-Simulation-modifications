import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os
import csv
import glob

#This code takes a box converts its contents into Io_mod
#Equation 10 of Colafrancesco, Emritte and Paolo

#This version picks out only redshifts in the spectrum file

k = 1.32*(10**(-23)) #Boltzman constant
h = 6.63*(10**(-34)) #Planck constant
c = 3.00*(10**(8)) #speed of light
To = 2.725

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
            print (str('%06.2f'%z))
            print ('redshift from filename is: ' + str(z))

            fi = float(1420/(1+z))
            f = float("{0:.2f}".format(fi))

            print ('frequency is:' + str(f))

            def Unperturbed(ndim):
                #create a onesarray
                shape = (ndim, ndim, ndim)
                Io_st = ((2*((k*To)**3))/((h*c)**2))*((((h*f)/(k*To))**3)/((np.exp(((h*f)/(k*To))))-1))
                array = np.ones(shape)
                Io_st_array = array * Io_st
                return Io_st_array
            fs = Unperturbed(256)

            #Take slice of the box

            def readslice(ndim):
                shape = (ndim,ndim,ndim)
                fd = open(filename, 'rb')
                print filename
                data = np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape)
                dI = ((2*k*(f**2))/(c**2)) * data
                fd.close()
                return dI
            ff =  readslice(256)


            def I_mod(ndim):
                    shape = (ndim,ndim,ndim)
                    fd1 = open(filename, 'rb')
                    print filename
                    data1 = np.fromfile(file=fd1, dtype= np.dtype('f4')).reshape(shape)
                    dI1 = ((2*k*(f**2))/(c**2)) * data1
                    Io_st1 = ((2*((k*To)**3))/((h*c)**2))*((((h*f)/(k*To))**3)/((np.exp(((h*f)/(k*To))))-1))
                    array1 = np.ones(shape)
                    Io_st_array1 = array1 * Io_st1
                    Io_mod= Io_st_array1 + dI1
                    Io_mod[118:138, 118:138, 118:138] = x
                    fd.close()
                    return I_mod

            Io_mod1 = fs + ff
            SZE = I_mod - Io_mod1


            #The rest is just imaging and naming

            v = ['f = ', str(f) , ' MHz ', ', z = ', str(z), ' with signal_T = ',' mK']
            v2= ['f_=_', str(f) , '_MHz_', '_and_z_=_', str(z),'_signal_x=_','.png']

            print 'Outputting Simulation Box'

            plt.subplot(2,2,1)

            plt.title('21cmFAST BOX slice at ' + ''.join(v) )
            plt.xlabel('Mpc')
            plt.ylabel('Mpc')
            plt.imshow(fs[118,:,:], cmap='jet')
            plt.colorbar()

            plt.subplot(2,2,2)
            plt.xlabel('Mpc')
            plt.ylabel('Mpc')
            plt.imshow(ff[118,:,:], cmap='jet')
            plt.colorbar()

            plt.subplot(2,2,3)
            plt.xlabel('Mpc')
            plt.ylabel('Mpc')
            plt.imshow(I0_mod1[125,:,:], cmap='jet')
            plt.colorbar()


            plt.subplot(2,2,4)
            #plt.title('Masked signal without foreground at'+ ''.join(v) )
            plt.xlabel('Mpc')
            plt.ylabel('Mpc')
            plt.imshow(SZE[118,:,:], cmap='jet')
            plt.colorbar()
            plt.savefig(location + 'All_plots_at_' +''.join(v2), bbox_inches='tight')
            #plt.show()
            plt.close()

            break
        print ' *****************Done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ******************'

                       
