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

#The comptonised spectrum



wd = '/Volumes/CHARLES/Boxes/theprog/hunds'
os.chdir(wd)

f = open('spectrum.csv', 'rb')
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

Lets work on the Slice and other signal
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

            for i in range (len(redshift)):
                #Working through redshifts
                if redshift[i] == z:
                    x = mod[i]
                    I_mod=((2*k*(f**2))/(c**2))*x
                    print ( 'x, the modification is:' +str(x))


                    def Unperturbed(ndim):
                        #create a onesarray
                        shape = (ndim, ndim, ndim)
                        Io_st = ((2*((k*To)**3))/((h*c)**2))*((((h*f)/(k*To))**3)/((np.exp(((h*f)/(k*To))))-1))
                        array = np.ones(shape)
                        Io_st_array = array * Io_st
                        return Io_st_array
                    fa = Unperturbed(256)


                    #Take slice of the box
                    def readslice(ndim):
                        shape = (ndim,ndim,ndim)
                        fd = open(filename, 'rb')
                        print filename
                        data = np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape)
                        dI = ((2*k*(f**2))/(c**2)) * data
                        fd.close()
                        return dI
                    fb =  readslice(256)

                    def Imod(ndim):
                        shape = (ndim,ndim,ndim)
                        a=np.zeros(shape)
                        a[118:138,118:138,118:138]=I_mod
                        return a
                    fc = Imod(256)

                    SZE = fc - (fa + fb)

                    a =np.mean(SZE[118:138,118:138,118:138])
                    b = ((c**2)/(2*k*(f**2)))*(10**3)*(a)

                    file = open(location +'log' +str(x) +'.txt', 'w')
                    file.write(str(z))
                    file.write(' '+str(f))
                    file.write(' '+str(a))
                    file.write(' '+str(b))
                    file.write(' '+str(x))
                    #file.write(' '+str(c))
                    #file.write(' '+str(d))
                    file.close()


                    #The rest is just imaging and naming


                    v = ['f = ', str(f) , ' MHz ', ', z = ', str(z), ' with signal_T = ', str(x), ' mK']
                    v2= ['f_=_', str(f) , '_MHz_', '_and_z_=_', str(z),'_signal_x=_', str(x),'.pdf']


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
                else:
                    x=0


                print 'Now I hope this is gonna work and we are done'

