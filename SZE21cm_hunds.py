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





f = open('paolocmbcompt_21cmdt6_10kev_1em3.csv', 'rb')
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

wd = '/Volumes/CHARLES/Boxes/theprog/hunds'
os.chdir(wd)

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
            v = float("{0:.2f}".format(fi))
            f = (v*(10**6))

            print ('frequency is:' + str(f))

            for i in range (len(redshift)):
                #Working through redshifts
                if redshift[i] == z:
                    x = mod[i]
                    I_mod= x
                    p = ((h*f)/(k*To))
                    print ( 'x, the modification is:' +str(x))


                    def Unperturbed(ndim):
                        #create a onesarray
                        shape = (ndim, ndim, ndim)
                        Io_st = ((2*(((k*To)**3)/((h*c)**2))*((p**3)/((np.exp(p))-1)))*(1/(10**(-26))))
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
                        dI = ((2*(k*(f**2)))/(c**2))*(1/(10**(-26)))*data
                        fd.close()
                        return dI
                    fb =  readslice(256)

                    def Imod(ndim):
                        shape = (ndim,ndim,ndim)
                        a=np.zeros(shape)
                        a[118:138,118:138,118:138]=I_mod
                        return a
                    fc = Imod(256)

                    SZEI = fc - (fa + fb)
                    SZET = SZEI*(10**(-26))

                    a =np.mean(SZET[118:138,118:138,118:138])

                    file = open(location +'SZEspec' +str(x) +'.txt', 'w')
                    file.write(str(z))
                    file.write(' '+str(f))
                    file.write(' '+str(a))
                    file.write(' '+str(x))
                    #file.write(' '+str(x))
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
                    plt.imshow(SZET[118,:,:], cmap='jet')

                    plt.colorbar()
                    plt.savefig(location + 'All_plots_at_' +''.join(v2), bbox_inches='tight')
                    #plt.show()
                    plt.close()

                    break
                else:
                    x=0


                print 'Now I hope this is gonna work and we are done'

