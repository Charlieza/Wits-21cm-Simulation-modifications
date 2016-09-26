import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os
import csv

#This code takes a spectrum file and making a signal out of it at various
#redshifts and adds the signal to a delta_T output box from 21cmFAST

#This version picks out only redshifts in the spectrum files

###### Lets start by bringing in the spectrum #######


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

wd = '../Boxes'
os.chdir(wd)

#Lets work on the Slice and other signal
for i in os.listdir(os.getcwd()):
        if i.startswith('delta_T_v3_no_halos'):

                filename = i
                print ' ****************** THE BEGINING OF WORKING ON NEW FILE***********************'
                #print filename
                print '*************creating new directory***************'
                location = '../only_'+ os.path.basename(os.path.normpath(os.getcwd()))+ 'pictures'

                if not os.path.exists(location):
                        os.makedirs(location)
                z = float(filename[21:-85])

                print ('redshift from filename is: ' + str(z))
                fi = float(1420/(1+z))
                f = float("{0:.2f}".format(fi))
                print ('frequency is:' + str(f))

                #Lets pick only the boxes at the redshifts defined in the spectrum file

                for i in range (len(redshift)):
                    #Working through redshifts
                    if redshift[i] == z:
                        x = mod[i]

                        print ( 'x, the modification is:' +str(x))

                        w = int(round(z))
                        g = float("{0:.2f}".format(w))
                        
                        print w
                        print g


                        #Lets build the signal to be added

                        

                        def another(ndim):
                            #create a zeroarray
                            shape = (ndim, ndim, ndim)
                            array = np.zeros(shape)

                            #choose region to mask, and the value to mask it by

                            array[100:150, 100:150, 100:150] = x
                            return array
                        fs = another(256)

                        #Take slice of the box

                        def readslice(ndim):
                            shape = (ndim,ndim,ndim)
                            fd = open(filename, 'rb')
                            print filename
                            data = np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape)
                            fd.close()
                            return data
                        ff =  readslice(256)

                        def whole(ndim):
                            shape = (ndim, ndim, ndim)
                            s = '../Boxes'
                            os.chdir(s)
                            for k in os.listdir(os.getcwd()):
                                if k.startswith('delta_T_v3_no_halos'+str(g)):
                                    myfile = k
                                    return myfile
                                print myfile
                            print myfile
                            fd = open(myfile, 'rb')
                            data = np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape)
                            fd.close()
                            return data
                        fg = whole(256)
                        

                        #Add the signal to the Box slice

                        fz = ff + fs

                        fk = ff - fg + fs

                        #The rest is just imaging
                        
                        v = ['f = ', str(f) , ' MHz ', ' and z = ', str(z), ' signal_T = ', str(x), ' mK']
                        v2= ['f_= ', str(f) , ' MHz ', ' and z = ', str(z),' signal_x= ', str(x),'.png']
                        print 'Outputting Box signal only'
                        plt.title('21cmFAST BOX slice at ' + ''.join(v) )
                        plt.xlabel('Mpc')
                        plt.ylabel('Mpc')
                        plt.imshow(ff[0,:,:],cmap = cm.jet)
                        plt.colorbar()
                        plt.savefig(location +'Origbox_at_' +''.join(v2), bbox_inches='tight')
                        #plt.show()
                        plt.close()
                        
                        print 'Outputting Just masked region'
                        plt.title('masked array signal at ' + ''.join(v))
                        plt.xlabel('Mpc')
                        plt.ylabel('Mpc')
                        plt.imshow(fs[100,:,:],cmap = cm.jet)
                        plt.colorbar()
                        plt.savefig(location + 'masked_array_at'+''.join(v2), bbox_inches='tight')
                        #plt.show()
                        plt.close()
                        
                        print 'Outputting Box with signal and masked region '
                        plt.title('21cmFAST Box Slice with added signal at ' + ''.join(v))
                        plt.xlabel('Mpc')
                        plt.ylabel('Mpc')
                        plt.imshow(fz[100,:,:],cmap = cm.jet)
                        plt.colorbar()
                        plt.savefig(location + 'box_plus_masked_at' +''.join(v2), bbox_inches='tight')
                        #plt.show()
                        plt.close()


                        print 'Outputting Box with signal and foreground '
                        plt.title('21cmFAST Box Slice foreground and signal at ' + ''.join(v))
                        plt.xlabel('Mpc')
                        plt.ylabel('Mpc')
                        plt.imshow(fk[100,:,:],cmap = cm.jet)
                        plt.colorbar()
                        plt.savefig(location + 'forground_box' +''.join(v2), bbox_inches='tight')
                        #plt.show()
                        plt.close()

                        break
                    else:
                        x=0
                        
                        
                print ( 'x, the modification is:' +str(x))





                print ' *****************Done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ******************'
