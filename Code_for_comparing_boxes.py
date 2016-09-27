import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import os
import csv

d = '/Volumes/CHARLES/from_z_=5_to_z=1100/Boxes'

os.chdir(d)



for i in os.listdir(os.getcwd()):
        if i.startswith('delta_'):
                filename = i
                e = '/Volumes/CHARLES/from_z_=5_to_z=1100/Boxes'
                os.chdir(e)
                for j in os.listdir(os.getcwd()):
                        if j.startswith('delta_'):
                                files = j
                                print ( 'I have file     ' +filename +'    and    ' +files)

                                z1 = float(filename[21:-85])
                                z2 = float(files[21:-85])
                                
                                if z1 == z2:
                                        print 'Hey      !!!!!!!!!!!!!! I got one              !!!!!!!!!!!!!!!!'
                                        
                                        def box_x(ndim):
                                                shape = (ndim,ndim,ndim)
                                                fd = open(filename, 'rb')
                                                print ('The one is     ' +filename)
                                                data = np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape)
                                                fd.close()
                                                return data
                                        fa = box_x(256)



                                        def box_y(ndim):
                                                shape = (ndim,ndim,ndim)
                                                fd = open(files, 'rb')
                                                print (' The second is      ' +files)
                                                data = np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape)
                                                fd.close()
                                                return data
                                        fb = box_y(256)

                                        ff = fa + fb
                                        plt.imshow(ff[0,:,:],cmap = cm.jet)
                                        plt.colorbar()
                                        #plt.savefig()
                                        plt.show()
                                        plt.close()

    
        
        
