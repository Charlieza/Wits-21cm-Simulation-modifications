import os
import re
import sys
import math
import string
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import mayavi.mlab as mlab
import time
import commands
import numpy
import Image
import numpy as np
import pyfits as pf
def readslice(ndim):
        shape = (ndim,ndim,ndim)
        fd = open("delta_T_v3_no_halos_z024.61_nf0.999804_useTs1_zetaX2.0e+56_alphaX1.2_TvirminX3.0e+04_aveTb-04.33_Pop2_256_300Mpc", 'rb')
        data = np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape)
        fd.close()
        return data
        
ff =  readslice(256)

ndim =256

use_color_flag=True
reverse_redshift_flag=True
print "i"
shape=(ndim,ndim,ndim)
print shape
mycolor='black-white'

if use_color_flag:
    mycolor='blue-red'
    mycolor='RdBu'
    #need to set color scale to be symmetric around zero
    
if reverse_redshift_flag:
    print 'want to be able to reseverse long axis ordering'
    ff=numpy.flipud(ff)
    #box_data=box_data[:,::-1,:]
    
print ff.min(), ff.max()
#clip_data
clip=30.0
vmax=clip
vmin= -clip
ff[ff>vmax]=vmax
ff[ff<vmin]=vmin
print ff.min(), ff.max()
c_black=(0.0,0.0,0.0)
c_white=(1.0,1.0,1.0)

#s=box_data
mlab.figure(1,bgcolor=c_white,fgcolor=c_black)
mlab.pipeline.volume(mlab.pipeline.scalar_field(ff))

mlab.draw()  #force update of figure
mlab.colorbar(orientation='vertical')
mlab.outline(color=c_black)

mlab.show()

