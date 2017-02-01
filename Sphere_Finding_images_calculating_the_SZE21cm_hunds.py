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


va,vb,vc,vd,ve,vf,vg,vh,vi = randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50)
vj,vk,vl,vm,vn,vo,vp,vq,vr = randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50)
vs,vt,vu,vv,vw,vx,vy,vz,xa = randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50)
xb,xc,xd,xe,xf,xg,xh,xi,xj = randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50)
xk,xl,xm,xn,xo,xp,xq,xr,xs = randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50)
xt,xu,xv,xw,xx,xy,xz= randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50),randint(0,50)
vvv = 25
print (va,vb,vc,vd,ve,vf,vg,vh,vi,vj,vk,vl,vm,vn,vo,vp,vq,vr,vs,vt,vu,vv,vw,vx,vy,vz,xa,xb,xc,xd,xe,xf,xg,xh,xi,xj,xk,xl,xm,xn,xo,xp,xq,xr,xs,xt,xu,xv,xw,xx,xy,xz)

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
                                    
                                    print 'masking box with cluster'
                                    
                                    def sector_mask(shape,centre,radius,angle1_range,angle2_range):
                                        """ Return a boolean mask for a circular sector. The start/stop angles in`angle_range` should be given in clockwise order."""
                                        x,y,z = np.ogrid[0:50,0:50,0:50]
                                        cx,cy,cz = centre
                                        tmin,tmax = np.deg2rad(angle1_range)
                                        smin,smax = np.deg2rad(angle2_range)
                                        # ensure stop angle > start angle
                                        if tmax < tmin:
                                            tmax += 2*np.pi

                                        if smax < smin:
                                            smax += np.pi

                                        # convert cartesian --> polar coordinates
                                        r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz)
                                        phi = np.arccos((z-cz)/(np.sqrt(r2))) - smin
                                        theta = np.arcsin((y-cy)/((np.sqrt(r2))*np.sin(phi))) - tmin

                                        # wrap angles between 0 and 2*pi
                                        phi %= (np.pi)
                                        theta %= (2*np.pi)

                                        # circular mask

                                        circmask = r2 <= radius*radius

                                        # angular mask
                                        anglemask1 = theta <= (tmax-tmin)
                                        anglemask2 = phi <= (smax-smin)
                                        return circmask*anglemask1*anglemask2

                                    shape =(50,50,50)

                                    matrix = bc
                                    mask = sector_mask(matrix.shape,(va,vb,vc),5,(0,360),(0,180))
                                    mask1 = sector_mask(matrix.shape,(vd,ve,vf),5,(0,360),(0,180))
                                    mask2 = sector_mask(matrix.shape,(vg,vh,vi),5,(0,360),(0,180))
                                    mask3 = sector_mask(matrix.shape,(vj,vk,vl),5,(0,360),(0,180))
                                    mask4 = sector_mask(matrix.shape,(vm,vn,vo),5,(0,360),(0,180))
                                    mask5 = sector_mask(matrix.shape,(vp,vq,vr),5,(0,360),(0,180))
                                    mask6 = sector_mask(matrix.shape,(vs,vt,vu),5,(0,360),(0,180))
                                    mask7 = sector_mask(matrix.shape,(vv,vw,vx),5,(0,360),(0,180))
                                    mask8 = sector_mask(matrix.shape,(vy,vz,xa),5,(0,360),(0,180))
                                    mask9 = sector_mask(matrix.shape,(xb,xc,xd),5,(0,360),(0,180))
                                    mask10 = sector_mask(matrix.shape,(xe,xf,xg),5,(0,360),(0,180))
                                    mask11 = sector_mask(matrix.shape,(xh,xi,xj),5,(0,360),(0,180))
                                    mask12 = sector_mask(matrix.shape,(xk,xl,xm),5,(0,360),(0,180))
                                    mask13 = sector_mask(matrix.shape,(xn,xo,xp),5,(0,360),(0,180))
                                    mask14 = sector_mask(matrix.shape,(xq,xr,xs),5,(0,360),(0,180))
                                    mask15 = sector_mask(matrix.shape,(xt,xu,xv),5,(0,360),(0,180))
                                    mask16 = sector_mask(matrix.shape,(xw,xx,xy),5,(0,360),(0,180))
                                    mask17 = sector_mask(matrix.shape,(vvv,vvv,vvv),5,(0,360),(0,180))
                                    matrix[mask] = I_mod
                                    matrix[mask1] = I_mod
                                    matrix[mask2] = I_mod
                                    matrix[mask3] = I_mod
                                    matrix[mask4] = I_mod
                                    matrix[mask5] = I_mod
                                    matrix[mask6] = I_mod
                                    matrix[mask7] = I_mod
                                    matrix[mask8] = I_mod
                                    matrix[mask9] = I_mod
                                    matrix[mask10] = I_mod
                                    matrix[mask11] = I_mod
                                    matrix[mask12] = I_mod
                                    matrix[mask13] = I_mod
                                    matrix[mask14] = I_mod
                                    matrix[mask15] = I_mod
                                    matrix[mask16] = I_mod
                                    matrix[mask17] = I_mod
                                    fd.close()
                                    return matrix
                                
                                box_cluster = box_with_cluster(50)

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
                                fa = Unperturbed_CMB(50)

                                print 'working on the background box'


                                #Take slice of the box
                                def readslice(ndim):
                                    shape = (ndim,ndim,ndim)
                                    fd = open('/home/charles/Documents/simulation_test_folder/data2/hunds/'+myfile, 'rb')
                                    print myfile
                                    dI = (((2*(k*(f**2))/(c**2)))*(10**-3)*(10**26)*(1.18*(10**-7))*(np.fromfile(file=fd, dtype= np.dtype('f4')).reshape(shape)))
                                    fd.close()
                                    return dI
                                fb =  readslice(50)

                                print ' now adding background with unperturbed CMB'

                                background = (fa + fb)

                                print 'Now performing differential procedure'

                                SZEI =  box_cluster - background

                                print 'Now getting you a good temperature box'
                                SZET = ((SZEI*(10**(-26)))*(c**2)*(1/((2*k*(f**2))))*(1.18*(10**7)))*(10**3)

                                print 'Calculating cluster temperature'
                                a =np.mean(SZET[(vvv+1):(vvv+2),(vvv+1):(vvv+2),(vvv+1):(vvv+2)])

                                print 'Calculating foreground temperature'

                                b = np.mean(SZET[xq:(xq+2),xr:(xr+2),xs:(xs+2)])

                                print 'Outputting values to file'

                                file = open(location +'SZE-21cmSpec' +str(x) +'.txt', 'w')
                                file.write(str(z))
                                file.write(' '+str(v))
                                file.write(' '+str(x))
                                file.write(' '+str(a))
                                file.write(' '+str(b))
                                #file.write(' '+str(vvv))
                                #file.write(' '+str(vn))
                                #file.write(' '+str(vo))
                                file.close()

                                #The rest is just imaging and naming
                                v1 = ['f = ', str(v) , ' MHz ', ', z = ', str(z), ' with signal_T = ', str(x), ' mK']
                                v2= ['f_=_', str(v) , '_MHz_', '_and_z_=_', str(z),'_signal_x=_', str(x),'.png']

                                print 'Outputting Simulation Box'

                                #plt.subplot(2,2,1)
                                plt.title('Background with Cluster')
                                plt.xlabel('Mpc')
                                plt.ylabel('Mpc')
                                plt.imshow(box_cluster[(vvv+1),:,:], cmap='plasma')
                                c_bar = plt.colorbar(format='%.e')
                                c_bar.set_label(r'${\rm \delta I [\mathrm{Jy}]}$', fontsize=24, rotation=-90, labelpad=32)
                                plt.savefig(location + 'cluster+background/plot_at_' +''.join(v2), )
                                plt.close()



                                #plt.subplot(2,2,2)
                                plt.title('Background')
                                plt.xlabel('Mpc')
                                plt.ylabel('Mpc')
                                plt.imshow(background[(vvv+1),:,:], cmap='plasma')
                                c_bar = plt.colorbar(format='%.e')
                                c_bar.set_label(r'${\rm \delta I [\mathrm{Jy}]}$', fontsize=24, rotation=-90, labelpad=32)
                                plt.savefig(location + 'background/plot_at_' +''.join(v2), )
                                plt.close()



                                #plt.subplot(2,2,3)
                                plt.title('The difference')
                                plt.xlabel('Mpc')
                                plt.ylabel('Mpc')
                                plt.imshow(SZEI[(vvv+1),:,:], cmap='plasma')
                                c_bar = plt.colorbar(format='%.e')
                                c_bar.set_label(r'${\rm \delta I [\mathrm{Jy}]}$', fontsize=24, rotation=-90, labelpad=32)
                                plt.savefig(location + 'diff/plot_at_' +''.join(v2), )
                                plt.close()

                                
                                #plt.subplot(2,2,4)
                                plt.title('The Difference in temperature values T(mK)')
                                plt.xlabel('Mpc')
                                plt.ylabel('Mpc')
                                plt.imshow(SZET[(vvv+1),:,:], cmap='plasma')
                                c_bar = plt.colorbar()
                                c_bar.set_label(r'${\rm \delta I [\mathrm{mK}]}$', fontsize=24, rotation=-90, labelpad=32)
                                plt.savefig(location + 'tempdiff/plot_at_' +''.join(v2), )
                                #plt.show()
                                plt.close()

                                break
                            else:
                                x=0
                                print 'Now I hope this is gonna work and we are done'




                                
