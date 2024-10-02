import numpy as np
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import hazel
import h5py
import scipy.io
import math
#from ipdb import set_trace as stop
#
fa = h5py.File('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/dat_f.h5','r')#('/home/mbenko/data2/data/data2_gris_deg.h5', 'r')
#
msk = fa['mask'][:]
msk[:,100::] = 0
n_cycle=3
#
startt=0#25#115
endd=440
n_lambda=endd-startt
#
theta=55.535#65.9521
mu = np.cos(theta*np.pi/180)
clv = hazel.util.i0_allen(10830,mu) / hazel.util.i0_allen(10830,1.0)
#print("CLV = {0}".format(clv))
#
sto=fa['data'][:,:,0:440,:]#.reshape(52*220,234,4)#fa['stokes'][:,:,:].reshape((50,232,fa['stokes'].shape[1],4))

fw = h5py.File('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/wav_f.h5', 'r')

#print(np.max(max_stokes))
wav = fw['wavelength'][0:440]

np.savetxt('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/observations/10830_spot.wavelength', wav, header='lambda')

fa.close()

f = open('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/observations/10830_spot.weights', 'w')
f.write('# WeightI WeightQ WeightU WeightV\n')
for i in range(n_lambda):
    f.write('12.0    6.0   6.0   12.0\n')
f.close()

n_pixel = 105*440#11600

a = sto.reshape(105*440,440,4)
stokes = a * clv


stokes_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)
sigma_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)
los_3d = np.zeros((n_pixel,3), dtype=np.float64)
boundary_3d = np.zeros((n_pixel,n_lambda,4), dtype=np.float64)

noise=np.ones((n_pixel, n_lambda,4))#
boundary = np.array([1.0,0.0,0.0,0.0])
for i in range(n_pixel):

    noise[i,:,0]=noise[i,:,0]*np.std(stokes[i,:,0], ddof=1)# 
    noise[i,:,1]=noise[i,:,1]*np.std(stokes[i,:,1], ddof=1)# 
    noise[i,:,2]=noise[i,:,2]*np.std(stokes[i,:,2], ddof=1)# 
    noise[i,:,3]=noise[i,:,3]*np.std(stokes[i,:,3], ddof=1)# 

#    noise = np.std(stokes[i,0:20,1])

    stokes_3d[i,:,:] = stokes[i,:,:]
    sigma_3d[i,:,:] = noise[i,...]#*np.ones((220,4))
    los_3d[i,:] = np.array([0.0,0.0,90.0])
    boundary_3d[i,:,:] = np.repeat(np.atleast_2d(boundary), n_lambda, axis=0)

print(noise.shape)

f = h5py.File('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/observations/10830_spot_stokes.h5', 'w')
db_stokes = f.create_dataset('stokes', stokes_3d.shape, dtype=np.float64)
db_sigma = f.create_dataset('sigma', sigma_3d.shape, dtype=np.float64)
db_los = f.create_dataset('LOS', los_3d.shape, dtype=np.float64)
db_boundary = f.create_dataset('boundary', boundary_3d.shape, dtype=np.float64)
db_stokes[:] = stokes_3d
db_sigma[:] = sigma_3d
db_los[:] = los_3d
db_boundary[:] = boundary_3d
f.close()


#r = h5py.File('/home/matobenko/Data/160723/dat_mskf.h5', 'r') 

#msk = r['mask1'][:]
#mask = np.zeros((n_pixel,), dtype=np.int8)
#b = np.average(abs(stokes_3d[:,:,3]), axis=1)
#b = b.reshape(50*232,4)
#mask[:] = msk.reshape(105*440)#1*(b > 0.002)
#mask[0:15] = 1
#fa = h5py.File('/home/matobenko/Data/160723/inversion_he01/observations/10830_mask.h5', 'w')
#db_mask = fa.create_dataset('mask', mask.shape, dtype=np.int8)
#db_mask[:] = mask
fa.close()



#h5py.File('/home/matobenko/Data/160723/dat_mskf.h5', 'r') 

