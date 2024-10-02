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
#a1, a2 = 64, 92
#b1, b2 = 28, 57
#bound = 4
#   
x=132 #+ 17#98/21,110/31
y=76#+ 20#
#
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
sto=fa['data'][:,:,:,:]#,:,:].reshape((52,220,endd-startt,4))
#
a1 = sto.shape[0]#/home/mbenko/hazel_work/Hazel_program/hazel_model/chromospheres/init_spot1.1dpe[0]
a2 = sto.shape[1]
print(sto.shape)
#sto *= clv
#r = h5py.File('/home/matobenko/Data/160723/dat_mskf.h5', 'r')

#msk = fa['mask'][:]
#msk[:,100::] = 0

#
stokes=np.ones((n_lambda,4))
stokes[:,0]=sto[y,x,startt:endd,0]#*msk[y,x]#
stokes[:,1]=sto[y,x,startt:endd,1]#*msk[y,x]
stokes[:,2]=sto[y,x,startt:endd,2]#*msk[y,x]
stokes[:,3]=sto[y,x,startt:endd,3]#*msk[y,x]
#
stokes *=clv
noise=np.ones((n_lambda,4))# 
noise[:,0]=noise[:,0]*np.std(stokes[:,0], ddof=1)# 
noise[:,1]=noise[:,1]*np.std(stokes[:,1], ddof=1)# 
noise[:,2]=noise[:,2]*np.std(stokes[:,2], ddof=1)# 
noise[:,3]=noise[:,3]*np.std(stokes[:,3], ddof=1)# 
#

#
stokes *= clv
#sto *= clv
noise *= clv 
#
max_stokes = np.max(stokes[0,0])
#
cont=sto[:,:,160,0]

fw = h5py.File('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/wav_f.h5', 'r')

#print(np.max(max_stokes))
wav = fw['wavelength'][0:440]
#
tmp = hazel.tools.File_observation(mode='single')
tmp.set_size(n_lambda=n_lambda, n_pixel=1)
#print(tmp.obs['stokes'].shape)
tmp.obs['wavelength'][:] = wav#fw['wavelength'][startt:endd]
tmp.obs['stokes'][0,:,:] = stokes
tmp.obs['sigma'][:] =  noise
tmp.obs['los'][0,:] = np.array([0.0,0.0,90.0])
tmp.obs['boundary'][:] = 0.0
tmp.obs['boundary'][0,:,0] = 1.0 *max_stokes#0.55
tmp.obs['weights'][:,0] = 12.0
tmp.obs['weights'][:,1] = 6.0
tmp.obs['weights'][:,2] = 6.0
tmp.obs['weights'][:,3] = 12.0
tmp.save('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/observations/10830_spot')
#
#tmp = hazel.tools.File_chromosphere(mode='single')
#tmp.set_default(n_pixel=1, default='disk')
#tmp.save('/home/mbenko/hazel_work/Hazel_program/hazel_model/chromospheres/init_spot')
#
mod = hazel.Model('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/configuration/conf_single_he.ini', working_mode='inversion', verbose=3)#, randomization=2)
mod.read_observation()
mod.open_output()
mod.invert()
mod.write_output()
mod.close_output()

f = h5py.File('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/output.h5', 'r')

#label = ['I', 'Q', 'U', 'V']
    #
    #pl.figure(figsize=(6, 4))
pl.figure(figsize=(15,10))
gs1 = gridspec.GridSpec(4, 2)#3,2
gs1.update(left=0.05, right=0.20, wspace=0.05)
ax1 = pl.subplot(gs1[1:3, :])
ax1.imshow(cont,cmap=pl.get_cmap('gray'), origin='lower')#imshow
#ax1.contour(msk,[0.5],colors='green')
ax1=pl.plot([x,x],[y-10,y+10], color='r') #y os
ax1=pl.plot([x-10,x+10],[y,y], color='r')
ax1=pl.grid(True)

#Bch = np.sqrt( f['ch1']['Bx'][0,0,1,0]**2+f['ch1']['By'][0,0,1,0]**2+f['ch1']['Bz'][0,0,1,0]**2)
#incch = 180*(np.arccos(f['ch1']['Bz'][0,0,1,0]/Bch))/np.pi
#azch=np.arctan2(f['ch1']['By'][0,0,1,0], f['ch1']['Bx'][0,0,1,0]) * 180 / np.pi
#Bph = np.sqrt(f['ph1']['Bx'][0,0,1,0]**2+f['ph1']['By'][0,0,1,0]**2+f['ph1']['Bz'][0,0,1,0]**2)
#incph = 180*(np.arccos(f['ph1']['Bz'][0,0,1,0]/Bph))/np.pi
#azph=np.arctan2(f['ph1']['By'][0,0,1,0], f['ph1']['Bx'][0,0,1,0]) * 180 / np.pi
#az=np.arctan2(by, bx) * 180 / np.pi

ax1=pl.text(-2,920,'Photosphere')
try:
    ax1=pl.text(-5,870,'Bx={:+.2f}'.format(f['ph1']['Bx'][0,0,0]))
    ax1=pl.text(-5,830,'By={:+.2f}'.format(f['ph1']['By'][0,0,0]))
    ax1=pl.text(-5,780,'Bz={:+.2f}'.format(f['ph1']['Bz'][0,0,0]))
    ax1=pl.text(-5,730,'v={:+.2f}'.format(f['ph1']['v'][0,0,0]))
    ax1=pl.text(-5,680,'vmic={:+.2f}'.format(f['ph1']['vmic'][0,0,0]))
except KeyError:
    pass

    #ax1 = pl.text(-2,470, '{:+.2f}'.format(Bph) )
#ax1 = pl.text(110,470, '{:+.2f}'.format(incph))
#ax1 = pl.text(170,470, '{:+.2f}'.format(azph))

ch = ['ch1','ch2','ch3']
pos_x = [-2,200,400]
par_ph = ['Bx', 'By', 'Bz', 'deltav', 'v', 'tau', 'beta']
y_end = 630
y_step = 50

for i in range(3):
    ax1 = pl.text(pos_x[i],y_end+y_step*0.5,ch[i])
    k=0
    try:
        for j in range(7):
            k+=1
            ax1=pl.text(pos_x[i],y_end-(k)*y_step,par_ph[j]+'={:+.2f}'.format(f[ch[i]][par_ph[j]][0,0,0]))
    except:
        pass
#
#
ax1=pl.text(-5,280,'x={:+.2f}'.format(x))
ax1=pl.text(200,180,'y={:+.2f}'.format(y))
ax1=pl.text(400,180,'chi2 ={:+.2f}'.format(f['spec1']['chi2'][0,0,0]))
y_int = 0.03

# Define gridspec for gs5
gs5 = gridspec.GridSpec(4, 2)
gs5.update(left=0.05, right=0.20, wspace=0.05)

# Subplot for gs5 (placing it below gs1)
ax5 = pl.subplot(gs5[3:, :])  # Now occupying the last row (3:), all columns
ax5.plot(np.arange(-3.6, 1.3, 0.1), f['ph1']['T'][0, 0, :])  # Plotting ph1['T'] data
ax5.set_xlabel(r'$\tau_{500}$', fontsize=12)  # x-axis label
ax5.set_ylabel(r'T [K]', fontsize=12)  # y-axis label

# Customize the x-axis ticks to have 6 evenly spaced ticks
ax5.set_xticks(np.linspace(-3.6, 1.2, 5))

#
#plot, I, Q, U, V
gs2 = gridspec.GridSpec(2, 3)
gs2.update(left=0.26, right=0.59, hspace=0.15)
gs3 = gridspec.GridSpec(2, 3)
gs3.update(left=0.665, right=0.995, hspace=0.15)
label_text=['I','Q','U','V']
y_int = 0.01
ylim_values = [(0, 1.2), (-y_int, y_int), (-y_int, y_int), (-y_int, y_int)]


for i in range(4):
    if i < 2:
        ax = pl.subplot(gs2[-2+i, :])
    else:
        ax = pl.subplot(gs3[-4+i, :])

    ax.plot(f['spec1']['wavelength'][:]-10830,stokes[:,i],':' ,label='observation') 
    ax.plot(f['spec1']['wavelength'][:]-10830,f['spec1']['stokes'][0,0,i,:],label='synthesize')

    ax.set_ylabel('{0}/Ic'.format(label_text[i]))
    ax.set_xlabel('Wavelength - 10830[$\AA$]')
    ax.set_ylim(ylim_values[i])
    
    if i ==0:
        ax.legend()
 
fig_name='invers{0}'.format(x) +'_{0}'.format(y) +'.png'
pl.savefig(fig_name, dpi=200) 

#print(max_stokes)
f.close()
#fa.close()
pl.show()

