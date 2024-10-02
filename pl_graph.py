import numpy as np
import matplotlib.pyplot as pl
import matplotlib.gridspec as gridspec
import hazel
import h5py
import scipy.io
import math

f = h5py.File('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/output3ph2.h5', 'r')
fa = h5py.File('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/observations/10830_spot_stokes.h5','r')
fm= h5py.File('/mnt/RAIDofsF/mbenko/results/160823/inversion_he01/observations/10830_masksp.h5', 'r')

label=['I','Q','U','V']
print('(npix,nrand,ncycle,nstokes,nlambda) -> {0}'.format(f['spec1']['stokes'].shape))

import matplotlib.pyplot as plt
cont = fa['stokes'][:,0,0].reshape(105,440)

for t in range(46000):
    print(t)
    if fm['mask'][t] == 1:  # Assuming fm['mask'] is indexed by k
       # fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
       # ax = ax.flatten()

        pl.figure(figsize=(15,10))
        gs1 = gridspec.GridSpec(4, 2)#3,2
        gs1.update(left=0.05, right=0.20, wspace=0.05)
        ax1 = pl.subplot(gs1[1:3, :])
        ax1.imshow(cont,cmap=pl.get_cmap('gray'), origin='lower')#imshow
        x = t%440
        y = t//440
        ax1=pl.plot([x,x],[y-10,y+10], color='r') #y os
        ax1=pl.plot([x-10,x+10],[y,y], color='r')
        ax1=pl.grid(True)

        ax1=pl.text(-2,920,'Photosphere')
        try:
            ax1=pl.text(-5,870,'Bx={:+.2f}'.format(f['ph1']['Bx'][t,0,0]))
            ax1=pl.text(-5,830,'By={:+.2f}'.format(f['ph1']['By'][t,0,0]))
            ax1=pl.text(-5,780,'Bz={:+.2f}'.format(f['ph1']['Bz'][t,0,0]))
            ax1=pl.text(-5,730,'v={:+.2f}'.format(f['ph1']['v'][t,0,0]))
            ax1=pl.text(-5,680,'vmic={:+.2f}'.format(f['ph1']['vmic'][t,0,0]))
        except KeyError:
            pass


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
                    ax1=pl.text(pos_x[i],y_end-(k)*y_step,par_ph[j]+'={:+.2f}'.format(f[ch[i]][par_ph[j]][t,0,0]))
            except:
                pass
    #
    #
        ax1=pl.text(-5,180,'x={:+.2f}'.format(x))
        ax1=pl.text(200,180,'y={:+.2f}'.format(y))
        ax1=pl.text(400,180,'chi2 ={:+.2f}'.format(f['spec1']['chi2'][t,0,0]))
        y_int = 0.03

# Define gridspec for gs5
    #    gs5 = gridspec.GridSpec(4, 2)
    #    gs5.update(left=0.05, right=0.20, wspace=0.05)

        # Subplot for gs5 (placing it below gs1)
    #    ax5 = pl.subplot(gs5[3:, :])  # Now occupying the last row (3:), all columns
    #    ax5.plot(np.arange(-3.6, 1.3, 0.1), f['ph1']['T'][t, 0, :])  # Plotting ph1['T'] data
    #    ax5.set_xlabel(r'$\tau_{500}$', fontsize=12)  # x-axis label
    #    ax5.set_ylabel(r'T [K]', fontsize=12)  # y-axis label

# Customize the x-axis ticks to have 6 evenly spaced ticks
    #    ax5.set_xticks(np.linspace(-3.6, 1.2, 5))

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

            ax.plot(f['spec1']['wavelength'][:]-10830,fa['stokes'][t,:,i],':' ,label='observation') 
            ax.plot(f['spec1']['wavelength'][:]-10830,f['spec1']['stokes'][t,0,i,:],label='synthesize')

            ax.set_ylabel('{0}/Ic'.format(label_text[i]))
            ax.set_xlabel('Wavelength - 10830[$\AA$]')
            ax.set_ylim(ylim_values[i])
    
            if i ==0:
                ax.legend()
        
        fig_name='invers{0}'.format(x) +'_{0}'.format(y) +'.png'
        pl.savefig(fig_name, dpi=200) 





        #for i in range(4):
        #    ax[i].plot(f['spec1']['wavelength'][:] - 10830, fa['stokes'][k, :, i])
        #    ax[i].plot(f['spec1']['wavelength'][:] - 10830, f['spec1']['stokes'][k, 0, i, :])

        #    ax[i].set_xlabel('Wavelength - 10830 [$\AA$]')
        #    ax[i].set_ylabel('{0}/Ic'.format(label[i]))

        #plt.tight_layout()
        #fig_name='invers{0}'.format(k)+'.png'
        #pl.savefig(fig_name, dpi=200)
        #plt.show()
    
f.close()
fa.close()
fm.close()
