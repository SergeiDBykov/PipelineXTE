#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:06:23 2020

@author: s.bykov
"""


from pipeline_core import *

ObsID='90089-11-03-02' # 90089-11-03-00G or 90089-11-03-01G 90089-11-03-02
xte_obs=ObservationXTE(ObsID)
os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')


binsize=5

'''

Channel 9:
        energy = 5.8421955+-0.20359278
    eff. area = 2646.0869228735896  cm^2
    energy width = 0.40718556
    1/(eff_area*en width) = 0.0009281187047645966


Channel 10:
        energy = 6.249672+-0.20388365
    eff. area = 2816.576150699146  cm^2
    energy width = 0.4077673
    1/(eff_area*en width) = 0.0008706951129767611


Channel 11:
        energy = 6.86219+-0.4086342
    eff. area = 4367.3066357014395  cm^2
    energy width = 0.8172684
    1/(eff_area*en width) = 0.0002801699853467043

Channel 12:
        energy = 7.680602+-0.40977788
    eff. area = 4491.555703678442  cm^2
    energy width = 0.81955576
    1/(eff_area*en width) = 0.00027165937860035115


Channel 13:
    energy = 8.501286+-0.4109063
    eff. area = 4563.682412291262  cm^2
    energy width = 0.8218126
    1/(eff_area*en width) = 0.00026663170529039713

'''



path_to_lc=f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/sa_data_lc_{binsize}sec'
os.chdir(f'products/sa_data_lc_{binsize}sec')

class TimeSeries():

    def __init__(self,lcname):
        self.fits=fits.open(path_to_lc+f'/{lcname}.lc_bary_orb_corr')
        self.time=self.fits[1].data['time']
        self.rate=self.fits[1].data['rate']
        self.error=self.fits[1].data['error']
        self.binsize=np.median(np.diff(self.time))
        self.fits.close()
    def divide(self,val):
        self.rate=self.rate/val
        self.error=self.error/val

#%%read lcurves and plot mean spe
en=np.array([5.8421955,6.249672,6.86219,7.680602,8.501286])
enerr=np.array([0.2034719,0.20388365,0.4086342,0.40977788,0.4109063])

factor=[0.0009281187047645966,0.0008706951129767611,0.0002801699853467043,
        0.00027165937860035115,0.00026663170529039713]

lc9=TimeSeries('ch9')
lc10=TimeSeries('ch10')
lc11=TimeSeries('ch11')
lc12=TimeSeries('ch12')
lc13=TimeSeries('ch13')

lclist=[lc9,lc10,lc11,lc12,lc13]

for k,lc in enumerate(lclist):
    lc=lc.divide(1/factor[k])


rate=np.array([lc.rate.mean() for lc in lclist])
rate_error=np.array([lc.rate.std() for lc in lclist])


plt.figure()

spe=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin_spe/conversion_factors/data_area_en.qdp',skip_header=3)



plt.errorbar(en,rate,rate_error,enerr,label='spectra_From_lc')
plt.plot(spe[:,0],spe[:,2],'k.')



a,b=np.polyfit(en[[0,3,4]],rate[[0,3,4]],1)

enaxis=np.linspace(en[0],en[-1],10)
plt.plot(enaxis,a*enaxis+b,label='linear best model ignoring 6-7 kev')

def best_line(x,N):
    return x*a+N




#plt.xscale('log')
plt.show()
plt.ylabel(' counts/s / keV')
plt.xlabel('energy')
plt.legend()






#%% plot a few points

N=lc10.rate.shape[0]

fig,ax=plt.subplots(3)
from random import randint
randi=[randint(0, N) for p in range(0, 3)]
#randi=[100,101,102]
for k,i in enumerate(randi):
    rate=np.array([x.rate[i] for x in lclist])
    rate_error=np.array([x.error[i] for x in lclist])
    ax[k].errorbar(en,rate,rate_error,enerr,label='spectra_From_lc')

    N_opt,N_opt_err=curve_fit(best_line,en[[0,3,4]],rate[[0,3,4]],
                              p0=b,sigma=rate_error[[0,3,4]],absolute_sigma=True)
    myline=lambda x: a*x+N_opt
    ax[k].plot(enaxis,myline(enaxis),'k:',alpha=0.3)



    iron_band_flux=np.average([rate[1],rate[2]],weights=[enerr[1],enerr[2]])
    iron_band_en=np.mean([en[1],en[2]])
    iron_band_flux_err=np.sqrt(rate_error[1]**2+rate_error[2]**2)

    ax[k].errorbar(iron_band_en,iron_band_flux,iron_band_flux_err,[en[2]-en[1]],color='r')

    fe_flux=iron_band_flux-myline(iron_band_en)
    fe_flux_err=iron_band_flux_err

    ax[k].vlines(iron_band_en,iron_band_flux,iron_band_flux-fe_flux,color='c',alpha=0.99)



plt.ylim(0.1,0.35)
plt.show()



#%% find residual flux

def find_fe_flux(lclist,frac=1):
    N=len(lclist[0].rate)

    diff=[]
    diff_err=[]
    print(f'{frac*100}% OF THE LINEAR APPROX IS USED AS CONTINUUM')

    for i in range(N):
        rate=np.array([x.rate[i] for x in lclist])
        rate_error=np.array([x.error[i] for x in lclist])


        N_opt,N_opt_err=curve_fit(best_line,en[[0,3,4]],rate[[0,3,4]],
                                  p0=b,sigma=rate_error[[0,3,4]],absolute_sigma=True)
        Norm_err=np.sqrt(np.diag(N_opt_err))[0]
        myline=lambda x: a*x+N_opt


        iron_band_flux=np.average([rate[1],rate[2]],weights=[enerr[1],enerr[2]])
        iron_band_en=np.mean([en[1],en[2]])

        iron_band_flux_err=np.sqrt(rate_error[1]**2+rate_error[2]**2)

        fe_flux=iron_band_flux-myline(iron_band_en)*frac
        fe_flux_err=np.sqrt(iron_band_flux_err**2+Norm_err**2)

        diff.append(fe_flux)
        diff_err.append(fe_flux_err)

    diff=np.asarray(diff)
    diff=np.reshape(diff,N)
    diff_err=np.asarray(diff_err)
    diff_err=np.reshape(diff_err,N)
    return diff,diff_err

frac=0.95
fe_flux,fe_flux_err=find_fe_flux(lclist,frac=frac)

figure()
plt.errorbar(lclist[0].time,fe_flux,fe_flux_err)
plt.xlabel('Time, s')
plt.ylabel('Iron line flux')
plt.show()

figure()
plt.hist(fe_flux,bins=100,label='Fe_flux')
plt.xlabel('Iron line flux ')
plt.show()



print(f'{frac*100}% OF THE LINEAR APPROX IS USED AS CONTINUUM')
print(f'Mean Flux: {np.mean(fe_flux)}')
print(f'Std Flux: {np.std(fe_flux)}')
print(f' mean/ std : {np.mean(fe_flux)/np.std(fe_flux)}')
print(f'Mean significance (mean flux/err): {np.mean(fe_flux/fe_flux_err)}')


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))
wmean,wstd=weighted_avg_and_std(fe_flux, fe_flux_err**(-2))
print(f'Weighted mean flux: {wmean}')
print(f'Weighted std flux: {wstd}')
print(f'Weighted mean/weighted std : {wmean/wstd}')


create_dir('fe_line')
os.system(f'cp ch11.lc_bary_orb_corr fe_line/python_lin_approx_fe_line_{frac}.lc_bary_orb_corr')
with fits.open(f'fe_line/python_lin_approx_fe_line_{frac}.lc_bary_orb_corr', mode='update') as hdul:
    hdul[1].data['rate']=fe_flux
    hdul[1].data['error']=fe_flux_err
    hdul.flush()  # changes are written back to original.fits


