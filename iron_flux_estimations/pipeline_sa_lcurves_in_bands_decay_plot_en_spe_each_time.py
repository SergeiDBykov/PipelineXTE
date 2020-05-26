#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:06:23 2020

@author: s.bykov
"""


from pipeline_core import *

ObsID='90427-01-03-02' #  90427-01-03-00 90427-01-03-02
xte_obs=ObservationXTE(ObsID)
os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')


binsize=1

'''

Channel 9:
    energy = 5.4322805+-0.20317984
    eff. area = 4091.5338513619663  cm^2
    energy width = 0.40635968
    1/(eff_area*en width) = 0.0006014551413650553


Channel 10:
    energy = 5.838932+-0.2034719
    eff. area = 2645.4159092628215  cm^2
    energy width = 0.4069438
    1/(eff_area*en width) = 0.0009289056464977036



Channel 11:
    energy = 6.246167+-0.20376277
    eff. area = 2816.1715476717145  cm^2
    energy width = 0.40752554
    1/(eff_area*en width) = 0.0008713368112906551


Channel 12:
        energy = 6.653981+-0.20405173
    eff. area = 4309.658590127702  cm^2
    energy width = 0.40810346
    1/(eff_area*en width) = 0.0005685738168610067

Channel 13:
        energy = 7.062372+-0.20433903
    eff. area = 4424.693113729067  cm^2
    energy width = 0.40867806
    1/(eff_area*en width) = 0.0005530132315100311


Channel 14:
        energy = 7.4713354+-0.20462418
    eff. area = 4477.967104653102  cm^2
    energy width = 0.40924836
    1/(eff_area*en width) = 0.0005456726080212225


Channel 15:
        energy = 7.8808675+-0.2049079
    eff. area = 4505.861407994961  cm^2
    energy width = 0.4098158
    1/(eff_area*en width) = 0.0005415436557273444


Channel 16:
        energy = 8.290965+-0.2051897
    eff. area = 4548.214051969805  cm^2
    energy width = 0.4103794
    1/(eff_area*en width) = 0.000535764029380296

'''



path_to_lc=f'/Users/s.bykov/work/xray_pulsars/rxte/results/out90427-01-03-02/products/sa_data_lc_{binsize}sec'
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
en=np.array([5.4322805,
             5.838932,
             6.246167,
             6.653981,
             7.4713354,
             7.8808675
])


enerr=np.array([0.203]*len(en))

factor=np.array([0.0006014551413650553,
                 0.0009289056464977036,
                 0.0008713368112906551,
                 0.0005685738168610067,
                 0.0005456726080212225,
                 0.0005415436557273444
])

lc9=TimeSeries('ch9')
lc10=TimeSeries('ch10')
lc11=TimeSeries('ch11')
lc12=TimeSeries('ch12')
#lc13=TimeSeries('ch13')
lc14=TimeSeries('ch14')
lc15=TimeSeries('ch15')

lclist=[lc9,lc10,lc11,lc12,lc14,lc15]

for k,lc in enumerate(lclist):
    lc=lc.divide(1/factor[k])


rate=np.array([lc.rate.mean() for lc in lclist])
rate_error=np.array([lc.rate.std() for lc in lclist])


plt.figure()

spe=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/rxte/results/out90427-01-03-02/products/fasebin_spe/conversion_factors/data_area_en.qdp',skip_header=3)



plt.errorbar(en,rate,rate_error,enerr,label='spectra_From_lc')
plt.plot(spe[:,0],spe[:,2],'k.')



a,b=np.polyfit(en[[0,1,4,5]],rate[[0,1,4,5]],1)

enaxis=np.linspace(4,9,100)
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

    N_opt,N_opt_err=curve_fit(best_line,en[[0,1,4,5]],rate[[0,1,4,5]],
                              p0=b,sigma=rate_error[[0,1,4,5]],absolute_sigma=True)
    myline=lambda x: a*x+N_opt
    ax[k].plot(enaxis,myline(enaxis),'k:',alpha=0.3)



    iron_band_flux=np.mean([rate[2],rate[3]])
    iron_band_en=np.mean([en[2],en[3]])
    iron_band_flux_err=np.sqrt(rate_error[2]**2+rate_error[3]**2)

    ax[k].errorbar(iron_band_en,iron_band_flux,iron_band_flux_err,[0.40924836],color='r')

    fe_flux=iron_band_flux-myline(iron_band_en)
    fe_flux_err=iron_band_flux_err

    ax[k].vlines(iron_band_en,iron_band_flux,iron_band_flux-fe_flux,color='c',alpha=0.99)



plt.ylim(0.1,0.25)
plt.show()





#%% find fe flux in time

def find_fe_flux(lclist,frac=1):
    N=len(lclist[0].rate)

    diff=[]
    diff_err=[]
    print(f'{frac*100}% OF THE LINEAR APPROX IS USED AS CONTINUUM')

    for i in range(N):
        rate=np.array([x.rate[i] for x in lclist])
        rate_error=np.array([x.error[i] for x in lclist])


        N_opt,N_opt_err=curve_fit(best_line,en[[0,1,4,5]],rate[[0,1,4,5]],
                                  p0=b,sigma=rate_error[[0,1,4,5]],absolute_sigma=True)

        myline=lambda x: a*x+N_opt


        iron_band_flux=np.mean([rate[2],rate[3]])
        iron_band_en=np.mean([en[2],en[3]])
        iron_band_flux_err=np.sqrt(rate_error[2]**2+rate_error[3]**2)

        fe_flux=iron_band_flux-myline(iron_band_en)*frac
        fe_flux_err=iron_band_flux_err
        #fe_flux_err=np.sqrt(ferr[1]**2)

        diff.append(fe_flux)
        diff_err.append(fe_flux_err)

    diff=np.asarray(diff)
    diff=np.reshape(diff,N)
    diff_err=np.asarray(diff_err)
    diff_err=np.reshape(diff_err,N)
    return diff,diff_err

frac=1
fe_flux,fe_flux_err=find_fe_flux(lclist,frac=frac)



#%% plot iron intensity
figure()
plt.errorbar(lclist[0].time,fe_flux,fe_flux_err)
plt.xlabel('Time, s')
plt.ylabel('Iron line flux')
plt.show()

figure()
plt.hist(fe_flux,bins=100,label='Fe_flux')
plt.xlabel('Iron line flux ')
plt.show()

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))



print(f'{frac*100}% OF THE LINEAR APPROX IS USED AS CONTINUUM')
print(f'Mean Flux: {np.mean(fe_flux)}')
print(f'Std Flux: {np.std(fe_flux)}')
wmean,wstd=weighted_avg_and_std(fe_flux, fe_flux_err**(-2))
print(f'Weighted mean flux: {wmean}')
print(f'Weighted std flux: {wstd}')
print(f'Weighted mean/weighted std : {wmean/wstd}')

print(f'Mean Error: {np.mean(fe_flux_err)}')
print(f'Mean significance (mean flux/err): {np.mean(fe_flux/fe_flux_err)}')



#%% save file
create_dir('fe_line')
os.system(f'cp ch13.lc_bary_orb_corr fe_line/python_lin_approx_fe_line_{frac}.lc_bary_orb_corr')
with fits.open(f'fe_line/python_lin_approx_fe_line_{frac}.lc_bary_orb_corr', mode='update') as hdul:
    hdul[1].data['rate']=fe_flux
    hdul[1].data['error']=fe_flux_err
    hdul.flush()  # changes are written back to original.fits



'''
XSPEC12>flux 6 7
 Model Flux    0.2104 photons (2.1877e-09 ergs/cm^2/s) range (6.0000 - 7.0000 keV)
 3    1   cflux      lg10Flux   cgs      -9.80041     +/-  3.05832E-02


'''


#%% mean of three lc

#rate_iron_band=(lc11.rate+lc12.rate+lc13.rate)/3
rate_iron_band=(lc11.rate+lc12.rate)/2

frac=0.6

delta_rate=lc14.rate-frac*rate_iron_band

print(f'''
     mean: {np.mean(delta_rate)}
     std: {np.std(delta_rate)}
     signif: {np.mean(delta_rate)/np.std(delta_rate)}
      ''')

plt.figure()
plt.errorbar(lc11.time,delta_rate,lc14.error)
plt.show()



os.chdir(path_to_lc)
os.system(f'cp ch11.lc_bary_orb_corr fe_line/python_fe_line_{frac}_no13.lc_bary_orb_corr')
with fits.open(f'fe_line/python_fe_line_{frac}_no13.lc_bary_orb_corr', mode='update') as hdul:
    hdul[1].data['rate']=delta_rate
    hdul[1].data['error']=lc14.error
    hdul.flush()  # changes are written back to original.fits

