#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 13:13:01 2020

@author: s.bykov
"""

from core_for_eqw import eqw_timeseries
from pipeline_core import *
ObsID='90089-11-03-01G' # 90427-01-03-11  90089-11-03-03  90089-11-03-01G ###bad or old: 90089-11-04-03  90089-11-03-01G
binsize='20'

stop
#%% no gauss

result=eqw_timeseries(ObsID,binsize,no_gauss=1)
time,lcrate,chi2_phase,xspec_rate=result.rate_check

time=time[:-1]
lcrate=lcrate[:-1]
chi2_phase=chi2_phase[:]
xspec_rate=xspec_rate[:]


fig,axs=plt.subplots(2,1)
fig.subplots_adjust(hspace=0.5)
delay,lag,ccf=my_crosscorr(time,lcrate,xspec_rate,axs[0],axs[1],
          subtract_mean=1,divide_by_mean=1,only_pos_delays=0,
          y1label='lc rate 3-12',y2label='xspec rate 3 12',my_only_pos_delays=0,my_ccf=0)
plt.show()
fig.savefig(f'rate_ccf_{result.dt}.png')
delay=delay[0]

fig,axs=plt.subplots()
axs.plot(time,lcrate)

roll=int(delay/result.dt)
new_rate=np.roll(xspec_rate,roll)

axs.plot(time,new_rate,'m:')

plt.show()
fig.savefig(f'rate_{result.dt}s.png')


#%% all
result=eqw_timeseries(ObsID,binsize,roll=int(delay/result.dt))




result.read_and_plot_data()
result.eqw_flux_ccf()
result.norm_flux_ccf()

#%% crosscorr heasoft
SSTOP
import os
os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90427-01-03-11/products/fasebin_timeser')

data=np.genfromtxt('cutoffpl/ph_res_cutoffpl.dat')
time=data[:,0]*20
norm=data[:,-3]
flux=data[:,7]

hdu=fits.PrimaryHDU()
c1=fits.Column(name='TIME',format='1D',unit='s',array=time)
c2=fits.Column(name='RATE',format='1D',unit='counts',array=flux)
c3=fits.Column(name='ERROR',format='1D',unit='counts',array=flux/20)
c4=fits.Column(name='FRACEXP',format='1D',unit='',array=flux/flux)

cols=fits.ColDefs([c1,c2])
datahdu=fits.BinTableHDU.from_columns(cols)
datahdu.name='RATE'
datahdu.header['TIME']=len(time)
datahdu.header['RATE']=len(flux)
datahdu.header['ERROR']=len(flux/20)
datahdu.header['FRACEXP']=len(flux/flux)
hdulist=fits.HDUList([hdu,datahdu])
hdulist.writeto('ph_res_flux.fits')
hdulist.close()