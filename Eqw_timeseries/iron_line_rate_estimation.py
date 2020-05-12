#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 15:59:38 2020

@author: s.bykov
"""


from PipelineXTE.pipeline_core import *

stop



#%% part 1: real effective stuff
os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90089-11-03-01G/products/fasebin_spe')


data_area=np.genfromtxt('lda_area_en.qdp',skip_header=3)
spe=data_area
data_noarea=np.genfromtxt('lda_noarea_en.qdp',skip_header=3)

data_area_ch=np.genfromtxt('lda_area_ch.qdp',skip_header=3)
data_noarea_ch=np.genfromtxt('lda_noarea_ch.qdp',skip_header=3)



plt.figure()
plt.plot(data_area[:,0],data_noarea[:,2]/data_area[:,2])
plt.xlabel('en, keV')
plt.ylabel('data_noarea/data_area, cm^2')
plt.show()


plt.figure()
plt.plot(data_area_ch[:,0],data_noarea_ch[:,2]/data_area_ch[:,2])
plt.xlabel('channel')
plt.ylabel('data_noarea/data_area, cm^2')
plt.show()



plt.figure()
plt.errorbar(data_area[:,0],data_area_ch[:,0],data_area_ch[:,1],data_area[:,1],'b+')
plt.xlabel('en, keV')
plt.ylabel('channel')

plt.show()


#%% mean rates

os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90089-11-03-01G/products/ironline_rate_estim')

class TimeSeries():

    def __init__(self,lcname):
        self.fits=fits.open(f'{lcname}.lc_bary')
        self.time=self.fits[1].data['time']
        self.rate=self.fits[1].data['rate']
        self.error=self.fits[1].data['error']
        self.binsize=np.median(np.diff(self.time))
        self.fits.close()
    def divide(self,val):
        self.rate=self.rate/val
        self.error=self.error/val

lc46=TimeSeries('ch1415')
lc67=TimeSeries('ch1619')
lc79=TimeSeries('ch2021')
#lc712=TimeSeries('lc712')
en=dE=np.array([6.045788+5.232003,7.270824+6.045788,8.09038+7.270824])/2

dE=np.array([6.045788-5.232003,7.270824-6.045788,8.09038-7.270824])
en_err=dE/2

area_fact=np.array([4000,2733,4332])

for k,lc in enumerate([lc46,lc67,lc79]):
    lc=lc.divide(dE[k]*area_fact[k])

rate=np.array([lc.rate.mean() for lc in [lc46,lc67,lc79]])
rate_error=np.array([lc.rate.std() for lc in [lc46,lc67,lc79]])

plt.errorbar(en,rate,rate_error,en_err)


plt.plot(spe[:,0],spe[:,2],'k.',label='data from xspec + 36%')

plt.xlim(4,9)
plt.show()
