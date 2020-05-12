#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 10 13:50:01 2020

@author: s.bykov
"""


from PipelineXTE.pipeline_core import *

os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90089-11-03-01G/products/std1_lc')

create_dir('./monte_carlo_period_error')

lcname='std1_0.1s_bary.lc_orb_corr'

lcfile=fits.open(lcname)
time=lcfile[1].data['time']
rate=lcfile[1].data['rate']
error=lcfile[1].data['error']

def gen_lc(rate,error):
    np.random.RandomState()
    newrate=rate+np.random.normal(loc=0,scale=error) #it adds to actual rate a random value assigned from normal distribution with zero meand and var as rate_error**2 (it is a vector, so error is different for each entry)
    return newrate


def gen_and_save_lc(rate,error,num=0):
    newrate=gen_lc(rate, error)
    os.system(f'cp {lcname} ./monte_carlo_period_error/gen{num}.fits')
    with fits.open(f'./monte_carlo_period_error/gen{num}.fits', mode='update') as hdul:
        hdul[1].data['rate']=newrate
        hdul.flush()  # changes are written back to original.fits

gen_and_save_lc(rate, error)



STOP
#%% test generator
time=lcfile[1].data['time'][0:11]
rate=lcfile[1].data['rate'][0:11]
error=lcfile[1].data['error'][0:11]

error[10]=error[10]*100/2

plt.figure()
plt.errorbar(time,rate,error,capsize=5)

#def geterate_lc(rate,error):
newrate=rate+np.random.normal(loc=0,scale=error)

plt.plot(time,newrate)


plt.show()