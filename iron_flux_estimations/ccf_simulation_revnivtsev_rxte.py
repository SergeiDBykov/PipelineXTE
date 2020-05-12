#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 20:36:28 2020

@author: s.bykov
"""

import os
import seaborn as sns
from PipelineNuSTAR.core import *


#%% simulate lc #https://stingray.readthedocs.io/en/latest/notebooks/Simulator/Simulator%20Tutorial.html
import stingray
from stingray import Lightcurve, Crossspectrum, Powerspectrum
from stingray.simulator import simulator
from Misc.TimeSeries import cross_correlation

os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90089-11-03-01G/products/lc_sa')



def simulate_lc712_from_powerspectrum(N=10000,mean=3500,rms=0.16,dt=0.1,
                                   plot_results=0):

    sim = simulator.Simulator(N=N, mean=mean ,rms=rms, dt=dt)

    w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
    #from xspec
    xspec_pars=[0.119682,
           0.994604,
            25.6419,
           0.228065,
         0.00100764,
            7.50057,
           0.456827,
        0.000865741,
            31.3886,
           0.685431,
        0.000683068,
            1.85199,
                  0,
            1.45852,
           0.912076,
            1.91496]

    #https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node191.html
    def lorentzian( x, x0, gam, a ):
        #return a * gam**2 / ( gam**2 + ( x - x0 )**2)
        return a*(gam/(2*np.pi))/( (x-x0)**2 + (gam/2)**2  )
    def po(x,gamma,N):
        return N*x**(-gamma)

    spectrum = lorentzian(w, *xspec_pars[0:3])+lorentzian(w, *xspec_pars[3:6])+ lorentzian(w, *xspec_pars[6:9])+lorentzian(w,*xspec_pars[9:12])+po(w,*xspec_pars[12:14])+po(w,*xspec_pars[14:16])


    lc = sim.simulate(spectrum)
    #tmp=np.random.normal(lc.counts-lc.counts,270)
    #lc.counts=lc.counts+np.random.normal(lc.counts-lc.counts,270)
    if plot_results:
        plt.figure()
        plt.plot(lc.time,lc.counts)
        plt.show()

        plt.figure()
        ps=Powerspectrum(lc)
        ps=ps.rebin_log(0.01)
        plt.loglog(ps.freq,ps.power)
        plt.show()
    return lc


def iron_band_model_lc(lc712,A,deltaT):
    #lc67(t)= A*lc712(t)+(1-A)*lc712(t+dt)
    #dt>0 : delay with respect to lc712
    lc712_cts=lc712.counts

    tmpcounts=np.roll(lc712_cts,int(deltaT/lc712.dt))
    lc67_cts=A*lc712_cts+(1-A)*tmpcounts
    lc67=Lightcurve(lc712.time, lc67_cts)
    return lc67

def find_ccf(lc712,lc67,plot=0,deltaT=0,A=0):
    CCF_obj_crosscorr=cross_correlation.CrossCorrelation(lc712.time, lc67.counts,lc712.counts,circular=0)
    CCF_obj_crosscorr.calc_ccf()
    #plt.close()

    if plot:
        fig,ax=plt.subplots(1,sharex='all')
        ax.plot(CCF_obj_crosscorr.lag,CCF_obj_crosscorr.ccf,'b:.',label=f'simulations dT={deltaT}s; A={A}')
        ax.grid()
        N=int(CCF_obj_crosscorr.lag.shape[0]/2)
        ax.set_xlabel('Delay, s')
        ax.set_ylabel('CCF')
        ax.legend()
        plt.xlim(-15,15)
        plt.ylim(-0.5,1)
        plt.show()
    return CCF_obj_crosscorr

ccfs=[]
for i in range(500):
    print(i)
    lc712=simulate_lc712_from_powerspectrum()
    A=0.5
    deltaT=1.5
    lc67=iron_band_model_lc(lc712, A, deltaT)
    ccf=find_ccf(lc712, lc67,deltaT=deltaT,A=A,plot=0)
    ccfs.append(ccf.ccf)
ccfs=np.asarray(ccfs)


#7-12 autocorr
ccfs_autocorr=[]
for i in range(500):
    print(i)
    lc712=simulate_lc712_from_powerspectrum()
    ccf=find_ccf(lc712, lc712,deltaT=deltaT,A=A,plot=0)
    ccfs_autocorr.append(ccf.ccf)
ccfs_autocorr=np.asarray(ccfs_autocorr)

#%%plot simulations
fig,ax_ccf=plt.subplots(figsize=(16,6))

ax_ccf.errorbar(ccf.lag,ccfs.mean(axis=0),ccfs.std(axis=0),label=f'simulations dT={deltaT}s; A={A}')
ax_ccf.errorbar(ccf.lag,ccfs_autocorr.mean(axis=0),ccfs_autocorr.std(axis=0),label=f'simulations of autocorr (7-12 keV)',color='c')
ax_ccf.set_xlim(-15,15)
ax_ccf.set_ylim(-0.5,1)



def plot_ccf(filepath,ax):
    ccf=np.genfromtxt(filepath,skip_header=3)
    N=int(ccf[:,0].shape[0]/2)
    ax.errorbar(ccf[:,0],ccf[:,2],ccf[:,3],drawstyle='steps-mid',label='data')

    #ax.errorbar(-ccf[N:,0],ccf[N:,2],ccf[N:,3],alpha=0.5,color='r',drawstyle='steps-mid')
    ax.errorbar(ccf[N:,0],ccf[0:N+1:,2][::-1],ccf[0:N+1:,3][::-1],alpha=0.5,color='m',drawstyle='steps-mid',label='data (neg delay)')
    ax.set_xlabel('Delay, s')

    ax.set_ylabel('CCF')
    fig.tight_layout()
    sns.despine(fig,top=1,right=0)
    plt.show()

plot_ccf('ccf_0.1sec.qdp',ax_ccf)

ax_ccf.set_xlabel('iron delay, s')
ax_ccf.legend()
plt.show()

#plt.savefig(f'simulations/ccf_dt{deltaT}s_A{A}.png')
#%%plot difference ccf vs autocorr
fig,ax_ccf=plt.subplots(figsize=(16,6))

ax_ccf.errorbar(ccf.lag,ccfs.mean(axis=0)-ccfs_autocorr.mean(axis=0),np.sqrt(ccfs.std(axis=0)**2+ccfs_autocorr.std(axis=0)**2),label=f'ccf-acf')

ax_ccf.set_xlim(-15,15)
ax_ccf.set_ylim(-0.5,1)
ax_ccf.set_xlabel('iron delay, s')
ax_ccf.legend()
plt.show()



#%% plot real diff between ccf and acf
ccf=np.genfromtxt('ccf_0.1sec.qdp',skip_header=3)
acf=np.genfromtxt('acf_2031_0.1sec.qdp',skip_header=3)

fig,ax_diff=plt.subplots(figsize=(16,6))


ax_diff.errorbar(ccf[:,0],ccf[:,2]-0.65*acf[:,2],np.sqrt(ccf[:,3]**2+acf[:,3]**2))
ax_diff.set_xlim(-5,5)
plt.show()


#%% plot real diff between ccf and acf, but ccf in 2031 ch and 1415 channel
ccf=np.genfromtxt('ccf_0.1sec_2031_vs_1415.qdp',skip_header=3)
acf=np.genfromtxt('acf_2031_0.1sec.qdp',skip_header=3)

fig,ax_diff=plt.subplots(figsize=(16,6))

ax_diff.errorbar(ccf[:,0],ccf[:,2]-acf[:,2],np.sqrt(ccf[:,3]**2+acf[:,3]**2))
ax_diff.set_xlim(-5,5)
plt.show()



#%% plot the hist of a few correlation coeffs

fig,ax=plt.subplots()
ax.hist(ccfs[:,19999//2+12],bins=50)
plt.show()

fig,ax=plt.subplots()
ax.hist(ccfs[:,19999//2+23],bins=50)
plt.show()



#%% save data (trash)

tmp=np.vstack((lc67.time,lc67.counts)).T
np.savetxt(f'simulations/lc67_simul_dt{deltaT}s_A{A}.txt',tmp,delimiter=' ')

tmp=np.vstack((lc712.time,lc712.counts)).T
np.savetxt(f'simulations/lc712_simul_dt{deltaT}s_A{A}.txt',tmp,delimiter=' ')

