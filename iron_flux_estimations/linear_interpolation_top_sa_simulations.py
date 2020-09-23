#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 20:36:28 2020

@author: s.bykov
"""

import os
import seaborn as sns
from PipelineNuSTAR.core import *



from pipeline_core import *

ObsID='90089-11-03-02' # 90089-11-03-00G or 90089-11-03-01G 90089-11-03-02
xte_obs=ObservationXTE(ObsID)
os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/sa_data_lc_4.4sec')


binsize=4.4



#%% simulate lc #https://stingray.readthedocs.io/en/latest/notebooks/Simulator/Simulator%20Tutorial.html
import stingray
from stingray import Lightcurve, Crossspectrum, Powerspectrum
from stingray.simulator import simulator
from Misc.TimeSeries import cross_correlation

binsize=4.4
def simulate_lc712_from_powerspectrum(mean=3839,rms=202.54/3839,dt=binsize,
                                   plot_results=0):

    N=int(5000/dt)
    sim = simulator.Simulator(N=N, mean=mean ,rms=rms, dt=dt)

    sim = simulator.Simulator(N=N, mean=mean ,rms=rms, dt=dt)

    w = np.fft.rfftfreq(sim.N, d=sim.dt)[1:]
    #from xspec
    #model  lorentz + lorentz +  powerlaw  +   powerlaw

    xspec_pars=[ 0,
0.095447704,
  0.048784494,
0,
0.00044436465]

    #https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/node191.html
    def lorentzian( x, x0, gam, a ):
        return a*(gam/(2*np.pi))/( (x-x0)**2 + (gam/2)**2  )
    def po(x,gamma,N):
        return N*x**(-gamma)


    spectrum = lorentzian(w, *xspec_pars[0:3])+po(w,*xspec_pars[3:5])

    lc = sim.simulate(spectrum)
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

#%%

simulate_lc712_from_powerspectrum(plot_results=1)

#%% simulate
ccfs=[]
#from scipy.ndimage import gaussian_filter1d
for i in range(100):
    print(i)
    lc712=simulate_lc712_from_powerspectrum()
    A=0.8
    deltaT=4.4*3
    lc67=iron_band_model_lc(lc712, A, deltaT)
    #lc67.counts=gaussian_filter1d(lc67.counts,sigma=0.1)
    ccf=find_ccf(lc712, lc67,deltaT=deltaT,A=A,plot=0)
    ccfs.append(ccf.ccf)
ccfs=np.asarray(ccfs)


# #7-12 autocorr
# ccfs_autocorr=[]
# for i in range(5):
#     print(i)
#     lc712=simulate_lc712_from_powerspectrum()
#     ccf=find_ccf(lc712, lc712,deltaT=deltaT,A=A,plot=0)
#     ccfs_autocorr.append(ccf.ccf)
# ccfs_autocorr=np.asarray(ccfs_autocorr)


#%%plot simulations

import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook',rc={"xtick.top" : True,'xtick.direction':'inout','ytick.direction':'inout','xtick.minor.visible':True,'ytick.minor.visible':True})


def read_and_plot_ccf(ax,ObsID,name='python_lin_approx_fe_line_1',
                      xlims=[-150,150],
                      binsize='4.4'):
    matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
    matplotlib.rcParams['figure.subplot.left']=0.15
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=0.85
    matplotlib.rcParams['figure.subplot.top']=0.9

    ccf_fits=fits.open('fe_line/'+name+'.ccf')

    delay=ccf_fits[1].data['DELAY']
    crosscorr=ccf_fits[1].data['CROSS_CORR']
    error=ccf_fits[1].data['ERROR']
    N=int(len(delay)/2)


    ax.errorbar(delay,crosscorr,error,drawstyle='steps-mid',label=name)
    #ax.errorbar(delay[N:],crosscorr[0:N+1:][::-1],error[0:N+1:][::-1],alpha=0.5,color='r',drawstyle='steps-mid',ls=':')


    ax.set_xlim(xlims[0],xlims[1])
    ax.set_xlabel('Delay, s')
    ax.set_ylabel('CCF')
    #ax.set_title(f'{name}: {ObsID}, binsize={binsize} s')
    ax.grid()



    #plt.savefig(f'fe_line/python_lin_approx_fe_line_{frac}_{ObsID}_ccf.pdf')
    #plt.close()


fig,ax_ccf=plt.subplots()

frac=1
read_and_plot_ccf(ax_ccf,ObsID,name=f'python_lin_approx_fe_line_{frac}')


ax_ccf.errorbar(ccf.lag,ccfs.mean(axis=0)*0.4,ccfs.std(axis=0)*0.4,color='k',label=f'simulations dT={deltaT}s; A={A}',alpha=0.8,zorder=-10)

plt.savefig(f'simulations_{A}_{deltaT}.pdf')
plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/ccf_{frac}_{ObsID}_simul.pdf')

plt.show()

# def plot_ccf(filepath,ax):
#     ccf=np.genfromtxt(filepath,skip_header=3)
#     N=int(ccf[:,0].shape[0]/2)
#     norm=1#np.max(ccf[:,2])
#     ax.errorbar(ccf[:,0],ccf[:,2]/norm,ccf[:,3]/norm,color='b',drawstyle='steps-mid',label='data',alpha=0.8)
#     #ax.errorbar(-ccf[N:,0],ccf[N:,2],ccf[N:,3],alpha=0.5,color='r',drawstyle='steps-mid')
#     #ax.errorbar(ccf[N:,0],ccf[0:N+1:,2][::-1]/norm,ccf[0:N+1:,3][::-1]/norm,alpha=0.5,drawstyle='steps-mid',label='data (neg delay)')
#     ax.set_xlabel('Delay, s')

#     ax.set_ylabel('CCF')
#     fig.tight_layout()
#     sns.despine(fig,top=1,right=0)
#     plt.show()

# plot_ccf('/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/sa_data_lc_4.4sec/ccf_0.95.qdp',ax_ccf)

# ax_ccf.set_xlabel('Delay, s')
# #ax_ccf.legend()
# plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/ccf_{frac}_{Obs}_simul.pdf')

# plt.show()

# #plt.savefig(f'simulations/ccf_dt{deltaT}s_A{A}.png')
