#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 16:08:29 2019

@author: s.bykov
"""

#%% imports and definitions
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from glob import glob
import os
from scipy.optimize import curve_fit
from PipelineXTE.pipeline_core import ObservationXTE
def gauss(t,t0,sigma,N):
    return N*np.exp(-(t-t0)**2/(2*sigma**2))



from Misc import  doppler_correction as doppler
from Misc.doppler_correction import  day2sec


matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
matplotlib.rcParams['figure.subplot.left']=0.10
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.9
matplotlib.rcParams['figure.subplot.top']=0.9




plt.ion()

def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)


def vals_and_errors(ObsParams,name,funct=lambda x: x):
    if isinstance(ObsParams,pd.core.frame.DataFrame):
        par,Min,Max=funct(ObsParams[name].values),funct(ObsParams[name+'_lo'].values),funct(ObsParams[name+'_hi'].values)
    elif isinstance(ObsParams,pd.core.series.Series):
        par,Min,Max=funct(ObsParams[name]),funct(ObsParams[name+'_lo']),funct(ObsParams[name+'_hi'])

    uplim_ind=Min==0

    low = par - Min
    hi  =  Max - par
    parr = par


    parr[uplim_ind]=hi[uplim_ind]
    low[uplim_ind]=parr[uplim_ind]
    hi[uplim_ind]=0


    err=np.vstack((low,hi))

    return parr,err




savepath='/Users/s.bykov/work/xray_pulsars/rxte/plots_results/'
results_path='/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pandas_data/'

filename='standard_pipeline' #standard_pipeline standard_pipeline_edge_en_free standard_pipeline_with_phabs
ObsParams=pd.read_pickle(results_path+f'{filename}.pkl')
ObsParams=ObsParams.sort_values('MJD_START')

ObsParams.period_orb_corr= ObsParams.period_orb_corr.replace(to_replace='None',value=np.nan)
ObsParams.period_orb_corr_err= ObsParams.period_orb_corr_err.replace(to_replace='None',value=np.nan)


ObsParams.loc[ObsParams.fasebin_cfg=='se','config']='E\_125us\_64M\_0\_1s'
ObsParams.loc[ObsParams.fasebin_cfg=='sa','config']='B\_16ms\_46M\_0\_49\_H'
ObsParams.loc[ObsParams.fasebin_cfg=='None','config']='-'


STOP



#%% compare spectral parameters
time=ObsParams.MJD_START

for par in ['norm_line',
            'eqw',
            'edgeTau']:

    fig,ax=plt.subplots(figsize=(9,6))

    for model,color in zip(['cutoffpl','edge_cutoffpl',],['k','r','m','g']):

        if model=='cutofpl' and par=='edgeTau':
            1
        else:
            par_val,par_err=vals_and_errors(ObsParams,model+'_'+par,lambda x: x)
            ecol='k'
            ax.plot(time,par_val,marker='d',mfc=color,mec=ecol,mew=1,ls='None',label=model+'_'+par,alpha=0.6)
            ax.errorbar(time,par_val,par_err,ecolor=ecol,fmt='none',alpha=0.5)
            ax.set_ylim(0.85*min(par_val),1.2*max(par_val))
    ax.set_ylabel(par)
    ax.set_xlabel('Time, MJD')
    ax.legend()
    ax.grid()
    fig.savefig(savepath+f'{par}.png')


#%% compare lines energies
fig,ax=plt.subplots(figsize=(6,4))
for model in ['cutoffpl','edge_cutoffpl','ionize_edge_cutoffpl']:
    eline=ObsParams[model+'_eline'].values
    eline_lo=ObsParams[model+'_eline_lo']
    eline_hi=ObsParams[model+'_eline_hi']
    eline_err=np.max((eline-eline_lo,eline_hi-eline),axis=0)

    ax.hist(eline,alpha=0.5,bins=25,label=model+f', mean={"%.2f" %np.mean(eline)} \n weighted mean {"%.2f"%np.average(eline,weights=eline_err**(-2))} ')
    ax.set_xlabel('Line energy, keV')
    plt.show()
ax.legend()
plt.savefig(savepath+f'eline_hist.png',dpi=250)





#%% compare chi2
fig,ax=plt.subplots()
for model in ['cutoffpl','edge_cutoffpl']:
    chi2=ObsParams[model+'_chi2']
    ax.hist(chi2,alpha=0.5,bins=25,label=model+f', mean={"%.2f" %np.mean(chi2)}')
    ax.set_xlabel('chi2_red')
    plt.show()
ax.legend()
plt.savefig(savepath+f'chi2_hist.png',dpi=250)



#%% plot cutoffpl residuals no gauss - finding K-edge

def plot_cutoffpl_ratio_phabs(ObsIDs,ax=None):
    model='cutoffpl_no_gauss'
    if ax==None:
        fig,ax = plt.subplots(figsize=(8, 6))
    else:
        pass

    for ObsID in ObsIDs:

        data1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2_top/{model}/mean_spe_rat.dat')


        label=f"Obs {ObsID}"

        ax.errorbar(data1[0],data1[1],data1[2],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    ax.set_xscale('log')

    ax.legend(loc='upper left',fontsize=8)
    ax.grid('y')
    ax.axhline(1,color='k',ls=':',zorder=-10,alpha=0.6)
    ax.set_ylabel('data / best fit phabs*cutoffpl ',fontsize=8)
    ax.set_xlabel('Energy, keV')
    plt.savefig(savepath+f'cutoffpl_ratio_phabs.png',dpi=250)



def plot_cutoffpl_ratio_no_phabs(ObsIDs,ax=None):
    model='cutoffpl_no_gauss'
    if ax==None:
        fig,ax = plt.subplots(figsize=(8, 6))
    else:
        pass

    for ObsID in ObsIDs:

        data1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2_top/{model}/mean_spe_rat.dat')


        label=f"Obs {ObsID}"

        ax.errorbar(data1[0],data1[1],data1[2],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    ax.set_xscale('log')

    ax.legend(loc='upper left',fontsize=8)
    ax.grid('y')
    ax.axhline(1,color='k',ls=':',zorder=-10,alpha=0.6)
    ax.set_ylabel('data / best fit cutoffpl ',fontsize=8)
    ax.set_xlabel('Energy, keV')
    plt.savefig(savepath+f'cutoffpl_ratio_no_phabs.png',dpi=250)



plot_cutoffpl_ratio_phabs(['90089-11-02-03','90089-11-03-01G','90089-22-01-01G','90014-01-05-04'])
plot_cutoffpl_ratio_no_phabs(['90089-11-02-03','90089-11-03-01G','90089-22-01-01G','90014-01-05-04'])




#%% plot residuals after gauss and edge no phabs

fig = plt.figure(figsize=(6.6, 6.6))
plt.subplots_adjust(hspace=0)
rows=4
cols=3

ax_eeufs=plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_del1 = plt.subplot2grid((rows,cols), (2, 0), rowspan=1, colspan=3,sharex=ax_eeufs)
ax_del2=plt.subplot2grid((rows,cols), (3, 0), rowspan=1, colspan=3,sharex=ax_eeufs)
#ax_del3=plt.subplot2grid((rows,cols), (4, 0), rowspan=1, colspan=3,sharex=ax_eeufs)

for ObsID,color in zip(['90089-11-05-08G','90089-11-04-02G'],['k','b']):
        label=ObsID
        eeufs=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2_top/cutoffpl/mean_spe_eeufs.dat')
        ax_eeufs.errorbar(eeufs[0],eeufs[1],eeufs[2],eeufs[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6,color=color)

        del1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2_top/cutoffpl/mean_spe_del.dat')
        ax_del1.errorbar(del1[0],del1[1],del1[2],del1[3],drawstyle='steps-mid',ls=':',alpha=0.6,color=color)

        del2=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2_top/edge_cutoffpl/mean_spe_del.dat')
        ax_del2.errorbar(del2[0],del2[1],del2[2],del2[3],drawstyle='steps-mid',ls=':',alpha=0.6,color=color)

ax_eeufs.legend()
ax_eeufs.set_ylabel('Ph/cm/s/keV')
ax_del1.set_ylabel('Sigma \n cutoffpl+gauss')
ax_del1.grid('y')
ax_del2.set_ylabel('Sigma \n cutoffpl*edge+gauss')
ax_del2.grid('y')
ax_del2.set_xlabel('Energy, keV')
ax_eeufs.set_xscale('log')

plt.savefig(savepath+f'edge_residuals.png',dpi=250)




#%% plot residuals after gauss and edge with phabs

fig = plt.figure(figsize=(6.6, 6.6))
plt.subplots_adjust(hspace=0)
rows=4
cols=3

ax_eeufs=plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_del1 = plt.subplot2grid((rows,cols), (2, 0), rowspan=1, colspan=3,sharex=ax_eeufs)
ax_del2=plt.subplot2grid((rows,cols), (3, 0), rowspan=1, colspan=3,sharex=ax_eeufs)
#ax_del3=plt.subplot2grid((rows,cols), (4, 0), rowspan=1, colspan=3,sharex=ax_eeufs)

for ObsID,color in zip(['90089-11-05-08G','90089-11-04-02G'],['k','b']):
        label=ObsID
        eeufs=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2_top/cutoffpl/mean_spe_eeufs.dat')
        ax_eeufs.errorbar(eeufs[0],eeufs[1],eeufs[2],eeufs[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6,color=color)

        del1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2_top/phabs_cutoffpl/mean_spe_del.dat')
        ax_del1.errorbar(del1[0],del1[1],del1[2],del1[3],drawstyle='steps-mid',ls=':',alpha=0.6,color=color)

        del2=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2_top/absorb_edge_cutoffpl/mean_spe_del.dat')
        ax_del2.errorbar(del2[0],del2[1],del2[2],del2[3],drawstyle='steps-mid',ls=':',alpha=0.6,color=color)

ax_eeufs.legend()
ax_eeufs.set_ylabel('Ph/cm/s/keV')
ax_del1.set_ylabel('Sigma \n (cutoffpl+gauss) * phabs')
ax_del1.grid('y')
ax_del2.set_ylabel('Sigma \n (cutoffpl*edge+gauss)*phabs')
ax_del2.grid('y')
ax_del2.set_xlabel('Energy, keV')
ax_eeufs.set_xscale('log')

plt.savefig(savepath+f'edge_residuals_phabs.png',dpi=250)





#%% plot dip decline factor for iron

fig,ax=plt.subplots(figsize=(9,6))


time=ObsParams.MJD_START

par_val,par_err=ObsParams.cutoffpl_en_fix_rel_dip_iron.values,ObsParams.cutoffpl_en_fix_rel_dip_iron_err.values
ecol='k'
color='red'
ax.plot(time[:-1],par_val[:-1],marker='d',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6)
ax.errorbar(time[:-1],par_val[:-1],par_err[:-1],ecolor=ecol,fmt='none',alpha=0.5)
ax.set_ylim(0.5,4)
ax.set_ylabel('I_max/I_min (Iron)')
ax.set_xlabel('Time, MJD')
ax.grid()
#fig.savefig(savepath+f'{par}.png')



#%% plot dip decline factor for flux

fig,ax=plt.subplots(figsize=(9,6))


time=ObsParams.MJD_START

par_val,par_err=ObsParams.cutoffpl_en_fix_rel_dip_flux.values,ObsParams.cutoffpl_en_fix_rel_dip_flux_err.values
time=time[~np.isnan(par_val)]
par_err=par_err[~np.isnan(par_val)]
par_val=par_val[~np.isnan(par_val)]
ecol='k'
color='red'
ax.plot(time,par_val,marker='d',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6)
ax.errorbar(time,par_val,par_err,ecolor=ecol,fmt='none',alpha=0.5)
ax.set_ylim(0.85*min(par_val),1.2*max(par_val))
ax.set_ylabel('F_Max/F_min')
ax.set_xlabel('Time, MJD')
#ax.grid()
#fig.savefig(savepath+f'{par}.png')

time=ObsParams.MJD_START
ax_twinx=ax.twinx()
par_val,par_err=vals_and_errors(ObsParams,'edge_cutoffpl'+'_'+'edgeTau',lambda x: x)
ecol='k'
ax_twinx.plot(time,par_val,marker='d',mfc='g',mec=ecol,mew=1,ls='None',alpha=0.6)
ax_twinx.errorbar(time,par_val,par_err,ecolor=ecol,fmt='none',alpha=0.5)
#ax_twinx.set_ylim(-0.001,0.1)
ax_twinx.set_ylabel('Edge_tau')
align_yaxis(ax, 1.15, ax_twinx, 0)
fig.savefig(savepath+f'edge_vs_dip.png')




#%% compare models


fig = plt.figure(figsize=(6.6*4, 4*6.6))
plt.subplots_adjust(hspace=0)
rows=8
cols=5

ax_norm=plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=5)
ax_eqw = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=5,sharex=ax_norm)
ax_eline=plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=5,sharex=ax_norm)

ax_edge=plt.subplot2grid((rows,cols), (6, 0), rowspan=1, colspan=5,sharex=ax_norm)

ax_po=plt.subplot2grid((rows,cols), (7, 0), rowspan=1, colspan=5,sharex=ax_norm)
ax_ecut=ax_po.twinx()

time=ObsParams.MJD_START




for model,color in zip(['edge_cutoffpl','cutoffpl'],['k','r']):


    norm_line,norm_line_err=vals_and_errors(ObsParams,model+'_norm_line',lambda x: x*1000)

    ecol='k'
    ax_norm.semilogx(time,norm_line,marker='d',mfc=color,mec=ecol,mew=1,ls='None',label=model+'_eqw',alpha=0.6)
    ax_norm.errorbar(time,norm_line,norm_line_err,ecolor=ecol,fmt='none',alpha=0.5)
    ax_norm.set_ylim(0.9*min(norm_line),1.1*max(norm_line))


    eqw,eqw_err=vals_and_errors(ObsParams,model+'_eqw',funct=lambda x: x*1000)

    ecol='k'
    ax_eqw.plot(time,eqw,marker='s',mfc=color,mec=ecol,mew=1,ls='None',label=model+'_eqw',alpha=0.6)
    ax_eqw.errorbar(time,eqw,eqw_err,ecolor=ecol,fmt='none',alpha=0.5)




    eline,eline_err=vals_and_errors(ObsParams,model+'_eline',funct=lambda x:x)

    ax_eline.errorbar(time,eline,eline_err,fmt='.',color=color,marker='d',ms=4,alpha=0.8)
    ax_eline.axhline(6.4,color='b',zorder=-10,ls=':',alpha=0.7)
    if model=='edge_cutoffpl':
        maxtau,maxtau_err=vals_and_errors(ObsParams,model+'_edgeTau',funct=lambda x:100*x)
        ax_edge.errorbar(time,maxtau,maxtau_err,fmt='.',color=color,marker='d',ms=4,alpha=0.8)
        ax_edge.set_ylim(0,5)

    ecut,ecut_err=vals_and_errors(ObsParams, model+'_ecut')
    po,po_err=vals_and_errors(ObsParams, model+'_po')

    ax_po.errorbar(time,po,po_err,fmt='.',color=color,marker='d',ms=4,alpha=0.8,label=model+'_ecut')
    ax_ecut.errorbar(time,ecut,ecut_err,marker='.',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6,label=model+'_po')



ax_norm.set_xlabel('Time, MJD')
ax_norm.set_ylabel('Gaussian norm, 1000 ph/s/cm^2')
ax_norm.legend(loc='best',fontsize=7)
ax_norm.grid()


ax_po.set_ylabel('Gamma')
ax_ecut.set_ylabel('Cutoff energy, keV')

ax_eline.set_ylabel('Iron line energy \n eV',color='k')
ax_eline.grid()
ax_eline.set_ylim(6.3,6.6)

ax_eqw.set_title(filename)
ax_eqw.set_xlabel('Time, MJD')
ax_eqw.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best',fontsize=7)
ax.grid()

ax_ecut.set_xlabel('Time, MJD')
ax_ecut.legend(loc='best',fontsize=7)


plt.show()







#%%plot flux continuum and iron

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

time=ObsParams.MJD_START

flux,flux_err=vals_and_errors(ObsParams,'cutoffpl_tot_flux',funct=lambda x: x/1e-8)

ax.errorbar(time,flux,flux_err,fmt='.',color='b',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b')
ax.set_xlabel('Time, MJD')


ax_iron=ax.twinx()
flux_iron,flux_err_iron=vals_and_errors(ObsParams,'cutoffpl_fe_flux',funct=lambda x: x/1e-10)

ax_iron.errorbar(time,flux_iron,flux_err_iron,fmt='.',color='g',marker='s',ms=4,alpha=0.6,zorder=-1)


ax_iron.set_ylabel('Iron line Flux, \n $10^{-10}$ cgs',color='g')


ax.axvspan(53340.29,53360.00,alpha=0.05, color='gray')
ax.axvspan(53384.36,53428.51,alpha=0.05, color='gray')
#ax.axvspan(53380,53380.8,alpha=0.05,color='gray')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'flux_cont_and_iron.pdf',dpi=500)



#%% plot eqw vs time vs distance

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw',funct=lambda x: 1e3*x)
ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line equivalent width \n eV',color='k')

ax.set_xlabel('Time, MJD')

ax.set_ylim(50,110)

ax.axvspan(53340.29,53360.00,alpha=0.05, color='gray')
ax.axvspan(53384.36,53428.51,alpha=0.05, color='gray')
#ax.axvspan(53380,53380.8,alpha=0.05,color='gray')

orbtime=np.linspace(53335,53420,200)
r,_,_,_,_=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)


ax2=ax.twinx()


ax2.plot(orbtime,r,'b-.',alpha=0.6,lw=0.5)


ax2.set_ylabel('Distance to the companion, \n lt-sec.')




fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw_vs_dist_in_time.pdf',dpi=500)
plt.show()



#%% plot iron line energy
fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))

time=ObsParams.MJD_START

el,el_err=vals_and_errors(ObsParams,'cutoffpl_eline')
ax.errorbar(time,el,el_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line energy \n eV',color='k')

ax.set_xlabel('Time, MJD')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eline.pdf',dpi=500)
plt.show()




#%% plot eqw vs time  with flux and distance

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(8,6))
ax_flux=ax.twinx()
time=ObsParams.MJD_START
eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw',funct=lambda x: 1e3*x)
flux,flux_err=vals_and_errors(ObsParams,'cutoffpl_tot_flux',funct=lambda x: x/1e-8)
ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8,zorder=10)

ax_flux.errorbar(time,flux,flux_err,fmt='.',color='b',marker='s',ms=4,alpha=0.3,zorder=-10)

ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b')
ax_flux.set_xlabel('Time, MJD')

ax.set_ylabel('Iron line equivalent width \n eV',color='c')
ax.set_xlabel('Time, MJD')
ax.set_ylim(40,120)



ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best')



orbtime=np.linspace(53335,53430,200)
r,_,_,_,_=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)


ax2=ax.twinx()


ax2.plot(orbtime,r,'r-.',alpha=0.6,lw=1,zorder=0)


ax2.set_ylabel('Distance to the companion, \n lt-sec.',color='r')

ax2.spines["right"].set_position(("axes", +1.15))



fig.tight_layout()
sns.despine(fig,top=1,right=0)
#plt.savefig(savepath+f'eqw_flux_and_dist.pdf',dpi=500)
plt.show()




#%% plot pulsed fractions vs time

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

time=ObsParams.MJD_START

eqw_PF,eqw_PF_err=ObsParams.eqw_PF,ObsParams.eqw_PF_err
ax.errorbar(time,eqw_PF,eqw_PF_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8,label='Eq. width')
#ax.plot(time,eqw_PF,color='c',marker='s',ms=4,alpha=0.8,label='Eq. width')

flux_gauss_PF,flux_gauss_PF_err=ObsParams.flux_gauss_PF,ObsParams.flux_gauss_PF_err
ax.errorbar(time,flux_gauss_PF,flux_gauss_PF_err,fmt='.',color='r',marker='s',ms=4,alpha=0.8,label='Fe flux')
#ax.plot(time,flux_gauss_PF,color='r',marker='s',ms=4,alpha=0.8,label='Fe flux')


flux_712_PF,flux_712_PF_err=ObsParams.flux_712_PF,ObsParams.flux_712_PF_err
ax.errorbar(time,flux_712_PF,flux_712_PF_err,fmt='.',color='g',marker='s',ms=4,alpha=0.8,label='7-12 Flux')


ax.set_ylabel('Pulsed fraction',color='k')
ax.set_xlabel('Time, MJD')
ax.set_ylim(0,1)
ax.legend()




# ax.axvspan(53340.29,53360.00,alpha=0.05, color='gray')
# ax.axvspan(53384.36,53428.51,alpha=0.05, color='gray')
# ax.axvspan(53380,53380.8,alpha=0.05,color='gray')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'pulsed_fractions.png',dpi=500)
plt.show()



#%% spe_pars_stuff


fig = plt.figure(figsize=(6.6*4, 4*6.6/2))
rows=5
cols=5
ax = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=5)

ax2=plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=5)

ax3=plt.subplot2grid((rows,cols), (4, 0), rowspan=1, colspan=2)


time=ObsParams.MJD_START

#for model,color in zip(['cutoffpl','pohi'],['r','g']):
for model,color in zip(['cutoffpl'],['c']):

    eqw,eqw_err=vals_and_errors(ObsParams,model+'_eqw')
    eqw=eqw*1000
    eqw_err=eqw_err*1000


    ecol='k'
    ax.plot(time,eqw,marker='s',mfc=color,mec=ecol,mew=1,ls='None',label=model+'_eqw',alpha=0.6)
    ax.errorbar(time,eqw,eqw_err,ecolor=ecol,fmt='none',alpha=0.5)

    #ss=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/rxte/plots_results/SSandAA.csv')
    #ax.errorbar(ss[:,0],ss[:,1]*1e3*0.8,ss[:,1]/ss[:,1]*3*0.8,color='r',fmt='none',label='SS paper * 80\% , approx errors 3 eV',alpha=0.4)


    efold=ObsParams[model+'_efold']
    ecut=ObsParams[model+'_ecut']
    chi2=ObsParams[model+'_chi2']
    gamma=ObsParams[model+'_po']
    ax2.plot(time,efold,marker='d',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6,label=model+'_efold')
    #ax2.plot(time,ecut,marker='.',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6,label=model+'_ecut')
    ax2_twin=ax2.twinx()
    ax2_twin.plot(time,gamma,marker='.',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6,label=model+'_po')

    chi2.hist(ax=ax3,color=color,alpha=0.5,label=model+'_chi2')
    #ObsParams.iloc[where(ObsParams.cutoffpl_chi2>2)[0]][['ObsID','MJD_START','cutoffpl_chi2']]
ax.set_title(filename)
ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best',fontsize=7)
ax.grid()

ax2.set_xlabel('Time, MJD')
ax2.set_ylabel('E_fold/E_cut, keV')
ax2.legend(loc='best',fontsize=7)
ax2.grid()

ax3.legend()

plt.show()





#%% plot eqw for cutoffpl

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw',funct=lambda x: 1e3*x)
cond=ObsParams.cutoffpl_chi2<2
ax.errorbar(time[cond],eqw[cond],eqw_err.T[cond].T,fmt='.',color='c',marker='s',ms=4,label='cutoffpl ',alpha=0.8,zorder=10)


ax.set_ylabel('Iron line equivalent width \n eV',color='k')
ax.set_xlabel('Time, MJD')
ax.set_ylim(40,120)

ax.legend(loc='upper left')
fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw_cutoff.png',dpi=500)
plt.show()





#%% plot eqw for gabs

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'gabs_eqw',funct=lambda x: 1e3*x)
cond=ObsParams.gabs_chi2<10
ax.errorbar(time[cond],eqw[cond],eqw_err.T[cond].T,fmt='.',color='c',marker='s',ms=4,label='gabs ',alpha=0.8,zorder=10)


ax.set_ylabel('Iron line equivalent width \n eV',color='k')
ax.set_xlabel('Time, MJD')
ax.set_ylim(40,120)

ax.legend(loc='upper left')
fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw_gabs.png',dpi=500)
plt.show()



#%% plot eqw for cyclabs

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'cyclabs_eqw',funct=lambda x: 1e3*x)
cond=ObsParams.cyclabs_chi2<2
ax.errorbar(time[cond],eqw[cond],eqw_err.T[cond].T,fmt='.',color='c',marker='s',ms=4,label='cyclabs ',alpha=0.8,zorder=10)


ax.set_ylabel('Iron line equivalent width \n eV',color='k')
ax.set_xlabel('Time, MJD')
#ax.set_ylim(40,120)



ax.legend(loc='upper left')
fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw_cyclabs.png',dpi=500)
plt.show()




#%% plot cycle line energy for GABS

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

time=ObsParams.MJD_START

eline=ObsParams['gabs_eline']
cond=ObsParams.gabs_chi2<10
ax.plot(time[cond],eline[cond],color='c',marker='s',ms=4,label='gabs model (left)',alpha=0.8,zorder=10)


ax.set_ylabel('Cycle line energy, keV',color='k')
ax.set_xlabel('Time, MJD')
#ax.set_ylim(40,120)

ax.legend(loc='upper left')
fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'ecycle_gabsmodel.png',dpi=500)

plt.show()



#%% plot cycle line energy for cycleabs

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

time=ObsParams.MJD_START

eline=ObsParams['cycleabs_eline']
cond=ObsParams.gabs_chi2<10
ax.plot(time[cond],eline[cond],color='c',marker='s',ms=4,label='cycleabs model (left)',alpha=0.8,zorder=10)


ax.set_ylabel('Cycle line energy, keV',color='k')
ax.set_xlabel('Time, MJD')
#ax.set_ylim(40,120)



ax.legend(loc='upper left')
fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'ecycle_cycleabsmodel.png',dpi=500)

plt.show()



#%% latex table
from Misc.TeX_Tables import pandas_to_tex
from Misc.TeX_Tables.pandas_to_tex import *

model='edge_cutoffpl_'
def tex_tbl():
    null=lambda x: x
    free_columns=['ObsID','MJD_START','EXPOSURE','config',model+'chi2']
    free_columns_functions=[null,null,null,null,null]
    free_columns_formats=[0,1,0,0,2]

    err_columns=['eqw','cutoffpl_flux',
                 ]
    err_functions=[lambda x: 1000*x, lambda x: x/1e-9,
                   ]
    err_formats=[0,2]

    err_columns=[model+item for item in err_columns]

    headers=['ObsID','Time, MJD','Exposure, s','Configuration','$\chi^2_{red}$',
                  'Eq. width', 'Flux (3-12 keV)'
                 ]

    transpose=0
    df_tex=make_latex_table(df,
                          free_columns=free_columns, free_columns_functions=free_columns_functions,
                          free_columns_formats=free_columns_formats,
                          err_columns=err_columns, err_functions=err_functions,
                          err_formats=err_formats)
    df_tex.columns=headers
    save_latex_table(df_tex, savepath=savepath+'/tex/allpars_'+name+'.tex',
                     columns_to_write='DEFAULT',
                     columns_names='DEFAULT',
                     transpose=transpose)
df=ObsParams
tex_tbl()




#%% cutoffpl with edge: line energy, etc


#%%eqw
fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_edge_eqw',funct=lambda x: 1e3*x)
ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line equivalent width \n eV',color='k')

ax.set_xlabel('Time, MJD')

ax.set_ylim(50,110)

ax.axvspan(53340.29,53360.00,alpha=0.05, color='gray')
ax.axvspan(53384.36,53428.51,alpha=0.05, color='gray')
#ax.axvspan(53380,53380.8,alpha=0.05,color='gray')

orbtime=np.linspace(53335,53420,200)
r,_,_,_,_=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)


ax2=ax.twinx()


ax2.plot(orbtime,r,'b-.',alpha=0.6,lw=0.5)


ax2.set_ylabel('Distance to the companion, \n lt-sec.')




fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.show()



#%% line energy

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))

time=ObsParams.MJD_START

eline,eline_err=vals_and_errors(ObsParams,'cutoffpl_edge_eline',funct=lambda x:x)
ax.errorbar(time,eline,eline_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line energy \n eV',color='k')

ax.set_xlabel('Time, MJD')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.show()



#%% absorbtion edge tau


fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))

time=ObsParams.MJD_START

eline,eline_err=vals_and_errors(ObsParams,'cutoffpl_edge_maxtau',funct=lambda x:x*100)
ax.errorbar(time,eline,eline_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('K-edge tau *100',color='k')

ax.set_xlabel('Time, MJD')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.show()

#%% spectral stuff


fig = plt.figure(figsize=(6.6*4, 4*6.6))
rows=6
cols=5

ax = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=5)
ax2=plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=5)
ax3=plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=5)


time=ObsParams.MJD_START

model='cutoffpl_edge'
color='k'

eqw,eqw_err=vals_and_errors(ObsParams,model+'_eqw')
eqw=eqw*1000
eqw_err=eqw_err*1000


ecol='k'
ax.plot(time,eqw,marker='s',mfc=color,mec=ecol,mew=1,ls='None',label=model+'_eqw',alpha=0.6)
ax.errorbar(time,eqw,eqw_err,ecolor=ecol,fmt='none',alpha=0.5)


efold=ObsParams[model+'_efold']

maxtau,maxtau_err=vals_and_errors(ObsParams,model+'_maxtau',funct=lambda x:100*x)
ax3.errorbar(time,maxtau,maxtau_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)
ax3.set_ylabel('max_tau*100')
eline,eline_err=vals_and_errors(ObsParams,'cutoffpl_edge_eline',funct=lambda x:x)
ax4=ax3.twinx()
ax4.errorbar(time,eline,eline_err,fmt='.',color='r',marker='d',ms=4,alpha=0.8)

ax4.set_ylabel('Iron line energy \n eV',color='k')
ax4.axhline(6.4,color='b',zorder=-10,ls=':',alpha=0.7)


gamma=ObsParams[model+'_po']
ax2.plot(time,efold,marker='d',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6,label=model+'_efold')

ax2_twin=ax2.twinx()
ax2_twin.plot(time,gamma,marker='.',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6,label=model+'_po')
ax2_twin.set_ylabel('Gamma')


ax.set_title(filename)
ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best',fontsize=7)
ax.grid()

ax2.set_xlabel('Time, MJD')
ax2.set_ylabel('E_fold, keV')
ax2.legend(loc='best',fontsize=7)
ax2.grid()

ax3.legend()

plt.show()






#%% compare fluxes



fig,ax = plt.subplots(figsize=(6.6*4, 4*6.6))
rows=4
cols=5
ax = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=5)

ax2=plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=5)

time=ObsParams.MJD_START




for model,color in zip(['cutoffpl_edge','cutoffpl'],['k','r']):

    flux_line,flux_line_err=vals_and_errors(ObsParams,model+'_fe_flux',lambda x: x/1e-10)

    ecol='k'
    ax.plot(time,flux_line,marker='s',mfc=color,mec=ecol,mew=1,ls='None',label=model+'_eqw',alpha=0.6)
    ax.errorbar(time,flux_line,flux_line_err,ecolor=ecol,fmt='none',alpha=0.5)

    norm_line,norm_line_err=vals_and_errors(ObsParams,model+'_norm_line',lambda x: x*1000)

    ecol='k'
    ax2.plot(time,norm_line,marker='d',mfc=color,mec=ecol,mew=1,ls='None',label=model+'_eqw',alpha=0.6)
    ax2.errorbar(time,norm_line,norm_line_err,ecolor=ecol,fmt='none',alpha=0.5)


ax.set_title(filename)
ax.set_xlabel('Time, MJD')
ax.set_ylabel('Gaussian flux, 1e-10 cgs')
ax.legend(loc='best',fontsize=7)
ax.grid()

ax2.set_xlabel('Time, MJD')
ax2.set_ylabel('Gaussian norm, 1000 ph/s/cm^2')
ax2.legend(loc='best',fontsize=7)
ax2.grid()

plt.show()







#%% plot ration or delchi for cutoffpl continua
from pipeline_core import xspec_scripts_path
fig,ax = plt.subplots(figsize=(8, 3))


for ObsID in ['90014-01-05-03']:#,'90089-11-01-04','90089-11-02-00','90089-11-03-01G','90089-11-05-08G','90427-01-03-00','90014-01-03-020', '90427-01-04-02']:

    xte_obs=ObservationXTE(ObsID)
    ser=xte_obs.pandas_series()
    model='cutoffpl_no_gauss'
    spe_path='/products/pcu2_top'
    os.chdir(xte_obs.out_path+spe_path)
    os.system('rm -rf cutoffpl_no_gauss')
    os.mkdir('cutoffpl_no_gauss')

    os.system(f'xspec - {xspec_scripts_path}{model}.txt')

    data=np.genfromtxt('cutoffpl_no_gauss/mean_spe_rat.dat')

    label=f"Obs {ObsID}, F_{{3-12 keV}}={'%.2f'%(ser.edge_cutoffpl_cutoffpl_flux/1e-8)} 10^{-8} csg"

    ax.errorbar(data[0],data[1],data[2],data[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    ax.set_xscale('log')



ax.legend(loc='upper left',fontsize=8)
ax.axhline(1,color='k',ls=':',zorder=-10,alpha=0.6)
ax.set_ylabel('data/ [cutoffpl model] OR sqrt(chi^2)',fontsize=8)

plt.show()




#%% JUNE TRASH
'''

#%% fun with correlations

corrdata=ObsParams[['cutoffpl_efold','cutoffpl_eline','cutoffpl_eqw',
                    'cutoffpl_fe_flux','cutoffpl_po','cutoffpl_tot_flux',
                   'deltat' ]]

sns.heatmap(corrdata.corr(),cmap='jet')



#%%plot  period stuff with approximation

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

time=ObsParams.MJD_START[~ObsParams.period_orb_corr.isna()]
ax.set_xlabel('Time, MJD')


ax_period=ax

period=ObsParams.period_orb_corr[~ObsParams.period_orb_corr.isna()]
period_err=ObsParams.period_orb_corr_err[~ObsParams.period_orb_corr.isna()]

factor=doppler.kepler_solution(time*day2sec, doppler.orb_params_v0332)[3]

period_bary=factor*period


ax_period.errorbar(time,period,period_err,color='g',marker='s',ls='',ms=4,alpha=0.8)
ax_period.plot(time,period_bary,color='m',marker='.',ls='',ms=2,alpha=0.6)

ax_period.set_ylabel('Period, s',color='g')
ax_period.set_ylim(4.373,4.377)

ax.axvspan(53340.29,53360.00,alpha=0.2, color='gray')
ax.axvspan(53384.36,53428.51,alpha=0.2, color='gray')

from scipy.optimize import curve_fit


t0=53340
def line(t,p0,pdot):
    return p0+(t-t0)*pdot

popt,pcov=curve_fit(line,time,period,p0=[4.373,4e-6])

ax.plot(time,line(time,*popt),'k:',label='')
print(f'All outburst: P(t)=P0+Pdot*(t-53340), P0={popt[0]} Pdot={popt[1]}+-{np.sqrt(np.diag(pcov)[1])}')


popt,pcov=curve_fit(line,time[(53340<time) & (time<53360)],period[(53340<time) & (time<53360)],p0=[4.373,4e-6],
                    sigma=period_err[(53340<time) & (time<53360)],absolute_sigma=1)
ax.plot(time[(53340<time) & (time<53360)],line(time[(53340<time) & (time<53360)],*popt),'r:',label='')
print(f'Rising of the outburst: P(t)=P0+Pdot*(t-53340), P0={popt[0]} Pdot={popt[1]}+-{np.sqrt(np.diag(pcov)[1])}')



popt,pcov=curve_fit(line,time[(53384.36<time) & (time<53428.51)],period[(53384.36<time) & (time<53428.51)],p0=[4.373,4e-6],
                    sigma=period_err[(53384.36<time) & (time<53428.51)],absolute_sigma=1)
ax.plot(time[(53384.36<time) & (time<53428.51)],line(time[(53384.36<time) & (time<53428.51)],*popt),'m:',label='')
print(f'Decay part of the outburst: P(t)=P0+Pdot*(t-53340), P0={popt[0]} Pdot={popt[1]}+-{np.sqrt(np.diag(pcov)[1])}')


popt,pcov=curve_fit(line,time[(53384.36<time) & (time<53400)],period[(53384.36<time) & (time<53400)],p0=[4.373,4e-6],
                    sigma=period_err[(53384.36<time) & (time<53400)],absolute_sigma=1)
ax.plot(time[(53384.36<time) & (time<53400)],line(time[(53384.36<time) & (time<53400)],*popt),'r:',label='')

print(f'Decay part of the outburst before 53400: P(t)=P0+Pdot*(t-53340), P0={popt[0]} Pdot={popt[1]}+-{np.sqrt(np.diag(pcov)[1])}')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'period_evol.png',dpi=500)

plt.show()




#%% plot eqw vs time with orbital convolution

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw',funct=lambda x: 1e3*x)

ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line equivalent width \n eV',color='k')
ax.set_xlabel('Time, MJD')
ax.set_ylim(40,120)


orbtime=np.linspace(53330,53420,200)
r,_,_,_,_=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)

#popt=np.array([5.33661218e+04, 1.58437049e+01, 3.58940948e+00])
#popt,pcov=curve_fit(gauss,time[time<53385],eqw[time<53385],p0=[53370,20,50])
#popt,pcov=curve_fit(gauss,time,eqw,sigma=eqw_err.max(axis=0),p0=[53370,20,50])
popt,pcov=curve_fit(gauss,time,eqw,p0=[53370,20,50])

gaussian_depend=gauss(orbtime,popt[0],popt[1],popt[2])
inv=r**(-2)
y=gaussian_depend+70000*(inv-np.mean(inv)) #2000 for d^(-1), 70000 for d^(-2)

ax.plot(orbtime,gaussian_depend,'k--',alpha=0.7)
ax.plot(orbtime,y,'k-',alpha=0.7)

ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw_dist_convolve.pdf',dpi=500)
plt.show()







#%% plot eqw vs time with  my orbital convolution

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw',funct=lambda x: 1e3*x)

flux,flux_err=vals_and_errors(ObsParams,'cutoffpl_tot_flux',funct=lambda x: x/1e-8)
flux_time=time


#asm light curve
#asm=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/rxte/plots_results/ASM_LC_V0332+53.txt',skip_header=5)
#flux_time=asm[:,0]
#flux=asm[:,7]



ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line equivalent width \n eV',color='k')
ax.set_xlabel('Time, MJD')
ax.set_ylim(40,120)


orbtime=np.linspace(53330,53420,200)
r,_,_,_,_=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)

star_eqw=1/r**2
#star_eqw=1-np.sqrt(r**2-9**2)/r
star_eqw /=np.max(star_eqw)
star_eqw=star_eqw*35

from scipy import interpolate
#flux_funct=interpolate.interp1d(time,flux,fill_value='extrapolate')(orbtime)
flux_funct=interpolate.interp1d(flux_time,flux,fill_value=0.5,bounds_error=0)(orbtime)

disk_eqw=77*(flux_funct/np.max(flux_funct))**(3/20)+10

ax.plot(orbtime,star_eqw+disk_eqw,'k-',alpha=0.7)

ax.plot(orbtime,disk_eqw,'k-.',alpha=0.7)
#ax.plot(orbtime,star_eqw,'k-.',alpha=0.7)


ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
#plt.savefig(savepath+f'eqw_dist_convolve_disk_and_star.pdf',dpi=500)
plt.show()



#%% plot eqw vs time with BOTH orb convolutions

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw',funct=lambda x: 1e3*x)

ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line equivalent width \n eV',color='k')
ax.set_xlabel('Time, MJD')
ax.set_ylim(40,120)


orbtime=np.linspace(53330,53420,200)
r,_,_,_,_=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)

#popt=np.array([5.33661218e+04, 1.58437049e+01, 3.58940948e+00])
#popt,pcov=curve_fit(gauss,time[time<53385],eqw[time<53385],p0=[53370,20,50])
#popt,pcov=curve_fit(gauss,time,eqw,sigma=eqw_err.max(axis=0),p0=[53370,20,50])
popt,pcov=curve_fit(gauss,time,eqw,p0=[53370,20,50])

gaussian_depend=gauss(orbtime,popt[0],popt[1],popt[2])
inv=r**(-2)
y=gaussian_depend+70000*(inv-np.mean(inv)) #2000 for d^(-1), 70000 for d^(-2)

ax.plot(orbtime,gaussian_depend,'g:',alpha=0.7)
ax.plot(orbtime,y,'g-',alpha=0.7)

ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best')



flux,flux_err=vals_and_errors(ObsParams,'cutoffpl_tot_flux',funct=lambda x: x/1e-8)
flux_time=time


orbtime=np.linspace(53330,53420,200)
r,_,_,_,_=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)

star_eqw=1/r**2
#star_eqw=1-np.sqrt(r**2-9**2)/r
star_eqw /=np.max(star_eqw)
star_eqw=star_eqw*35

from scipy import interpolate
#flux_funct=interpolate.interp1d(time,flux,fill_value='extrapolate')(orbtime)
flux_funct=interpolate.interp1d(flux_time,flux,fill_value=0.5,bounds_error=0)(orbtime)

disk_eqw=77*(flux_funct/np.max(flux_funct))**(3/20)

ax.plot(orbtime,star_eqw+disk_eqw,'r-',alpha=0.7)

ax.plot(orbtime,disk_eqw,'r:',alpha=0.7)
#ax.plot(orbtime,star_eqw,'k-.',alpha=0.7)


ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw_dist_convolve_both.pdf',dpi=500)
plt.show()





#%%plot flux and period stuff

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

time=ObsParams.MJD_START

flux,flux_err=vals_and_errors(ObsParams,'cutoffpl_tot_flux',funct=lambda x: x/1e-8)

ax.errorbar(time,flux,flux_err,fmt='.',color='b',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('RXTE/PCA Flux (3-12 keV), \n 10^(-8) cgs',color='b')
ax.set_xlabel('Time, MJD')


ax_period=ax.twinx()

period=ObsParams.period_orb_corr

factor=doppler.kepler_solution(time*day2sec, doppler.orb_params_v0332)[3]
period_bary=factor*ObsParams.period_orb_corr


ax_period.plot(time,period,color='g',marker='s',ls='',ms=4,alpha=0.8)
ax_period.plot(time,period_bary,color='m',marker='.',ls='',ms=2,alpha=0.6)

ax_period.set_ylabel('Period, s',color='g')
ax_period.set_ylim(4.373,4.377)


ax.axvspan(53340.29,53360.00,alpha=0.05, color='gray')
ax.axvspan(53384.36,53428.51,alpha=0.05, color='gray')
ax.axvspan(53380,53380.8,alpha=0.05,color='gray')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'flux_and_period.png',dpi=500)

from scipy.optimize import curve_fit

popt,pcov=curve_fit(gauss,time[time<53385],flux[time<53385],p0=[53370,20,3])

ax.plot(time,gauss(time,*popt),'k:',label='')

plt.show()




#%% plot eqw vs distance

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

time=ObsParams.MJD_START
orbtime=time
r,_,_,_,_=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)

eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw',funct=lambda x: 1e3*x)
flux,flux_err=vals_and_errors(ObsParams,'cutoffpl_tot_flux',funct=lambda x: x/1e-8)

ax.errorbar(r,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line equivalent width \n eV',color='k')
ax.set_xlabel('Distance to the companion, \n lt-sec.')
ax.set_ylim(40,120)

fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw_vs_dist.png',dpi=500)
plt.show()




#%% plot eqw vs time without models but with flux

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6,6))
ax_flux=ax.twinx()
time=ObsParams.MJD_START
eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw',funct=lambda x: 1e3*x)
flux,flux_err=vals_and_errors(ObsParams,'cutoffpl_tot_flux',funct=lambda x: x/1e-8)
ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax_flux.errorbar(time,flux,flux_err,fmt='.',color='gray',marker='o',ms=4,alpha=0.3)

ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='gray')
ax_flux.set_xlabel('Time, MJD')

ax.set_ylabel('Iron line equivalent width \n eV',color='c')
ax.set_xlabel('Time, MJD')
ax.set_ylim(40,120)



ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw.pdf',dpi=500)
plt.show()


#%% find observations after 53400 with small exposure in phase-resolved

expos_ph_res=[]
for ObsID in ObsParams.ObsID:
    try:
        ff=fits.open(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/fasebin_sys.pha')
        exposure=ff[1].data['exposure'].mean()
        expos_ph_res.append(exposure)
    except:
        expos_ph_res.append(0)

expos_ph_res=np.array(expos_ph_res)

for ex,obs,mjd, conf in zip(expos_ph_res,ObsParams.ObsID,ObsParams.MJD_START,ObsParams.fasebin_cfg):
    print(f'{obs}: expo/phase bin {round(ex)}s (MJD {mjd}, conf {conf})')

plt.figure()
ObsParams.plot(x='MJD_START',y='EXPOSURE')

plt.plot(ObsParams.MJD_START,expos_ph_res)

plt.show()
#bad periods or very small expo after 53400 mjd
badpers=['90014-01-04-02','90014-01-04-03','90014-01-05-00']

small_expo=['90014-01-05-03','90014-01-05-06','90014-01-06-02',
            '90014-01-06-03','90014-01-07-01','90014-01-07-03',
            '90014-01-07-02','90014-01-07-04','90014-01-07-00',
            '90014-01-08-00','90014-01-08-01']

#ignore_fasebin - '90014-01-05-03'


#%% plot delay stuff vs time vs distance

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

df=ObsParams[~ObsParams.ObsID.isin(ignored_obs)]

deltat,deltat_err=df.deltat,df.deltat_err


time=df.MJD_START

deltat,deltat_err=df.deltat,df.deltat_err
ax.errorbar(time,deltat,deltat_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line equivalent width delay, sec',color='k')
ax.set_xlabel('Time, MJD')


# ax.axvspan(53340.29,53360.00,alpha=0.05, color='gray')
# ax.axvspan(53384.36,53428.51,alpha=0.05, color='gray')
# ax.axvspan(53380,53380.8,alpha=0.05,color='gray')

orbtime=np.linspace(53335,53420,200)
r,_,_,_,z=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)


ax2=ax.twinx()


ax2.plot(orbtime,r,'b-.',alpha=0.6,lw=0.5)


ax2.set_ylabel('Distance to the companion, \n lt-sec.')

fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'delay.pdf',dpi=500)
plt.show()





#%% plot delay (pos and absolute) stuff vs time vs distance

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

df=ObsParams[~ObsParams.ObsID.isin(ignored_obs)]

deltat_pos,deltat_pos_err=df.deltat_pos,df.deltat_pos_err
time=df.MJD_START

ax.errorbar(time,deltat_pos,deltat_pos_err,fmt='.',color='b',marker='s',ms=4,alpha=0.8)


deltat_abs,deltat_abs_err=df.deltat_absolute,df.deltat_pos_err
ax.errorbar(time,deltat_abs,deltat_abs_err,fmt='.',color='k',marker='s',ms=3,alpha=0.5)

#deltat_neg,deltat_neg_err=df.deltat_neg,df.deltat_neg_err
#ax.errorbar(time,deltat_neg,deltat_neg_err,fmt='.',color='r',marker='s',ms=3,alpha=0.5)



ax.set_ylabel('Iron line equivalent width delay, sec',color='k')
ax.set_xlabel('Time, MJD')


# orbtime=np.linspace(53335,53420,200)
# r,_,_,_,z=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)
# ax2=ax.twinx()
# ax2.plot(orbtime,r,'b-.',alpha=0.6,lw=0.5)
#ax2.set_ylabel('Distance to the companion, \n lt-sec.')

fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'delays.pdf',dpi=500)
plt.show()



#%% plot normalization delay stuff vs time vs distance

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

df=ObsParams[~ObsParams.ObsID.isin(ignored_obs)]

deltat,deltat_err=df.norm_deltat,df.norm_deltat_err
time=df.MJD_START

ax.errorbar(time,deltat,deltat_err,fmt='.',color='r',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line normalization width delay, sec',color='k')
ax.set_xlabel('Time, MJD')


# ax.axvspan(53340.29,53360.00,alpha=0.05, color='gray')
# ax.axvspan(53384.36,53428.51,alpha=0.05, color='gray')
# ax.axvspan(53380,53380.8,alpha=0.05,color='gray')

orbtime=np.linspace(53335,53420,200)
r,_,_,_,z=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)


ax2=ax.twinx()


ax2.plot(orbtime,r,'b-.',alpha=0.6,lw=0.5)


ax2.set_ylabel('Distance to the companion, \n lt-sec.')

fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'norm_delay.pdf',dpi=500)
plt.show()


#%% plot delay stuff from my monte carlo estimation

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

df=ObsParams[~ObsParams.ObsID.isin(ignored_obs)]
df=df[~np.isnan(df.delay_lo)]


time=df.MJD_START

ax.fill_between(time,df.delay_lo,df.delay_hi,color='c',alpha=0.6)

ax.set_ylabel('Iron line equivalent width delay, sec',color='k')
ax.set_xlabel('Time, MJD')


# ax.axvspan(53340.29,53360.00,alpha=0.05, color='gray')
# ax.axvspan(53384.36,53428.51,alpha=0.05, color='gray')
# ax.axvspan(53380,53380.8,alpha=0.05,color='gray')

orbtime=np.linspace(53335,53420,200)
r,_,_,_,_=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)


ax2=ax.twinx()


ax2.plot(orbtime,r,'b-.',alpha=0.6,lw=0.5)


ax2.set_ylabel('Distance to the companion, \n lt-sec.')

fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'mc_delay_err.png',dpi=500)
plt.show()



#%% plot delay stuff vs flux

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

df=ObsParams[(ObsParams.MJD_START<53400) | (ObsParams.EXPOSURE>1500)]
deltat,deltat_err=df.deltat,df.deltat_err


flux,flux_err=vals_and_errors(df,'cutoffpl_tot_flux',funct=lambda x: x/1e-8)

nu_412_flux,nu_412_flux_max_err=0.83,0.02#1.220422,1.220422-1.219378
nu_deltat,nu_deltat_err=1.3,0.3

ax.errorbar(flux,deltat,deltat_err,flux_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.errorbar(nu_412_flux,nu_deltat,nu_deltat_err,nu_412_flux_max_err,fmt='.',color='r',marker='s',ms=4,alpha=0.8,label='NuSTAR')


ax.set_ylabel('Iron line equivalent width delay, sec',color='k')
ax.set_xlabel('Flux 7-12')
ax.legend()


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'delay_vs_flux.png',dpi=500)
plt.show()


'''



#%% OLD


'''
#%% spe_pars_stuff


fig = plt.figure(figsize=(6.6*4, 4*6.6/2))
rows=5
cols=5
ax = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=5)

ax2=plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=5)

ax3=plt.subplot2grid((rows,cols), (4, 0), rowspan=1, colspan=2)


time=ObsParams.MJD_START

#for model,color in zip(['cutoffpl','pohi'],['r','g']):
for model,color in zip(['cutoffpl'],['c']):

    eqw,eqw_err=vals_and_errors(ObsParams,model+'_eqw')
    eqw=eqw*1000
    eqw_err=eqw_err*1000


    ecol='k'
    ax.plot(time,eqw,marker='s',mfc=color,mec=ecol,mew=1,ls='None',label=model+'_eqw',alpha=0.6)
    ax.errorbar(time,eqw,eqw_err,ecolor=ecol,fmt='none',alpha=0.5)

    #ss=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/rxte/plots_results/SSandAA.csv')
    #ax.errorbar(ss[:,0],ss[:,1]*1e3*0.8,ss[:,1]/ss[:,1]*3*0.8,color='r',fmt='none',label='SS paper * 80\% , approx errors 3 eV',alpha=0.4)


    efold=ObsParams[model+'_efold']
    ecut=ObsParams[model+'_ecut']
    chi2=ObsParams[model+'_chi2']
    gamma=ObsParams[model+'_po']
    ax2.plot(time,efold,marker='d',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6,label=model+'_efold')
    ax2.plot(time,ecut,marker='.',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6,label=model+'_ecut')
    ax2_twin=ax2.twinx()
    ax2_twin.time,gamma,marker='d',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6,label=model+'_po')

    chi2.hist(ax=ax3,color=color,alpha=0.5,label=model+'_chi2')
    #ObsParams.iloc[where(ObsParams.cutoffpl_chi2>2)[0]][['ObsID','MJD_START','cutoffpl_chi2']]
ax.set_title(filename)
ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best',fontsize=7)
ax.grid()

ax2.set_xlabel('Time, MJD')
ax2.set_ylabel('E_fold/E_cut, keV')
ax2.legend(loc='best',fontsize=7)
ax2.grid()

ax3.legend()

plt.show()
#%%

ax.axvspan(53340.29,53360.00,alpha=0.1, color='blue',label='Phase spectroscopy available')
ax.axvspan(53384.36,53428.51,alpha=0.1, color='blue')
ax.axvspan(53380,53380.8,alpha=0.1,color='blue')


ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best',fontsize=7)
ax.grid()
plt.show()
#fig.savefig(results_path+f'{filename}_eqw.png')




#
##%% plot eqw stuff
#matplotlib.rcParams['figure.figsize'] = 2*6.6, 6.6
#fig,ax=plt.subplots()
#time=ObsParams['MJD_START']
#eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw')
#eqw=eqw*1000
#eqw_err=eqw_err*1000
#
#
##======= my data =======
#
#ecol='k'
#ax.plot(time,eqw,marker='s',mfc='green',mec=ecol,mew=1,ls='None',label=filename)
#ax.errorbar(time,eqw,eqw_err,ecolor=ecol,fmt='none',alpha=0.7)
#
#
#
##======= their data =====
#
#ss=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/rxte/plots_results/SSandAA.csv')
#ax.errorbar(ss[:,0],ss[:,1]*1e3,ss[:,1]/ss[:,1]*3,color='r',fmt='none',label='SS paper , approx errors 3 eV',alpha=0.4)
#
#
#og=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/rxte/plots_results/Orkhan.csv')
#ax.errorbar(og[:,0],og[:,1]*1e3,og[:,1]/og[:,1]*5,color='k',fmt='none',label='OG arrpox errors 5 eV',alpha=0.4)
#
#
#
#ax.axvspan(53340.29,53360.00,alpha=0.1, color='blue',label='Phase spectroscopy available')
#ax.axvspan(53384.36,53428.51,alpha=0.1, color='blue')
#ax.axvspan(53380,53380.8,alpha=0.1,color='blue')
#
#
#ax.set_xlabel('Time, MJD')
#ax.set_ylabel('Equivalent width of iron line, eV')
#ax.legend(loc='best',fontsize=7)
#ax.grid()
#plt.show()
##fig.savefig(results_path+f'{filename}_eqw.png')



#%% plot flux stuff

fig,ax=plt.subplots()
time=ObsParams['MJD_START']
flux,flux_err=vals_and_errors(ObsParams,'cutoffpl_tot_flux')


#======= my data =======

ecol='k'
ax.plot(time,flux,marker='s',mfc='blue',mec=ecol,mew=1,ls='None',label=filename)
ax.errorbar(time,flux,flux_err,ecolor=ecol,fmt='none',alpha=0.7)

# ===== fit with gauss ===


#def gauss(t,t0,sigma,N):
#    return N*np.exp(-(t-t0)**2/(2*sigma**2))
#
#
#from scipy.optimize import curve_fit
#
#popt,pcov=curve_fit(gauss,time[~np.isnan(flux)],flux[~np.isnan(flux)],p0=[53370,20,3e-8])
#
#ax.plot(time,gauss(time,*popt),'k:',label='')



ax.axvspan(53340.29,53360.00,alpha=0.1, color='blue',label='Phase spectroscopy available')
ax.axvspan(53384.36,53428.51,alpha=0.1, color='blue')
ax.axvspan(53380,53380.8,alpha=0.1,color='blue')


ax.set_xlabel('Time, MJD')
ax.set_ylabel('Flux (3-12 keV), erg/s/cm^2')
ax.legend(loc='best')
#ax.grid()
plt.show()
fig.savefig(results_path+f'{filename}_flux.png')




#%% plot period stuff

msv_res=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/rxte/plots_results/msv_res.txt')
msv_time=msv_res[:,0]
msv_per=msv_res[:,4]
msv_fl=msv_res[:,6]


ObsParams.plot(x='MJD_START',y='period',yerr='period_err',kind='scatter',label='my')
plt.plot(msv_time,msv_per,'k+',label='MSV')


plt.axvspan(53340.29,53360.00,alpha=0.1, color='blue',label='Phase spectroscopy available')
plt.axvspan(53384.36,53428.51,alpha=0.1, color='blue')
plt.axvspan(53380,53380.8,alpha=0.1,color='blue')


my_per=ObsParams.period
mytime=ObsParams.MJD_START
_,_,_,r,_=doppler.kepler_solution(mytime*day2sec,doppler.orb_params_v0332)

plt.plot(mytime,my_per/r,'m.',label='corrected per')
plt.show()
plt.legend()



#%% plot period stuff orb_corr

msv_res=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/rxte/plots_results/msv_res.txt')
msv_time=msv_res[:,0]
msv_per=msv_res[:,4]
msv_fl=msv_res[:,6]


ObsParams.plot(x='MJD_START',y='period_orb_corr',kind='scatter',label='my')
plt.plot(msv_time,msv_per,'k+',label='MSV')


plt.axvspan(53340.29,53360.00,alpha=0.1, color='blue',label='Phase spectroscopy available')
plt.axvspan(53384.36,53428.51,alpha=0.1, color='blue')
plt.axvspan(53380,53380.8,alpha=0.1,color='blue')


my_per=ObsParams.period_orb_corr
mytime=ObsParams.MJD_START
_,_,_,r,_=doppler.kepler_solution(mytime*day2sec,doppler.orb_params_v0332)

plt.plot(mytime,my_per*r,'m.',label='corrected per')
plt.show()
plt.legend()


#%% plot eqw and inv distance

fig,ax=plt.subplots()
ax2=ax.twinx()
model='cutoffpl'
time=ObsParams.MJD_START
eqw,eqw_err=vals_and_errors(ObsParams,model+'_eqw')
eqw=eqw*1000
eqw_err=eqw_err*1000
ecol='k'
ax.plot(time,eqw,marker='s',mfc='green',mec=ecol,mew=1,ls='None',label=filename)
ax.errorbar(time,eqw,eqw_err,ecolor=ecol,fmt='none',alpha=0.7)



orbtime=np.linspace(53335,53440,200)
_,_,_,r=apply_doppler(orbtime*day2sec,orbtime,orbtime,P,T_p,e,w,asini)


inv=r

ax2.plot(orbtime,inv,'b-.')


ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax2.set_ylabel('Distance from the NS to the companion, lt-sec.')
ax.legend(loc='best')
ax.grid()
plt.show()
fig.savefig(results_path+f'{filename}_eqw_and_dist.png')




#%% plot delat and inv distance

fig,ax=plt.subplots()
ax2=ax.twinx()
model='cutoffpl'
time=ObsParams.MJD_START
delay=ObsParams.wrap_deltat
delay_err=ObsParams.deltat_err
ecol='k'
ax.plot(time,delay,marker='s',mfc='green',mec=ecol,mew=1,ls='None',label=filename)
ax.errorbar(time,delay,delay_err,ecolor=ecol,fmt='none',alpha=0.7)



orbtime=np.linspace(53335,53440,200)
_,_,_,r=apply_doppler(orbtime*day2sec,orbtime,orbtime,P,T_p,e,w,asini)


inv=r

ax2.plot(orbtime,inv,'b-.')


ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax2.set_ylabel('Distance from the NS to the companion, lt-sec.')
ax.legend(loc='best')
ax.grid()
plt.show()
fig.savefig(results_path+f'{filename}_eqw_and_dist.png')


#%% plot eqw convolved with  inv distance

fig,ax=plt.subplots()



ecol='k'
ax.plot(time,eqw,marker='s',mfc='green',mec=ecol,mew=1,ls='None',label=filename)
ax.errorbar(time,eqw,eqw_err,ecolor=ecol,fmt='none',alpha=0.7)



orbtime=np.linspace(53320,53440,200)
_,_,factor,r=apply_doppler(orbtime*day2sec,orbtime,orbtime,P,T_p,e,w,asini)


gaussian_depend=gauss(orbtime,popt[0],popt[1]*2,N=90)

#gaussian_depend_conv=gaussian_depend*factor

#inv=r**(-1)*np.mean(r)
#y=gaussian_depend+30*(inv-mean(inv))

ax.plot(orbtime,gaussian_depend,'b-.')
ax.plot(orbtime,gaussian_depend_conv,'b:')

#ax.plot(orbtime,y,'b:')

ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width of iron line, eV')
ax.legend(loc='best')
ax.grid()
plt.show()
#fig.savefig(results_path+f'{filename}_eqw_and_dist_conv.png')







#%% make latex table

def error(val,low,hi,funct=lambda x:x ,k=2):
    low=funct(val-low)
    hi=funct(hi-val)
    val=funct(val)

    strval=f"%.{k}f" % val
    strlow=f"%.{k}f" % (low)
    strhi=f"%.{k}f" % (hi)

    y=('$'+strval+'_{-'+strlow+'}^{+'+strhi+'}$')


    if strlow==strhi:
        y=('$'+strval+'\pm'+strlow+'$')
        if low==0 and hi==0 :
            y=('$'+strval+'$')

    if y=='$0.0\pm0.0$':
        return '-'
    if y=='$0.0$':
        return '-'
    return y

def latex_error(row,par,funct=lambda x:x,k=3):
    err=error(row[par],row[par+'_lo'],row[par+'_hi'],funct=funct,k=k)
    return err

#df['latex_error_f']=df.apply(lambda x: latex_error(x,'fe_flux',funct=lambda x: x/1e-8),axis=1)
#df['latex_flux']=df.apply(lambda x: error(x['fe_flux'],x['fe_flux_lo'],x['fe_flux_hi'],funct=lambda x: x/1e-9),axis=1)

#def unique_values(myList):
#    output = []
#    for x in myList:
#        if x not in output:
#            output.append(x)
#    return output

def frm2(x):
    return f"%.2f" % x

def frm3(x):
    return f"%.3f" % x


temp=ObsParams.sort_values('MJD_START')
#for colname in ['SE_config','SA_config']:
#    temp[colname]=temp[colname].apply(lambda x: list(set(x)))
#    temp[colname]=temp[colname].apply(lambda x: 'None' if len(x)==0 else x[0])
temp['cutoffpl_tot_flux']=temp['cutoffpl_tot_flux'].apply(lambda x: x/1e-8)

temp.to_latex(buf=results_path+f'latex_{filename}.txt',
              columns=['ObsID','MJD_START','SA_config', 'SE_config','cutoffpl_tot_flux','cutoffpl_chi2'],
              formatters=[lambda x: x ,frm2 ,lambda x: x,lambda x: x,frm2,frm2],
              index=False,escape=1)


'''
