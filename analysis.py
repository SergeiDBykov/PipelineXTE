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

from Miscellaneous import  doppler_correction as doppler
from Miscellaneous.doppler_correction import  day2sec


font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 11}


matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['xtick.top']=1
matplotlib.rcParams['ytick.right']=1

matplotlib.rcParams['figure.figsize'] = 6.6*4, 6.6*2
matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.85
matplotlib.rcParams['figure.subplot.top']=0.90
#matplotlib.rcParams['text.usetex'] = True
#matplotlib.rcParams['text.latex.unicode'] = True

plt.ion()

def vals_and_errors(ObsParams,name,funct=lambda x: x):
    val=ObsParams[name].values
    hi=ObsParams[name+'_hi'].values-val
    try:
        low=val-ObsParams[name+'_lo'].values
    except: #just stuoed error in var names in spe_info.txt
        low=val-ObsParams[name+'_low'].values

    err=np.vstack((low,hi))
    return val,err




results_path='/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pandas_data/'

filename='standard_pipeline_orb_corr_periods' # standard_pipeline    standard_pipeline_free_sigma
# standard_pipeline standard_pipeline_025sigma standard_pipeline_ign_12 standard_pipeline_sts_period
#standard_pipeline_orb_corr_periods
ObsParams=pd.read_pickle(results_path+f'{filename}.pkl')
ObsParams=ObsParams.sort_values('MJD_START')

ObsParams.period_orb_corr= ObsParams.period_orb_corr.replace(to_replace='None',value=np.nan)
ObsParams.period_orb_corr_err= ObsParams.period_orb_corr_err.replace(to_replace='None',value=np.nan)


STOP

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

    ss=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/rxte/plots_results/SSandAA.csv')
    ax.errorbar(ss[:,0],ss[:,1]*1e3*0.8,ss[:,1]/ss[:,1]*3*0.8,color='r',fmt='none',label='SS paper * 80\% , approx errors 3 eV',alpha=0.4)


    efold=ObsParams[model+'_efold']
    ecut=ObsParams[model+'_ecut']
    chi2=ObsParams[model+'_chi2']
    ax2.plot(time,efold,marker='d',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6,label=model+'_efold')
    ax2.plot(time,ecut,marker='.',mfc=color,mec=ecol,mew=1,ls='None',alpha=0.6,label=model+'_ecut')


    chi2.hist(ax=ax3,color=color,alpha=0.5,label=model+'_chi2')
    #ObsParams.iloc[where(ObsParams.cutoffpl_chi2>2)[0]][['ObsID','MJD_START','cutoffpl_chi2']]
ax.set_title(filename)
ax.set_xlabel('Time, MJD')
ax.set_ylabel('Equivalent width, eV')
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
ax.set_ylabel('Equivalent width, eV')
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
#ax.set_ylabel('Equivalent width, eV')
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
ax.set_ylabel('Equivalent width, eV')
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
ax.set_ylabel('Equivalent width, eV')
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
ax.set_ylabel('Equivalent width, eV')
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

