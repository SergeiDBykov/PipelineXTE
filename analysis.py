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

from Misc import  doppler_correction as doppler
from Misc.doppler_correction import  day2sec


matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.85
matplotlib.rcParams['figure.subplot.top']=0.95


import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')

plt.ion()

def vals_and_errors(ObsParams,name,funct=lambda x: x):
    if isinstance(ObsParams,pd.core.frame.DataFrame):
        par,Min,Max=funct(ObsParams[name].values),funct(ObsParams[name+'_lo'].values),funct(ObsParams[name+'_hi'].values)
    elif isinstance(ObsParams,pd.core.series.Series):
        par,Min,Max=funct(ObsParams[name]),funct(ObsParams[name+'_lo']),funct(ObsParams[name+'_hi'])

    low = par - Min
    hi  =  Max - par
    parr = par

    err=np.vstack((low,hi))
    #return np.array((parr,err))
    return parr,err




savepath='/Users/s.bykov/work/xray_pulsars/rxte/plots_results/'
results_path='/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pandas_data/'

filename='standard_pipeline'
ObsParams=pd.read_pickle(results_path+f'{filename}.pkl')
ObsParams=ObsParams.sort_values('MJD_START')

ObsParams.period_orb_corr= ObsParams.period_orb_corr.replace(to_replace='None',value=np.nan)
ObsParams.period_orb_corr_err= ObsParams.period_orb_corr_err.replace(to_replace='None',value=np.nan)


ObsParams.loc[ObsParams.fasebin_cfg=='se','config']='E\_125us\_64M\_0\_1s'
ObsParams.loc[ObsParams.fasebin_cfg=='sa','config']='B\_16ms\_46M\_0\_49\_H'
ObsParams.loc[ObsParams.fasebin_cfg=='None','config']='-'
#ObsParams[ObsParams.fasebin_cfg=='se']['config']='E\_125us\_64M\_0\_1s'
#ObsParams['config']='E\_125us\_64M\_0\_1s' if ObsParams.fasebin_cfg=='se'
#B\_16ms\_46M\_0\_49\_H
STOP

#%% fun with correlations

corrdata=ObsParams[['cutoffpl_efold','cutoffpl_eline','cutoffpl_eqw',
                    'cutoffpl_fe_flux','cutoffpl_po','cutoffpl_tot_flux',
                   'deltat' ]]

sns.heatmap(corrdata.corr(),cmap='jet')

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
def gauss(t,t0,sigma,N):
    return N*np.exp(-(t-t0)**2/(2*sigma**2))

from scipy.optimize import curve_fit

popt,pcov=curve_fit(gauss,time[time<53385],flux[time<53385],p0=[53370,20,3])

ax.plot(time,gauss(time,*popt),'k:',label='')

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





#%% plot eqw vs time stuff

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw',funct=lambda x: 1e3*x)

ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,label='Sigma=0.3 keV',alpha=0.8)


ax.set_ylabel('Iron line equivalent width \n eV',color='k')
ax.set_xlabel('Time, MJD')
ax.set_ylim(40,120)




from Misc.doppler_correction import orb_params_v0332
T_p=orb_params_v0332['T_p']/86400
P=orb_params_v0332['P']/86400


for i in range(1,4):
    ax.axvline(T_p+(i-114)*P,color='k',ls=':',alpha=0.3)




ax.axvspan(53340.29,53360.00,alpha=0.05, color='gray')
ax.axvspan(53384.36,53428.51,alpha=0.05, color='gray')
ax.axvspan(53380,53380.8,alpha=0.05,color='gray')


fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw.png',dpi=500)
plt.show()


#%% plot eqw vs flux

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})


eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw',funct=lambda x: 1e3*x)
flux,flux_err=vals_and_errors(ObsParams,'cutoffpl_tot_flux',funct=lambda x: x/1e-8)

ax.errorbar(flux,eqw,eqw_err,flux_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line equivalent width \n eV',color='k')
ax.set_xlabel('RXTE/PCA Flux (3-12 keV), \n 10^(-8) cgs')
ax.set_ylim(40,120)




fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw_vs_flux.png',dpi=500)
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

#%% plot eqw vs time vs distance

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

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

#%% plot eqw vs time vs flux

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

time=ObsParams.MJD_START

eqw,eqw_err=vals_and_errors(ObsParams,'cutoffpl_eqw',funct=lambda x: 1e3*x)
ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line equivalent width \n eV',color='k')

ax.set_xlabel('Time, MJD')

ax.set_ylim(50,110)

ax.axvspan(53340.29,53360.00,alpha=0.05, color='gray')
ax.axvspan(53384.36,53428.51,alpha=0.05, color='gray')
#ax.axvspan(53380,53380.8,alpha=0.05,color='gray')



ax2=ax.twinx()

flux,flux_err=vals_and_errors(ObsParams,'cutoffpl_tot_flux',funct=lambda x: x/1e-8)

ax2.errorbar(time,flux,flux_err,fmt='.',color='g',marker='s',ms=4,alpha=0.8)

ax2.set_ylabel('RXTE/PCA Flux (3-12 keV), \n 10^(-8) cgs',color='b')




fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'eqw_vs_flux_in_time.png',dpi=500)
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
ignored_obs=[
'90014-01-05-03',
'90014-01-05-06',
'90014-01-06-00',
'90014-01-07-01',
'90014-01-07-03',
'90014-01-07-02',
'90014-01-07-00',
'90014-01-08-00',
'90014-01-08-01',
'90014-01-06-03']

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
r,_,_,_,_=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)


ax2=ax.twinx()


ax2.plot(orbtime,r,'b-.',alpha=0.6,lw=0.5)


ax2.set_ylabel('Distance to the companion, \n lt-sec.')

fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'delay.pdf',dpi=500)
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



#%% latex table
from Misc.TeX_Tables import pandas_to_tex
from Misc.TeX_Tables.pandas_to_tex import *

model='cutoffpl_'
def tex_tbl():
    null=lambda x: x
    free_columns=['ObsID','MJD_START','EXPOSURE','config',model+'chi2']
    free_columns_functions=[null,null,null,null,null]
    free_columns_formats=[0,1,0,0,2]

    err_columns=['eqw','tot_flux',
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
