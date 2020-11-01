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




savepath='/Users/s.bykov/work/xray_pulsars/rxte/plots_results/mean_spe/'
results_path='/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pandas_data/'

filename='standard_pipeline' #standard_pipeline standard_pipeline_edge_en_free standard_pipeline_with_phabs
ObsParams=pd.read_pickle(results_path+f'{filename}.pkl')
ObsParams=ObsParams.sort_values('MJD_START')

ObsParams.period_orb_corr= ObsParams.period_orb_corr.replace(to_replace='None',value=np.nan)
ObsParams.period_orb_corr_err= ObsParams.period_orb_corr_err.replace(to_replace='None',value=np.nan)


ObsParams.loc[ObsParams.fasebin_cfg=='se','config']='E\_125us\_64M\_0\_1s'
ObsParams.loc[ObsParams.fasebin_cfg=='sa','config']='B\_16ms\_46M\_0\_49\_H'
ObsParams.loc[ObsParams.fasebin_cfg=='None','config']='-'

#%% assign groups
gr1=['90089-11-03-00G','90089-11-03-01G','90089-11-03-02']
gr2=['90427-01-03-00','90427-01-03-01','90427-01-03-02']
gr3=['90427-01-03-14G','90014-01-02-00','90427-01-03-05']
gr4=['90427-01-03-06','90427-01-03-07','90014-01-02-08']
gr5=['90427-01-03-09','90427-01-03-11','90427-01-03-12']
gr6=['90014-01-02-13','90014-01-03-00','90014-01-03-01']
gr7=['90014-01-03-020','90014-01-03-02','90014-01-03-03']


def set_group(row):
    if row['ObsID'] in gr1:
        return 'i'
    elif row['ObsID'] in gr2:
        return 'ii'
    elif row['ObsID'] in gr3:
        return 'iii'
    elif row['ObsID'] in gr4:
        return 'iv'
    elif row['ObsID'] in gr5:
        return 'v'
    elif row['ObsID'] in gr6:
        return 'vi'
    elif row['ObsID'] in gr7:
        return 'vii'
    else:
        return '-'



ObsParams['Gr']=ObsParams.apply(set_group,axis=1)
# ObsParams['Gr']=['i' if x in gr1 else pass for x in ObsParams['ObsID']]
#  #['i' if x in gr1 else '-' for x in ObsParams['ObsID']]
# ObsParams['Gr']=['ii' if x in gr2   else pass  for x in ObsParams['ObsID']]
# ObsParams['Gr']=['iii' if x in gr3  else pass  for x in ObsParams['ObsID']]
# ObsParams['Gr']=['iv' if x in gr4  else pass  for x in ObsParams['ObsID']]
# ObsParams['Gr']=['v' if x in gr5  else pass  for x in ObsParams['ObsID']]
# ObsParams['Gr']=['vi' if x in gr6  else pass for x in ObsParams['ObsID']]
# ObsParams['Gr']=['vii' if x in gr7  else pass for x in ObsParams['ObsID']]


MJD_REF=53000

STOP


#%% plot line cutoffpl flux, line norm

for model in ['edge_cutoffpl','cutoffpl']:

    fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

    time=ObsParams.MJD_START-MJD_REF

    flux,flux_err=vals_and_errors(ObsParams,model+'_cutoffpl_flux',funct=lambda x: x/1e-8)

    ax.errorbar(time,flux,flux_err,fmt='.',color='b',marker='o',ms=4,alpha=0.8)

    ax.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b')
    ax.set_xlabel(f'Time, MJD-{MJD_REF}')
    ax.set_yscale('log')

    ax.set_ylim(0.1,5)

    ax_iron=ax.twinx()
    flux_iron,flux_err_iron=vals_and_errors(ObsParams,model+'_norm_line',funct=lambda x: x/1e-3)

    ax_iron.errorbar(time,flux_iron,flux_err_iron,fmt='.',color='g',marker='o',ms=4,alpha=0.6,zorder=-1)


    ax_iron.set_ylabel('Iron line norm, \n $10^{3}$ ph $cm^{-2}$ $s^{-1}$',color='g')
    ax_iron.set_yscale('log')


    ax.axvspan(53340.29-MJD_REF,53360.00-MJD_REF,alpha=0.05, color='gray')
    ax.axvspan(53384.36-MJD_REF,53428.51-MJD_REF,alpha=0.05, color='gray')
    #ax.axvspan(53380,53380.8,alpha=0.05,color='gray')


    fig.tight_layout()
    plt.savefig(savepath+f'flux_cont_and_iron_{model}.pdf',dpi=500)
    plt.savefig(savepath+f'flux_cont_and_iron_{model}.png',dpi=500)




fig, ax_iron = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

# flux,flux_err=vals_and_errors(ObsParams,'edge_cutoffpl'+'_cutoffpl_flux',funct=lambda x: x/1e-8)
# ax.errorbar(time,flux,flux_err,fmt='.',color='g',marker='s',ms=4,alpha=0.8)
# ax.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='g')
# ax.set_xlabel(f'Time, MJD-{MJD_REF}')
# ax.set_yscale('log')

# ax.set_ylim(0.1,5)
# ax_iron=ax.twinx()
ax_iron.set_ylabel('Iron line norm, \n $10^{3}$ ph $cm^{-2}$ $s^{-1}$',color='k')
ax_iron.set_yscale('log')
ax.axvspan(53340.29-MJD_REF,53360.00-MJD_REF,alpha=0.05, color='gray')
ax.axvspan(53384.36-MJD_REF,53428.51-MJD_REF,alpha=0.05, color='gray')

for model in ['edge_cutoffpl','cutoffpl']:

    flux_iron,flux_err_iron=vals_and_errors(ObsParams,model+'_norm_line',funct=lambda x: x/1e-3)
    ax_iron.errorbar(time,flux_iron,flux_err_iron,fmt='.',marker='s',ms=4,alpha=0.6,zorder=-1,label=model)
    fig.tight_layout()

ax_iron.legend()
plt.savefig(savepath+f'flux_cont_and_iron.pdf',dpi=500)
plt.savefig(savepath+f'flux_cont_and_iron.png',dpi=500)





#%% plot line cutoffpl flux, eqw, separation


for model in ['edge_cutoffpl','cutoffpl']:

    fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=[8, 8/2])

    time=ObsParams.MJD_START-MJD_REF

    eqw,eqw_err=vals_and_errors(ObsParams,model+'_eqw',funct=lambda x: x*1000)

    ax.errorbar(time,eqw,eqw_err,fmt='.',color='turquoise',marker='o',ms=4,alpha=0.8)

    ax.set_ylabel('Iron line \n eq. width, eV',color='turquoise')
    ax.set_xlabel(f'Time, MJD-{MJD_REF}')
    #ax.set_yscale('log')

    ax.set_ylim(40,110)

    ax_flux=ax.twinx()
    flux,flux_err=vals_and_errors(ObsParams,model+'_cutoffpl_flux',funct=lambda x: x/1e-8)

    ax_flux.errorbar(time,flux,flux_err,fmt='.',color='gray',marker='o',ms=4,alpha=0.4)


    ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='gray')
    ax_flux.set_yscale('log')


    orbtime=np.linspace(53335,53430,200)
    r,_,_,_,_=doppler.kepler_solution(orbtime*day2sec,doppler.orb_params_v0332)

    ax2=ax.twinx()
    ax2.plot(orbtime-MJD_REF,r,'r-.',alpha=0.6,lw=1,zorder=0)
    ax2.set_ylabel('Distance to the companion, \n lt-sec.',color='r')
    ax2.spines["right"].set_position(("axes", +1.15))

    fig.tight_layout()
    plt.savefig(savepath+f'eqw_{model}.pdf',dpi=500)
    plt.savefig(savepath+f'eqw_{model}.png',dpi=500)





fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})


ax.set_ylabel('Iron line \n eq. width, eV',color='k')
ax.set_xlabel(f'Time, MJD-{MJD_REF}')

for model in ['edge_cutoffpl','cutoffpl']:
    eqw,eqw_err=vals_and_errors(ObsParams,model+'_eqw',funct=lambda x: x*1000)
    ax.errorbar(time,eqw,eqw_err,fmt='.',marker='o',ms=4,alpha=0.8,label=model)
    fig.tight_layout()

ax.legend()
plt.savefig(savepath+f'eqw.pdf',dpi=500)
plt.savefig(savepath+f'eqw.png',dpi=500)




#%% plot line K edge depth


fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=[8, 8/2])

time=ObsParams.MJD_START-MJD_REF

tau,tau_err=vals_and_errors(ObsParams,'edge_cutoffpl_edgeTau',funct=lambda x: x*100)

ax.errorbar(time,tau,tau_err,fmt='.',color='r',marker='o',ms=4,alpha=0.8,uplims=tau_err[1]==0)

ax.set_ylabel('K-edge $\\tau$, 10$^{-2}$',color='r')
ax.set_xlabel(f'Time, MJD-{MJD_REF}')
#ax.set_yscale('log')

ax.set_ylim(-0.2,5)

ax_flux=ax.twinx()
flux,flux_err=vals_and_errors(ObsParams,model+'_cutoffpl_flux',funct=lambda x: x/1e-8)

ax_flux.errorbar(time,flux,flux_err,fmt='.',color='gray',marker='o',ms=4,alpha=0.4)


ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='gray')
ax_flux.set_yscale('log')


fig.tight_layout()
plt.savefig(savepath+f'tau_edge_cutoffpl.pdf',dpi=500)
plt.savefig(savepath+f'tau_edge_cutoffpl.png',dpi=500)






#%% plot pulsed emission
fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=[8, 8/2])

time=ObsParams.MJD_START-MJD_REF

iron_puls,iron_puls_err=ObsParams.cutoffpl_rel_dip_iron.values,ObsParams.cutoffpl_rel_dip_iron_err.values
iron_puls -=1
iron_puls[-4]=np.nan
iron_puls_err[-4]=np.nan

ax.plot(time,iron_puls,'mo',ls='None',alpha=0.6,label='Iron')
ax.errorbar(time,iron_puls,iron_puls_err,ecolor='m',fmt='none',alpha=0.5)



cont_puls,cont_puls_err=ObsParams.cutoffpl_rel_dip_flux.values,ObsParams.cutoffpl_rel_dip_flux_err.values
cont_puls -=1

ax.plot(time,cont_puls,'co',ls='None',alpha=0.6,label='Cont. (3-12 keV)')
ax.errorbar(time,cont_puls,cont_puls_err,ecolor='c',fmt='none',alpha=0.5)


#ax.set_ylim(0.9-1,4-1)
ax.set_yscale('log')
ax.set_ylabel('$\\frac{I_{max}}{I_{min}}$-1 \n cutoffpl+gauss model',fontsize=12)
ax.set_xlabel(f'Time, MJD-{MJD_REF}')
#ax.grid()
ax.legend()

fig.savefig(savepath+f'cutoffpl_pulsed_fractions.png')







#%%plot spectra stuff


for parname in ['po','ecut','chi2']:
    fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})
    for model in ['edge_cutoffpl','cutoffpl']:
        time=ObsParams.MJD_START-MJD_REF
        if parname!='chi2':
            par,par_err=vals_and_errors(ObsParams,model+'_'+parname)
        else:
            par=ObsParams[f'{model}_chi2']
            par_err=1e-2
        ax.errorbar(time,par,par_err,fmt='.',marker='o',ms=6,alpha=0.8,label=f'{parname}, mo: {model}')

        ax.set_ylabel(parname,color='k')
        ax.set_xlabel(f'Time, MJD-{MJD_REF}')

    ax.legend()
    fig.tight_layout()
    plt.savefig(savepath+f'{parname}.png',dpi=500)



fig,ax=plt.subplots()
for model in ['edge_cutoffpl','cutoffpl']:
    chi2=ObsParams[model+'_chi2']
    ax.hist(chi2,alpha=0.5,bins=25,label=model+f', mean={"%.2f" %np.mean(chi2)}')
    ax.set_xlabel('chi2_red')
    plt.show()
ax.legend()
plt.savefig(savepath+f'chi2_hist.png',dpi=250)



#%% plot residuals before and after edge (and gauss),  delchi

print(ObsParams.cutoffpl_chi2-ObsParams.edge_cutoffpl_chi2>0.5)

fig = plt.figure(figsize=(6.6, 6.6))
plt.subplots_adjust(hspace=0)
rows=4
cols=3

ax_del_no_gauss=plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_del_no_edge = plt.subplot2grid((rows,cols), (2, 0), rowspan=1, colspan=3,sharex=ax_del_no_gauss)
ax_del_edge=plt.subplot2grid((rows,cols), (3, 0), rowspan=1, colspan=3,sharex=ax_del_no_gauss)

for ObsID,color in zip(['90089-11-05-08G','90089-11-04-02G','90089-22-01-00G','90427-01-03-14G'],['b','g','r','m']):
        label=ObsID
        del_no_gauss=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2_top/cutoffpl_no_gauss/mean_spe_del.dat')
        ax_del_no_gauss.errorbar(del_no_gauss[0],del_no_gauss[1],del_no_gauss[2],del_no_gauss[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6,color=color)

        del1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2_top/cutoffpl/mean_spe_del.dat')
        ax_del_no_edge.errorbar(del1[0],del1[1],del1[2],del1[3],drawstyle='steps-mid',ls=':',alpha=0.6,color=color)

        del2=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2_top/edge_cutoffpl/mean_spe_del.dat')
        ax_del_edge.errorbar(del2[0],del2[1],del2[2],del2[3],drawstyle='steps-mid',ls=':',alpha=0.6,color=color)

ax_del_no_gauss.legend()
ax_del_no_gauss.set_ylabel('$\Delta \chi$ \n cutoffpl',fontsize=8)
ax_del_no_edge.set_ylabel('$\Delta \chi$ \n cutoffpl+gauss',fontsize=8)
ax_del_no_edge.grid('y')
ax_del_edge.set_ylabel('$\Delta \chi$ \n (cutoffpl+gauss)*edge',fontsize=8)
ax_del_edge.grid('y')
ax_del_edge.set_xlabel('Energy, keV')
ax_del_no_gauss.set_xscale('log')
#fig.tight_layout()
plt.savefig(savepath+f'edge_residuals.png',dpi=250)





#%% latex table
from Misc.TeX_Tables import pandas_to_tex
from Misc.TeX_Tables.pandas_to_tex import *

for model in ['cutoffpl_','edge_cutoffpl_']:
    def tex_tbl():
        null=lambda x: x
        free_columns=['Gr','ObsID','MJD_START','EXPOSURE','config',model+'chi2']
        free_columns_functions=[null,null,null,null,null,null]
        free_columns_formats=[0,0,1,0,0,2]

        err_columns=['eqw','cutoffpl_flux',
                     ]
        err_functions=[lambda x: 1000*x, lambda x: x/1e-9,
                       ]
        err_formats=[0,2]

        err_columns=[model+item for item in err_columns]

        headers=['Gr.','ObsID','Time, MJD','Exposure, s','Configuration','$\chi^2_{red}$',
                      'Eq. width', 'Flux (3-12 keV)'
                     ]

        transpose=0
        df_tex=make_latex_table(df,
                              free_columns=free_columns, free_columns_functions=free_columns_functions,
                              free_columns_formats=free_columns_formats,
                              err_columns=err_columns, err_functions=err_functions,
                              err_formats=err_formats)
        df_tex.columns=headers
        save_latex_table(df_tex, savepath=savepath+'/tex/allpars_'+model+'.tex',
                         columns_to_write='DEFAULT',
                         columns_names='DEFAULT',
                         transpose=transpose)
    df=ObsParams
    tex_tbl()




#%% astropy units: calculating tau
import astropy
from astropy import units as u
from astropy.units import cds
cds.enable()
k_e=0.34*(u.cm)**2/u.g

def Mdot(L38,R6=1,Msun=1.4):
    tmp=((L38*10**38*(u.erg/u.s))*(R6*10**6*u.cm))/((Msun*u.M_sun)*(u.cds.G)).cgs
    return tmp.to(u.g/u.s)


def R_M(L38,Msun=1.4,R6=1,B12=2.5,L=0.5):
    #tmp=L*(Msun*u.M_sun)**(1/7)*(R6*10**6*u.cm)**(10/7)*(B12*10**12*u.G)**(4/7)*(0.1*L38*10**38*(u.erg/u.s))**(-2/7)
    tmp=7*10**(7)*u.cm*L*(Msun)**(1/7)*(R6)**(10/7)*(B12)**(4/7)*(0.1*L38)**(-2/7)
    return tmp

def vff(R,Msun=1.4):
    tmp=(2*u.cds.G*Msun*u.Msun/R)**(1/2)
    return tmp.cgs


def tau_fe_exp(L38,R6=1,Msun=1.4,B12=2.5,dphi=1/3,L=0.5,Z_fe=1):
    rm=R_M(Msun,R6,B12,L38,L=L)
    tmp=(Mdot(L38,Msun,Msun)*k_e)/(4*np.pi*dphi * (vff(rm,Msun) * rm ) )
    return tmp*2*Z_fe


