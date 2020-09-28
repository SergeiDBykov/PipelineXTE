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
from PipelineXTE.pipeline_core import *
def gauss(t,t0,sigma,N):
    return N*np.exp(-(t-t0)**2/(2*sigma**2))



from Misc import  doppler_correction as doppler
from Misc.doppler_correction import  day2sec


matplotlib.rcParams['figure.figsize'] = 24, 28
matplotlib.rcParams['figure.subplot.left']=0.07
matplotlib.rcParams['figure.subplot.bottom']=0.07
matplotlib.rcParams['figure.subplot.right']=0.94
matplotlib.rcParams['figure.subplot.top']=0.93


import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook',rc={"xtick.top" : True,'xtick.direction':'inout','ytick.direction':'inout','xtick.minor.visible':True,'ytick.minor.visible':True})

plt.ioff()

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
pulse_profile_save_path='/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/'

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
'90014-01-06-03',
'90427-01-04-03']


ObsList_RISE=['90089-11-02-00',
 '90089-11-02-01',
 '90089-11-02-02',
 '90089-11-02-03',
 '90089-11-02-03G',
 '90089-11-02-04',
 '90089-11-02-05',
 '90089-11-02-06',
 '90089-11-02-07',
 '90089-11-02-08',
 '90089-11-02-09',
 '90089-11-02-10',
 '90089-11-01-02',
 '90089-11-01-03',
 '90089-11-01-04']

ObsList_TOP=['90089-11-04-04','90089-11-04-03','90089-11-04-02G','90089-11-04-01',
             '90089-11-04-00G','90089-11-03-05','90089-11-03-04','90089-11-03-03',
             '90089-11-03-02','90089-11-03-01G','90089-11-03-00G']


ObsList_DECAY_1=['90427-01-03-00',
 '90427-01-03-01',
 '90427-01-03-02',
 '90427-01-03-05',
 '90427-01-03-06',
 '90427-01-03-07',
 '90427-01-03-09',
 '90427-01-03-11',
 '90427-01-03-12',
 '90427-01-03-14G',
 '90014-01-02-00',
 '90014-01-02-03',
 '90014-01-02-08',
 '90014-01-02-10',
 '90014-01-02-13',
 '90014-01-02-15',
 '90014-01-03-00',
 '90014-01-03-01',
 '90014-01-03-02',
 '90014-01-03-020',
 '90014-01-03-03',
 '90014-01-04-00',
 '90014-01-04-01',
 '90014-01-04-02',
 '90014-01-04-03']


ObsList_DECAY_2=[ '90014-01-05-00',
 '90014-01-05-01',
 '90014-01-05-02',
 '90014-01-05-04',
 '90014-01-05-05',
 '90014-01-06-01',
 '90014-01-06-02',
 '90014-01-06-03',
 '90014-01-07-04',
 '90427-01-04-00',
 '90427-01-04-01',
 '90427-01-04-02',
 '90427-01-04-04',
 '90427-01-04-03',
 '90427-01-04-05']



#%% read data of any observation

def read_ph_res_results_data(ObsID='90089-11-03-01G',model='cutoffpl_en_fix',
                             roll=0):


    xte_obs=ObservationXTE(ObsID)
    ser=xte_obs.pandas_series()

    mjd=ser['MJD_START']
    expo=ser['EXPOSURE']/1000
    expo=np.round(expo,1)
    tot_flux=ser['cutoffpl_cutoffpl_flux']/1e-8
    tot_flux=np.round(tot_flux,2)

    datalist=ser['fasebin_cfg']
    if datalist=='se':
        datamode='E_125us_64M_0_1s'
    elif datalist=='sa':
        datamode='B_16ms_46M_0_49_H/B_16ms_64M_0_249'



    fasebin_file=fits.open(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/fasebin_orig.pha')
    cts=fasebin_file[1].data['counts'].sum()
    fasebin_file.close()


    data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{model}/ph_res_{model}.dat')

    if roll:
        tmp=data[:,4]
        tmp_arg=np.argmin(tmp)
        data=np.roll(data,-tmp_arg,axis=0)
    else:
        data=data


    N_sp=(data[0,1]-1)/2
    spe_num=data[:,0]

    data=np.vstack((data,data))
    nph=data[0,1]
    data[:,0]=np.arange(1,nph) #create new spe_num
    spe_num=data[:,0]
    phase=((spe_num-1)/(N_sp))
    chi2=data[:,2]

    #parameters


    flux312=data[:,4]
    flux312_low=flux312-data[:,5]
    flux312_hi=data[:,6]-flux312
    flux312_err=np.vstack((flux312_low, flux312_hi))#.max(axis=0)


    po=data[:,7]
    po_lo=po-data[:,8]
    po_hi=data[:,9]-po
    po_err=np.vstack((po_lo,po_hi))


    ecut=data[:,10]
    ecut_lo=ecut-data[:,11]
    ecut_hi=data[:,12]-ecut
    ecut_err=np.vstack((ecut_lo,ecut_hi))


    eline=data[:,13]
    eline_lo=eline-data[:,14]
    eline_hi=data[:,15]-eline
    eline_err=np.vstack((eline_lo,eline_hi))

    norm_line=data[:,16]
    norm_line_lo=norm_line-data[:,17]
    norm_line_hi=data[:,18]-norm_line
    norm_line_err=np.vstack((norm_line_lo,norm_line_hi))

    edgeE=data[:,19]
    edgeE_lo=edgeE-data[:,20]
    edgeE_hi=data[:,21]-edgeE
    edgeE_err=np.vstack((edgeE_lo,edgeE_hi))

    edgeTau=data[:,22]
    edgeTau_lo=edgeTau-data[:,23]
    edgeTau_hi=data[:,24]-edgeTau
    edgeTau_err=np.vstack((edgeTau_lo,edgeTau_hi))

    sigma=data[:,25]

    eqw=data[:,26]
    eqw_lo=eqw-data[:,27]
    eqw_hi=data[:,28]-eqw
    eqw_err=np.vstack((eqw_lo,eqw_hi))

    return [model,phase,
            mjd,
            tot_flux,
            datamode,
            cts,
            flux312,flux312_err,
            po,po_err,
            ecut,ecut_err,
            eline,eline_err,
            norm_line,norm_line_err,
            edgeE,edgeE_err,
            edgeTau,edgeTau_err,
            chi2,
            sigma,
            eqw,eqw_err]




def plot_ph_res_results(ObsID='90089-11-03-01G',model='cutoffpl_en_fix',roll=0,plot_evol=1,phases=None):
    [model,phase,
            mjd,
            tot_flux,
            datamode,
            cts,
            flux312,flux312_err,
            po,po_err,
            ecut,ecut_err,
            eline,eline_err,
            norm_line,norm_line_err,
            edgeE,edgeE_err,
            edgeTau,edgeTau_err,
            chi2,
            sigma,
            eqw,eqw_err]=read_ph_res_results_data(ObsID,model,roll)
    sigma=np.mean(sigma)

    fig = plt.figure(figsize=(12,16))
    plt.subplots_adjust(hspace=0)
    plt.subplots_adjust(wspace=1)
    rows=9
    cols=11
    #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
    ax_norm = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=7)
    ax_eline = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=7)
    ax_edgeTau = plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=7)
    ax_po = plt.subplot2grid((rows,cols), (6, 0), rowspan=2, colspan=7)
    ax_ecut = ax_po.twinx()
    ax_chi2 = plt.subplot2grid((rows,cols), (8, 0), rowspan=1, colspan=7)

    if phases==None:
        ind_flux_min=np.argmin(flux312)
        ind_flux_max=np.argmax(flux312)
    else:
        ind_flux_min=phases[0]-1
        ind_flux_max=phases[1]-1

    for ax in [ax_norm,ax_eline,ax_edgeTau]:
        ax_efold=ax.twinx()
        ax_efold.errorbar(phase,flux312/1e-8,flux312_err/1e-8,color='k',drawstyle='steps-mid',ls=':',alpha=0.6)
        ax_efold.set_ylabel('Flux 3-12 keV')
        ax_efold.axvline(phase[ind_flux_min], color='k',ls='-.', alpha=0.4,lw=2,zorder=-10)
        ax_efold.axvline(phase[ind_flux_max], color='k', ls='-.', alpha=0.4,lw=2,zorder=-10)


    ax_norm.errorbar(phase,norm_line*1e3,norm_line_err*1e3,color='c',drawstyle='steps-mid',alpha=0.6)
    mo_norm,_,chi2_red,_=fit_const_chi_square(norm_line*1000,1000*norm_line_err.max(axis=0))
    ax_norm.axhline(mo_norm,alpha=0.3,ls=':')
    ax_norm.set_ylabel(f'Iron line \n norm, (1000 ph/cm^2/s) \n $chi2={"%.1f" % chi2_red}',fontsize=12)
    ax_norm.set_title(ObsID+f'({datamode});$F_{{3-12}}$={round(tot_flux,2)} $10^{{-8}}cgs;\n {round(cts/1e6,1)} Mcts $, Model {model}, sigma={sigma}',fontsize=10)


    ax_eline.errorbar(phase,eline,eline_err,color='r',drawstyle='steps-mid',alpha=0.6)
    mo_eline,_,chi2_red,_=fit_const_chi_square(eline,eline_err.max(axis=0))
    ax_eline.axhline(mo_eline,alpha=0.3,ls=':')
    ax_eline.set_ylabel(f'Line energy, keV \n $chi2={"%.1f" % chi2_red}',fontsize=12)
    ax_eline.axhline(6.4,color='k',alpha=0.6,zorder=-10)


    ax_edgeTau.errorbar(phase,edgeTau*1e2,edgeTau_err*1e2,color='g',drawstyle='steps-mid',alpha=0.6)
    mo_eline,_,chi2_red,_=fit_const_chi_square(edgeTau,edgeTau_err.max(axis=0))
    ratioEdge=edgeE[0]/eline[0]
    ax_edgeTau.set_ylabel(f'Iron K-edge $\tau$*100, \n $chi2={"%.1f" % chi2_red} \n Fe K_Edge={"%.2f"%ratioEdge}*Fe K_alpha ',fontsize=9,color='g')
    ax_edgeTau.set_ylim(0.9*np.min(edgeTau*1e2),1.1*np.max(edgeTau*1e2))

    # ax_edgeE=ax_edgeTau.twinx()
    # ax_edgeE.errorbar(phase,edgeE,edgeE_err*0,color='c',drawstyle='steps-mid',alpha=0.6)
    # ax_edgeE.set_ylabel(f'Iron K-edge energy \n (=const*Iron line energy)',fontsize=12,color='c')

    ax_po.errorbar(phase,po,po_err,color='b',drawstyle='steps-mid',alpha=0.6)
    ax_po.set_ylabel('Phot. index')

    # ax_ecut.errorbar(phase,ecut,ecut_err,color='m',drawstyle='steps-mid',alpha=0.6)
    # ax_ecut.set_ylabel('E_cut, keV')

    ax_ecut.errorbar(phase,eqw*1e3,eqw_err*1e3,color='m',drawstyle='steps-mid',alpha=0.6)
    ax_ecut.set_ylabel('Iron line EQW, eV')



    ax_chi2.errorbar(phase,chi2,0,color='k',drawstyle='steps-mid',alpha=0.6)
    ax_chi2.set_ylabel('chi2')
    ax_chi2.set_xlabel('Phase',fontsize=12)




    ax_flux= plt.subplot2grid((rows,cols), (0, 8), rowspan=2, colspan=3)
    time=ObsParams.MJD_START
    flux=ObsParams.cutoffpl_cutoffpl_flux/1e-8
    ax_flux.semilogy(time,flux,color='b',marker='s',lw=0,ms=4,alpha=0.8)
    ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b',fontsize=8)
    ax_flux.set_xlabel('Time, MJD')
    ax_flux.yaxis.set_label_position("right")
    ax_flux.yaxis.tick_right()
    ax_flux.axvline(mjd,color='r',ls='-.')


    #fig.tight_layout()

    ax_ratio=plt.subplot2grid((rows,cols), (6, 8), rowspan=3, colspan=3)
    plot_spe_ratio(ObsID, ind_flux_max+1, ind_flux_min+1,ax=ax_ratio)

    ax_edge=plt.subplot2grid((rows,cols), (3, 8), rowspan=3, colspan=3)
    plot_cutoffpl_ratio(ObsID, phases=[ind_flux_max+1, ind_flux_min+1],ax=ax_edge)

    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/'+f'MJD_{round(mjd,3)}_{ObsID}_{model}.png',dpi=500)


    plt.show()



def plot_ph_res_results_no_eline(ObsID='90089-11-03-01G',model='cutoffpl_en_fix',roll=0,plot_evol=1,phases=None):
    [model,phase,
            mjd,
            tot_flux,
            datamode,
            cts,
            flux312,flux312_err,
            po,po_err,
            ecut,ecut_err,
            eline,eline_err,
            norm_line,norm_line_err,
            edgeE,edgeE_err,
            edgeTau,edgeTau_err,
            chi2,
            sigma,
            eqw,eqw_err]=read_ph_res_results_data(ObsID,model,roll)
    sigma=np.mean(sigma)

    fig = plt.figure(figsize=(12,16))
    plt.subplots_adjust(hspace=0)
    plt.subplots_adjust(wspace=1)
    rows=9
    cols=11
    #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
    ax_norm = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=7)
    ax_eqw = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=7)
    ax_edgeTau = plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=7)
    ax_chi2 = plt.subplot2grid((rows,cols), (8, 0), rowspan=1, colspan=7)

    if phases==None:
        ind_flux_min=np.argmin(flux312)
        ind_flux_max=np.argmax(flux312)
    else:
        ind_flux_min=phases[0]-1
        ind_flux_max=phases[1]-1

    for ax in [ax_norm,ax_eqw,ax_edgeTau]:
        ax_efold=ax.twinx()
        ax_efold.errorbar(phase,flux312/1e-8,flux312_err/1e-8,color='k',drawstyle='steps-mid',ls=':',alpha=0.6)
        ax_efold.set_ylabel('Flux 3-12 keV')
        ax_efold.axvline(phase[ind_flux_min], color='k',ls='-.', alpha=0.4,lw=2,zorder=-10)
        ax_efold.axvline(phase[ind_flux_max], color='k', ls='-.', alpha=0.4,lw=2,zorder=-10)


    ax_norm.errorbar(phase,norm_line*1e3,norm_line_err*1e3,color='c',drawstyle='steps-mid',alpha=0.6)
    mo_norm,_,chi2_red,_=fit_const_chi_square(norm_line*1000,1000*norm_line_err.max(axis=0))
    ax_norm.axhline(mo_norm,alpha=0.3,ls=':')
    ax_norm.set_ylabel(f'Iron line \n norm, (1000 ph/cm^2/s) \n $chi2={"%.1f" % chi2_red}',fontsize=12)
    ax_norm.set_title(ObsID+f'({datamode});$F_{{3-12}}$={round(tot_flux,2)} $10^{{-8}}cgs;\n {round(cts/1e6,1)} Mcts $, Model {model}, sigma={sigma}',fontsize=10)


    ax_eqw.errorbar(phase,eqw*1e3,eqw_err*1e3,color='r',drawstyle='steps-mid',alpha=0.6)
    ax_eqw.set_ylabel(f'Line Eq. width, eV \n $chi2={"%.1f" % chi2_red}',fontsize=12)


    uplim_ind=edgeTau==0
    edgeTau[uplim_ind]=edgeTau_err[1][uplim_ind]
    edgeTau_err[1][uplim_ind]=0
    edgeTau_err[0][uplim_ind]=edgeTau[uplim_ind]

    ax_edgeTau.errorbar(phase,edgeTau*1e2,edgeTau_err*1e2,color='g',drawstyle='steps-mid',alpha=0.6,uplims=uplim_ind)
    ratioEdge=edgeE[0]/eline[0]
    ax_edgeTau.set_ylabel(f'Iron K-edge $\tau$*100, \n Fe K_Edge={"%.2f"%ratioEdge}*Fe K_alpha ',fontsize=9,color='g')
    ax_edgeTau.set_ylim(-0.4,1.1*np.max(edgeTau*1e2))

    # ax_edgeE=ax_edgeTau.twinx()
    # ax_edgeE.errorbar(phase,edgeE,edgeE_err*0,color='c',drawstyle='steps-mid',alpha=0.6)
    # ax_edgeE.set_ylabel(f'Iron K-edge energy \n (=const*Iron line energy)',fontsize=12,color='c')


    ax_chi2.errorbar(phase,chi2,0,color='k',drawstyle='steps-mid',alpha=0.6)
    ax_chi2.set_ylabel('chi2')
    ax_chi2.set_xlabel('Phase',fontsize=12)




    ax_flux= plt.subplot2grid((rows,cols), (0, 8), rowspan=2, colspan=3)
    time=ObsParams.MJD_START
    flux=ObsParams.cutoffpl_cutoffpl_flux/1e-8
    ax_flux.semilogy(time,flux,color='b',marker='s',lw=0,ms=4,alpha=0.8)
    ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b',fontsize=8)
    ax_flux.set_xlabel('Time, MJD')
    ax_flux.yaxis.set_label_position("right")
    ax_flux.yaxis.tick_right()
    ax_flux.axvline(mjd,color='r',ls='-.')


    #fig.tight_layout()

    ax_ratio=plt.subplot2grid((rows,cols), (6, 8), rowspan=3, colspan=3)
    plot_spe_ratio(ObsID, ind_flux_max+1, ind_flux_min+1,ax=ax_ratio)

    ax_edge=plt.subplot2grid((rows,cols), (3, 8), rowspan=3, colspan=3)
    plot_cutoffpl_ratio(ObsID, phases=[ind_flux_max+1, ind_flux_min+1],ax=ax_edge)

    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/'+f'MJD_{round(mjd,3)}_{ObsID}_{model}.png',dpi=500)


    plt.show()

def plot_spe_ratio(ObsID,phase1,phase2,ax=None):
    model='cutoffpl_no_gauss'
    data1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{model}/spe_plots/ph_spe_no_gauss_lda_{phase1}.dat')

    data2=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{model}/spe_plots/ph_spe_no_gauss_lda_{phase2}.dat')
    if ax==None:
        fig,ax = plt.subplots(figsize=(8, 6))
    else:
        pass
    label=f"Obs {ObsID}, phase {phase1} / phase {phase2}, \n {model}"
    def ratio_error(a,b,da,db):
        f=a/b
        sigma=np.abs(f)*np.sqrt( (da/a)**2 + (db/b)**2  )
        return f, sigma

    ax.errorbar(data1[0],data1[1]/data2[1],ratio_error(data1[1],data2[1],data1[2],data2[2])[1],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    ax.set_xscale('log')

    ax.legend(loc='upper left',fontsize=8)
    ax.grid('y')
    ax.set_ylabel('spectral ratio (lda)',fontsize=8)

    #plt.show()


# def plot_spe_difference(ObsID,phase1,phase2,ax=None):
#     model='cutoffpl_no_gauss'
#     data1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{model}/spe_plots/ph_spe_no_gauss_lda_{phase1}.dat')

#     data2=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{model}/spe_plots/ph_spe_no_gauss_lda_{phase2}.dat')
#     if ax==None:
#         fig,ax = plt.subplots(figsize=(8, 6))
#     else:
#         pass
#     label=f"Obs {ObsID}, phase {phase1} / phase {phase2}, \n {model}"

#     ax.errorbar(data1[0],data1[1]-data2[1],np.sqrt(data1[2]**2+data2[2]**2),data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
#     ax.set_xscale('log')

#     ax.legend(loc='upper left',fontsize=8)
#     ax.grid('y')
#     ax.set_ylabel('spectral difference (lda)',fontsize=8)

#     #plt.show()



def plot_cutoffpl_ratio(ObsID,phases,ax=None):
    model='cutoffpl_no_gauss'
    if ax==None:
        fig,ax = plt.subplots(figsize=(8, 6))
    else:
        pass

    for phase in phases:

        data1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{model}/spe_plots/ph_spe_no_gauss_ra_{phase}.dat')


        label=f"Obs {ObsID}, phase {phase} / cutoffpl best fit model"

        ax.errorbar(data1[0],data1[1],data1[2],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
        ax.set_xscale('log')

    ax.legend(loc='upper left',fontsize=8)
    ax.grid('y')
    ax.axhline(1,color='k',ls=':',zorder=-10,alpha=0.6)
    ax.set_ylabel(' ratio ',fontsize=8)

    #plt.show()


STOP
#%% run this

models=['cutoffpl_en_fix_edge_fix']#,'cutoffpl_en_fix']
#'90089-11-03-01G' '90427-01-03-00'
for ObsID in ['90427-01-03-00','90089-11-03-01G']:
    for model in models:
        plot_ph_res_results_no_eline(ObsID=ObsID,model=model)



#%% bulk run

ObsList_SA=['90427-01-03-00',
 '90427-01-03-01',
 '90427-01-03-02',
 '90427-01-03-05',
 '90427-01-03-06',
 '90427-01-03-07',
 '90427-01-03-09',
 '90427-01-03-11',
 '90427-01-03-12',
 '90427-01-03-14G',
 '90014-01-02-00',
 '90014-01-02-03',
 '90014-01-02-08',
 '90014-01-02-10',
 '90014-01-02-13',
 '90014-01-02-15',
 '90014-01-03-00',
 '90014-01-03-01',
 '90014-01-03-02',
 '90014-01-03-020',
 '90014-01-03-03',
 '90014-01-04-00',
 '90014-01-04-01',
 '90014-01-04-02',
 '90014-01-04-03']+['90089-11-03-02','90089-11-03-01G','90089-11-03-00G']


for ObsID in ObsList_SA:

    plot_ph_res_results_no_eline(ObsID=ObsID,model='cutoffpl_en_fix_edge_fix')
    plt.close('all')




#%% run all observations with cutoffpl
ObsList=[
'90089-11-04-04','90089-11-04-03','90089-11-04-02G','90089-11-04-01',
             '90089-11-04-00G','90089-11-03-05','90089-11-03-04','90089-11-03-03',
    '90089-11-03-02','90089-11-03-01G','90089-11-03-00G',
    '90089-11-02-00',
 '90089-11-02-01',
 '90089-11-02-02',
 '90089-11-02-03',
 '90089-11-02-03G',
 '90089-11-02-04',
 '90089-11-02-05',
 '90089-11-02-06',
 '90089-11-02-07',
 '90089-11-02-08',
 '90089-11-02-09',
 '90089-11-02-10',
 '90089-11-01-02',
 '90089-11-01-03',
 '90089-11-01-04',
 '90427-01-03-00',
 '90427-01-03-01',
 '90427-01-03-02',
 '90427-01-03-05',
 '90427-01-03-06',
 '90427-01-03-07',
 '90427-01-03-09',
 '90427-01-03-11',
 '90427-01-03-12',
 '90427-01-03-14G',
 '90014-01-02-00',
 '90014-01-02-03',
 '90014-01-02-08',
 '90014-01-02-10',
 '90014-01-02-13',
 '90014-01-02-15',
 '90014-01-03-00',
 '90014-01-03-01',
 '90014-01-03-02',
 '90014-01-03-020',
 '90014-01-03-03',
 '90014-01-04-00',
 '90014-01-04-01',
 '90014-01-04-02',
 '90014-01-04-03',
 '90014-01-05-00',
 '90014-01-05-01',
 '90014-01-05-02',
 '90014-01-05-03',
 '90014-01-05-04',
 '90014-01-05-05',
 '90014-01-05-06',
 '90014-01-06-00',
 '90014-01-06-01',
 '90014-01-06-02',
 '90014-01-06-03',
 '90014-01-07-00',
 '90014-01-07-01',
 '90014-01-07-02',
 '90014-01-07-03',
 '90014-01-07-04',
 '90014-01-08-00',
 '90014-01-08-01',
 '90427-01-04-00',
 '90427-01-04-01',
 '90427-01-04-02',
 '90427-01-04-03',
 '90427-01-04-04',
 '90427-01-04-05' ]


err=[]
msg=[]
plt.ioff()
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            plot_ph_res_results_no_eline(ObsID=ObsID,model='cutoffpl_en_fix')
            plt.close('all')

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg







'''
xspec - ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res_cutoffpl_no_gauss.txt

xspec - ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res_cutoffpl_en_free.txt

xspec - ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res_cutoffpl_en_fix.txt

xspec - ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res_cutoffpl_en_fix_edge_fix.txt

xspec - ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res_cutoffpl_en_free_edge_free.txt

xspec - ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res_cutoffpl_en_free_edge_free_ionize.txt


'''