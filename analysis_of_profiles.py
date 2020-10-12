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


from pipeline_core import ObsList_SA,ObsList_SE,ObsList_all

from Misc import  doppler_correction as doppler
from Misc.doppler_correction import  day2sec


matplotlib.rcParams['figure.figsize'] = 24, 28
matplotlib.rcParams['figure.subplot.left']=0.1
matplotlib.rcParams['figure.subplot.bottom']=0.07
matplotlib.rcParams['figure.subplot.right']=0.9
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


# read data of any observation

def read_ph_res_results_data(ObsID,model,
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

    nh=data[:,29]

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
            eqw,eqw_err,
            nh]




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
            eqw,eqw_err,nh]=read_ph_res_results_data(ObsID,model,roll)
    sigma=np.mean(sigma)
    nh=np.mean(nh)

    fig = plt.figure(figsize=(12,10))
    plt.subplots_adjust(hspace=0)
    plt.subplots_adjust(wspace=1)
    rows=9
    cols=11
    #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
    ax_norm = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=7)
    ax_eqw = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=7)
    ax_edgeTau = plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=7)
    ax_po = plt.subplot2grid((rows,cols), (6, 0), rowspan=2, colspan=7)
    ax_ecut = ax_po.twinx()
    ax_chi2 = plt.subplot2grid((rows,cols), (8, 0), rowspan=1, colspan=7)

    if phases==None:
        ind_flux_min=np.argmin(norm_line)
        ind_flux_max=np.argmax(norm_line)
    else:
        ind_flux_min=phases[0]-1
        ind_flux_max=phases[1]-1
    ind_cont_flux_max=np.argmin(flux312)

    for ax in [ax_eqw,ax_edgeTau,ax_norm]:
        ax_efold=ax.twinx()
        ax_efold.errorbar(phase,flux312/1e-8,flux312_err/1e-8,color='k',drawstyle='steps-mid',ls=':',alpha=0.6)
        ax_efold.set_ylabel('Flux 3-12 keV')
        ax_efold.axvline(phase[ind_flux_min], color='r',ls='-.', alpha=0.2,lw=2,zorder=-10)
        ax_efold.axvline(phase[ind_flux_max], color='blue', ls='-.', alpha=0.2,lw=2,zorder=-10)
        ax_efold.axvline(phase[ind_cont_flux_max], color='g', ls='-.', alpha=0.2,lw=2,zorder=-10)

    for i in range(int(len(phase)/2)):
        ax_efold.text(phase[i],min(flux312/1e-8),str(i+1),fontsize=8)



    ax_norm.errorbar(phase,norm_line*1e3,norm_line_err*1e3,color='c',drawstyle='steps-mid',alpha=0.6)
    mo_norm,_,chi2_red,_=fit_const_chi_square(norm_line*1000,1000*norm_line_err.max(axis=0))
    ax_norm.axhline(mo_norm,alpha=0.3,ls=':')
    ax_norm.set_ylabel(f'Iron line \n norm, (1000 ph/cm^2/s) \n $chi2={"%.1f" % chi2_red}',fontsize=12)

    ax_norm.set_title(ObsID+f'({datamode});\n {round(cts/1e6,1)} Mcts $, Model {model}, Fe_sigma={sigma}, H_N={nh}',fontsize=10)


    ax_eqw.errorbar(phase,eqw*1e3,eqw_err*1e3,color='r',drawstyle='steps-mid',alpha=0.6)
    ax_eqw.set_ylabel(f'Line Eq. width, keV ',fontsize=12)


    uplim_ind=edgeTau==0
    edgeTau[uplim_ind]=edgeTau_err[1][uplim_ind]
    edgeTau_err[1][uplim_ind]=0
    edgeTau_err[0][uplim_ind]=edgeTau[uplim_ind]

    ax_edgeTau.errorbar(phase,edgeTau*1e2,edgeTau_err*1e2,color='g',drawstyle='steps-mid',alpha=0.6,uplims=uplim_ind)
    ax_edgeTau.set_ylabel(f'Iron K-edge tau*100',fontsize=9,color='g')
    ax_edgeTau.set_ylim(0.9*np.min(edgeTau*1e2),1.1*np.max(edgeTau*1e2))

    # ax_edgeE=ax_edgeTau.twinx()
    # ax_edgeE.errorbar(phase,edgeE,edgeE_err*0,color='c',drawstyle='steps-mid',alpha=0.6)
    # ax_edgeE.set_ylabel(f'Iron K-edge energy \n (=const*Iron line energy)',fontsize=12,color='c')

    ax_po.errorbar(phase,po,po_err,color='b',drawstyle='steps-mid',alpha=0.6)
    ax_po.set_ylabel('Phot. index')

    # ax_ecut.errorbar(phase,ecut,ecut_err,color='m',drawstyle='steps-mid',alpha=0.6)
    # ax_ecut.set_ylabel('E_cut, keV')

    ax_ecut.errorbar(phase,ecut,ecut_err,color='m',drawstyle='steps-mid',alpha=0.6)
    ax_ecut.set_ylabel('E_cut, keV')



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
    plot_spe_ratio(ObsID, ind_flux_max+1, ind_cont_flux_max+1,ax=ax_ratio)

    ax_ratio.yaxis.set_label_position("right")
    ax_ratio.yaxis.tick_right()


    ax_edge=plt.subplot2grid((rows,cols), (3, 8), rowspan=3, colspan=3)
    plot_cutoffpl_ratio(ObsID, phases=[ind_flux_max+1, ind_flux_min+1,ind_cont_flux_max+1],ax=ax_edge)
    ax_edge.yaxis.set_label_position("right")
    ax_edge.yaxis.tick_right()

    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/{model}/'+f'MJD_{round(mjd,3)}_{ObsID}_{model}.png',dpi=100)

    fig.tight_layout()
    plt.show()


def plot_ph_res_results_with_eline(ObsID='90089-11-03-01G',model='cutoffpl_en_fix',roll=0,plot_evol=1,phases=None):
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
            eqw,eqw_err,nh]=read_ph_res_results_data(ObsID,model,roll)
    sigma=np.mean(sigma)
    nh=np.mean(nh)

    fig = plt.figure(figsize=(8,10))
    plt.subplots_adjust(hspace=0)
    plt.subplots_adjust(wspace=1)
    rows=9
    cols=11
    #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
    ax_norm = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=7)
    ax_eqw = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=7)
    ax_edgeTau = plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=7)
    ax_eline = plt.subplot2grid((rows,cols), (6, 0), rowspan=2, colspan=7)
    ax_chi2 = plt.subplot2grid((rows,cols), (8, 0), rowspan=1, colspan=7)

    if phases==None:
        ind_flux_min=np.argmin(norm_line)
        ind_flux_max=np.argmax(norm_line)
    else:
        ind_flux_min=phases[0]-1
        ind_flux_max=phases[1]-1

    for ax in [ax_norm,ax_eqw,ax_edgeTau]:
        ax_efold=ax.twinx()
        ax_efold.errorbar(phase,flux312/1e-8,flux312_err/1e-8,color='k',drawstyle='steps-mid',ls=':',alpha=0.6)
        ax_efold.set_ylabel('Flux 3-12 keV')
        ax_efold.axvline(phase[ind_flux_min], color='r',ls='-.', alpha=0.2,lw=2,zorder=-10)
        ax_efold.axvline(phase[ind_flux_max], color='blue', ls='-.', alpha=0.2,lw=2,zorder=-10)


    ax_norm.errorbar(phase,norm_line*1e3,norm_line_err*1e3,color='c',drawstyle='steps-mid',alpha=0.6)
    mo_norm,_,chi2_red,_=fit_const_chi_square(norm_line*1000,1000*norm_line_err.max(axis=0))
    ax_norm.axhline(mo_norm,alpha=0.3,ls=':')
    ax_norm.set_ylabel(f'Iron line \n norm, (1000 ph/cm^2/s) \n $chi2={"%.1f" % chi2_red}',fontsize=12)

    ax_norm.set_title(ObsID+f'({datamode});\n {round(cts/1e6,1)} Mcts $, Model {model}, Fe_sigma={sigma}',fontsize=10)


    ax_eqw.errorbar(phase,eqw*1e3,eqw_err*1e3,color='r',drawstyle='steps-mid',alpha=0.6)
    ax_eqw.set_ylabel(f'Line Eq. width, keV ',fontsize=12)


    uplim_ind=edgeTau==0
    edgeTau[uplim_ind]=edgeTau_err[1][uplim_ind]
    edgeTau_err[1][uplim_ind]=0
    edgeTau_err[0][uplim_ind]=edgeTau[uplim_ind]

    ax_edgeTau.errorbar(phase,edgeTau*1e2,edgeTau_err*1e2,color='g',drawstyle='steps-mid',alpha=0.6,uplims=uplim_ind)
    ax_edgeTau.set_ylabel(f'Iron K-edge tau*100',fontsize=9,color='g')
    ax_edgeTau.set_ylim(0.9*np.min(edgeTau*1e2),1.1*np.max(edgeTau*1e2))

    # ax_edgeE=ax_edgeTau.twinx()
    # ax_edgeE.errorbar(phase,edgeE,edgeE_err*0,color='c',drawstyle='steps-mid',alpha=0.6)
    # ax_edgeE.set_ylabel(f'Iron K-edge energy \n (=const*Iron line energy)',fontsize=12,color='c')

    ax_eline.errorbar(phase,eline,eline_err,color='b',drawstyle='steps-mid',alpha=0.6)
    ax_eline.axhline(6.4,color='k',alpha=0.4,zorder=-10)
    ax_eline.set_ylabel('Iron line energy')


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
    ax_ratio.yaxis.set_label_position("right")
    ax_ratio.yaxis.tick_right()


    ax_edge=plt.subplot2grid((rows,cols), (3, 8), rowspan=3, colspan=3)
    plot_cutoffpl_ratio(ObsID, phases=[ind_flux_max+1, ind_flux_min+1],ax=ax_edge)
    ax_edge.yaxis.set_label_position("right")
    ax_edge.yaxis.tick_right()

    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/{model}/'+f'MJD_{round(mjd,3)}_{ObsID}_{model}.png',dpi=100)

    fig.tight_layout()
    plt.show()



def plot_spe_ratio(ObsID,phase1,phase2,ax=None,po_pars=None):
    model='cutoffpl_no_gauss'
    data1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{model}/spe_plots/ph_spe_no_gauss_lda_{phase1}.dat')

    data2=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{model}/spe_plots/ph_spe_no_gauss_lda_{phase2}.dat')
    if ax==None:
        fig,ax = plt.subplots(figsize=(8, 6))
    else:
        pass
    label=f"phase {phase1} / phase {phase2}"
    def ratio_error(a,b,da,db):
        f=a/b
        sigma=np.abs(f)*np.sqrt( (da/a)**2 + (db/b)**2  )
        return f, sigma

    ax.errorbar(data1[0],data1[1]/data2[1],ratio_error(data1[1],data2[1],data1[2],data2[2])[1],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    ax.set_xscale('log')

    ax.legend(loc='upper left',fontsize=8)
    ax.grid('y')
    ax.set_ylabel('spectral ratio (lda)',fontsize=8)

    if po_pars is not None:
        en=np.linspace(data1[0][0],data1[0][-1],100)
        sigma=0.3
        eline=6.4
        po1=po_pars[0]*en**(-po_pars[1])*np.exp(-en/po_pars[2])+po_pars[3]*np.exp(-(en-eline)**2/(2*sigma**2))/np.sqrt(sigma**2*2*np.pi)
        po2=po_pars[4]*en**(-po_pars[5])*np.exp(-en/po_pars[6])+po_pars[7]*np.exp(-(en-eline)**2/(2*sigma**2))/np.sqrt(sigma**2*2*np.pi)
        ax.plot(en,po1/po2,'k-.')

# plot_spe_ratio('90427-01-03-00',7,14,po_pars=[0.179063,-0.671173,5.8954,0.018753,0.250967,-0.431976,6.23211,0.0167662])
# plt.show()

# plot_spe_ratio('90089-11-03-00G',8,11,po_pars=[0.196072,-0.704540,5.47549,1.93215E-02,0.329151,-0.430498,7.44940,2.02309E-02])
# plt.show()


def plot_cutoffpl_ratio(ObsID,phases,ax=None):
    model='cutoffpl_no_gauss'
    if ax==None:
        fig,ax = plt.subplots(figsize=(8, 6))
    else:
        pass

    for phase in phases:

        data1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{model}/spe_plots/ph_spe_no_gauss_ra_{phase}.dat')


        label=f"phase {phase} / cutoffpl "

        ax.errorbar(data1[0],data1[1],data1[2],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
        ax.set_xscale('log')

    ax.legend(loc='upper left',fontsize=8)
    ax.grid('y')
    ax.axhline(1,color='k',ls=':',zorder=-10,alpha=0.6)
    ax.set_ylabel(' ratio ',fontsize=8)



def plot_cutoffpl_del(ObsID,phases,ax=None):
    model='cutoffpl_no_gauss'
    if ax==None:
        fig,ax = plt.subplots(figsize=(8, 6))
    else:
        pass

    for phase in phases:

        data1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{model}/spe_plots/ph_spe_no_gauss_del_{phase}.dat')


        label=f"phase {phase} / cutoffpl "

        ax.errorbar(data1[0],data1[1],data1[2],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
        ax.set_xscale('log')

    ax.legend(loc='upper left',fontsize=8)
    ax.grid('y')
    ax.axhline(1,color='k',ls=':',zorder=-10,alpha=0.6)
    ax.set_ylabel(' ratio ',fontsize=8)






def calculate_parameters_data(ObsID,model):
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
            eqw,eqw_err]=read_ph_res_results_data(ObsID,model,0)
    sigma=np.mean(sigma)

    fl=flux312
    fl_e=flux312_err.max(axis=0)

    flux_argmax=np.argmax(fl)
    flux_argmin=np.argmin(fl)

    rel_dip_flux,rel_dip_flux_err=ratio_error(fl[flux_argmax],fl[flux_argmin],fl_e[flux_argmax],fl_e[flux_argmin])

    fl=norm_line
    fl_e=norm_line_err.max(axis=0)

    flux_argmax=np.argmax(fl)
    flux_argmin=np.argmin(fl)

    rel_dip_iron,rel_dip_iron_err=ratio_error(fl[flux_argmax],fl[flux_argmin],fl_e[flux_argmax],fl_e[flux_argmin])

    xte_obs=ObservationXTE(ObsID)
    xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,model+'_rel_dip_flux',rel_dip_flux)
    xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,model+'_rel_dip_flux_err',rel_dip_flux_err)

    xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,model+'_rel_dip_iron',rel_dip_iron)
    xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,model+'_rel_dip_iron_err',rel_dip_iron_err)



STOP
#%% run this


#for ObsID in ['90427-01-03-00','90089-11-03-01G']:
plot_ph_res_results(ObsID='90089-11-03-01G',model='cutoffpl_en_fix')

plot_ph_res_results(ObsID='90089-11-03-00G',model='cutoffpl_en_fix_edge_fix')

plot_ph_res_results_with_eline(ObsID='90089-11-03-00G',model='cutoffpl_en_free_edge_fix')


#%% run sa observations with edge

err=[]
msg=[]
plt.ioff()
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_SA
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            plot_ph_res_results(ObsID=ObsID,model='cutoffpl_en_fix_edge_fix')
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







#%% run all observations with cutoffpl


err=[]
msg=[]
plt.ioff()
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_SA+ObsList_SE
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            plot_ph_res_results(ObsID=ObsID,model='cutoffpl_en_fix')
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




#%% calculate fasebin_stuff for observation
err=[]
msg=[]
plt.ioff()
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_SA+ObsList_SE
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            calculate_parameters_data(ObsID=ObsID,model='cutoffpl_en_fix')


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

#%% stack observations

def roll_fasebin_files(ObsList=['90089-11-01-02', '90089-11-01-03', '90089-11-01-04'],
                       Roll_list=[0,-3,-1],folder_name='fasebin_combine_rising_phase'):

    fig,ax=plt.subplots(3,figsize=(12,12))

    os.system(f'rm -rf /Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}')
    os.system(f'mkdir /Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}')

    all_cts=[]

    for ObsID,roll in zip(ObsList,Roll_list):

        os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin')
        os.system('rm -f fasebin_sys_roll.pha')
        os.system('cp fasebin_sys.pha fasebin_sys_roll.pha ')


        fasebin_roll=fits.open('fasebin_sys_roll.pha',mode='update')
        cts=fasebin_roll[1].data['counts'].sum(axis=1)
        #ax[0].plot(cts,'+-',label=f'{ObsID} - original')
        ax[0].plot(cts/fasebin_roll[1].data['exposure'][0],'+-',label=f'{ObsID} - original')
        ax[0].set_ylabel('counts/s/ phase_bin')

        tmp_data=fasebin_roll[1].data

        Indeces=np.arange(0,len(tmp_data))
        fasebin_roll[1].data=tmp_data[np.roll(Indeces,roll)]
        fasebin_roll.close()
        fasebin_check=fits.open('fasebin_sys_roll.pha')
        cts=fasebin_check[1].data['counts'].sum(axis=1)
        fasebin_check.close()
        all_cts.append(cts)

        # ax[1].plot(cts-np.mean(cts),'+-',label=f'{ObsID} - rolled')
        # ax[1].set_ylabel('counts-mean')

        ax[1].plot(cts/fasebin_roll[1].data['exposure'][0],'+-',label=f'{ObsID} - rolled ')
        ax[1].set_ylabel('counts/s')


    all_cts=np.array(all_cts)
    all_cts=np.sum(all_cts,axis=0)
    ax[2].plot(all_cts,'d-',label='counts sum (rolled)')


    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}')
    path=[f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/fasebin_sys_roll.pha' for ObsID in ObsList]
    path='\n'.join(path)
    print(path)
    os.system(f"echo '{path}' > alldays_pha.list")

    os.system("fbadd infile='@alldays_pha.list' outfile='fasebin_sys.pha'")

    for ObsID in ObsList:
        os.system(f"cp /Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/response.rsp ./response_{ObsID}.rsp")
        os.system(f"cp /Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/fasebin_sys_roll.pha  ./fasebin_sys_roll_{ObsID}.pha")


    fasebin_joined=fits.open('fasebin_sys.pha')
    cts=fasebin_joined[1].data['counts'].sum(axis=1)
    ax[2].plot(cts,'o',label='fasebin_sys.pha - counts (all observations summed)')


    ax[0].legend()
    ax[1].legend()
    ax[2].legend()


    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}/pulses.png')
    plt.show()

#roll_fasebin_files(ObsList=['90089-11-01-02', '90089-11-01-03', '90089-11-01-04'], Roll_list=[0,-3,-1])

roll_fasebin_files(ObsList=['90089-11-03-00G','90089-11-03-01G','90089-11-03-02'],Roll_list=[0,-5,+6],
                   folder_name='fasebin_combine_top_phase')


