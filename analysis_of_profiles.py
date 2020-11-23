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


# matplotlib.rcParams['figure.figsize'] = 24, 28
matplotlib.rcParams['figure.subplot.left']=0.05
matplotlib.rcParams['figure.subplot.bottom']=0.07
matplotlib.rcParams['figure.subplot.right']=0.9
matplotlib.rcParams['figure.subplot.top']=0.93


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


    po_norm=data[:,13]
    po_norm_lo=po_norm-data[:,14]
    po_norm_hi=data[:,15]-po_norm
    po_norm_err=np.vstack((po_norm_lo,po_norm_hi))

    eline=data[:,16]
    eline_lo=eline-data[:,17]
    eline_hi=data[:,18]-eline
    eline_err=np.vstack((eline_lo,eline_hi))

    norm_line=data[:,19]
    norm_line_lo=norm_line-data[:,20]
    norm_line_hi=data[:,21]-norm_line
    norm_line_err=np.vstack((norm_line_lo,norm_line_hi))

    edgeE=data[:,22]
    edgeE_lo=edgeE-data[:,23]
    edgeE_hi=data[:,24]-edgeE
    edgeE_err=np.vstack((edgeE_lo,edgeE_hi))

    edgeTau=data[:,25]
    edgeTau_lo=edgeTau-data[:,26]
    edgeTau_hi=data[:,27]-edgeTau
    edgeTau_err=np.vstack((edgeTau_lo,edgeTau_hi))

    sigma=data[:,28]

    eqw=data[:,29]
    eqw_lo=eqw-data[:,30]
    eqw_hi=data[:,31]-eqw
    eqw_err=np.vstack((eqw_lo,eqw_hi))

    nh=data[:,32]

    sys_err=data[:,33]


    if ObsID=='group_1':
        ObsID='90089-11-03-01G'
    elif ObsID=='group_2':
        ObsID='90427-01-03-00'
    elif ObsID=='group_3':
        ObsID='90427-01-03-14G'
    elif ObsID=='group_4':
         ObsID='90014-01-02-10'
    elif ObsID=='group_5':
        ObsID='90427-01-03-09'
    elif ObsID=='group_6':
        ObsID='90014-01-02-13'
    elif ObsID=='group_7':
        ObsID='90014-01-03-020'
    else:
        pass
    xte_obs=ObservationXTE(ObsID)
    ser=xte_obs.pandas_series()

    mjd=ser['MJD_START']
    expo=ser['EXPOSURE']/1000
    expo=np.round(expo,1)
    tot_flux=ser['cutoffpl_cutoffpl_flux']/1e-8
    tot_flux=np.round(tot_flux,2)

    datalist=ser['fasebin_cfg']
    if datalist=='se':
        datamode='E_125us'
    elif datalist=='sa':
        datamode='B_16ms'


    return [model,phase,
            mjd,
            tot_flux,
            datamode,
            cts,
            flux312,flux312_err,
            po,po_err,
            ecut,ecut_err,
            po_norm,po_norm_err,
            eline,eline_err,
            norm_line,norm_line_err,
            edgeE,edgeE_err,
            edgeTau,edgeTau_err,
            chi2,
            sigma,
            eqw,eqw_err,
            nh,
            sys_err]




def plot_ph_res_results(ObsID='90089-11-03-01G',model='cutoffpl',roll=0,plot_evol=1,phases=None):
    [model,phase,
            mjd,
            tot_flux,
            datamode,
            cts,
            flux312,flux312_err,
            po,po_err,
            ecut,ecut_err,
            po_norm,po_norm_err,
            eline,eline_err,
            norm_line,norm_line_err,
            edgeE,edgeE_err,
            edgeTau,edgeTau_err,
            chi2,
            sigma,
            eqw,eqw_err,nh,sys_err]=read_ph_res_results_data(ObsID,model,roll)
    sigma=np.mean(sigma)
    nh=np.mean(nh)
    sys_err=np.mean(sys_err)

    fig = plt.figure(figsize=(18,10))
    plt.subplots_adjust(hspace=0)
    plt.subplots_adjust(wspace=1)
    rows=9
    cols=14
    #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
    ax_norm = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=7)
    ax_eqw = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=7)
    ax_edgeTau = plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=7)
    ax_po = plt.subplot2grid((rows,cols), (6, 0), rowspan=2, colspan=7)
    ax_po_norm = ax_po.twinx()
    ax_chi2 = plt.subplot2grid((rows,cols), (8, 0), rowspan=1, colspan=7)

    if phases==None:
        ind_flux_min=np.argmin(norm_line)
        ind_flux_max=np.argmax(norm_line)
    else:
        ind_flux_min=phases[0]-1
        ind_flux_max=phases[1]-1
    ind_cont_flux_min=np.argmin(flux312)
    ind_cont_flux_max=np.argmax(flux312)

    for ax in [ax_eqw,ax_edgeTau,ax_norm]:
        ax_efold=ax.twinx()
        ax_efold.errorbar(phase,flux312/1e-8,flux312_err/1e-8,color='k',drawstyle='steps-mid',ls=':',alpha=0.6)
        ax_efold.set_ylabel('Flux 3-12 keV')
        ax_efold.axvline(phase[ind_flux_min], color='r',ls='-.', alpha=0.2,lw=2,zorder=-10)
        ax_efold.axvline(phase[ind_flux_max], color='blue', ls='-.', alpha=0.2,lw=2,zorder=-10)
        ax_efold.axvline(phase[ind_cont_flux_min], color='g', ls='-.', alpha=0.2,lw=2,zorder=-10)
        ax_efold.axvline(phase[ind_cont_flux_max], color='m', ls='-.', alpha=0.2,lw=2,zorder=-10)

    for i in range(int(len(phase)/2)):

        ax_efold.text(phase[i],min(flux312/1e-8),str(i+1),fontsize=8)



    ax_norm.errorbar(phase,norm_line*1e3,norm_line_err*1e3,color='turquoise',drawstyle='steps-mid',alpha=0.6)
    #mo_norm,_,chi2_red,_=fit_const_chi_square(norm_line*1000,1000*norm_line_err.max(axis=0))
    #ax_norm.axhline(mo_norm,alpha=0.3,ls=':')
    ax_norm.set_ylabel(f'Iron line norm, \n (1000 $ph cm^{-2}s^{-1}$)',fontsize=12)

    ax_norm.set_title(ObsID+f'({datamode});\n {round(cts/1e6,1)} Mcts $, Model {model}, $\sigma_{{Fe}}$={sigma}, $H_N$={nh}, sys_err={sys_err*100} %',fontsize=10)


    ax_eqw.errorbar(phase,eqw*1e3,eqw_err*1e3,color='r',drawstyle='steps-mid',alpha=0.6)
    ax_eqw.set_ylabel(f'Line Eq. width, eV ',fontsize=12)


    uplim_ind=edgeTau_err[0]==edgeTau
    edgeTau[uplim_ind]=edgeTau_err[1][uplim_ind]+edgeTau[uplim_ind]
    edgeTau_err[1][uplim_ind]=0
    edgeTau_err[0][uplim_ind]=edgeTau[uplim_ind]

    ax_edgeTau.errorbar(phase,edgeTau*1e2,edgeTau_err*1e2,color='g',drawstyle='steps-mid',alpha=0.6,uplims=uplim_ind)
    ax_edgeTau.set_ylabel(f'Iron K-edge $\\tau$*100',fontsize=9,color='g')
    ax_edgeTau.set_ylim(-0.2,8)

    # ax_edgeE=ax_edgeTau.twinx()
    # ax_edgeE.errorbar(phase,edgeE,edgeE_err*0,color='c',drawstyle='steps-mid',alpha=0.6)
    # ax_edgeE.set_ylabel(f'Iron K-edge energy \n (=const*Iron line energy)',fontsize=12,color='c')

    ax_po.errorbar(phase,po,po_err,color='b',drawstyle='steps-mid',alpha=0.6)
    ax_po.set_ylabel('Phot. index',color='b')


    ax_po_norm.errorbar(phase,po_norm,po_norm_err,color='m',drawstyle='steps-mid',alpha=0.6)
    ax_po_norm.set_ylabel('cutoffpl norm.', color='m')

    ax_ecut=ax_po.twinx()
    ax_ecut.errorbar(phase,ecut,ecut_err,color='g',drawstyle='steps-mid',alpha=0.6)
    ax_ecut.set_ylabel('$E_{cut}$, keV',color='g')
    ax_ecut.spines["right"].set_position(("axes", 1.1))


    ax_chi2.errorbar(phase,chi2,0,color='k',drawstyle='steps-mid',alpha=0.6)
    ax_chi2.set_ylabel('$\chi^2_{red}$')
    ax_chi2.set_xlabel('Phase',fontsize=12)




    ax_flux= plt.subplot2grid((rows,cols), (0, 8), rowspan=1, colspan=6)
    time=ObsParams.MJD_START
    flux=ObsParams.cutoffpl_cutoffpl_flux/1e-8
    ax_flux.semilogy(time,flux,color='b',marker='s',lw=0,ms=4,alpha=0.8)
    ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b',fontsize=8)
    ax_flux.set_xlabel('Time, MJD')
    ax_flux.yaxis.set_label_position("right")
    ax_flux.yaxis.tick_right()
    ax_flux.axvline(mjd,color='r',ls='-.')


    #fig.tight_layout()

    ax_ratio_1=plt.subplot2grid((rows,cols), (6, 8), rowspan=3, colspan=6)
    plot_spe_ratio(ObsID,  ind_flux_min+1,ind_flux_max+1,ax=ax_ratio_1,color='r')
    plot_spe_ratio(ObsID, ind_cont_flux_min+1,ind_flux_max+1,ax=ax_ratio_1,color='g')
    plot_spe_ratio(ObsID, ind_cont_flux_max+1,ind_flux_max+1,ax=ax_ratio_1,color='m')


    ax_ratio_1.legend(loc='lower right')


    ax_edge=plt.subplot2grid((rows,cols), (2, 8), rowspan=2, colspan=6)
    plot_cutoffpl_ratio(ObsID, phases=[ind_flux_max+1, ind_flux_min+1,ind_cont_flux_min+1,ind_cont_flux_max+1],ax=ax_edge)
    ax_edge.yaxis.set_label_position("right")
    ax_edge.yaxis.tick_right()
    ax_edge.xaxis.set_label_position("top")
    ax_edge.xaxis.tick_top()

    for ax in [ax_ratio_1]:#,ax_ratio_2,ax_edge]:
        ax.axvline(6.4,color='k',lw=0.3,ls=':')
        ax.axvline(7.1,color='k',lw=0.3,ls=':')



    ax_delchi=plt.subplot2grid((rows,cols), (4, 8), rowspan=1, colspan=6,sharex=ax_edge)
    phase_bad_chi=np.argmax(chi2)+1
    os.chdir(f'products/fasebin/{model}/spe_plots')
    delchi_data=np.genfromtxt(f'ph_spe_{phase_bad_chi}_del.dat')
    ax_delchi.errorbar(delchi_data[0],delchi_data[1],1,label=f'phase {phase_bad_chi}')
    ax_delchi.set_xscale('log')
    ax_delchi.grid()
    ax_delchi.yaxis.set_label_position("right")
    ax_delchi.yaxis.tick_right()
    ax_delchi.set_ylabel('(data-model)/error ')


    phase_cont_flux_min=ind_cont_flux_min+1
    delchi_data=np.genfromtxt(f'ph_spe_{phase_cont_flux_min}_del.dat')
    ax_delchi.errorbar(delchi_data[0],delchi_data[1],1,label=f'phase {phase_cont_flux_min}')
    ax_delchi.legend()



    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/{model}/'+f'MJD_{"%.3f"%mjd}_{ObsID}_{model}.png',dpi=100)

    #fig.tight_layout()
    plt.show()
    return [fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]


def plot_spe_ratio(ObsID,phase1,phase2,color='k',ax=None,po_pars=None):
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

    rat,delta=ratio_error(data1[1],data2[1],data1[2],data2[2])
    energy=data1[0]

    ax.errorbar(energy,rat,delta,data1[3],color=color,label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    ax.set_xscale('log')

    ax.legend(loc='upper left',fontsize=8)
    ax.grid('y')
    ax.set_ylabel('spectral ratio (lda)',fontsize=8)



    #from scipy.optimize import curve_fit
    #def cutoffpl(e,gamma,ecut):
    #    return e**(-gamma)*np.exp(-e/ecut)

    #popt,pcov=curve_fit(cutoffpl,energy,rat,sigma=delta,absolute_sigma=1,p0=[0,40])
    #print('cutoffpl: delta G, ecut:',popt)
    #enaxis=np.linspace(energy[0],energy[-1],100)
    #plt.plot(enaxis,cutoffpl(enaxis,*popt),'k:')


    if po_pars is not None:
        en=np.linspace(data1[0][0],data1[0][-1],100)
        sigma=0.3
        eline=6.4
        po1=po_pars[0]*en**(-po_pars[1])*np.exp(-en/po_pars[2])+po_pars[3]*np.exp(-(en-eline)**2/(2*sigma**2))/np.sqrt(sigma**2*2*np.pi)
        po2=po_pars[4]*en**(-po_pars[5])*np.exp(-en/po_pars[6])+po_pars[7]*np.exp(-(en-eline)**2/(2*sigma**2))/np.sqrt(sigma**2*2*np.pi)
        ax.plot(en,po1/po2,'k-.')



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




def calculate_parameters_data(ObsID,model):
    [model,phase,
            mjd,
            tot_flux,
            datamode,
            cts,
            flux312,flux312_err,
            po,po_err,
            ecut,ecut_err,
            po_norm,po_norm_err,
            eline,eline_err,
            norm_line,norm_line_err,
            edgeE,edgeE_err,
            edgeTau,edgeTau_err,
            chi2,
            sigma,
            eqw,eqw_err,nh,sys_err]=read_ph_res_results_data(ObsID,model,0)
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



def roll_fasebin_files(ObsList=['90089-11-03-00G','90089-11-03-01G','90089-11-03-02'],
                       Roll_list=[0,0,-11],folder_name='90089-11-03-00G_group'):
    #if any(np.array(Roll_list)>0):
    #    raise Exception('Only negative roll due to the last spectra equal to the first')
    fig,ax=plt.subplots(3,figsize=(12,12))

    os.system(f'rm -rf /Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}')
    os.system(f'mkdir /Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}')
    os.system(f'mkdir /Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}/products')
    os.system(f'mkdir /Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}/products/fasebin')


    all_cts=[]
    obs_cts=[]
    obs_num_pcu_on=[]
    for ObsID,roll in zip(ObsList,Roll_list):

        os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin')
        os.system('rm -f fasebin_orig_roll.pha')
        os.system('cp fasebin_orig.pha fasebin_orig_roll.pha ')


        fasebin_roll=fits.open('fasebin_orig_roll.pha',mode='update')
        cts=fasebin_roll[1].data['counts'].sum(axis=1)
        obs_cts.append(fasebin_roll[1].data['counts'].sum())
        #ax[0].plot(cts,'+-',label=f'{ObsID} - original')
        ax[0].plot(cts/fasebin_roll[1].data['exposure'][0],'+-',label=f'{ObsID} - original')
        ax[0].set_ylabel('counts/s/ phase_bin')

        tmp_data=fasebin_roll[1].data

        Indeces=np.arange(0,len(tmp_data))
        fasebin_roll[1].data=tmp_data[np.roll(Indeces,roll)]
        fasebin_roll.close()
        fasebin_check=fits.open('fasebin_orig_roll.pha')
        cts=fasebin_check[1].data['counts'].sum(axis=1)
        fasebin_check.close()
        all_cts.append(cts)

        xte_obs=ObservationXTE(ObsID)

        # ax[1].plot(cts-np.mean(cts),'+-',label=f'{ObsID} - rolled')
        # ax[1].set_ylabel('counts-mean')

        ax[1].plot(cts/fasebin_roll[1].data['exposure'][0],'+-',label=f'{ObsID} - rolled ')
        try:
            where_flat=np.where(np.diff(cts)==0)[0]
            ax[1].axvline(where_flat,color='k',lw=1,ls=':')
        except:
            1
        #if any(np.diff(cts)==0):
        #    print(np.diff(cts))
        #    raise Exception(f'Diff=0 somewhere!!, ObsID={ObsID}')
        ax[1].set_ylabel('counts/s')


    all_cts=np.array(all_cts)
    all_cts=np.sum(all_cts,axis=0)
    ax[2].plot(all_cts,'d-',label='counts sum (rolled)')


    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}/products/fasebin')
    path=[f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/fasebin_orig_roll.pha' for ObsID in ObsList]
    path='\n'.join(path)
    print(path)
    os.system(f"echo '{path}' > alldays_pha.list")

    os.system("fbadd infile='@alldays_pha.list' outfile='fasebin_orig.pha'")

    for ObsID in ObsList:
        os.system(f"cp /Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/response.rsp ./response_{ObsID}.rsp")
        os.system(f"cp /Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/fasebin_orig_roll.pha  ./fasebin_orig_roll_{ObsID}.pha")
        os.system(f"cp /Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/dets.png  ./dets_{ObsID}.png")



    #obs_cts=np.array(obs_cts)
    #weights=obs_cts/obs_cts.sum()

    #addrmf_str=''.join(['response_'+x+'.rsp,' for x in ObsList]).rstrip(',')+' '+ ''.join([str(x)+',' for x in weights]).rstrip(',')+' response.rsp'
    #print('creating response: '+'addrmf '+addrmf_str)
    #os.system('addrmf '+addrmf_str)
    os.system(f'cp response_{ObsID}.rsp response.rsp')

    fasebin_joined=fits.open('fasebin_orig.pha')
    cts=fasebin_joined[1].data['counts'].sum(axis=1)
    ax[2].plot(cts,'o',label='fasebin_orig.pha - counts (all observations summed)')


    ax[0].legend()
    ax[1].legend()
    ax[2].legend()


    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}/products/fasebin/pulses.png')
    plt.show()




def xspec_command_edge(dataname):
    return f'''
    query yes
    data {dataname}.pha
    response ../response.rsp
    ign **-3. 12.-**
    system 0.004


    mo (gauss+cutoffpl)*edge
    6.4 -1
    0.3 -1
    1e-3
    -0.5
    6
    0.3
    7.1 -1
    0

    #0 works fine, 6e-2 for the second group suffuce

    renorm
    fit 1000
    fit





    parallel error 10

    error 1. maximum 100 1-8


    set eline [tcloutr par 1]
    set eline_err [tcloutr error 1]

    set sigma [tcloutr par 2]

    set norm_line [tcloutr par 3]
    set norm_line_err [tcloutr error 3]

    set po [tcloutr par 4]
    set po_err [tcloutr error 4]

    set ecut [tcloutr par 5]
    set ecut_err [tcloutr error 5]

    set norm_po [tcloutr par 6]
    set norm_po_err [tcloutr error 6]


    set edgeE [tcloutr par 7]
    set edgeE_err [tcloutr error 7]


    set edgeTau [tcloutr par 8]
    set edgeTau_err [tcloutr error 8]

    if [lindex $edgeTau_err 1]==0 {{
            newpar 8 1e-2
            fit 1000
    error maximum 100 2.71 8
        set edgeTau [tcloutr par 8]
        set edgeTau_err [tcloutr error 8]
    }}


if [lindex $edgeTau_err 0]==0 {{
error maximum 100 2.71 8
    set edgeTau [tcloutr par 8]
    set edgeTau_err [tcloutr error 8]
}}




    try {{
    eqw 1 err 100 68
    set eqw [tcloutr eqw 1]

    }} on error {{result options}} {{
    echo 'eqw err fail'
    eqw 1
    tclout eqwidth 1
    set eqw $xspec_tclout
    set eqw [list [lindex $eqw 0] [lindex $eqw 0] [lindex $eqw 0]]
    }}

        steppar log 8 1e-4 7e-2 250

        tcloutr plot contour x
        tcloutr plot contour


        set f [open ./{dataname}_edge_steppar.dat w]
        setpl en
        puts $f [tcloutr plot contour  x ]
        puts $f [tcloutr plot contour y ]

        close $f



        set f [open ./{dataname}_lda.dat w]
        setpl en
        puts $f [tcloutr plot lda  x ]
        puts $f [tcloutr plot lda y ]
        puts $f [tcloutr plot lda yerr ]
        puts $f [tcloutr plot lda xerr ]

        close $f


    #######
    ### additional data
    #######

    tclout dof
    set dof $xspec_tclout

    tclout stat
    set stat $xspec_tclout

    set dof [lindex $dof 0]

    set chi_r [format "%.3f" [expr {{$stat/$dof}}]]



    fit




    flux 3 12
    tclout flux 1
    set fluxx $xspec_tclout
    set fluxx [lindex $fluxx 0]
    set flux_init [expr {{log10($fluxx)}}]

    freeze 3,6



    add 2 cflux
    3
    12
    -10



    newpar 6 $flux_init
    fit

    error  6





    ##############
    ### FLUX CUTOFF 3-12
    ##############
    set parname cutoff312
    set parind 6



    #function for flux calculation for some component
    tclout par $parind
    set logF $xspec_tclout
    set logF [lindex $logF 0]

    set mean [expr {{pow(10,$logF)}}]

    tclout error $parind
    set logFerr $xspec_tclout

    set lo [lindex $logFerr 0]
    set lo [expr {{pow(10, $lo)}}]

    set hi [lindex $logFerr 1]
    set hi [expr {{pow(10,$hi)}}]

    set err_lo [expr {{$mean-$lo}}]
    set err_hi [expr {{$hi-$mean}}]

    #for error propagation
    set e${{parname}} [expr {{max($err_lo,$err_hi)}}]

    #for puts in file
    set ${{parname}}_mean $mean
    set ${{parname}}_lo $lo
    set ${{parname}}_hi $hi





    set fileid [open ./{dataname}_edge_all.dat a]

    puts $fileid "0 0 $chi_r $dof   [lindex $cutoff312_mean 0] [lindex $cutoff312_lo 0]  [lindex $cutoff312_hi 0]  [lindex $po 0] [lindex $po_err 0] [lindex $po_err 1] [lindex $ecut 0] [lindex $ecut_err 0] [lindex $ecut_err 1] [lindex $norm_po 0] [lindex $norm_po_err 0] [lindex $norm_po_err 1] [lindex $eline 0] [lindex $eline_err 0] [lindex $eline_err 1] [lindex $norm_line 0] [lindex $norm_line_err 0] [lindex $norm_line_err 1] [lindex $edgeE 0] [lindex $edgeE_err 0] [lindex $edgeE_err 1] [lindex $edgeTau 0] [lindex $edgeTau_err 0] [lindex $edgeTau_err 1] [lindex $sigma 0] [lindex $eqw 0] [lindex $eqw 1] [lindex $eqw 2] 0 0.04"

    close $fileid


    setpl en
    cpd {dataname}_edge.gif/VGIF
    plot eeufs del ra
    cpd none
    rm {dataname}_edge.gif
    mv {dataname}_edge.gif_2 {dataname}_edge.gif
    rm -f {dataname}_edge.ps


    set fileid [open ./{dataname}_edge.dat a]

    puts $fileid "$chi_r $dof  [lindex $edgeTau 0] [lindex $edgeTau_err 0] [lindex $edgeTau_err 1] "
    exit
'''


def xspec_command_smedge(dataname):
    return f'''
    data {dataname}.pha
    response ../response.rsp
    ign **-3. 12.-**
    system 0.003

    mo gauss+cutoffpl*smedge
    6.4 -1
    0.3 -1
    1e-3
    -0.5
    6
    0.3
    7.1 -1
    0
    -2.67000  -1
    0.3 -1


    renorm
    fit 1000

    error maxumum 100 3.84 8

    set edgeTau [tcloutr par 8]
    set edgeTau_err [tcloutr error 8]

    tclout dof
    set dof $xspec_tclout

    tclout stat
    set stat $xspec_tclout

    set dof [lindex $dof 0]

    set chi_r [format "%.3f" [expr {{$stat/$dof}}]]


    setpl en
    cpd {dataname}_edge.gif/VGIF
    plot eeufs del ra
    cpd none
    rm {dataname}_edge.gif
    mv {dataname}_edge.gif_2 {dataname}_edge.gif
    rm -f {dataname}_edge.ps


    set fileid [open ./{dataname}_edge.dat a]

    puts $fileid "$chi_r $dof  [lindex $edgeTau 0] [lindex $edgeTau_err 0] [lindex $edgeTau_err 1] "
    exit
'''



def xspec_command_no_edge(dataname):
    return f'''
    data {dataname}.pha
    response ../response.rsp
    ign **-3. 12.-**
    system 0.003

    mo gauss+cutoffpl
    6.4 -1
    0.3 -1
    1e-3
    -0.5
    6
    0.3


    renorm
    fit 1000

    tclout dof
    set dof $xspec_tclout

    tclout stat
    set stat $xspec_tclout

    set dof [lindex $dof 0]

    set chi_r [format "%.3f" [expr {{$stat/$dof}}]]


    setpl en
    cpd {dataname}_no_edge.gif/VGIF
    plot eeufs del ra
    cpd none
    rm {dataname}_no_edge.gif
    mv {dataname}_no_edge.gif_2 {dataname}_no_edge.gif
    rm -f {dataname}_no_edge.ps

    set fileid [open ./{dataname}_no_edge.dat a]

    puts $fileid "$chi_r $dof  0 0 0 "


    del 1
    fit 1000

    set f [open ./{dataname}_no_gauss_ra.dat w]
    setpl en
    puts $f [tcloutr plot ra  x ]
    puts $f [tcloutr plot ra y ]
    puts $f [tcloutr plot ra yerr ]
    puts $f [tcloutr plot ra xerr ]

    close $f





    exit

'''

def average_phase_bins_addspec(ObsID,indeces,
                       names,colors,
                       folder='average_bins'):
    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin')
    os.system(f'rm -rf {folder}')
    os.mkdir(f'{folder}')
    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{folder}')

    fasebin_file=fits.open('../fasebin_orig.pha')

    fasebin_spectra=fasebin_file[1].data['counts']
    fasebin_counts=fasebin_file[1].data['counts'].sum(axis=1)/1e7

    fig,ax=plt.subplots(figsize=(10,10))

    rows=5
    cols=1
    #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
    ax = plt.subplot2grid((rows,cols), (0, 0), rowspan=3, colspan=1)
    ax_steppar = plt.subplot2grid((rows,cols), (3, 0), rowspan=2, colspan=1)



    bin_num=np.arange(1,len(fasebin_counts)+1)
    ax.plot(bin_num,fasebin_counts,drawstyle='steps-mid',marker='o',color='k')
    ax.set_xlabel('bin number ')
    ax.set_ylabel('counts in bin, 1e7 cts')

    if ObsID=='90089-11-03-00G_group':
        x=bin_num[np.array([1,2,3,4,5,10,11,12,13,18,19,20,21])-1]+0.
        y=fasebin_counts[np.array([1,2,3,4,5,10,11,12,13,18,19,20,21])-1]+0.
        def cosine(x,A,B,C,D):
            return A+B*np.cos((x-C)/D)
        popt,pcov=curve_fit(cosine,x,y,p0=[1.2,0.1,0,3])
        xaxis=np.linspace(0,32,50)
        ax.plot(xaxis,cosine(xaxis,*popt),'m:')
        ax.plot(x,y,'r.')
        plt.show()


    for i in range(int(len(fasebin_counts)/2)):
        ax.text(bin_num[i],min(fasebin_counts)*0.9,str(i+1),fontsize=8)
        ax.axvline(bin_num[i]+0.5,ls=':',lw=0.4,color='k')

    ax.set_ylim(min(fasebin_counts)*0.87,max(fasebin_counts)*1.05)

    for ind,name,color in zip(indeces,names,colors):
        ind=np.array(ind)-1
        for k,i in enumerate(ind):
            if k==0:
                ax.axvspan(bin_num[i]-0.5,bin_num[i]+0.5, alpha=0.2,zorder=-10,color=color,label=name)
            else:
                ax.axvspan(bin_num[i]-0.5,bin_num[i]+0.5, alpha=0.2,zorder=-10,color=color)

    ax_edge=ax.twinx()

    for ind,name,color in zip(indeces,names,colors):
        print(ind,name,color)
        print(indeces,names,colors)
        fig,ax=plt.subplots(figsize=(10,6))

        sum_spectra=[]
        for phase_num in ind:
            phase_num=phase_num
            cmp=f'cmppha infile=../fasebin_orig.pha outfile={name}_phase_{phase_num}.tmp_pha row={phase_num} cmpmode=expand'
            print(cmp)
            os.system(cmp)

            spe_file=fits.open(f'{name}_phase_{phase_num}.tmp_pha')
            spectra=spe_file[1].data['counts']
            ax.plot(spectra,label='phase spectra: '+name+', phase '+str(phase_num))
            sum_spectra.append(spectra)

        sum_spectra=np.array(sum_spectra).sum(axis=0)
        ax.loglog(sum_spectra,'ko-',label='sum of spectra')


        os.system(f'ls {name}*.tmp_pha > list_{name}.txt')

        #ls_of_spectra=glob('hump*.tmp_pha')

        add_spec_humps=f'addspec infil=list_{name}.txt  outfil=all_{name} qaddrmf=no qsubback=no errmeth=Gauss'
        os.system(add_spec_humps)

        addspe_file=fits.open(f'all_{name}.pha')
        addspe_spectra=addspe_file[1].data['counts']
        ax.plot(addspe_spectra,label='addspe_spectra')

        assert all(addspe_spectra==sum_spectra), 'addspe spectra and sum of phases are not equal'

        ax.legend()
        ax.set_xlabel('channel')
        ax.set_ylabel('counts')
        plt.tight_layout()
        plt.show()
        plt.savefig(f'{name}_counts_check.png')
        plt.close(fig)
        xsp_cmd=xspec_command_edge(f'all_{name}')
        print(xsp_cmd)
        file = open(f"xspec_script_edge_{name}.txt", "w")
        file.write(xsp_cmd)
        file.close()
        os.system(f'xspec - xspec_script_edge_{name}.txt')
        xsp_cmd=xspec_command_no_edge(f'all_{name}')
        file = open(f"xspec_script_no_edge_{name}.txt", "w")
        file.write(xsp_cmd)
        file.close()
        os.system(f'xspec - xspec_script_no_edge_{name}.txt')

        xspec_data_edge=np.genfromtxt(f'all_{name}_edge.dat')
        chi2_r,_,edge,edge_lo,edge_hi=xspec_data_edge
        uplim=edge_lo==0
        edge_lo=edge-edge_lo
        edge_hi=edge_hi-edge

        edge*=100
        edge_lo *=100
        edge_hi *=100

        print('+++++++++++')
        print(f'PHASE SPAN: {name} \n  EDGE TAU * 10000={100*edge}(+{100*edge_hi}, - {100*edge_lo})')
        print('+++++++++++')

        import more_itertools as mit
        for group in mit.consecutive_groups(ind):
            group=list(group)
            group_err=(group[-1]-group[0])/2+0.5


            group=np.mean(group)

            if uplim:
                _,_,edge,edge_lo,edge_hi=xspec_data_edge
                edge=edge_hi
                edge_lo=edge_hi
                edge_hi=0
                edge*=100
                edge_lo *=100
                edge_hi *=100
                edge_err=np.vstack((edge_lo,edge_hi))
                print(edge,edge_err)
            else:
                edge_err=np.vstack((edge_lo,edge_hi))



            ax_edge.plot(group,edge,color=color,marker='d',)
            if not uplim:
                ax_edge.text(group,edge*1.1,f'{name}\n chi2_red={chi2_r}, \n 68% c.l. \n {"%.1f"%edge}\n -{"%.1f"%edge_err[0]}+{"%.1f"%edge_err[1]}',fontsize=8)
            else:
                ax_edge.text(group,edge*1.1,f'{name}\n chi2_red={chi2_r}, \n 90% c.l. \n {"%.1f"%edge}\n -{"%.1f"%edge_err[0]}+{"%.1f"%edge_err[1]}',fontsize=8)
            ax_edge.errorbar(group,edge,edge_err,group_err,fmt='none',ecolor=color,marker='d',uplims=uplim)

    ax_edge.set_ylabel('Iron K edge tau * 100')
    ax_edge.set_ylim(-0.1,8)
    ax.legend()
    plt.title(f'{ObsID}, {folder}')
    #plt.show()
    #plt.savefig('phase_regions.png')


    #fig,ax_steppar=plt.subplots(figsize=(8,4))

    for name in names:
        steppar_data=np.genfromtxt(f'all_{name}_edge_steppar.dat').T
        tau_dips=steppar_data[:,0]*100
        dchi2_dips=steppar_data[:,1]-min(steppar_data[:,1])
        ax_steppar.semilogx(tau_dips,dchi2_dips,'-',label=name)

    for cl,sigm in zip([0.682689,0.9,0.95,0.9973,0.999936],['(1 sigma)','','','(3 sigma)','(4 sigma)']):
        dc=stats.chi2.ppf(cl,df=1)
        ax_steppar.axhline(dc,color='k',ls=':',lw=1)
        cl=cl*100
        ax_steppar.text(0.3,dc,f'{"%.4f"%cl} % cl {sigm}')

    ax_steppar.set_ylim(0,20)
    ax_steppar.set_xlim(0.1,10)
    ax_steppar.grid()
    ax_steppar.legend()
    ax_steppar.set_ylabel('delta chi2')
    ax_steppar.set_xlabel('Edge Tau * 100')
    #plt.title(f'{ObsID}, {folder}')
    plt.tight_layout()
    plt.show()
    plt.savefig('confidence_levels.png')
    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{ObsID}_cutoffpl_edge_average_phases.png',dpi=100)


    return


def add_average_analysis_on_plots(ObsID,folder,phases,names,colors,axes,n_ph=16):
    for name,color,phase in zip(names,colors,phases):
        os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{folder}')
        phase=(np.array(phase)-1)/n_ph
        phase_err=len(phase)/(n_ph)/2
        phase=np.mean(phase)
        data=np.genfromtxt(f'all_{name}_edge_all.dat')
        data=data.reshape(1,len(data))
        nph=data[0,1]
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


        po_norm=data[:,13]
        po_norm_lo=po_norm-data[:,14]
        po_norm_hi=data[:,15]-po_norm
        po_norm_err=np.vstack((po_norm_lo,po_norm_hi))

        eline=data[:,16]
        eline_lo=eline-data[:,17]
        eline_hi=data[:,18]-eline
        eline_err=np.vstack((eline_lo,eline_hi))

        norm_line=data[:,19]
        norm_line_lo=norm_line-data[:,20]
        norm_line_hi=data[:,21]-norm_line
        norm_line_err=np.vstack((norm_line_lo,norm_line_hi))

        edgeE=data[:,22]
        edgeE_lo=edgeE-data[:,23]
        edgeE_hi=data[:,24]-edgeE
        edgeE_err=np.vstack((edgeE_lo,edgeE_hi))

        edgeTau=data[:,25]
        edgeTau_lo=edgeTau-data[:,26]
        edgeTau_hi=data[:,27]-edgeTau
        edgeTau_err=np.vstack((edgeTau_lo,edgeTau_hi))

        sigma=data[:,28]

        eqw=data[:,29]
        eqw_lo=eqw-data[:,30]
        eqw_hi=data[:,31]-eqw
        eqw_err=np.vstack((eqw_lo,eqw_hi))

        nh=data[:,32]

        sys_err=data[:,33]


        # uplim_ind=edgeTau==0
        # edgeTau[uplim_ind]=edgeTau_err[1][uplim_ind]
        # edgeTau_err[1][uplim_ind]=0
        # edgeTau_err[0][uplim_ind]=edgeTau[uplim_ind]

        uplim_ind=edgeTau_err[0]==edgeTau
        edgeTau[uplim_ind]=edgeTau_err[1][uplim_ind]+edgeTau[uplim_ind]
        edgeTau_err[1][uplim_ind]=0
        edgeTau_err[0][uplim_ind]=edgeTau[uplim_ind]


        axes[0].errorbar(phase,edgeTau*1e2,edgeTau_err*1e2,phase_err,color=color,drawstyle='steps-mid',alpha=0.6,lw=4,uplims=uplim_ind)
        axes[1].errorbar(phase,eqw*1e3,eqw_err*1000,phase_err,color=color,drawstyle='steps-mid',alpha=0.6,lw=4)

        axes[2].errorbar(phase,norm_line*1e3,norm_line_err*1000,phase_err,color=color,drawstyle='steps-mid',alpha=0.6,lw=4)

        axes[4].errorbar(phase,flux312/1e-8,flux312_err/1e-8,phase_err,color='gray',drawstyle='steps-mid',alpha=0.3,lw=2)

        ax_ratio=axes[3]

    data1=np.genfromtxt(f'all_dips_lda.dat')

    data2=np.genfromtxt(f'all_off_dip_lda.dat')
    label=f" dip / off-dip"
    def ratio_error(a,b,da,db):
        f=a/b
        sigma=np.abs(f)*np.sqrt( (da/a)**2 + (db/b)**2  )
        return f, sigma

    # rat,delta=ratio_error(data1[1],data2[1],data1[2],data2[2])
    # energy=data1[0]

    # ax_ratio.errorbar(energy,rat,delta,data1[3],color='k',label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    # ax_ratio.legend()



def plot_cutoffpl_ratio_for_group(ObsID, phases=['dips','off_dip'],ax=None):
    model='average_phase_all_one_dip_all_off_dip'
    if ax==None:
        fig,ax = plt.subplots(figsize=(8, 6))
    else:
        pass

    for phase in phases:

        data1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{model}/all_{phase}_no_gauss_ra.dat')

        label=f" {phase} / cutoffpl "

        ax.errorbar(data1[0],data1[1],data1[2],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
        ax.set_xscale('log')

    ax.legend(loc='upper left',fontsize=8)
    ax.grid('y')
    ax.axhline(1,color='k',ls=':',zorder=-10,alpha=0.6)
    ax.set_ylabel(' ratio ',fontsize=8)




STOP

#%% run sa observations with edge

err=[]
msg=[]
plt.ioff()
matplotlib.use('Agg')

if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_SA
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            plot_ph_res_results(ObsID=ObsID,model='cutoffpl_edge')
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
matplotlib.use('Agg')
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_SA+ObsList_SE
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:

            plot_ph_res_results(ObsID=ObsID,model='cutoffpl')
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
#plt.ioff()
#matplotlib.use('Agg')

if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_SA+ObsList_SE
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            calculate_parameters_data(ObsID=ObsID,model='cutoffpl')


        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg




#%% fit spectra for observation in groups 1-4
ObsList=['90089-11-03-02','90089-11-03-01G','90089-11-03-00G']+['90427-01-03-00','90427-01-03-01','90427-01-03-02']+['90427-01-03-14G','90014-01-02-00','90427-01-03-05']+['90427-01-03-06','90427-01-03-07','90014-01-02-08','90014-01-02-10']

msg=[]
for k,ObsID in enumerate(ObsList):
    print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
    xte_obs=ObservationXTE(ObsID)
    #xte_obs.fit_ph_res(model='cutoffpl_no_gauss',error=0)
    #xte_obs.fit_ph_res(model='cutoffpl',chmin=6,chmax=8,error=0.00)
    try:
        #xte_obs.fit_ph_res(model='cutoffpl_edge',error=0)
        plot_ph_res_results(ObsID,'cutoffpl_edge')
    except:
        msg.append(f'Fail, {ObsID}')
print(msg)


#%% for top-phase observations group 1, 16 bins
groupid='group_1'

# roll_fasebin_files(ObsList=['90089-11-03-00G','90089-11-03-01G','90089-11-03-02'],Roll_list=[0,-5,-10],
#                     folder_name='out'+groupid)
# os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{groupid}/products/fasebin')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_no_gauss.txt ')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl.txt ')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_edge.txt ')
[fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]=plot_ph_res_results(groupid,'cutoffpl_edge')

plot_ph_res_results(groupid,'cutoffpl')


avg_bins_res=average_phase_bins_addspec(ObsID=groupid,
                            indeces=[[2,3,4,5,10,11,12,13],[7,8,15,16]],
                            names=['off_dip','dips'],
                            colors=['b','r'],
                            folder='average_phase_all_no_borders',
                            )

add_average_analysis_on_plots(groupid,'average_phase_all_no_borders',
                              [[7,8],[15,16],[2,3,4,5],[10,11,12,13]],
                              ['dips','dips','off_dip','off_dip'],
                              ['r','r','b','b'], [ax_edgeTau,ax_eqw,ax_norm,ax_ratio_1,ax_efold])
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{groupid}_cutoffpl_edge_average_phases_with_results.png',dpi=100)

#%% for decline phase observations: group 2, 16 bins
groupid='group_2'

# roll_fasebin_files(ObsList=['90427-01-03-00','90427-01-03-01','90427-01-03-02'],Roll_list=[0,0,-8],
#                   folder_name='out'+groupid)

# os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{groupid}/products/fasebin')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_no_gauss.txt ')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl.txt ')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_edge.txt ')
[fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]=plot_ph_res_results(groupid,'cutoffpl_edge')




average_phase_bins_addspec(ObsID=groupid,
                            indeces=[[1,2,3,4,5,6,7,8,9,10,11],[13,14,15]],
                            names=['off_dip','dips'],
                            colors=['b','r'],
                            folder='average_phase_all_one_dip_all_off_dip')

add_average_analysis_on_plots(groupid,'average_phase_all_one_dip_all_off_dip',
                              [[13,14,15],[1,2,3,4,5,6,7,8,9,10,11]],
                              ['dips','off_dip'],
                              ['r','b','b'], [ax_edgeTau,ax_eqw,ax_norm,ax_ratio_1,ax_efold])
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{groupid}_cutoffpl_edge_average_phases_with_results.png',dpi=100)



#%% for decline phase observations: group 3, 16 bins
groupid='group_3'
# roll_fasebin_files(ObsList=['90427-01-03-14G','90014-01-02-00','90427-01-03-05'],Roll_list=[0,-10,-12],
#                   folder_name='out'+groupid)

# os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{groupid}/products/fasebin')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_no_gauss.txt ')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl.txt ')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_edge.txt ')
[fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]=plot_ph_res_results(groupid,'cutoffpl_edge')


average_phase_bins_addspec(groupid,
                            indeces=[[1,2,3,4,5,6,7,8],[10,11,12]],
                            names=['off_dip','dips'],
                            colors=['b','r'],
                            folder='average_phase_all_one_dip_all_off_dip')

add_average_analysis_on_plots(groupid,'average_phase_all_one_dip_all_off_dip',
                              [[10,11,12],[1,2,3,4,5,6,7,8]],
                              ['dips','off_dip','off_dip'],
                              ['r','b','b'], [ax_edgeTau,ax_eqw,ax_norm,ax_ratio_1,ax_efold])
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{groupid}_cutoffpl_edge_average_phases_with_results.png',dpi=100)






#%%for decline phase observations: group 4, 16 bins

groupid='group_4'

roll_fasebin_files(ObsList=['90427-01-03-06','90427-01-03-07','90014-01-02-08'],
                    Roll_list=[0,-8,-9,],folder_name='out'+groupid)


os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{groupid}/products/fasebin')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_no_gauss.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_edge.txt ')
[fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]=plot_ph_res_results(groupid,'cutoffpl_edge')


average_phase_bins_addspec(groupid,
                            indeces=[[6,7,8,9,10,11,12,13,14,15],[2,3,4]],
                            names=['off_dip','dips'],
                            colors=['b','r'],
                            folder='average_phase_all_one_dip_all_off_dip')

add_average_analysis_on_plots(groupid,'average_phase_all_one_dip_all_off_dip',
                              [[2,3,4],[6,7,8,9,10,11,12,13,14,15]],
                              ['dips','off_dip','off_dip'],
                              ['r','b','b'], [ax_edgeTau,ax_eqw,ax_norm,ax_ratio_1,ax_efold])
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{groupid}_cutoffpl_edge_average_phases_with_results.png',dpi=100)







#%%for decline phase observations: group 5

groupid='group_5'
roll_fasebin_files(ObsList=['90427-01-03-09','90427-01-03-11','90427-01-03-12'],
                     Roll_list=[0,-11,-6,-4],folder_name='out'+groupid)

os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{groupid}/products/fasebin')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_no_gauss.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_edge.txt ')
[fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]=plot_ph_res_results(groupid,'cutoffpl_edge')


average_phase_bins_addspec(groupid,
                            indeces=[[1,2,3,4,5,6,7,8,9],[11,12,13,14]],
                            names=['off_dip','dips'],
                            colors=['b','r'],
                            folder='average_phase_all_one_dip_all_off_dip')

add_average_analysis_on_plots(groupid,'average_phase_all_one_dip_all_off_dip',
                              [[11,12,13,14],[1,2,3,4,5,6,7,8,9]],
                              ['dips','off_dip'],
                              ['r','b'], [ax_edgeTau,ax_eqw,ax_norm,ax_ratio_1,ax_efold])
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{groupid}_cutoffpl_edge_average_phases_with_results.png',dpi=100)











#%% for decline phase observations: group 6



groupid='group_6'
roll_fasebin_files(ObsList=['90014-01-02-13','90014-01-03-00','90014-01-03-01'],Roll_list=[0,-4,-2],folder_name='out'+groupid)

os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{groupid}/products/fasebin')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_no_gauss.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_edge.txt ')
[fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]=plot_ph_res_results(groupid,'cutoffpl_edge')




average_phase_bins_addspec(groupid,
                            indeces=[[10,11,12,13,14,15,16],[5,6,7,8]],
                            names=['off_dip','dips'],
                            colors=['b','r'],
                            folder='average_phase_all_one_dip_all_off_dip')

add_average_analysis_on_plots(groupid,'average_phase_all_one_dip_all_off_dip',
                              [[5,6,7,8],[10,11,12,13,14,15,16]],
                              ['dips','off_dip'],
                              ['r','b'], [ax_edgeTau,ax_eqw,ax_norm,ax_ratio_1,ax_efold])

fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{groupid}_cutoffpl_edge_average_phases_with_results.png',dpi=100)






#%% for decline phase observations: group 7
groupid='group_7'
roll_fasebin_files(ObsList=['90014-01-03-020','90014-01-03-02','90014-01-03-03'],Roll_list=[0,-5,-4],folder_name='out'+groupid)

os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{groupid}/products/fasebin')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_no_gauss.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_edge.txt ')
[fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]=plot_ph_res_results(groupid,'cutoffpl_edge')




average_phase_bins_addspec(groupid,
                            indeces=[[1,2,3,4,5,6,12,13,14,15,16],[8,9,10]],
                            names=['off_dip','dips'],
                            colors=['b','r'],
                            folder='average_phase_all_one_dip_all_off_dip')

add_average_analysis_on_plots(groupid,'average_phase_all_one_dip_all_off_dip',
                              [[8,9,10],[1,2,3,4,5,6],[12,13,14,15,16]],
                              ['dips','off_dip','off_dip'],
                              ['r','b','b'], [ax_edgeTau,ax_eqw,ax_norm,ax_ratio_1,ax_efold])

fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{groupid}_cutoffpl_edge_average_phases_with_results.png',dpi=100)


#%% 32 bins stuff


#%% make 32 phase fasebin for group 1-4
ObsList=['90089-11-03-02','90089-11-03-01G','90089-11-03-00G']+['90427-01-03-00','90427-01-03-01','90427-01-03-02']+['90427-01-03-14G','90014-01-02-00','90427-01-03-05']+['90427-01-03-06','90427-01-03-07','90014-01-02-08','90014-01-02-10']

for k,ObsID in enumerate(ObsList):
    print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
    xte_obs=ObservationXTE(ObsID)
    xte_obs.make_fasebin(nph=16)
    # xte_obs.fit_ph_res(model='cutoffpl_no_gauss',error=0)
    # xte_obs.fit_ph_res(model='cutoffpl',chmin=6,chmax=8,error=0.00)
    # xte_obs.fit_ph_res(model='cutoffpl_edge',error=0)

    # plot_ph_res_results(ObsID,'cutoffpl_edge')

#%% for top-phase observations group 1,  32 bins

groupid='group_1'

#roll_fasebin_files(ObsList=['90089-11-03-00G','90089-11-03-01G','90089-11-03-02'],Roll_list=[0,-10,-20],
#                  folder_name='out'+groupid)
# os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{groupid}/products/fasebin')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_no_gauss.txt ')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl.txt ')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_edge.txt ')
[fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]=plot_ph_res_results(groupid,'cutoffpl_edge')

# avg_bins_res=average_phase_bins_addspec(ObsID=groupid,
#                             indeces=[[1,2,3,4,5,6,7,8,9,10,11,17,18,19,20,21,22,23,24],[13,14,15,26,27,28,29,30,31,32]],
#                             names=['off_dip','dips'],
#                             colors=['b','r'],
#                             folder='average_phase_all_no_borders',
#                             )

add_average_analysis_on_plots(groupid,'average_phase_all_no_borders',
                              [[13,14,15],[26,27,28,29,30,31,32],[1,2,3,4,5,6,7,8,9,10,11],[17,18,19,20,21,22,23,24]],
                              ['dips','dips','off_dip','off_dip'],
                              ['r','r','b','b'], [ax_edgeTau,ax_eqw,ax_norm,ax_ratio_1,ax_efold],n_ph=32)
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{groupid}_cutoffpl_edge_average_phases_with_results.png',dpi=100)

# average_phase_bins_addspec(ObsID='90089-11-03-00G_group',
#                             indeces=[[2,10,11,12,13],[7,8,15,16]],
#                             names=['humps','dips'],
#                             colors=['b','r'],
#                             folder='average_phase_all_no_borders_skip_badchi2')


#%% for decline phase observations: group 2, 32 bins

groupid='group_2'

roll_fasebin_files(ObsList=['90427-01-03-00','90427-01-03-01','90427-01-03-02'],Roll_list=[0,-1,-18],
                  folder_name='out'+groupid)

os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{groupid}/products/fasebin')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_no_gauss.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_edge.txt ')
[fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]=plot_ph_res_results(groupid,'cutoffpl_edge')




average_phase_bins_addspec(ObsID=groupid,
                            indeces=[np.arange(1,21),np.arange(23,30)],
                            names=['off_dip','dips'],
                            colors=['b','r'],
                            folder='average_phase_all_one_dip_all_off_dip')

add_average_analysis_on_plots(groupid,'average_phase_all_one_dip_all_off_dip',
                              [np.arange(23,30),np.arange(1,21)],
                              ['dips','off_dip'],
                              ['r','b','b'], [ax_edgeTau,ax_eqw,ax_norm,ax_ratio_1,ax_efold],n_ph=32)
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{groupid}_cutoffpl_edge_average_phases_with_results.png',dpi=100)



#%% for decline phase observations: group 3, 32 bins
groupid='group_3'
#roll_fasebin_files(ObsList=['90427-01-03-14G','90014-01-02-00','90427-01-03-05'],Roll_list=[0,-20,-24],
#                  folder_name='out'+groupid)

# os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{groupid}/products/fasebin')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_no_gauss.txt ')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl.txt ')
# os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_edge.txt ')
[fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]=plot_ph_res_results(groupid,'cutoffpl_edge')


average_phase_bins_addspec(groupid,
                            indeces=[np.arange(1,17),np.arange(18,26)],
                            names=['off_dip','dips'],
                            colors=['b','r'],
                            folder='average_phase_all_one_dip_all_off_dip')

add_average_analysis_on_plots(groupid,'average_phase_all_one_dip_all_off_dip',
                              [np.arange(18,26),np.arange(1,17)],
                              ['dips','off_dip','off_dip'],
                              ['r','b','b'], [ax_edgeTau,ax_eqw,ax_norm,ax_ratio_1,ax_efold],n_ph=32)
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{groupid}_cutoffpl_edge_average_phases_with_results.png',dpi=100)



#%%for decline phase observations: group 4, 32 bins

groupid='group_4'

#roll_fasebin_files(ObsList=['90427-01-03-06','90427-01-03-07','90014-01-02-08'],
#                    Roll_list=[0,-15,-17,],folder_name='out'+groupid)


os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{groupid}/products/fasebin')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_no_gauss.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl.txt ')
os.system(f'xspec -  ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res/ph_res_cutoffpl_edge.txt ')
[fig,ax_norm,ax_eqw,ax_edgeTau,ax_chi2,ax_ratio_1,ax_efold]=plot_ph_res_results(groupid,'cutoffpl_edge')


average_phase_bins_addspec(groupid,
                            indeces=[np.arange(12,32),np.arange(1,10)],
                            names=['off_dip','dips'],
                            colors=['b','r'],
                            folder='average_phase_all_one_dip_all_off_dip')

add_average_analysis_on_plots(groupid,'average_phase_all_one_dip_all_off_dip',
                              [np.arange(1,10),np.arange(12,32)],
                              ['dips','off_dip','off_dip'],
                              ['r','b','b'], [ax_edgeTau,ax_eqw,ax_norm,ax_ratio_1,ax_efold],n_ph=32)
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/groups/'+f'{groupid}_cutoffpl_edge_average_phases_with_results.png',dpi=100)




#%% plot for paper: groups and observations

def load_and_plot_data_for_paper(ObsID,ax_eqw,ax_norm,ax_edgeTau,model='cutoffpl_edge',phases=None):
    fontsizes=9

    [model,phase,
            mjd,
            tot_flux,
            datamode,
            cts,
            flux312,flux312_err,
            po,po_err,
            ecut,ecut_err,
            po_norm,po_norm_err,
            eline,eline_err,
            norm_line,norm_line_err,
            edgeE,edgeE_err,
            edgeTau,edgeTau_err,
            chi2,
            sigma,
            eqw,eqw_err,nh,sys_err]=read_ph_res_results_data(ObsID,model,0)

    for ax in [ax_eqw,ax_edgeTau,ax_norm]:
        ax_efold=ax.twinx()
        ax_efold.errorbar(phase,flux312/1e-8,flux312_err/1e-8,color='k',drawstyle='steps-mid',ls=':',alpha=1)
        ax_efold.set_ylabel('Flux 3-12 keV, 10$^{-8} $cgs',fontsize=9)


    ax_norm.errorbar(phase,norm_line*1e3,norm_line_err*1e3,color='green',drawstyle='steps-mid',alpha=0.6)
    ax_norm.set_ylabel(f'Iron line norm, \n ($10^{-3}$ $ph cm^{-2}s^{-1}$)',fontsize=12)
    #ax_norm.fill_between(phase, (norm_line-norm_line_err[0])*1e3,(norm_line+norm_line_err[1])*1e3, color='g',alpha=0.1)
    title_ID=ObsID
    if title_ID=='group_1':
        title_ID='group i'
    elif title_ID=='group_2':
        title_ID='group ii'
    elif title_ID=='group_3':
        title_ID='group iii'
    elif title_ID=='group_4':
        title_ID='group iv'
    elif title_ID=='group_5':
        title_ID='group v'
    elif title_ID=='group_6':
        title_ID='group vi'


    ax_norm.set_title(title_ID,fontsize=15)


    ax_eqw.errorbar(phase,eqw*1e3,eqw_err*1e3,color='b',drawstyle='steps-mid',alpha=0.6)
    ax_eqw.set_ylabel(f'Line Eq. width, eV ',fontsize=12)


    uplim_ind=edgeTau_err[0]==edgeTau
    edgeTau[uplim_ind]=edgeTau_err[1][uplim_ind]+edgeTau[uplim_ind]
    edgeTau_err[1][uplim_ind]=0
    edgeTau_err[0][uplim_ind]=edgeTau[uplim_ind]

    ax_edgeTau.errorbar(phase,edgeTau*1e2,edgeTau_err*1e2,color='r',drawstyle='steps-mid',alpha=0.6,uplims=uplim_ind)
    ax_edgeTau.set_ylabel(f'Iron K-edge $\\tau$*100',fontsize=12,color='k')
    #ax_edgeTau.set_ylim(-0.2,8)
    ax_edgeTau.set_xlabel('Phase',fontsize=12)

    clrs=['r','b']
    if phases!=None:
        for clr,part in zip(clrs,phases):
            for ph in part:
                ph=ph-1
                ax_edgeTau.axvspan(phase[ph]-0.5/16,phase[ph]+0.5/16, alpha=0.2,zorder=-10,color=clr)




def add_average_analysis_on_plots_for_paper(ObsID,folder,phases,names,colors,axes,n_ph=16):
    for name,color,phase in zip(names,colors,phases):
        os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{folder}')
        phase=(np.array(phase)-1)/n_ph
        phase_err=len(phase)/(n_ph)/2
        phase=np.mean(phase)
        data=np.genfromtxt(f'all_{name}_edge_all.dat')
        data=data.reshape(1,len(data))

        edgeTau=data[:,25]
        edgeTau_lo=edgeTau-data[:,26]
        edgeTau_hi=data[:,27]-edgeTau
        edgeTau_err=np.vstack((edgeTau_lo,edgeTau_hi))

        uplim_ind=edgeTau_err[0]==edgeTau
        edgeTau[uplim_ind]=edgeTau_err[1][uplim_ind]+edgeTau[uplim_ind]
        edgeTau_err[1][uplim_ind]=0
        edgeTau_err[0][uplim_ind]=edgeTau[uplim_ind]


        axes[0].errorbar(phase,edgeTau*1e2,edgeTau_err*1e2,phase_err,ecolor='k',drawstyle='steps-mid',alpha=0.6,lw=2,uplims=uplim_ind,zorder=10)


def plot_cl(ObsID,folder,ax_steppar,names=['dips','off_dip']):

    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/{folder}')
    clrs=['r','b']
    for clr,name in zip(clrs,names):
        steppar_data=np.genfromtxt(f'all_{name}_edge_steppar.dat').T
        tau_dips=steppar_data[:,0]*100
        dchi2_dips=steppar_data[:,1]-min(steppar_data[:,1])
        ax_steppar.semilogx(tau_dips,dchi2_dips,'-',label=name,color=clr)

    for cl,sigm in zip([0.682689,0.9,0.95,0.9973,0.999936],['(1 sigma)','','','(3 sigma)','(4 sigma)']):
        dc=stats.chi2.ppf(cl,df=1)
        ax_steppar.axhline(dc,color='k',ls=':',lw=1)
        cl=cl*100
        ax_steppar.text(0.3,dc,f'{"%.4f"%cl} % cl {sigm}')


#%% group 1,2
fig = plt.figure(figsize=(16,8))

plt.subplots_adjust(wspace=0.05)
plt.subplots_adjust(hspace=0.0)
dh=0.05
matplotlib.rcParams['figure.subplot.left']=dh
matplotlib.rcParams['figure.subplot.bottom']=dh
matplotlib.rcParams['figure.subplot.right']=1-dh
matplotlib.rcParams['figure.subplot.top']=1-dh


rows=6
cols=7

ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.05, 'A', transform=ax_eqw.transAxes,
            size=20, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
ax_edgeTau = plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=3)
load_and_plot_data_for_paper('group_1', ax_fe_norm,ax_eqw, ax_edgeTau)




add_average_analysis_on_plots_for_paper('group_1','average_phase_all_no_borders',
                              [[13,14,15],[26,27,28,29,30,31,32],[1,2,3,4,5,6,7,8,9,10,11],[17,18,19,20,21,22,23,24]],
                              ['dips','dips','off_dip','off_dip'],
                              ['r','r','b','b'], [ax_edgeTau,ax_fe_norm,ax_eqw],n_ph=32)




ax_eqw = plt.subplot2grid((rows,cols), (0, 4), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.05, 'B', transform=ax_eqw.transAxes,
            size=20, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 4), rowspan=2, colspan=3)
ax_edgeTau = plt.subplot2grid((rows,cols), (4, 4), rowspan=2, colspan=3)
load_and_plot_data_for_paper('group_2', ax_fe_norm,ax_eqw, ax_edgeTau)



add_average_analysis_on_plots_for_paper('group_2','average_phase_all_one_dip_all_off_dip',
                              [np.arange(23,30),np.arange(1,21)],
                              ['dips','off_dip'],
                              ['r','b','b'], [ax_edgeTau,ax_fe_norm,ax_eqw],n_ph=32)

plt.show()






groupid='group_12'
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/paper/'+f'{groupid}_cutoffpl_edge.png',dpi=100)





#%% group 2,3
fig = plt.figure(figsize=(16,8))

plt.subplots_adjust(wspace=0.05)
plt.subplots_adjust(hspace=0.0)
dh=0.05
matplotlib.rcParams['figure.subplot.left']=dh
matplotlib.rcParams['figure.subplot.bottom']=dh
matplotlib.rcParams['figure.subplot.right']=1-dh
matplotlib.rcParams['figure.subplot.top']=1-dh


rows=6
cols=7

ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.05, 'A', transform=ax_eqw.transAxes,
            size=20, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
ax_edgeTau = plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=3)
load_and_plot_data_for_paper('group_3', ax_fe_norm,ax_eqw, ax_edgeTau)




add_average_analysis_on_plots_for_paper('group_3','average_phase_all_one_dip_all_off_dip',
                              [np.arange(18,26),np.arange(1,17)],
                              ['dips','off_dip'],
                              ['r','b','b'], [ax_edgeTau,ax_fe_norm,ax_eqw],n_ph=32)




ax_eqw = plt.subplot2grid((rows,cols), (0, 4), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.05, 'B', transform=ax_eqw.transAxes,
            size=20, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 4), rowspan=2, colspan=3)
ax_edgeTau = plt.subplot2grid((rows,cols), (4, 4), rowspan=2, colspan=3)
load_and_plot_data_for_paper('group_4', ax_fe_norm,ax_eqw, ax_edgeTau)



add_average_analysis_on_plots_for_paper('group_4','average_phase_all_one_dip_all_off_dip',
                              [np.arange(1,10),np.arange(12,32)],
                              ['dips','off_dip'],
                              ['r','b','b'], [ax_edgeTau,ax_fe_norm,ax_eqw],n_ph=32)

plt.show()






groupid='group_23'
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/paper/'+f'{groupid}_cutoffpl_edge.png',dpi=100)





#%% rising phase

figx,ax_edgeTau=plt.subplots()

fig = plt.figure(figsize=(16,6))

plt.subplots_adjust(wspace=0.00)
plt.subplots_adjust(hspace=0.0)
dh=0.05
matplotlib.rcParams['figure.subplot.left']=dh
matplotlib.rcParams['figure.subplot.bottom']=0.1
matplotlib.rcParams['figure.subplot.right']=1-dh
matplotlib.rcParams['figure.subplot.top']=1-dh


rows=4
cols=7

ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.05, 'A', transform=ax_eqw.transAxes,
            size=15, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
load_and_plot_data_for_paper('90089-11-02-04', ax_fe_norm,ax_eqw, ax_edgeTau,model='cutoffpl')
ax_fe_norm.set_xlabel('Phase')



ax_eqw = plt.subplot2grid((rows,cols), (0, 4), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.05, 'B', transform=ax_eqw.transAxes,
            size=15, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 4), rowspan=2, colspan=3)
load_and_plot_data_for_paper('90089-11-02-06', ax_fe_norm,ax_eqw, ax_edgeTau,model='cutoffpl')
ax_fe_norm.set_xlabel('Phase')



plt.show()


groupid='rising'
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/paper/'+f'{groupid}_cutoffpl.pdf',dpi=300)




#%% rising phase, only flux

fig = plt.figure(figsize=(16,6/2))

plt.subplots_adjust(wspace=0.00)
plt.subplots_adjust(hspace=0.0)
dh=0.2
matplotlib.rcParams['figure.subplot.left']=0.05
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=1-0.02
matplotlib.rcParams['figure.subplot.top']=0.9


rows=2
cols=7

ax_efold = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_efold.text(-0.1, 1.05, 'A', transform=ax_efold.transAxes,
            size=15, weight='bold')
ax_efold.set_xlabel('Phase')
ax_efold.set_title('90089-11-02-04')
ax_efold.set_ylabel('Flux 3-12 keV, $10^{-8}$ cgs')
q=read_ph_res_results_data('90089-11-02-04','cutoffpl',0)
ax_efold.errorbar(q[1],q[6]/1e-8,q[7]/1e-8,color='k',drawstyle='steps-mid',ls='-.',alpha=1)




ax_efold = plt.subplot2grid((rows,cols), (0, 4), rowspan=2, colspan=3)
ax_efold.text(-0.1, 1.05, 'B', transform=ax_efold.transAxes,
            size=15, weight='bold')
ax_efold.set_xlabel('Phase')
ax_efold.set_title('90089-11-02-06')
ax_efold.set_ylabel('Flux 3-12 keV, $10^{-8}$ cgs')
q=read_ph_res_results_data('90089-11-02-06','cutoffpl',0)
ax_efold.errorbar(q[1],q[6]/1e-8,q[7]/1e-8,color='k',drawstyle='steps-mid',ls='-.',alpha=1)

plt.show()


groupid='rising'
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/paper/'+f'{groupid}_cutoffpl.pdf',dpi=300)






#%% top phase

#figx,ax_edgeTau_dumm=plt.subplots()

fig = plt.figure(figsize=(16/2,6))

plt.subplots_adjust(wspace=0.00)
plt.subplots_adjust(hspace=0.0)
dh=0.1
matplotlib.rcParams['figure.subplot.left']=dh
matplotlib.rcParams['figure.subplot.bottom']=0.1
matplotlib.rcParams['figure.subplot.right']=1-dh
matplotlib.rcParams['figure.subplot.top']=1-dh


rows=6
cols=3

ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.05, 'A', transform=ax_eqw.transAxes,
            size=15, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
ax_edgeTau = plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=3)

load_and_plot_data_for_paper('group_1', ax_fe_norm,ax_eqw, ax_edgeTau,model='cutoffpl_edge',phases=[[7,8,15,16],[2,3,4,5,10,11,12,13]])
#ax_fe_norm.set_xlabel('Phase')



# ax_eqw = plt.subplot2grid((rows,cols), (0, 4), rowspan=2, colspan=3)
# ax_eqw.text(-0.1, 1.05, 'B', transform=ax_eqw.transAxes,
#             size=15, weight='bold')
# ax_fe_norm = plt.subplot2grid((rows,cols), (2, 4), rowspan=2, colspan=3)
# load_and_plot_data_for_paper('90089-11-04-00G', ax_fe_norm,ax_eqw, ax_edgeTau_dumm,model='cutoffpl')
# ax_fe_norm.set_xlabel('Phase')



plt.show()


groupid='top'
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/paper/'+f'{groupid}_cutoffpl.pdf',dpi=300)




#%% decline phase

fig = plt.figure(figsize=(16,6))

plt.subplots_adjust(wspace=0.00)
plt.subplots_adjust(hspace=0.0)
dh=0.05
matplotlib.rcParams['figure.subplot.left']=dh
matplotlib.rcParams['figure.subplot.bottom']=0.1
matplotlib.rcParams['figure.subplot.right']=1-dh
matplotlib.rcParams['figure.subplot.top']=1-dh

rows=6
cols=7

ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.05, 'A', transform=ax_eqw.transAxes,
            size=15, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
ax_edgeTau = plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=3)

load_and_plot_data_for_paper('group_2', ax_fe_norm,ax_eqw, ax_edgeTau,model='cutoffpl_edge',
                             phases=[[13,14,15],[1,2,3,4,5,6,7,8,9,10,11]])
#ax_fe_norm.set_xlabel('Phase')
#q=read_ph_res_results_data('group_2','cutoffpl',0)
#ax_eqw.errorbar(q[1],q[16]*1e3,q[17]*1e3,color='darkgreen',drawstyle='steps-mid',ls='-.',alpha=0.3)



ax_eqw = plt.subplot2grid((rows,cols), (0, 4), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.05, 'B', transform=ax_eqw.transAxes,
            size=15, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 4), rowspan=2, colspan=3)
ax_edgeTau = plt.subplot2grid((rows,cols), (4, 4), rowspan=2, colspan=3)

load_and_plot_data_for_paper('group_4', ax_fe_norm,ax_eqw, ax_edgeTau,model='cutoffpl_edge',
                             phases=[[2,3,4],[6,7,8,9,10,11,12,13,14,15]])
ax_fe_norm.set_xlabel('Phase')
#q=read_ph_res_results_data('group_4','cutoffpl',0)
#ax_eqw.errorbar(q[1],q[16]*1e3,q[17]*1e3,color='darkgreen',drawstyle='steps-mid',ls='-.',alpha=0.6)



plt.show()


groupid='decline_1'
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/paper/'+f'{groupid}_cutoffpl.pdf',dpi=300)





fig = plt.figure(figsize=(16,6))

plt.subplots_adjust(wspace=0.00)
plt.subplots_adjust(hspace=0.0)
dh=0.05
matplotlib.rcParams['figure.subplot.left']=dh
matplotlib.rcParams['figure.subplot.bottom']=0.1
matplotlib.rcParams['figure.subplot.right']=1-dh
matplotlib.rcParams['figure.subplot.top']=1-dh


rows=6
cols=7

ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.05, 'A', transform=ax_eqw.transAxes,
            size=15, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
ax_edgeTau = plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=3)

load_and_plot_data_for_paper('group_5', ax_fe_norm,ax_eqw, ax_edgeTau,model='cutoffpl_edge',
                             phases=[[11,12,13,14],[1,2,3,4,5,6,7,8,9]])
#q=read_ph_res_results_data('group_5','cutoffpl',0)
#ax_eqw.errorbar(q[1],q[16]*1e3,q[17]*1e3,color='darkgreen',drawstyle='steps-mid',ls='-.',alpha=0.6)



ax_eqw = plt.subplot2grid((rows,cols), (0, 4), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.05, 'B', transform=ax_eqw.transAxes,
            size=15, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 4), rowspan=2, colspan=3)
ax_edgeTau = plt.subplot2grid((rows,cols), (4, 4), rowspan=2, colspan=3)

load_and_plot_data_for_paper('group_6', ax_fe_norm,ax_eqw, ax_edgeTau,model='cutoffpl_edge',
                             phases=[[5,6,7,8],[10,11,12,13,14,15,16]])
ax_fe_norm.set_xlabel('Phase')
#q=read_ph_res_results_data('group_6','cutoffpl',0)
#ax_eqw.errorbar(q[1],q[16]*1e3,q[17]*1e3,color='darkgreen',drawstyle='steps-mid',ls='-.',alpha=0.6)



plt.show()


groupid='decline_2'
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/paper/'+f'{groupid}_cutoffpl.pdf',dpi=300)




fig = plt.figure(figsize=(16,6/2))

plt.subplots_adjust(wspace=0.00)
plt.subplots_adjust(hspace=0.0)
dh=0.2
matplotlib.rcParams['figure.subplot.left']=0.05
matplotlib.rcParams['figure.subplot.bottom']=0.1
matplotlib.rcParams['figure.subplot.right']=1-0.05
matplotlib.rcParams['figure.subplot.top']=0.9


rows=2
cols=7

ax_efold = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_efold.text(-0.1, 1.05, 'A', transform=ax_efold.transAxes,
            size=15, weight='bold')
ax_efold.set_xlabel('Phase')
ax_efold.set_title('90014-01-05-00')
ax_efold.set_ylabel('Flux 3-12 keV, $10^{-8}$ cgs')
q=read_ph_res_results_data('90014-01-05-00','cutoffpl',0)
ax_efold.errorbar(q[1],q[6]/1e-8,q[7]/1e-8,color='k',drawstyle='steps-mid',ls='-.',alpha=1)




ax_efold = plt.subplot2grid((rows,cols), (0, 4), rowspan=2, colspan=3)
ax_efold.text(-0.1, 1.05, 'B', transform=ax_efold.transAxes,
            size=15, weight='bold')
ax_efold.set_xlabel('Phase')
ax_efold.set_title('90014-01-05-02')
ax_efold.set_ylabel('Flux 3-12 keV, $10^{-8}$ cgs')
q=read_ph_res_results_data('90014-01-05-02','cutoffpl',0)
ax_efold.errorbar(q[1],q[6]/1e-8,q[7]/1e-8,color='k',drawstyle='steps-mid',ls='-.',alpha=1)

plt.show()


groupid='decline_3'
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/paper/'+f'{groupid}_cutoffpl.pdf',dpi=300)





#%% confidence intervals


fig,ax_steppar = plt.subplots(figsize=(4,2))

groupid='group_2'
plot_cl(groupid,'average_phase_all_one_dip_all_off_dip',ax_steppar)

ax_steppar.set_ylim(0,20)
ax_steppar.set_xlim(0.1,10)
ax_steppar.grid()
#ax_steppar.legend()
ax_steppar.set_ylabel('$\Delta\chi^2 $')
ax_steppar.set_xlabel('$\\tau\\times100$')
ax_steppar.set_title('group ii')
#plt.title(f'{ObsID}, {folder}')
plt.tight_layout()
plt.show()
fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/paper/'+f'{groupid}conf_int.pdf',dpi=300)



#%% plots for github images


def plot_ph_res_results_github(ObsID,model):
    [model,phase,
            mjd,
            tot_flux,
            datamode,
            cts,
            flux312,flux312_err,
            po,po_err,
            ecut,ecut_err,
            po_norm,po_norm_err,
            eline,eline_err,
            norm_line,norm_line_err,
            edgeE,edgeE_err,
            edgeTau,edgeTau_err,
            chi2,
            sigma,
            eqw,eqw_err,nh,sys_err]=read_ph_res_results_data(ObsID,model,roll)
    sigma=np.mean(sigma)
    nh=np.mean(nh)
    sys_err=np.mean(sys_err)

    fig = plt.figure(figsize=(18,10))
    plt.subplots_adjust(hspace=0)
    plt.subplots_adjust(wspace=1)
    rows=9
    cols=14
    #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
    ax_norm = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=7)
    ax_eqw = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=7)
    ax_edgeTau = plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=7)
    ax_po = plt.subplot2grid((rows,cols), (6, 0), rowspan=2, colspan=7)
    ax_po_norm = ax_po.twinx()
    ax_chi2 = plt.subplot2grid((rows,cols), (8, 0), rowspan=1, colspan=7)


    for ax in [ax_eqw,ax_edgeTau,ax_norm]:
        ax_efold=ax.twinx()
        ax_efold.errorbar(phase,flux312/1e-8,flux312_err/1e-8,color='k',drawstyle='steps-mid',ls=':',alpha=0.6)
        ax_efold.set_ylabel('Flux 3-12 keV cgs')

    for i in range(int(len(phase)/2)):
        ax_efold.text(phase[i],min(flux312/1e-8),str(i+1),fontsize=8)


    ax_norm.errorbar(phase,norm_line*1e3,norm_line_err*1e3,color='green',drawstyle='steps-mid',alpha=0.6)
    #mo_norm,_,chi2_red,_=fit_const_chi_square(norm_line*1000,1000*norm_line_err.max(axis=0))
    #ax_norm.axhline(mo_norm,alpha=0.3,ls=':')
    ax_norm.set_ylabel(f'Iron line norm, \n ($10^{-3}$ $ph cm^{-2}s^{-1}$)',fontsize=12,color='g')

    ax_norm.set_title(ObsID+f" (MJD {'%.1f'%mjd})",fontsize=10)


    ax_eqw.errorbar(phase,eqw*1e3,eqw_err*1e3,color='blue',drawstyle='steps-mid',alpha=0.6)
    ax_eqw.set_ylabel(f'Line Eq. width, eV ',fontsize=12,color='b')


    uplim_ind=edgeTau_err[0]==edgeTau
    edgeTau[uplim_ind]=edgeTau_err[1][uplim_ind]+edgeTau[uplim_ind]
    edgeTau_err[1][uplim_ind]=0
    edgeTau_err[0][uplim_ind]=edgeTau[uplim_ind]

    ax_edgeTau.errorbar(phase,edgeTau*1e2,edgeTau_err*1e2,color='r',drawstyle='steps-mid',alpha=0.6,uplims=uplim_ind)
    ax_edgeTau.set_ylabel(f'Iron K-edge $\\tau$*100',fontsize=9,color='r')
    ax_edgeTau.set_ylim(-0.2,8)

    # ax_edgeE=ax_edgeTau.twinx()
    # ax_edgeE.errorbar(phase,edgeE,edgeE_err*0,color='c',drawstyle='steps-mid',alpha=0.6)
    # ax_edgeE.set_ylabel(f'Iron K-edge energy \n (=const*Iron line energy)',fontsize=12,color='c')

    ax_po.errorbar(phase,po,po_err,color='b',drawstyle='steps-mid',alpha=0.6)
    ax_po.set_ylabel('Phot. index',color='b')


    ax_po_norm.errorbar(phase,po_norm,po_norm_err,color='m',drawstyle='steps-mid',alpha=0.6)
    ax_po_norm.set_ylabel('cutoffpl norm.', color='m')

    ax_ecut=ax_po.twinx()
    ax_ecut.errorbar(phase,ecut,ecut_err,color='g',drawstyle='steps-mid',alpha=0.6)
    ax_ecut.set_ylabel('$E_{cut}$, keV',color='g')
    ax_ecut.spines["right"].set_position(("axes", 1.1))


    ax_chi2.errorbar(phase,chi2,0,color='k',drawstyle='steps-mid',alpha=0.6)
    ax_chi2.set_ylabel('$\chi^2_{red}$')
    ax_chi2.set_xlabel('Phase',fontsize=12)




    ax_flux= plt.subplot2grid((rows,cols), (0, 8), rowspan=3, colspan=6)
    time=ObsParams.MJD_START
    flux=ObsParams.cutoffpl_cutoffpl_flux/1e-8
    ax_flux.semilogy(time-53000,flux,color='b',marker='s',lw=0,ms=4,alpha=0.8)
    ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b',fontsize=8)
    ax_flux.set_xlabel('Time, MJD-53000')
    ax_flux.yaxis.set_label_position("right")
    ax_flux.yaxis.tick_right()
    ax_flux.axvline(mjd-53000,color='r',ls='-.')


    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/paper/github/'+f'MJD_{"%.3f"%mjd}_{ObsID}.png',dpi=100)

    #fig.tight_layout()
    plt.show()
    return None


# ObsList_SA=ObsList_SA+['group_1','group_2','group_3','group_4','group_5','group_6']
# for ObsID in ObsList_SA:
#     plot_ph_res_results_github(ObsID,'cutoffpl_edge')
#     plt.close('all')

def plot_ph_res_results_github_pulse_profiles(ObsID,model='cutoffpl'):
    [model,phase,
            mjd,
            tot_flux,
            datamode,
            cts,
            flux312,flux312_err,
            po,po_err,
            ecut,ecut_err,
            po_norm,po_norm_err,
            eline,eline_err,
            norm_line,norm_line_err,
            edgeE,edgeE_err,
            edgeTau,edgeTau_err,
            chi2,
            sigma,
            eqw,eqw_err,nh,sys_err]=read_ph_res_results_data(ObsID,model,roll)
    sigma=np.mean(sigma)
    nh=np.mean(nh)
    sys_err=np.mean(sys_err)

    fig,[ax_efold,ax_flux] = plt.subplots(2,figsize=(10,5))

    plt.subplots_adjust(wspace=0.00)
    plt.subplots_adjust(hspace=0.3)
    matplotlib.rcParams['figure.subplot.left']=0.1
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=1-0.02
    matplotlib.rcParams['figure.subplot.top']=0.9

    ax_efold.errorbar(phase,flux312/1e-8,flux312_err/1e-8,color='k',drawstyle='steps-mid',ls=':',alpha=0.6)
    ax_efold.set_ylabel('Flux 3-12 keV cgs')
    ax_efold.set_xlabel('Phase')

    for i in range(int(len(phase)/2)):
        ax_efold.text(phase[i],min(flux312/1e-8),str(i+1),fontsize=8)


    time=ObsParams.MJD_START
    flux=ObsParams.cutoffpl_cutoffpl_flux/1e-8
    ax_flux.semilogy(time-53000,flux,color='b',marker='s',lw=0,ms=4,alpha=0.8)
    ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b',fontsize=8)
    ax_flux.set_xlabel('Time, MJD-53000')
    ax_flux.axvline(mjd-53000,color='r',ls='-.')


    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pulse_profiles/paper/github/pulse_profiles/'+f'MJD_{"%.3f"%mjd}_{ObsID}.png',dpi=100)

    #fig.tight_layout()
    plt.show()
    return None


ObsList=ObsList_SA+['group_1','group_2','group_3','group_4','group_5','group_6']+ObsList_SE
for ObsID in ObsList:
    plot_ph_res_results_github_pulse_profiles(ObsID)
    plt.close('all')
