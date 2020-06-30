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


matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.85
matplotlib.rcParams['figure.subplot.top']=0.95


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


#%% pulse profiles analysis


def load_and_plot_data(ObsID,plot_evol=1):
    chi2_crit=1.2
    xte_obs=ObservationXTE(ObsID)



    ser=xte_obs.pandas_series()

    period=ser['period_orb_corr']
    mjd=ser['MJD_START']
    expo=ser['EXPOSURE']/1000
    expo=np.round(expo,1)
    tot_flux=ser['cutoffpl_tot_flux']/1e-8
    tot_flux=np.round(tot_flux,2)

    datalist=ser['fasebin_cfg']
    if datalist=='se':
        datamode='E_125us_64M_0_1s'
    elif datalist=='sa':
        datamode='B_16ms_46M_0_49_H'


    data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cutoffpl/ph_res_cutoffpl.dat')
    tmp=data[:,7]
    tmp_arg=np.argmin(tmp)

    data=np.roll(data,-tmp_arg,axis=0)



    N_sp=(data[0,1]-1)/2
    spe_num=data[:,0]

    data=np.vstack((data,data))
    nph=data[0,1]
    data[:,0]=np.arange(1,nph) #create new spe_num
    spe_num=data[:,0]
    phase=((spe_num-1)/(N_sp))

    eqw=data[:,4]
    eqw_low=eqw-data[:,5]
    eqw_hi=data[:,6]-eqw

    eqw=eqw*1e3
    eqw_low=eqw_low*1e3
    eqw_hi=eqw_hi*1e3
    eqw_err=np.vstack((eqw_low, eqw_hi)).max(axis=0)

    flux712=data[:,7]
    flux712_low=flux712-data[:,8]
    flux712_hi=data[:,9]-flux712

    flux712=flux712/1e-8
    flux712_hi=flux712_hi/1e-8
    flux712_low=flux712_low/1e-8

    norm_line=data[:,14]*1000
    norm_line_low=norm_line-data[:,15]*1000
    norm_line_hi=data[:,16]*1000-norm_line
    norm_line_err=np.vstack((norm_line_low, norm_line_hi)).max(axis=0)


    flux712_err=np.vstack((flux712_low, flux712_hi)).max(axis=0)


    if plot_evol:

        matplotlib.rcParams['figure.figsize'] = 12, 12/3
        matplotlib.rcParams['figure.subplot.left']=0.1
        matplotlib.rcParams['figure.subplot.bottom']=0.15
        matplotlib.rcParams['figure.subplot.right']=0.95
        matplotlib.rcParams['figure.subplot.top']=0.9
        plt.subplots_adjust(wspace=2)
        plt.subplots_adjust(hspace=1)

        fig = plt.figure()
        rows=7
        cols=7
        #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
        ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
        ax_efold=ax_eqw.twinx()

        #ax_fe_flux = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
        #ax_efold_2=ax_fe_flux.twinx()

        ax_fe_norm = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
        ax_efold_3=ax_fe_norm.twinx()
        ax_ccf= plt.subplot2grid((rows,cols), (5, 0), rowspan=2, colspan=3)






        ax_efold.errorbar(phase,flux712,flux712_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
        ax_eqw.errorbar(phase,eqw,eqw_err,color='c',drawstyle='steps-mid',alpha=0.6)

        #ax_eqw.set_ylabel('Iron line \n Eq. width, eV',color='c',fontsize=8)

        _,_,chi2_red,pval=fit_const_chi_square(eqw,eqw_err)
        no_puls=0
        if chi2_red<chi2_crit:
            no_puls=1

        chi2_red='%.2f'%chi2_red
        logp=-np.log10(pval)
        logp='%.1f'%logp

        _,_,chi2_red_norm,_=fit_const_chi_square(norm_line,norm_line_err)
        chi2_red_norm='%.2f'%chi2_red_norm

        #ax_eqw.set_ylabel(f'Iron line \n Eq. width, eV \n $\chi^2_{{red}}$={chi2_red} \n -log(p)={logp}',color='c',fontsize=8)
        ax_eqw.set_ylabel(f'Iron line \n Eq. width, eV \n $\chi^2_{{red}}$={chi2_red}',color='c',fontsize=8)
        ax_fe_norm.set_ylabel(f'Iron line \n Norm. (a. u.)  \n $\chi^2_{{red}}$={chi2_red_norm}',color='r',fontsize=8)



        ax_efold.set_ylabel('Flux (7-12 keV) \n $10^{-8}$ cgs',fontsize=6)
        ax_eqw.set_xlabel('Phase',fontsize=8)
        ax_eqw.set_title(ObsID+f'({datamode});$F_{{3-12}}$={tot_flux} $10^{{-8}}cgs; exp. {expo}ks $')


        ax_efold_3.errorbar(phase,flux712,flux712_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
        ax_fe_norm.errorbar(phase,norm_line,norm_line_err,color='r',drawstyle='steps-mid',alpha=0.6)
        ax_efold_3.set_ylabel('Flux (7-12 keV) \n $10^{-8}$ cgs',fontsize=6)
        ax_fe_norm.set_xlabel('Phase',fontsize=8)


        CCF=cross_correlation.CrossCorrelation(phase*period,eqw,flux712,circular=1)
        lag,ccf=CCF.calc_ccf()

        ax_ccf.plot(lag,ccf,color='b',alpha=0.6)
        ax_ccf.set_xlim(-period,+period)
        ax_ccf.set_xlabel('Delay, sec',fontsize=8)
        ax_ccf.set_ylabel('Pearson r',fontsize=8)
        if no_puls:
            ax_ccf.text(0,0,f'$\chi^2_{{red}}<{chi2_crit}$',size=20,bbox=dict(boxstyle="square",
                   ec=(1., 0.5, 0.5),
                   fc=(1., 0.8, 0.8),alpha=0.6
                   ))

        if not no_puls:

            mc_lag,mc_ccfs=CCF.mc_errors_ccfs(y1_err=eqw_err,y2_err=flux712_err,N_trials=500,
                      divide_by_mean=1, subtract_mean=1)

            first_max_arr=[]
            second_max_arr=[]
            first_neg_max_arr=[]
            second_neg_max_arr=[]
            for i in range(len(mc_ccfs)):
                trial_lag,trial_ccf=mc_lag,mc_ccfs[i]

                cond=(trial_lag>=0) & (trial_lag<=2.5)
                ccf_tmp=trial_ccf[cond]
                lag_tmp=trial_lag[cond]
                first_max=lag_tmp[np.argmax(ccf_tmp)]
                first_max_arr.append(first_max)


                cond=(trial_lag>2.5) & (trial_lag<=4.4)
                ccf_tmp=trial_ccf[cond]
                lag_tmp=trial_lag[cond]
                second_max=lag_tmp[np.argmax(ccf_tmp)]
                second_max_arr.append(second_max)


                cond=(trial_lag>-2.5) & (trial_lag<0)
                ccf_tmp=trial_ccf[cond]
                lag_tmp=trial_lag[cond]
                first_neg_max=lag_tmp[np.argmax(ccf_tmp)]
                first_neg_max_arr.append(first_neg_max)


                cond=(trial_lag>=-4.4) & (trial_lag<-2.5)
                ccf_tmp=trial_ccf[cond]
                lag_tmp=trial_lag[cond]
                second_neg_max=lag_tmp[np.argmax(ccf_tmp)]
                second_neg_max_arr.append(second_neg_max)



            first_max_arr=np.array(first_max_arr)
            second_max_arr=np.array(second_max_arr)
            first_neg_max_arr=np.array(first_neg_max_arr)
            second_neg_max_arr=np.array(second_neg_max_arr)

            ax_lag_distr=ax_ccf.twinx()
            ax_lag_distr.hist(first_max_arr,color='k',alpha=0.3,lw=0.8,histtype='step',bins=25)
            ax_lag_distr.hist(second_max_arr,color='k',alpha=0.3,lw=0.8,histtype='step',bins=25)
            ax_lag_distr.hist(first_neg_max_arr,color='k',alpha=0.3,lw=0.8,histtype='step',bins=25)
            ax_lag_distr.hist(second_neg_max_arr,color='k',alpha=0.3,lw=0.8,histtype='step',bins=25)
            ax_lag_distr.set_ylabel('peak distribution', color='k')

            #quant=np.quantile(firs_max_arr,q=[0.16,0.84])
            #ax_lag_distr.axvline(quant[0],color='r',alpha=0.5)
            #ax_lag_distr.axvline(quant[1],color='r',alpha=0.5)


            xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'first_max',first_max_arr.mean())
            xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'first_max_err',first_max_arr.std())

            xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'second_max',second_max_arr.mean())
            xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'second_max_err',second_max_arr.std())


            xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'first_neg_max',first_neg_max_arr.mean())
            xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'first_neg_max_err',first_neg_max_arr.std())

            xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'second_neg_max',second_neg_max_arr.mean())
            xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'second_neg_max_err',second_neg_max_arr.std())

            xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'ph_res_binsize',period/nph)


        #     cond=(lag>=0) & (lag<=2.5)
        #     ccf_tmp=ccf[cond]
        #     lag_tmp=lag[cond]
        #     first_max=lag_tmp[np.argmax(ccf_tmp)]

        #     ax_ccf.axvline(first_max,ls=':',color='b',alpha=0.5)
        #     xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'first_max',first_max)
        #     xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'first_max_err',period/nph)


        #     cond=(lag>2.5) & (lag<=4.4)
        #     ccf_tmp=ccf[cond]
        #     lag_tmp=lag[cond]
        #     second_max=lag_tmp[np.argmax(ccf_tmp)]

        #     ax_ccf.axvline(second_max,ls=':',color='b',alpha=0.5)
        #     xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'second_max',second_max)


        #     cond=(lag>-2.5) & (lag<0)
        #     ccf_tmp=ccf[cond]
        #     lag_tmp=lag[cond]
        #     first_neg_max=lag_tmp[np.argmax(ccf_tmp)]

        #     ax_ccf.axvline(first_neg_max,ls=':',color='b',alpha=0.5)
        #     xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'first_neg_max',first_neg_max)


        #     cond=(lag>-4.4) & (lag<-2.5)
        #     ccf_tmp=ccf[cond]
        #     lag_tmp=lag[cond]
        #     second_neg_max=lag_tmp[np.argmax(ccf_tmp)]

        #     ax_ccf.axvline(second_neg_max,ls=':',color='b',alpha=0.5)
        #     xte_obs.write_to_obs_info(xte_obs.fasebin_info_file,'second_neg_max',second_neg_max)



        ax_flux= plt.subplot2grid((rows,cols), (0, 4), rowspan=7, colspan=4)
        time=ObsParams.MJD_START
        flux=ObsParams.cutoffpl_tot_flux/1e-8

        ax_flux.plot(time,flux,color='b',marker='s',lw=0,ms=4,alpha=0.8)
        ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b',fontsize=8)
        ax_flux.set_xlabel('Time, MJD')
        ax_flux.axvline(mjd,color='r',ls='-.')


        fig.tight_layout()
        #sns.despine(fig,top=1,right=0)
        #sns.set(font_scale=0.5)
        fig.savefig(f'{pulse_profile_save_path}/Day{np.round(mjd,3)}_{ObsID}_evol.png',dpi=500)
        if ObsID in ObsList_RISE:
            folder='rising_part'
        elif ObsID in ObsList_TOP:
            folder='top_part'
        elif ObsID in ObsList_DECAY_1:
            folder='decay_1_part'
        elif ObsID in ObsList_DECAY_2:
            folder='decay_2_part'


        fig.savefig(f'{pulse_profile_save_path}/{folder}/Day{np.round(mjd,3)}_{ObsID}_evol.png',dpi=500)



        plt.close(fig)
        plt.close('all')





#%% test
STOP
load_and_plot_data('90089-11-04-01')





#%% run all obs

ObsList=ObsList_RISE+ObsList_TOP+ObsList_DECAY_1+ObsList_DECAY_2



err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            load_and_plot_data(ObsID)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)





#%% read pandas series
err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    ObsParams=pd.DataFrame()
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)

            ser=xte_obs.pandas_series()
            ObsParams=ObsParams.append(ser)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

name='standard_pipeline'
pd.to_pickle(ObsParams,f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pandas_data/{name}.pkl')
ObsParams.to_csv(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pandas_data/{name}.csv',index=0)



#%% plot delays (easy error):

delay_list=['first_max','second_max','first_neg_max','second_neg_max']

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

df=ObsParams
time=df.MJD_START

for delay_name in delay_list:
    deltat,deltat_err=df[delay_name],df['first_max_err']

    ax.errorbar(time,deltat,deltat_err,fmt='.',marker='s',ms=4,alpha=0.8,label=delay_name)

ax.set_ylabel('Iron line equivalent width delay, sec',color='k')
ax.set_xlabel('Time, MJD')
ax.legend()
ax.grid()

fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'delay_all.pdf',dpi=500)
plt.show()




#%% plot delays (mc error):


delay_list=['first_max','second_max','first_neg_max','second_neg_max']

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

df=ObsParams[~ObsParams.ObsID.isin(ignored_obs)]
time=df.MJD_START

for delay_name in delay_list:
    deltat,deltat_err=df[delay_name],df[delay_name+'_err']

    ax.errorbar(time,deltat,deltat_err,fmt='.',marker='s',ms=4,alpha=0.8,label=delay_name)

ax.set_ylabel('Argument of loc. max of CCF, sec',color='k')
ax.set_xlabel('Time, MJD')
#ax.legend([])
ax.grid()

fig.tight_layout()
#sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'delay_all_mc_err.pdf',dpi=500)
plt.show()





#%% landscape plot



def load_and_plot_data_for_paper(ObsID,ax_eqw,ax_norm,ax_ccf):
    fontsizes=9

    chi2_crit=1.2
    xte_obs=ObservationXTE(ObsID)

    ser=xte_obs.pandas_series()

    period=ser['period_orb_corr']
    mjd=ser['MJD_START']
    expo=ser['EXPOSURE']/1000
    expo=np.round(expo,1)
    tot_flux=ser['cutoffpl_tot_flux']/1e-8
    tot_flux=np.round(tot_flux,2)


    data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cutoffpl/ph_res_cutoffpl.dat')
    tmp=data[:,7]
    tmp_arg=np.argmin(tmp)

    data=np.roll(data,-tmp_arg,axis=0)



    N_sp=(data[0,1]-1)/2
    spe_num=data[:,0]

    data=np.vstack((data,data))
    nph=data[0,1]
    data[:,0]=np.arange(1,nph) #create new spe_num
    spe_num=data[:,0]
    phase=((spe_num-1)/(N_sp))

    eqw=data[:,4]
    eqw_low=eqw-data[:,5]
    eqw_hi=data[:,6]-eqw

    eqw=eqw*1e3
    eqw_low=eqw_low*1e3
    eqw_hi=eqw_hi*1e3
    eqw_err=np.vstack((eqw_low, eqw_hi)).max(axis=0)

    flux712=data[:,7]
    flux712_low=flux712-data[:,8]
    flux712_hi=data[:,9]-flux712

    flux712=flux712/1e-8
    flux712_hi=flux712_hi/1e-8
    flux712_low=flux712_low/1e-8
    flux712_err=np.vstack((flux712_low, flux712_hi)).max(axis=0)

    norm_line=data[:,14]*1000
    norm_line_low=norm_line-data[:,15]*1000
    norm_line_hi=data[:,16]*1000-norm_line
    norm_line_err=np.vstack((norm_line_low, norm_line_hi)).max(axis=0)


    ax_efold=ax_eqw.twinx()
    ax_efold_3=ax_norm.twinx()


    ax_efold.errorbar(phase,flux712,flux712_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
    ax_eqw.errorbar(phase,eqw,eqw_err,color='c',drawstyle='steps-mid',alpha=0.6)

    _,_,chi2_red,pval=fit_const_chi_square(eqw,eqw_err)
    no_puls=0
    if chi2_red<chi2_crit:
        no_puls=1

    chi2_red='%.2f'%chi2_red
    logp=-np.log10(pval)
    logp='%.1f'%logp

    _,_,chi2_red_norm,_=fit_const_chi_square(norm_line,norm_line_err)
    chi2_red_norm='%.2f'%chi2_red_norm


    ax_eqw.set_ylabel(f'Iron line \n Eq. width, eV \n $\chi^2_{{red}}$={chi2_red}',color='c',fontsize=fontsizes)
    ax_fe_norm.set_ylabel(f'Iron line \n Norm. (a. u.)  \n $\chi^2_{{red}}$={chi2_red_norm}',color='r',fontsize=fontsizes)



    ax_efold.set_ylabel('Flux (7-12 keV) \n $10^{-8}$ cgs',fontsize=fontsizes)
    ax_eqw.set_xlabel('Phase',fontsize=fontsizes)
    ax_eqw.set_title('ObsID '+ObsID+f'; \t\t $F_{{3-12}}$={tot_flux} $10^{{-8}}cgs$',
                     fontsize=fontsizes,loc='left')


    ax_efold_3.errorbar(phase,flux712,flux712_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
    ax_fe_norm.errorbar(phase,norm_line,norm_line_err,color='r',drawstyle='steps-mid',alpha=0.6)
    ax_efold_3.set_ylabel('Flux (7-12 keV) \n $10^{-8}$ cgs',fontsize=fontsizes)
    ax_fe_norm.set_xlabel('Phase',fontsize=fontsizes)


    CCF=cross_correlation.CrossCorrelation(phase*period,eqw,flux712,circular=1)
    lag,ccf=CCF.calc_ccf()

    ax_ccf.plot(lag,ccf,color='b',alpha=0.6)
    ax_ccf.set_xlim(-period,+period)
    ax_ccf.set_xlabel('Delay, sec',fontsize=fontsizes)
    ax_ccf.set_ylabel('Pearson r',fontsize=fontsizes)




#%% make figure


fig = plt.figure(figsize=(14,10))


plt.subplots_adjust(wspace=0.3)
plt.subplots_adjust(hspace=0.01)
dh=0.1
matplotlib.rcParams['figure.subplot.left']=dh
matplotlib.rcParams['figure.subplot.bottom']=dh
matplotlib.rcParams['figure.subplot.right']=1-dh
matplotlib.rcParams['figure.subplot.top']=1-dh


rows=16
cols=7
#(rows,cols), (y,x) <- those are coordinates of an axis in subplots


rise_phase=['90089-11-01-03','90089-11-02-00','90089-11-02-06','90089-11-02-04']
top_phase=['90089-11-03-01G','90089-11-03-02','90089-11-04-00G','90089-11-04-04']
decay_phase_1=['90427-01-03-01','90427-01-03-06','90014-01-03-00','90014-01-03-00']
decay_phase_saw=['90014-01-03-01','90014-01-03-020','90014-01-03-03','90014-01-04-00']
decay_phase_2=['90014-01-04-01','90014-01-04-03','90014-01-05-04','90014-01-05-02']


PlotList=decay_phase_2
PlotName='decay_phase_2'


ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.1, 'A', transform=ax_eqw.transAxes,
            size=20, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)

ax_ccf = plt.subplot2grid((rows,cols), (5, 0), rowspan=2, colspan=3)

load_and_plot_data_for_paper(PlotList[0], ax_eqw, ax_fe_norm, ax_ccf)






ax_eqw = plt.subplot2grid((rows,cols), (9, 0), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.1, 'B', transform=ax_eqw.transAxes,
            size=20, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (11, 0), rowspan=2, colspan=3)
ax_ccf = plt.subplot2grid((rows,cols), (14, 0), rowspan=2, colspan=3)
load_and_plot_data_for_paper(PlotList[1], ax_eqw, ax_fe_norm, ax_ccf)





ax_eqw = plt.subplot2grid((rows,cols), (0, 4), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.1, 'C', transform=ax_eqw.transAxes,
            size=20, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 4), rowspan=2, colspan=3)
ax_ccf = plt.subplot2grid((rows,cols), (5, 4), rowspan=2, colspan=3)

load_and_plot_data_for_paper(PlotList[2], ax_eqw, ax_fe_norm, ax_ccf)




ax_eqw = plt.subplot2grid((rows,cols), (9, 4), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.1, 'D', transform=ax_eqw.transAxes,
            size=20, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (11, 4), rowspan=2, colspan=3)
ax_ccf = plt.subplot2grid((rows,cols), (14, 4), rowspan=2, colspan=3)

load_and_plot_data_for_paper(PlotList[3], ax_eqw, ax_fe_norm, ax_ccf)


plt.show()

fig.savefig(f'{pulse_profile_save_path}/landscape_figs/{PlotName}.pdf',dpi=500)







#%% make figure (TWO OBS INSTEAD OF 4)


fig = plt.figure(figsize=(14,5))


plt.subplots_adjust(wspace=0.7)
plt.subplots_adjust(hspace=0.01)
dh=0.1
matplotlib.rcParams['figure.subplot.left']=dh
matplotlib.rcParams['figure.subplot.bottom']=dh
matplotlib.rcParams['figure.subplot.right']=1-dh
matplotlib.rcParams['figure.subplot.top']=1-dh


rows=7
cols=7
#(rows,cols), (y,x) <- those are coordinates of an axis in subplots


rise_phase=['90089-11-01-03','90089-11-02-04']
top_phase=['90089-11-03-01G','90089-11-04-00G']
decay_phase_1=['90427-01-03-01','90014-01-03-00']
decay_phase_saw=['90014-01-03-01','90014-01-04-00']
decay_phase_2=['90014-01-04-01','90014-01-05-02']


PlotList=rise_phase
PlotName='rise_phase'


ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.1, 'A', transform=ax_eqw.transAxes,
            size=20, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)

ax_ccf = plt.subplot2grid((rows,cols), (5, 0), rowspan=2, colspan=3)

load_and_plot_data_for_paper(PlotList[0], ax_eqw, ax_fe_norm, ax_ccf)





ax_eqw = plt.subplot2grid((rows,cols), (0, 4), rowspan=2, colspan=3)
ax_eqw.text(-0.1, 1.1, 'B', transform=ax_eqw.transAxes,
            size=20, weight='bold')
ax_fe_norm = plt.subplot2grid((rows,cols), (2, 4), rowspan=2, colspan=3)
ax_ccf = plt.subplot2grid((rows,cols), (5, 4), rowspan=2, colspan=3)

load_and_plot_data_for_paper(PlotList[1], ax_eqw, ax_fe_norm, ax_ccf)




plt.show()

fig.savefig(f'{pulse_profile_save_path}/landscape_figs/two_obs/{PlotName}.pdf',dpi=500)

