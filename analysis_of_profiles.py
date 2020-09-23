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
        datamode='B_16ms_46M_0_49_H/B_16ms_64M_0_249'


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

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6.6,6.6))

df=ObsParams[~ObsParams.ObsID.isin(ignored_obs)]
time=df.MJD_START

for delay_name in delay_list:
    deltat,deltat_err=df[delay_name],df[delay_name+'_err']

    ax.errorbar(time,deltat,deltat_err,fmt='.',marker='s',ms=4,alpha=0.8,label=delay_name)

ax.set_ylabel('Argument of loc. max of CCF, sec',color='k')
ax.set_xlabel('Time, MJD')
ax.yaxis.set_minor_formatter(FormatStrFormatter("%.2f"))
#ax.legend([])
ax.grid()

fig.tight_layout()
#sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'delay_all_mc_err.pdf',dpi=500)
plt.show()




#%% plot delays (mc error)- only +-1 sec:


delay_list=['first_max','first_neg_max']

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0},figsize=(6.6,6.6/2))

df=ObsParams[~ObsParams.ObsID.isin(ignored_obs)]
time=df.MJD_START

for delay_name in delay_list:
    deltat,deltat_err=df[delay_name],df[delay_name+'_err']

    ax.errorbar(time,deltat,deltat_err,fmt='.',marker='s',ms=4,alpha=0.8,label=delay_name)

ax.set_ylabel('Argument of loc. max of CCF, sec',color='k')
ax.set_xlabel('Time, MJD')
ax.yaxis.set_minor_formatter(FormatStrFormatter("%.2f"))

#ax.legend([])
ax.grid()

fig.tight_layout()
#sns.despine(fig,top=1,right=0)
plt.savefig(savepath+f'delay_all_mc_err_1_sec.pdf',dpi=500)
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



#%% results for K-edge analysis




def read_ph_res_results_data():





def load_and_plot_data_edge(ObsID,plot_evol=1):
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
        datamode='B_16ms_46M_0_49_H/B_16ms_64M_0_249'


    data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cutoffpl_edge/ph_res_cutoffpl_edge.dat')
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
        #fig.savefig(f'{pulse_profile_save_path}/Day{np.round(mjd,3)}_{ObsID}_evol.png',dpi=500)
        if ObsID in ObsList_RISE:
            folder='rising_part'
        elif ObsID in ObsList_TOP:
            folder='top_part'
        elif ObsID in ObsList_DECAY_1:
            folder='decay_1_part'
        elif ObsID in ObsList_DECAY_2:
            folder='decay_2_part'


        #fig.savefig(f'{pulse_profile_save_path}/{folder}/Day{np.round(mjd,3)}_{ObsID}_evol.png',dpi=500)


    plt.show()
#        plt.close(fig)
#        plt.close('all')


load_and_plot_data_edge('90089-11-03-01G')


load_and_plot_data_edge('90427-01-03-00')



#%% free line energy analysis



def load_and_plot_data_line_en(ObsID,plot_evol=1):
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
        datamode='B_16ms_46M_0_49_H/B_16ms_64M_0_249'


    data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cutoffpl_en_free/ph_res_cutoffpl_en_free.dat')
    #tmp=data[:,7]
    #tmp_arg=np.argmin(tmp)

    #data=np.roll(data,-tmp_arg,axis=0)



    N_sp=(data[0,1]-1)/2
    spe_num=data[:,0]

    data=np.vstack((data,data))
    nph=data[0,1]
    data[:,0]=np.arange(1,nph) #create new spe_num
    spe_num=data[:,0]
    phase=((spe_num-1)/(N_sp))

    flux312=data[:,4]
    flux312_low=flux312-data[:,5]
    flux312_hi=data[:,6]-flux312

    flux312=flux312/1e-8
    flux312_hi=flux312_hi/1e-8
    flux312_low=flux312_low/1e-8
    flux312_err=np.vstack((flux312_low, flux312_hi))


    norm_line=data[:,12]*1000
    norm_line_low=norm_line-data[:,13]*1000
    norm_line_hi=data[:,14]*1000-norm_line
    norm_line_err=np.vstack((norm_line_low, norm_line_hi))

    gamma=data[:,7]
    efold=data[:,8]

    eline=data[:,9]
    eline_lo=eline-data[:,10]
    eline_hi=data[:,11]-eline
    eline_err=np.vstack((eline_lo, eline_hi))

    chi2_red=data[:,2]


    matplotlib.rcParams['figure.figsize'] = 8, 8/1.5
    matplotlib.rcParams['figure.subplot.left']=0.1
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=0.95
    matplotlib.rcParams['figure.subplot.top']=0.9
    plt.subplots_adjust(wspace=2)
    plt.subplots_adjust(hspace=1)

    fig = plt.figure()
    rows=6
    cols=5
    #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
    ax_norm = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
    ax_efold=ax_norm.twinx()

    #ax_fe_flux = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
    #ax_efold_2=ax_fe_flux.twinx()

    ax_eline = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
    ax_efold_3=ax_eline.twinx()
    ax_po= plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=3)


    ax_efold.errorbar(phase,flux312,flux312_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
    ax_norm.errorbar(phase,norm_line,norm_line_err,color='c',drawstyle='steps-mid',alpha=0.6)



    #ax_eqw.set_ylabel(f'Iron line \n Eq. width, eV \n $\chi^2_{{red}}$={chi2_red} \n -log(p)={logp}',color='c',fontsize=8)
    ax_norm.set_ylabel(f'Iron line \n norm, (a.u.)',color='c',fontsize=8)

    ax_efold.set_ylabel('Flux (3-12 keV) \n $10^{-8}$ cgs',fontsize=6)
    ax_norm.set_xlabel('Phase',fontsize=8)
    ax_norm.set_title(ObsID+f'({datamode});$F_{{3-12}}$={tot_flux} $10^{{-8}}cgs; exp. {expo}ks $')


    ax_efold_3.errorbar(phase,flux312,flux312_err,color='k',label='line energy',drawstyle='steps-mid',ls=':',alpha=0.6)
    ax_efold_3.set_ylabel('Flux (3-12 keV) \n $10^{-8}$ cgs',fontsize=6)
    ax_eline.errorbar(phase,eline,eline_err,color='r',drawstyle='steps-mid',alpha=0.6)
    ax_eline.set_ylabel('Line energy, keV \n cutoffpl',fontsize=6)
    ax_norm.set_xlabel('Phase',fontsize=8)
    ax_eline.axhline(6.4,color='k',alpha=0.6,zorder=-10)


    ax_po.plot(phase,gamma,color='b',alpha=0.6)
    ax_po.set_xlabel('Phase',fontsize=8)
    ax_po.set_ylabel('Gamma',fontsize=8,color='b')
    ax_chi2=ax_po.twinx()
    ax_chi2.plot(phase,chi2_red,color='r')
    ax_chi2.set_ylabel('$\chi^2_{red}$',fontsize=8,color='r')




    ax_flux= plt.subplot2grid((rows,cols), (0, 3), rowspan=6, colspan=2)
    time=ObsParams.MJD_START
    flux=ObsParams.cutoffpl_tot_flux/1e-8
    ax_flux.semilogy(time,flux,color='b',marker='s',lw=0,ms=4,alpha=0.8)
    ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b',fontsize=8)
    ax_flux.set_xlabel('Time, MJD')
    ax_flux.yaxis.set_label_position("right")
    ax_flux.yaxis.tick_right()
    ax_flux.axvline(mjd,color='r',ls='-.')




    #fig.tight_layout()
    plt.savefig(pulse_profile_save_path+'energy_variations/'+f'MJD_{mjd}_{ObsID}_cutoffpl_en_free.png',dpi=500)
    plt.close()
    #plt.show()




def load_and_plot_data_line_en_edge(ObsID,plot_evol=1):
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
        datamode='B_16ms_46M_0_49_H/B_16ms_64M_0_249'


    data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cutoffpl_en_free_edge/ph_res_cutoffpl_en_free_edge.dat')
    #tmp=data[:,7]
    #tmp_arg=np.argmin(tmp)

    #data=np.roll(data,-tmp_arg,axis=0)



    N_sp=(data[0,1]-1)/2
    spe_num=data[:,0]

    data=np.vstack((data,data))
    nph=data[0,1]
    data[:,0]=np.arange(1,nph) #create new spe_num
    spe_num=data[:,0]
    phase=((spe_num-1)/(N_sp))

    flux312=data[:,4]
    flux312_low=flux312-data[:,5]
    flux312_hi=data[:,6]-flux312

    flux312=flux312/1e-8
    flux312_hi=flux312_hi/1e-8
    flux312_low=flux312_low/1e-8
    flux312_err=np.vstack((flux312_low, flux312_hi))


    norm_line=data[:,12]*1000
    norm_line_low=norm_line-data[:,13]*1000
    norm_line_hi=data[:,14]*1000-norm_line
    norm_line_err=np.vstack((norm_line_low, norm_line_hi))

    gamma=data[:,7]
    efold=data[:,8]

    eline=data[:,9]
    eline_lo=eline-data[:,10]
    eline_hi=data[:,11]-eline
    eline_err=np.vstack((eline_lo, eline_hi))

    chi2_red=data[:,2]


    matplotlib.rcParams['figure.figsize'] = 8, 8/1.5
    matplotlib.rcParams['figure.subplot.left']=0.1
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=0.95
    matplotlib.rcParams['figure.subplot.top']=0.9
    plt.subplots_adjust(wspace=2)
    plt.subplots_adjust(hspace=1)

    fig = plt.figure()
    rows=6
    cols=5
    #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
    ax_norm = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
    ax_efold=ax_norm.twinx()

    #ax_fe_flux = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
    #ax_efold_2=ax_fe_flux.twinx()

    ax_eline = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
    ax_efold_3=ax_eline.twinx()
    ax_po= plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=3)


    ax_efold.errorbar(phase,flux312,flux312_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
    ax_norm.errorbar(phase,norm_line,norm_line_err,color='c',drawstyle='steps-mid',alpha=0.6)



    #ax_eqw.set_ylabel(f'Iron line \n Eq. width, eV \n $\chi^2_{{red}}$={chi2_red} \n -log(p)={logp}',color='c',fontsize=8)
    ax_norm.set_ylabel(f'Iron line \n norm, (a.u.)',color='c',fontsize=8)

    ax_efold.set_ylabel('Flux (3-12 keV) \n $10^{-8}$ cgs',fontsize=6)
    ax_norm.set_xlabel('Phase',fontsize=8)
    ax_norm.set_title(ObsID+f'({datamode});$F_{{3-12}}$={tot_flux} $10^{{-8}}cgs; exp. {expo}ks $')


    ax_efold_3.errorbar(phase,flux312,flux312_err,color='k',label='line energy',drawstyle='steps-mid',ls=':',alpha=0.6)
    ax_efold_3.set_ylabel('Flux (3-12 keV) \n $10^{-8}$ cgs',fontsize=6)
    ax_eline.errorbar(phase,eline,eline_err,color='r',drawstyle='steps-mid',alpha=0.6)
    ax_eline.set_ylabel('Line energy, keV \n cutoffpl*edge',fontsize=6)
    ax_norm.set_xlabel('Phase',fontsize=8)
    ax_eline.axhline(6.4,color='k',alpha=0.6,zorder=-10)


    ax_po.plot(phase,gamma,color='b',alpha=0.6)
    ax_po.set_xlabel('Phase',fontsize=8)
    ax_po.set_ylabel('Gamma',fontsize=8,color='b')
    ax_chi2=ax_po.twinx()
    ax_chi2.plot(phase,chi2_red,color='r')
    ax_chi2.set_ylabel('$\chi^2_{red}$',fontsize=8,color='r')


    ax_flux= plt.subplot2grid((rows,cols), (0, 3), rowspan=6, colspan=2)
    time=ObsParams.MJD_START
    flux=ObsParams.cutoffpl_tot_flux/1e-8
    ax_flux.semilogy(time,flux,color='b',marker='s',lw=0,ms=4,alpha=0.8)
    ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b',fontsize=8)
    ax_flux.set_xlabel('Time, MJD')
    ax_flux.yaxis.set_label_position("right")
    ax_flux.yaxis.tick_right()
    ax_flux.axvline(mjd,color='r',ls='-.')


    #fig.tight_layout()
    plt.savefig(pulse_profile_save_path+'energy_variations/'+f'MJD_{mjd}_{ObsID}_cutoffpl_en_free_edge.png',dpi=500)

    plt.close()
    #plt.show()


def plot_spe_ratio_en(ObsID,phase1,phase2):
    xte_obs=ObservationXTE(ObsID)

    data1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cutoffpl_en_free/spe_plots/ph_spe_{phase1}.dat')


    data2=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cutoffpl_en_free/spe_plots/ph_spe_{phase2}.dat')



    fig,[ax,ax2] = plt.subplots(2,figsize=(8, 3))


    label=f"Obs {ObsID}, phase {phase1} / phase {phase2}, \n cutoffpl model"

    def ratio_error(a,b,da,db):
        f=a/b
        sigma=np.abs(f)*np.sqrt( (da/a)**2 + (db/b)**2  )
        return f, sigma



    ax.errorbar(data1[0],data1[1]/data2[1],ratio_error(data1[1],data2[1],data1[2],data2[2])[1],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    ax.set_xscale('log')



    ax.legend(loc='upper left',fontsize=8)
    ax.axhline(1,color='k',ls=':',zorder=-10,alpha=0.6)
    ax.set_ylabel('eeufs1/eeufs2',fontsize=8)

    # ax2.errorbar(data1[0],data1[1],data1[2],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    # ax2.set_xscale('log')
    # ax2.set_yscale('log')
    # ax2.errorbar(data2[0],data2[1],data2[2],data2[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    # ax2.set_ylabel('eeufs1, eeufs 2')


    ax2.errorbar(data1[0],data1[1]-data2[1],np.sqrt(data1[2]**2+data2[2]**2),data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    ax2.set_xscale('log')
    ax2.set_ylabel('eeufs1 - eeufs 2')
    ax2.axhline(0,color='k',ls=':',zorder=-10,alpha=0.6)


    plt.show()



def plot_spe_ratio_en_edge(ObsID,phase1,phase2):
    xte_obs=ObservationXTE(ObsID)

    data1=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cutoffpl_en_free_edge/spe_plots/ph_spe_{phase1}.dat')


    data2=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cutoffpl_en_free_edge/spe_plots/ph_spe_{phase2}.dat')



    fig,[ax,ax2] = plt.subplots(2,figsize=(8, 6))


    label=f"Obs {ObsID}, phase {phase1} / phase {phase2}, \n cutoffpl*edge model"

    def ratio_error(a,b,da,db):
        f=a/b
        sigma=np.abs(f)*np.sqrt( (da/a)**2 + (db/b)**2  )
        return f, sigma



    ax.errorbar(data1[0],data1[1]/data2[1],ratio_error(data1[1],data2[1],data1[2],data2[2])[1],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    ax.set_xscale('log')



    ax.legend(loc='upper left',fontsize=8)
    ax.axhline(1,color='k',ls=':',zorder=-10,alpha=0.6)
    ax.set_ylabel('eeufs1/eeufs2',fontsize=8)

    # ax2.errorbar(data1[0],data1[1],data1[2],data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    # ax2.set_xscale('log')
    # ax2.set_yscale('log')
    # ax2.errorbar(data2[0],data2[1],data2[2],data2[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    # ax2.set_ylabel('eeufs1, eeufs 2')


    ax2.errorbar(data1[0],data1[1]-data2[1],np.sqrt(data1[2]**2+data2[2]**2),data1[3],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
    ax2.set_xscale('log')
    ax2.set_ylabel('eeufs1 - eeufs 2')
    ax2.axhline(0,color='k',ls=':',zorder=-10,alpha=0.6)


    plt.show()


ObsList_SA=['90089-11-03-02','90089-11-03-01G','90089-11-03-00G',
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
 '90014-01-04-03']



# ObsID='90427-01-03-00'
# ph1=1
# ph2=5


# load_and_plot_data_line_en(ObsID)
# #plot_spe_ratio_en(ObsID, ph1,ph2)


# load_and_plot_data_line_en_edge(ObsID)
# #plot_spe_ratio_en_edge(ObsID, ph1,ph2)


for ObsID in ObsList_SA:

    load_and_plot_data_line_en(ObsID)

    load_and_plot_data_line_en_edge(ObsID)


# #%% reflection models




# def load_and_plot_data_pexrav(ObsID,plot_evol=1):
#     xte_obs=ObservationXTE(ObsID)



#     ser=xte_obs.pandas_series()

#     period=ser['period_orb_corr']
#     mjd=ser['MJD_START']
#     expo=ser['EXPOSURE']/1000
#     expo=np.round(expo,1)
#     tot_flux=ser['cutoffpl_tot_flux']/1e-8
#     tot_flux=np.round(tot_flux,2)

#     datalist=ser['fasebin_cfg']
#     if datalist=='se':
#         datamode='E_125us_64M_0_1s'
#     elif datalist=='sa':
#         datamode='B_16ms_46M_0_49_H/B_16ms_64M_0_249'


#     data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/pexrav/ph_res_cutoffpl_en_free.dat')
#     #tmp=data[:,7]
#     #tmp_arg=np.argmin(tmp)

#     #data=np.roll(data,-tmp_arg,axis=0)



#     N_sp=(data[0,1]-1)/2
#     spe_num=data[:,0]

#     data=np.vstack((data,data))
#     nph=data[0,1]
#     data[:,0]=np.arange(1,nph) #create new spe_num
#     spe_num=data[:,0]
#     phase=((spe_num-1)/(N_sp))

#     flux312=data[:,4]
#     flux312_low=flux312-data[:,5]
#     flux312_hi=data[:,6]-flux312

#     flux312=flux312/1e-8
#     flux312_hi=flux312_hi/1e-8
#     flux312_low=flux312_low/1e-8
#     flux312_err=np.vstack((flux312_low, flux312_hi))


#     norm_line=data[:,12]*1000
#     norm_line_low=norm_line-data[:,13]*1000
#     norm_line_hi=data[:,14]*1000-norm_line
#     norm_line_err=np.vstack((norm_line_low, norm_line_hi))

#     gamma=data[:,7]
#     efold=data[:,8]

#     eline=data[:,9]
#     eline_lo=eline-data[:,10]
#     eline_hi=data[:,11]-eline
#     eline_err=np.vstack((eline_lo, eline_hi))

#     chi2_red=data[:,2]


#     matplotlib.rcParams['figure.figsize'] = 8, 8/1.5
#     matplotlib.rcParams['figure.subplot.left']=0.1
#     matplotlib.rcParams['figure.subplot.bottom']=0.15
#     matplotlib.rcParams['figure.subplot.right']=0.95
#     matplotlib.rcParams['figure.subplot.top']=0.9
#     plt.subplots_adjust(wspace=2)
#     plt.subplots_adjust(hspace=1)

#     fig = plt.figure()
#     rows=6
#     cols=3
#     #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
#     ax_norm = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=3)
#     ax_efold=ax_norm.twinx()

#     #ax_fe_flux = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
#     #ax_efold_2=ax_fe_flux.twinx()

#     ax_eline = plt.subplot2grid((rows,cols), (2, 0), rowspan=2, colspan=3)
#     ax_efold_3=ax_eline.twinx()
#     ax_po= plt.subplot2grid((rows,cols), (4, 0), rowspan=2, colspan=3)


#     ax_efold.errorbar(phase,flux312,flux312_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
#     ax_norm.errorbar(phase,norm_line,norm_line_err,color='c',drawstyle='steps-mid',alpha=0.6)



#     #ax_eqw.set_ylabel(f'Iron line \n Eq. width, eV \n $\chi^2_{{red}}$={chi2_red} \n -log(p)={logp}',color='c',fontsize=8)
#     ax_norm.set_ylabel(f'Iron line \n norm, (a.u.)',color='c',fontsize=8)

#     ax_efold.set_ylabel('Flux (3-12 keV) \n $10^{-8}$ cgs',fontsize=6)
#     ax_norm.set_xlabel('Phase',fontsize=8)
#     ax_norm.set_title(ObsID+f'({datamode});$F_{{3-12}}$={tot_flux} $10^{{-8}}cgs; exp. {expo}ks $')


#     ax_efold_3.errorbar(phase,flux312,flux312_err,color='k',label='line energy',drawstyle='steps-mid',ls=':',alpha=0.6)
#     ax_efold_3.set_ylabel('Flux (3-12 keV) \n $10^{-8}$ cgs',fontsize=6)
#     ax_eline.errorbar(phase,6.4*(1-eline/(300000)),6.4*eline_err/300000,color='r',drawstyle='steps-mid',alpha=0.6)
#     ax_eline.set_ylabel('6.4 keV *(1-[Shift, km/s]/c]',fontsize=6)




#     ax_eline.axhline(6.4,color='k',alpha=0.6,zorder=-10)
#     ax_norm.set_xlabel('Phase',fontsize=8)


#     ax_po.plot(phase,gamma,color='b',alpha=0.6)
#     ax_po.set_xlabel('Phase',fontsize=8)
#     ax_po.set_ylabel('Gamma',fontsize=8,color='b')
#     ax_chi2=ax_po.twinx()
#     ax_chi2.plot(phase,chi2_red,color='r')
#     ax_chi2.set_ylabel('$\chi^2_{red}$',fontsize=8,color='r')

#     fig.tight_layout()

#     plt.show()



# #%% load data
# load_and_plot_data_pexrav('90089-11-03-01G')
# load_and_plot_data_line_en('90089-11-03-01G')


# load_and_plot_data_pexrav('90427-01-03-00')
# load_and_plot_data_line_en('90427-01-03-00')


# load_and_plot_data_pexrav('90427-01-03-06')
# load_and_plot_data_line_en('90427-01-03-06')


# load_and_plot_data_pexrav('90014-01-03-020')
# load_and_plot_data_line_en('90014-01-03-020')


#%% fixed energy and width, without K-edge


def load_and_plot_data_line_en_fix(ObsID,shift_phase=1,plot_gamma=0):
    xte_obs=ObservationXTE(ObsID)

    ser=xte_obs.pandas_series()

    mjd=ser['MJD_START']
    expo=ser['EXPOSURE']/1000
    expo=np.round(expo,1)
    tot_flux=ser['cutoffpl_tot_flux']/1e-8
    tot_flux=np.round(tot_flux,2)

    datalist=ser['fasebin_cfg']
    if datalist=='se':
        datamode='E_125us_64M_0_1s'
    elif datalist=='sa':
        datamode='B_16ms_46M_0_49_H/B_16ms_64M_0_249'


    data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cutoffpl_en_fix/ph_res_cutoffpl_en_fix.dat')
    if shift_phase==1:
        tmp=data[:,7]
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

    flux312=data[:,4]
    flux312_low=flux312-data[:,5]
    flux312_hi=data[:,6]-flux312

    flux312=flux312/1e-8
    flux312_hi=flux312_hi/1e-8
    flux312_low=flux312_low/1e-8
    flux312_err=np.vstack((flux312_low, flux312_hi))


    norm_line=data[:,12]*1000
    norm_line_low=norm_line-data[:,13]*1000
    norm_line_hi=data[:,14]*1000-norm_line
    norm_line_err=np.vstack((norm_line_low, norm_line_hi))

    gamma=data[:,7]
    efold=data[:,8]

    eline=data[:,9]
    eline_lo=eline-data[:,10]
    eline_hi=data[:,11]-eline
    eline_err=np.vstack((eline_lo, eline_hi))

    chi2_red=data[:,2]


    matplotlib.rcParams['figure.figsize'] = 13, 4
    matplotlib.rcParams['figure.subplot.left']=0.07
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=0.99
    matplotlib.rcParams['figure.subplot.top']=0.9
    plt.subplots_adjust(wspace=0)
    plt.subplots_adjust(hspace=0)

    fig = plt.figure()
    rows=2
    cols=7
    #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
    ax_norm = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=4)
    ax_efold=ax_norm.twinx()


    ax_efold.errorbar(phase,flux312,flux312_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
    ax_norm.errorbar(phase,norm_line,norm_line_err,color='c',drawstyle='steps-mid',alpha=0.6)

    ax_norm.set_ylabel(f'Iron line \n norm, (1000 ph/cm^2/s)',fontsize=12)

    ax_efold.set_ylabel('Flux (3-12 keV) \n $10^{-8}$ cgs',fontsize=12)
    ax_norm.set_xlabel('Phase',fontsize=12)
    ax_norm.set_title(ObsID+f'({datamode});$F_{{3-12}}$={tot_flux} $10^{{-8}}cgs; exp. {expo}ks $',fontsize=10)

    if plot_gamma:
        ax_norm.cla()
        ax_norm.plot(phase,gamma)


    ax_flux= plt.subplot2grid((rows,cols), (0, 5), rowspan=2, colspan=2)
    time=ObsParams.MJD_START
    flux=ObsParams.cutoffpl_tot_flux/1e-8
    ax_flux.semilogy(time,flux,color='b',marker='s',lw=0,ms=4,alpha=0.8)
    ax_flux.set_ylabel('Flux (3-12 keV), \n $10^{-8}$ cgs',color='b',fontsize=8)
    ax_flux.set_xlabel('Time, MJD')
    #ax_flux.yaxis.set_label_position("right")
    #ax_flux.yaxis.tick_right()
    ax_flux.axvline(mjd,color='r',ls='-.')


    #fig.tight_layout()
    plt.savefig(pulse_profile_save_path+'iron_line_intensity_en_fix/'+f'MJD_{round(mjd,3)}_{ObsID}_cutoffpl_en_fix_no_edge.png',dpi=500)

    plt.close()
    #plt.show()


#%% run stuff
ObsList=ObsList_RISE+ObsList_TOP+ObsList_DECAY_1+ObsList_DECAY_2

for ObsID in ObsList:
    load_and_plot_data_line_en_fix(ObsID)





ObsList_one_peaked=[
 '90089-11-01-02',
 '90089-11-01-03',
 '90089-11-01-04']
for ObsID in ObsList_one_peaked:
    load_and_plot_data_line_en_fix(ObsID,shift_phase=0,plot_gamma=1)




#%%stack observations


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
        ax[0].set_ylabel('counts/ phase_bin')

        tmp_data=fasebin_roll[1].data

        Indeces=np.arange(0,len(tmp_data))
        fasebin_roll[1].data=tmp_data[np.roll(Indeces,roll)]
        fasebin_roll.close()
        fasebin_check=fits.open('fasebin_sys_roll.pha')
        cts=fasebin_check[1].data['counts'].sum(axis=1)
        fasebin_check.close()
        all_cts.append(cts)

        ax[1].plot(cts-np.mean(cts),'+-',label=f'{ObsID} - rolled')
        ax[1].set_ylabel('counts-mean')

    all_cts=np.array(all_cts)
    all_cts=np.sum(all_cts,axis=0)
    ax[2].plot(all_cts,'d-',label='counts sum (rolled)')

    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}')
    path=[f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/fasebin_sys_roll.pha' for ObsID in ObsList]
    path='\n'.join(path)
    print(path)
    os.system(f"echo '{path}' > alldays_pha.list")

    os.system("fbadd infile='@alldays_pha.list' outfile='fasebin_sys.pha'")

    os.system(f"cp /Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/response.rsp response.rsp")


    fasebin_joined=fits.open('fasebin_sys.pha')
    cts=fasebin_joined[1].data['counts'].sum(axis=1)
    ax[2].plot(cts,'o',label='fasebin_sys.pha - counts')


    ax[0].legend()
    ax[1].legend()
    ax[2].legend()


    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}/pulses.png')
    plt.show()

roll_fasebin_files(ObsList=['90089-11-01-02', '90089-11-01-03', '90089-11-01-04'],
                       Roll_list=[0,-3,-1])

#run something like this:  xspec - ~/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res_cutoffpl_en_fix_ign11.txt


#%% plot results
def load_and_plot_data_line_en_fix_stacked(folder_name,shift_phase=1,plot_gamma=0):

    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}')


    data=np.genfromtxt('./cutoffpl_en_fix/ph_res_cutoffpl_en_fix.dat')
    if shift_phase==1:
        tmp=data[:,7]
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

    flux312=data[:,4]
    flux312_low=flux312-data[:,5]
    flux312_hi=data[:,6]-flux312

    flux312=flux312/1e-8
    flux312_hi=flux312_hi/1e-8
    flux312_low=flux312_low/1e-8
    flux312_err=np.vstack((flux312_low, flux312_hi))


    norm_line=data[:,12]*1000
    norm_line_low=norm_line-data[:,13]*1000
    norm_line_hi=data[:,14]*1000-norm_line
    norm_line_err=np.vstack((norm_line_low, norm_line_hi))

    gamma=data[:,7]
    efold=data[:,8]

    eline=data[:,9]
    eline_lo=eline-data[:,10]
    eline_hi=data[:,11]-eline
    eline_err=np.vstack((eline_lo, eline_hi))

    chi2_red=data[:,2]


    matplotlib.rcParams['figure.figsize'] = 13, 4
    matplotlib.rcParams['figure.subplot.left']=0.07
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=0.99
    matplotlib.rcParams['figure.subplot.top']=0.9
    plt.subplots_adjust(wspace=0)
    plt.subplots_adjust(hspace=0)

    fig = plt.figure()
    rows=2
    cols=5
    #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
    ax_norm = plt.subplot2grid((rows,cols), (0, 0), rowspan=2, colspan=4)
    ax_efold=ax_norm.twinx()


    ax_efold.errorbar(phase,flux312,flux312_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
    ax_norm.errorbar(phase,norm_line,norm_line_err,color='c',drawstyle='steps-mid',alpha=0.6)

    ax_norm.set_ylabel('Iron line \n norm, (1000 ph/cm^2/s)',fontsize=12)

    ax_efold.set_ylabel('Flux (3-12 keV) \n $10^{-8}$ cgs',fontsize=12)
    ax_norm.set_xlabel('Phase',fontsize=12)
    ax_norm.set_title(f'{folder_name}')

    if plot_gamma:
        ax_norm.cla()
        ax_norm.plot(phase,gamma)


    #fig.tight_layout()
    plt.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/results/{folder_name}/cutoffpl_en_fix_no_edge.png',dpi=500)

    plt.close()
    #plt.show()


load_and_plot_data_line_en_fix_stacked(folder_name='fasebin_combine_rising_phase',
                                       plot_gamma=0)


