#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 13:19:32 2020

@author: s.bykov
"""

from PipelineXTE.pipeline_core import *


#xspec - /Users/s.bykov/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res_cutoffpl.txt
#os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90089-11-03-01G/products/fasebin_test_diff_T0/')
os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90089-11-02-06/products/fasebin_test_diff_T0/')




def ph_res_results(folder='fasebin_t01'):
    os.chdir(folder)
    model='cutoffpl'
    os.chdir(model)
    filename=f'ph_res_{model}.dat'
    period=1/0.22854999208145163

    data=np.genfromtxt(filename)

    N_sp=(data[0,1]-1)/2
    spe_num=data[:,0]


    data=np.vstack((data,data))
    nph=data[0,1]
    data[:,0]=np.arange(1,nph) #create new spe_num
    spe_num=data[:,0]
    phase=((spe_num-1)/(N_sp))

    chi2_red=data[:,2]


    eqw=data[:,4]
    eqw_low=eqw-data[:,5]
    eqw_hi=data[:,6]-eqw

    eqw=eqw*1e3
    eqw_low=eqw_low*1e3
    eqw_hi=eqw_hi*1e3
    eqw_err=np.vstack((eqw_low, eqw_hi))


    po=data[:,10]
    efold=data[:,11]
    ecut=data[:,12]
    eline=data[:,13]
    norm_line=data[:,14]

    norm_line_low=norm_line-data[:,15]
    norm_line_hi=data[:,16]-norm_line
    norm_line_err=np.vstack((norm_line_low, norm_line_hi))



    flux712=data[:,7]
    flux712_low=flux712-data[:,8]
    flux712_hi=data[:,9]-flux712

    flux712=flux712/1e-8
    flux712_hi=flux712_hi/1e-8
    flux712_low=flux712_low/1e-8

    flux712_err=np.vstack((flux712_low, flux712_hi))


    flux_gauss=data[:,17]
    flux_gauss_low=flux_gauss-data[:,18]
    flux_gauss_hi=data[:,19]-flux_gauss

    flux_gauss=flux_gauss/1e-10
    flux_gauss_hi=flux_gauss_hi/1e-10
    flux_gauss_low=flux_gauss_low/1e-10

    flux_gauss_err=np.vstack((flux_gauss_low, flux_gauss_hi))



    matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
    matplotlib.rcParams['figure.subplot.left']=0.15
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=0.85
    matplotlib.rcParams['figure.subplot.top']=0.9
    plt.subplots_adjust(wspace=2)
    plt.subplots_adjust(hspace=1)

    fig = plt.figure()
    rows=7
    cols=3
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

    ax_eqw.set_ylabel('Iron line \n Eq. width, eV',color='c',fontsize=8)
    ax_efold.set_ylabel('Flux (7-12 keV) \n $10^{-8}$ cgs',fontsize=6)
    ax_eqw.set_xlabel('Phase',fontsize=8)
    #ax_eqw.set_title(self.ObsID+f' ({datamode})')



    #ax_efold_2.errorbar(phase,flux712,flux712_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
    #ax_fe_flux.errorbar(phase,flux_gauss,flux_gauss_err,color='g',drawstyle='steps-mid',alpha=0.6)

    #ax_fe_flux.set_ylabel('Fe Ka flux ',color='g')
    #ax_efold_2.set_ylabel('Flux 7-12, 1e-8 cgs')
    #ax_fe_flux.set_xlabel('Phase')


    ax_efold_3.errorbar(phase,flux712,flux712_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
    ax_fe_norm.errorbar(phase,norm_line,norm_line_err,color='r',drawstyle='steps-mid',alpha=0.6)
    ax_fe_norm.set_ylabel('Iron line \n Norm. ',color='r',fontsize=8)
    ax_efold_3.set_ylabel('Flux (7-12 keV) \n $10^{-8}$ cgs',fontsize=6)
    ax_fe_norm.set_xlabel('Phase',fontsize=8)


    CCF=cross_correlation.CrossCorrelation(phase*period,eqw,flux712,circular=True)
    #CCF=cross_correlation.CrossCorrelation(phase*period,norm_line,flux712,circular=True)
    lag,ccf=CCF.calc_ccf()
    peaks,_,_=CCF.find_max()
    delay=min(peaks[peaks>0])
    #ax_ccf.axvline(delay,ls=':',color='g',alpha=0.5)

    ax_ccf.plot(lag,ccf,color='b',alpha=0.6)
    #ax_ccf.set_title(f'Flux lags <--- 0 ---> Eqw lags',fontsize=8)
    ax_ccf.set_xlim(0,2*period)
    #ax_ccf.set_xlabel('Eqw Delay, sec')
    ax_ccf.set_xlabel('Eqw Delay, sec',fontsize=8)
    ax_ccf.set_ylabel('Pearson r',fontsize=8)

    fig.tight_layout()
    sns.despine(fig,top=1,right=0)
    #sns.set(font_scale=0.5)
    fig.savefig(f'ph_res_report.pdf',dpi=500)

ph_res_results('fasebin_t04')
