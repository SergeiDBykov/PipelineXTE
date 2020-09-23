#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 12:53:20 2020

@author: s.bykov
"""
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from glob import glob
import os
from Misc.TimeSeries.cross_corr import my_crosscorr



#xspec scripts

#/Users/s.bykov/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res_cutoffpl_eqw_ts.txt

#/Users/s.bykov/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res_cutoffpl_no_gauss_eqw_ts.txt



class eqw_timeseries:
    def __init__(self,ObsID,dt=16,no_gauss=0,roll=0):
        self.ObsID=ObsID
        self.dt=float(dt)
        self.roll=roll

        os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{self.ObsID}/products/fasebin_timeser')


        filename=f'cutoffpl_no_gauss/ph_res_cutoffpl_no_gauss.dat'
        data_no_g=np.genfromtxt(filename)
        self.chi2_r_no_gauss=data_no_g[:,2]
        self.rate=data_no_g[:,4]
        self.chi2_phase=data_no_g[:,0]*self.dt
        if no_gauss:
            fig,ax=plt.subplots()
            ax.plot(self.chi2_phase,self.chi2_r_no_gauss)
            ax.set_label(f'chi2 no gauss vs time \n mean chi2 {np.mean(self.chi2_r_no_gauss)}')
            print(f'chi2 no gauss vs time \n mean chi2 {np.mean(self.chi2_r_no_gauss)}')
            fig.savefig(f'chi2_{dt}s.png')
            fig,axs=plt.subplots(2)
            axs[0].plot(self.chi2_phase,self.rate)
            os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{self.ObsID}/products/fasebin_timeser')

            lc=fits.open('lc312.lc')
            rate=lc[1].data['rate']
            time=lc[1].data['time']
            time=time-time[0]
            axs[1].plot(time,rate)

            plt.show()



            self.rate_check=[time,rate,self.chi2_phase,self.rate]




    def read_and_plot_data(self,no_gauss=0):
        model='cutoffpl'
        os.chdir(model)
        filename=f'ph_res_{model}.dat'
        tmp=np.genfromtxt(filename)
        tmp=np.roll(tmp,self.roll,axis=0)
        tmp[:,0]=np.arange(1,tmp.shape[0]+1)
        self.data=tmp
        #self.data=np.genfromtxt(filename)


        data=self.data
        N_sp=(data[0,1]-1)/2
        spe_num=data[:,0]

        if len(spe_num)!=N_sp:
            raise Exception('Not all spectra were fitted')

        spe_num=data[:,0]
        self.phase=spe_num*self.dt
        phase=self.phase
        chi2_red=data[:,2]


        eqw=data[:,4]
        eqw_low=eqw-data[:,5]
        eqw_hi=data[:,6]-eqw

        eqw=eqw*1e3
        eqw_low=eqw_low*1e3
        eqw_hi=eqw_hi*1e3
        eqw_err=np.vstack((eqw_low, eqw_hi))
        self.eqw=eqw

        po=data[:,10]
        efold=data[:,11]
        ecut=data[:,12]
        eline=data[:,13]
        norm_line=data[:,14]
        self.norm_line=norm_line
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
        self.flux712=flux712

        fig,[ax_eqw,ax_chi,ax_ecutpars,ax_po_and_norm]=plt.subplots(4,figsize=(15,6),sharex=True)
        ax_efold=ax_eqw.twinx()
        ax_po_and_norm_twin=ax_po_and_norm.twinx()




        ax_efold.errorbar(phase,flux712,flux712_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
        ax_eqw.errorbar(phase,eqw,eqw_err,color='r',drawstyle='steps-mid',alpha=0.8)
        ax_eqw.tick_params(axis='y', colors='red')
        ax_eqw.spines['left'].set_color('red')

        ax_eqw.set_ylabel('Fe Ka Eq. width, eV',color='r')
        ax_efold.set_ylabel('Flux 7-12, 1e-8 cgs')
        ax_eqw.set_xlabel('Time, s')



        ax_chi.step(phase,chi2_red,'b-',where='mid')
        ax_chi.axhline(1,color='b')
        ax_chi.set_ylabel('xi2',color='b')
        ax_eqw.set_xlabel('Spectrum #')
        ax_chi2_no_gauss=ax_chi#.twinx()
        ax_chi2_no_gauss.step(self.chi2_phase,self.chi2_r_no_gauss,color='r',where='mid')
        ax_chi.axhline(2,color='r')
        ax_chi.grid(1,'both')


        ax_ecutpars.plot(phase,efold,color='g',label='efold (left)')
        ax_ecutpars.legend(loc='upper right')
        ax_ecutpars.plot(phase,eline,color='y',label='eline (left)')
        ax_ecutpars.legend(loc='upper left')
        ax_po_and_norm_twin=ax_po_and_norm.twinx()
        ax_po_and_norm_twin.plot(phase,po,color='k',label='po gamma (right)')
        ax_po_and_norm.errorbar(phase,norm_line,norm_line_err,color='m',label='iron line norm (left)')
        ax_po_and_norm_twin.legend(loc='upper right')
        ax_po_and_norm.legend(loc='upper left')

        plt.subplots_adjust(wspace=0.7)
        plt.subplots_adjust(hspace=0.3)

        plt.show()
        fig.savefig(f'report_{self.dt}s.png')


    def eqw_flux_ccf(self):

        fig,axs=plt.subplots(2,1)
        fig.subplots_adjust(hspace=0.5)
        my_crosscorr(self.phase,self.eqw,self.flux712,axs[0],axs[1],
                  subtract_mean=1,divide_by_mean=1,only_pos_delays=0,
                  y1label='eqw',y2label='F(7-12)',my_only_pos_delays=0,my_ccf=0)
        plt.show()
        fig.savefig(f'report_ccf_{self.dt}s.png')


    def norm_flux_ccf(self):

        fig,axs=plt.subplots(2,1)
        fig.subplots_adjust(hspace=0.5)
        my_crosscorr(self.phase,self.norm_line,self.flux712,axs[0],axs[1],
                  subtract_mean=1,divide_by_mean=1,only_pos_delays=0,
                  y1label='norm line',y2label='F(7-12)',my_only_pos_delays=0,my_ccf=0)
        plt.show()
        fig.savefig(f'report_ccf_norm_{self.dt}s.png')


    def lc_and_flux(self):
        fig,axs=plt.subplots(2)

        axs[0].plot(self.chi2_phase,self.rate)
        os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{self.ObsID}/products/fasebin_timeser_{int(self.dt)}')

        lc=fits.open('lc312.lc')
        rate=lc[1].data['rate']
        time=lc[1].data['time']
        time=time-time[0]
        axs[1].plot(time,rate)

        plt.show()


        # fig,axs=plt.subplots(2,1)
        # fig.subplots_adjust(hspace=0.5)
        # delay,lag,ccf=my_crosscorr(self.phase,rate,self.flux712,axs[0],axs[1],
        #           subtract_mean=1,divide_by_mean=1,only_pos_delays=0,
        #           y1label='rate 7-12',y2label='F(7-12)',my_only_pos_delays=0,my_ccf=0)
        # plt.show()
        # delay=delay[0]/self.dt


        # fig,axs=plt.subplots()
        # axs.plot(self.phase,self.flux712)


        # new_rate=np.roll(rate,-int(delay))

        # axs.twinx().plot(time,new_rate,'m:')

        # plt.show()
