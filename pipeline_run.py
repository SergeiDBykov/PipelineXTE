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
import shutil
from Misc.TimeSeries import cross_correlation
from Misc.plot_spectra import Spectra
from Misc.doppler_correction import correct_times
from subprocess import call
from scipy.optimize import  curve_fit
from PipelineXTE.pipeline_core import ObservationXTE
from pipeline_core import ObsList_SA,ObsList_SE,ObsList_all

plt.ioff()

RXTE_path='/Users/s.bykov/work/xray_pulsars/rxte/'



#%% make configs and files

do_not_start_this_chunk

err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.filter_and_files()

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)

for e,m in zip(err,msg):
    print(e,m)

errors_init=err
msg_init=msg



#%% make spe_files

do_not_start_this_chunk


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)

            xte_obs.make_spectrum_and_rsp_std2()

            xte_obs.calc_bkg_files()
            xte_obs.calc_bkg_spe()
            xte_obs.calc_deadtime()

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

errors_spe=err
msg_spe=msg




#%% fit spectra with model cutoffpl 3-12 keV

do_not_start_this_chunk

err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.fit_std2_spe(model='cutoffpl',error=0)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

errors_spe_fit=err
msg_spe_fit=msg




#%% fit spectra with model cutoffpl and phabs 3-12 keV


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.fit_std2_spe(model='phabs_cutoffpl',error=0)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

errors_spe_fit=err
msg_spe_fit=msg



#%% fit spectra with model cutoffpl 3-12 keV with  edge

do_not_start_this_chunk

err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.fit_std2_spe(model='edge_cutoffpl',error=0)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

errors_spe_fit=err
msg_spe_fit=msg




#%% fit spectra with model cutoffpl 3-12 keV with  edge and phabs

do_not_start_this_chunk

err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.fit_std2_spe(model='absorb_edge_cutoffpl',error=0)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

errors_spe_fit=err
msg_spe_fit=msg





#%% fit spectra with model cutoffpl 3-12 keV with  ionized edge

do_not_start_this_chunk

err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.fit_std2_spe(model='ionize_edge_cutoffpl',error=0)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

errors_spe_fit=err
msg_spe_fit=msg





#%% fit spectra with model cutoffpl 3-12 keV without gauss with phabs

do_not_start_this_chunk

err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.fit_std2_spe(model='phabs_cutoffpl_no_gauss',error=0)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

errors_spe_fit=err
msg_spe_fit=msg





#%% fit spectra with model cutoffpl 3-12 keV without gauss without phabs

do_not_start_this_chunk

err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.fit_std2_spe(model='cutoffpl_no_gauss',error=0)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

errors_spe_fit=err
msg_spe_fit=msg





#%% fit spectra with model cutoffpl and gabs


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)

            xte_obs.fit_std2_spe(model='gabs')

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

errors_spe_fit=err
msg_spe_fit=msg




#%% fit spectra with model cutoffpl and cyclabs


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)

            xte_obs.fit_std2_spe(model='cyclabs',error=0)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

errors_spe_fit=err
msg_spe_fit=msg




#%% pd dataframe

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

name='standard_pipeline' #standard_pipeline_small_width standard_pipeline
pd.to_pickle(ObsParams,f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pandas_data/{name}.pkl')
ObsParams.to_csv(f'/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pandas_data/{name}.csv',index=0)


#%% gen std1 data


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.make_std1_lc()

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_lc=err
msg_lc=msg


#%% find period


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.make_efsearch(nper='64')
            xte_obs.find_period()

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_per=err
msg_per=msg



#%% find period with orbital correction


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.std1_lc_orb_corr()
            xte_obs.make_efsearch(nper='64')
            xte_obs.find_period()

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_per=err
msg_per=msg



#%% make responses for fasebin files


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.make_spe()

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_spe=err
msg_make_spe=msg




#%% fits fasebin spectra - all (cutoffpl)
err=[]
msg=[]

if input('start  calculation from the beginning?')=='y':
    for k,ObsID in enumerate(ObsList_SA+ObsList_SE):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList_SA+ObsList_SE))))
        try:
            xte_obs=ObservationXTE(ObsID)
            #xte_obs.make_fasebin(nph=16)
            xte_obs.fit_ph_res(model='cutoffpl_no_gauss',error=0)
            xte_obs.fit_ph_res(model='cutoffpl',chmin=6,chmax=8,error=0.00)
            #xte_obs.fit_ph_res(model='cutoffpl_en_fix_edge_fix',error=0)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg

#%% fits fasebin spectra - SA (edge cutoffpl)
if input('start  calculation from the beginning?')=='y':
    for k,ObsID in enumerate(ObsList_SA):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList_SA))))
        try:
            xte_obs=ObservationXTE(ObsID)
            #xte_obs.make_fasebin(nph=16)
            #xte_obs.fit_ph_res(model='cutoffpl_no_gauss',error=0)
            xte_obs.fit_ph_res(model='cutoffpl_edge',error=0)

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg




#%% SE fasebin - top of a flare:
ObsList_TOP=['90089-11-04-04','90089-11-04-03','90089-11-04-02G','90089-11-04-01',
             '90089-11-04-00G','90089-11-03-05','90089-11-03-04','90089-11-03-03']

my_ph_res_model='cutoffpl_no_gauss'
err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_TOP
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            #xte_obs.make_fasebin(nph=16)
            xte_obs.fit_ph_res(model='cutoffpl_no_gauss',error=0)
            xte_obs.fit_ph_res(model='cutoffpl_en_fix',chmin=6,chmax=8,error=0.00)
            #xte_obs.fit_ph_res(model='cutoffpl_en_fix_edge_fix',error=0)


        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg



#%% SA fasebin - top of a flare -new observations:
ObsList_TOP=['90089-11-03-02','90089-11-03-01G','90089-11-03-00G']

err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_TOP
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            #xte_obs.make_fasebin(nph=16)
            #xte_obs.fit_ph_res(model='cutoffpl_no_gauss',error=0)
            #xte_obs.fit_ph_res(model='cutoffpl_en_fix',chmin=6,chmax=8,error=0.00)
            xte_obs.fit_ph_res(model='cutoffpl_en_fix_edge_fix',error=0)


        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg



#%% SE fasebin - rising part of a flare: all 12 bins
ObsList_TOP=['90089-11-02-00',
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


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_TOP
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            #xte_obs.make_fasebin(nph=12)
            #xte_obs.fit_ph_res(model='cutoffpl_no_gauss',error=0)
            xte_obs.fit_ph_res(model='cutoffpl_en_fix',chmin=6,chmax=8,error=0.00)
            #xte_obs.fit_ph_res(model='cutoffpl_en_fix_edge_fix',error=0)


        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg


#%% SE fasebin - rising part of a flare: smaller number of bins

ObsList_TOP_small=[
 '90089-11-02-01',
 '90089-11-02-02',
 '90089-11-02-03',
 '90089-11-02-05',
 '90089-11-02-07',
 '90089-11-02-08',
 '90089-11-02-09',
 '90089-11-02-10',
 '90089-11-01-02',
 '90089-11-01-03',
 '90089-11-01-04']


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_TOP_small
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.make_fasebin(nph=8)
            xte_obs.fit_ph_res(model='cutoffpl_en_fix',chmin=6,chmax=8,error=0.00)
            #xte_obs.ph_res_results()

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg






#%% SA fasebin - fading part
ObsList_TOP=['90427-01-03-00',
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


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_TOP
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            #xte_obs.make_fasebin(nph=16)
            #xte_obs.fit_ph_res(model='cutoffpl_no_gauss',error=0)
            #xte_obs.fit_ph_res(model='cutoffpl_en_fix',chmin=6,chmax=8,error=0.00)
            xte_obs.fit_ph_res(model='cutoffpl_en_fix_edge_fix',error=0)


            #xte_obs.ph_res_results()

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg



#%% SE fasebin - fading part
ObsList_TOP=[ '90014-01-05-00',
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
 '90427-01-04-05']



err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_TOP
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            #xte_obs.make_fasebin(nph=8)
            #xte_obs.fit_ph_res(model='cutoffpl_no_gauss',error=0)
            xte_obs.fit_ph_res(model='cutoffpl_en_fix',chmin=6,chmax=8,error=0.00)
            #xte_obs.fit_ph_res(model='cutoffpl_en_fix_edge_fix',error=0)


        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg


#%% SE fasebin - rising part of a flare: one-peaked profile stack

ObsList_TOP_small=[
 '90089-11-01-02',
 '90089-11-01-03',
 '90089-11-01-04']


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList_TOP_small
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.make_fasebin(nph=8)
            xte_obs.fit_ph_res(model='cutoffpl_en_fix',chmin=6,chmax=8,error=0.00)
            #xte_obs.ph_res_results()

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg


#%% FASEBIN - TIME DELAYS, PF, etc
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

#for paper plots
#ObsList=['90089-11-02-06','90089-11-02-04',
#         '90089-11-03-01G','90089-11-04-00G',
#         '90427-01-03-02','90014-01-02-10']


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    ObsList=ObsList
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.ph_res_results()

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg



#%% fasebin without gauss


err=[]
msg=[]
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            #xte_obs.make_fasebin(nph=12)
            xte_obs.fit_ph_res(model='cutoffpl_no_gauss',chmin=6,chmax=8,error=0.00)
            xte_obs.ph_res_results(model='cutoffpl_no_gauss')

        except Exception as e:
            print(e)
            print('ERROR OCCURED WITH', ObsID)
            err.append(ObsID)
            msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg





#%%% cementery




#%% make pcu2 only SE data and create spectra



errors=[]
errors_msg=[]

from pipeline_pcu2_core import recreate_se_xdf
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(ObsList):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)

            os.chdir(xte_obs.out_path)
            recreate_se_xdf('se.xdf')
            xte_obs.make_spe()

        except Exception as e:
            print(e)
            errors.append(ObsID)
            errors_msg.append(e)

#%% inhomogenities
errors=[]
errors_msg=[]

from pipeline_pcu2_core import recreate_se_xdf
if input('start  calculation from the beginning?')=='y':
    os.chdir(RXTE_path+'/data/AO9')
    ObsList=glob('*')
    for k,ObsID in enumerate(['90014-01-04-02','90014-01-05-06','90089-11-01-03','90089-11-02-08']):
        print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
        try:
            xte_obs=ObservationXTE(ObsID)
            xte_obs.get_configs()
            #os.chdir(xte_obs.out_path)
            #recreate_se_xdf('se.xdf')
            #xte_obs.make_spe()

        except Exception as e:
            print(e)
            errors.append(ObsID)
            errors_msg.append(e)

#%% check channel 10 counts

fits_files=glob('/Users/s.bykov/work/xray_pulsars/rxte/results/*/products/fasebin_spe/se.pha')
for fpath in fits_files:
    ff=fits.open(fpath)
    ch=ff[1].data['counts'][10]
    if ch!=0:
        print('NOT zero at', fpath)






