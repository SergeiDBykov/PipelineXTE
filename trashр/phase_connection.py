#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 15:46:05 2020

@author: s.bykov
"""
#%% imports and defs
from PipelineXTE.pipeline_core import *
from Misc import doppler_correction as doppler
from Misc.doppler_correction import day2sec

results_path='/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pandas_data/'

filename='standard_pipeline'
ObsParams=pd.read_pickle(results_path+f'{filename}.pkl')
ObsParams=ObsParams.sort_values('MJD_START')

ObsParams.period_orb_corr= ObsParams.period_orb_corr.replace(to_replace='None',value=np.nan)
ObsParams.period_orb_corr_err= ObsParams.period_orb_corr_err.replace(to_replace='None',value=np.nan)

t0=53384.36
def line(t,p0,pdot):
    return p0+(t-t0)*pdot


time=ObsParams.MJD_START[~ObsParams.period_orb_corr.isna()]

period=ObsParams.period_orb_corr[~ObsParams.period_orb_corr.isna()]
period_err=ObsParams.period_orb_corr_err[~ObsParams.period_orb_corr.isna()]

popt,pcov=curve_fit(line,time[(53384.36<time) & (time<53400)],period[(53384.36<time) & (time<53400)],p0=[4.373,4e-6],
                    sigma=period_err[(53384.36<time) & (time<53400)],absolute_sigma=1)

print(f'Decay part of the outburst before 53400: P(t)=P0+Pdot*(t-{t0}), P0={popt[0]} Pdot={popt[1]}+-{np.sqrt(np.diag(pcov)[1])}')
P0=popt[0]
Pdot=popt[1]/86400

def get_mjd_tdb_and_dopplec_corr(ObsID):
    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin')
    with fits.open('fasebin_orig.pha') as f:
        mjd_tdb=f[1].header['TSTART']
        f.close()
    dc=doppler.kepler_solution(mjd_tdb*day2sec,doppler.orb_params_v0332)
    print(f'###\n ObsID {ObsID} \n MJD_TSTART (TDB): {mjd_tdb} \n romer delay (z/c): {dc[-1]} \n ####')
    return mjd_tdb,dc[0][-1]


def deltaT_if_pdot(N,P0,Pdot):
    '''
    t_i=t0+i*P0+1/2 i^2 P Pdot - the time of arrival of the i-th pulse
    '''
    return 0.5*N**2*P0*Pdot

#%% load obs
ObsID1='90427-01-03-00'

t1_tdb,z1=get_mjd_tdb_and_dopplec_corr(ObsID1)
t1=t1_tdb+z1/day2sec

ObsID2='90427-01-03-14G'  #ObsID2='90427-01-03-14G', 90014-01-03-03
t2_tdb,z2=get_mjd_tdb_and_dopplec_corr(ObsID2)
t2=t2_tdb+z2/day2sec

xte_obs=ObservationXTE(ObsID1)
ser=xte_obs.pandas_series()
period=ser['period_orb_corr']

N_pulses=(t2-t1)*day2sec/period
print('N pulses: ',N_pulses)

deltaT=deltaT_if_pdot(N_pulses, P0, Pdot)
print('deltaT (s)= ',deltaT)

mjd_null_obs1=t1
mjd_null_obs2=t1-deltaT/day2sec

print(f'MJD_starts: \n  {ObsID1}: {mjd_null_obs1} \n  {ObsID2}: {mjd_null_obs2}')


#%% check lcurves

def get_profile(ObsID):
    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin')
    with fits.open('fasebin_orig_newt0.pha') as f:
        phase=f[1].data['phase']
        counts=f[1].data['counts'].sum(axis=1)
        f.close()
    return phase,counts/np.mean(counts)

plt.figure()
pp1=get_profile(ObsID1)
plt.step(pp1[0],pp1[1],where='mid')

pp2=get_profile(ObsID2)
plt.step(pp2[0],pp2[1],where='mid')

plt.show()
