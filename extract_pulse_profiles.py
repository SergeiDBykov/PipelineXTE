#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 12:47:28 2020

@author: s.bykov
"""
#imports and defs
from PipelineXTE.pipeline_core import *
#28 SA observations

ObsList_top=[
'90089-11-03-00G',
'90089-11-03-01G',
'90089-11-03-02']

ObsList_decay=[
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
 '90014-01-04-03'
]

ObsList=ObsList_top+ObsList_decay

'''
TOP PART:
        relative- absolute
        ch 3 - 8 (3.2-3.6 keV)
channel 7 - 12 (4.825992-5.232003 keV)
channel 8 - 13 (5.232003 - 5.638603 keV)
channel 9 - 14 (5.638603-6.045788 keV)
channel 10 -  15 (6.045788-6.453556 keV)
channel 11 - 16-17 (6.453556-7.270824 keV)
channel 12 - 18-19 (7.270824 - 8.09038 keV)
channel 13 - 20-21  (8.09038-8.912192 keV)
channel 12-17 - 18-29 (7.270824- 12.22 keV)


DECAY PART:
        relative- absolute
        channel 4 (3.2-3.6 keV) - 8
channel 11 (6.042514-6.450047 keV)  - 15
channel 12 (6.450047-6.858159 keV) - 16
channel 13 (6.858159-7.266844) - 17
channel 14 (7.266844-7.676101) - 18
channel 15 (7.676101-8.085925) - 19
channel 16 (8.085925-8.496311) -20
channel 17 (8.496311-8.907259) - 21
channel 14-25 18-29 (7.26-12.21 keV)
'''





def gen_lc(ch_rel='9',ch_abs='14',binsize='0.1',outname='iron_band'):
    cmd=extract=f'saextrct infile=@sa.xdf gtiorfile=APPLY gtiandfile=maketime_file.fits outroot=products/sa_data_lc_{binsize}sec/{outname} timecol=TIME columns=GOOD accumulate=ONE binsz={binsize} printmode=lightcurve lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint=INDEF  chmin = INDEF chmax=INDEF chint={ch_abs} chbin=INDEF'
    run_command(extract,out_path=0,dull=0)


def make_light_curves(ObsID,binsize='0.05'):
    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')
    create_dir(f'products/sa_data_lc_{binsize}sec')
    xte_obs=ObservationXTE(ObsID)
    if ObsID in ObsList_top:
        gen_lc(ch_rel='3_17',ch_abs='8-29',outname='soft_band',binsize=binsize) #3-12.2 keV
    elif ObsID in ObsList_decay:
        gen_lc(ch_rel='4_25',ch_abs='8-29',outname='soft_band',binsize=binsize) #3-12 keV


    os.chdir(f'products/sa_data_lc_{binsize}sec')

    orbitfile=glob(xte_obs._orbit_path_+'/*')
    if len(orbitfile)==1:
        orbitfile=orbitfile[0]
    else:
        raise Exception('more than one orbit files found!!!')


    for lcname in glob('*lc'):
        #bcorr=f"barycorr infile={lcname} outfile={lcname}_bary orbitfiles={orbitfile} \n"
        #run_command(bcorr,out_path='./',dull=0)
        #print(bcorr)
        #os.system(bcorr)
        faxbary=f'faxbary infile={lcname} outfile={lcname}_bary orbitfiles={orbitfile} ra=5.37500000E+01  dec=5.31730003E+01 barytime=no'
        os.system(faxbary)

        import Misc
        for channel in glob('*lc_bary'):
            correct_times(f'{channel}',Misc.doppler_correction.orb_params_v0332)


STOP
#%% make individual

make_light_curves('90089-11-03-00G',binsize='0.02')
make_light_curves('90089-11-03-01G',binsize='0.02')
make_light_curves('90089-11-03-02',binsize='0.02')
make_light_curves('90427-01-03-00',binsize='0.02')





#%% make pulses in different bands

def make_light_curves_for_energy(ObsID='90427-01-03-00',binsize='0.1'):
    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')
    create_dir(f'products/sa_data_lc_{binsize}sec')
    xte_obs=ObservationXTE(ObsID)
    if ObsID in ObsList_top:
        gen_lc(ch_rel='3_17',ch_abs='8-29',outname='soft_band_5_12_keV',binsize=binsize)
        gen_lc(ch_rel='18_25',ch_abs='30-46',outname='hard_band_12_20_keV',binsize=binsize)
        gen_lc(ch_rel='26_43',ch_abs='46-116',outname='vhard_band_20_50_keV',binsize=binsize)
        gen_lc(ch_rel='44-63',ch_abs='117-249',outname='vvhard_band_50_110_keV',binsize=binsize)

    elif ObsID in ObsList_decay:
        gen_lc(ch_rel='4_25',ch_abs='8-29',outname='soft_band_5_12_keV',binsize=binsize)
        gen_lc(ch_rel='25_45',ch_abs='29-49',outname='hard_band_12_20_keV',binsize=binsize)


    os.chdir(f'products/sa_data_lc_{binsize}sec')

    orbitfile=glob(xte_obs._orbit_path_+'/*')
    if len(orbitfile)==1:
        orbitfile=orbitfile[0]
    else:
        raise Exception('more than one orbit files found!!!')


    for lcname in glob('*lc'):
        #bcorr=f"barycorr infile={lcname} outfile={lcname}_bary orbitfiles={orbitfile} \n"
        #run_command(bcorr,out_path='./',dull=0)
        #print(bcorr)
        #os.system(bcorr)
        faxbary=f'faxbary infile={lcname} outfile={lcname}_bary orbitfiles={orbitfile} ra=5.37500000E+01  dec=5.31730003E+01 barytime=no'
        os.system(faxbary)

    import Misc
    for channel in glob('*lc_bary'):
        correct_times(f'{channel}',Misc.doppler_correction.orb_params_v0332)


make_light_curves_for_energy(ObsID='90427-01-03-00')
make_light_curves_for_energy(ObsID='90089-11-03-01G')



#%% old stuff

# #%% run loop
# binsize='0.02'

# err=[]
# msg=[]
# for k,ObsID in enumerate(ObsList_top):
#     print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
#     try:
#         os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')

#         make_light_curves(ObsID,binsize=binsize)
#         os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/sa_data_lc_{binsize}sec')
#         #make_ccf(second_name='iron_band_1',binsize=binsize)
#         #make_ccf(second_name='iron_band_2',binsize=binsize)

#         read_and_plot_ccf(ObsID,second_name='iron_band_1',binsize=binsize)
#         read_and_plot_ccf(ObsID,second_name='iron_band_2',binsize=binsize)



#     except Exception as e:
#         print(e)
#         print('ERROR OCCURED WITH', ObsID)
#         err.append(ObsID)
#         msg.append(e)
# for e,m in zip(err,msg):
#     print(e,m)

# err_make_fb=err
# msg_make_fb=msg


