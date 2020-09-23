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
'90089-11-03-02',
'90427-01-03-00']

ObsList_decay=[
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


def make_light_curves(ObsID,binsize='0.1'):
    os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')
    create_dir(f'products/sa_data_lc_{binsize}sec')
    xte_obs=ObservationXTE(ObsID)
    if ObsID in ObsList_top:
        gen_lc(ch_rel='10',ch_abs='15',outname='iron_band_1',binsize=binsize) #6.045788-6.453556 keV
        gen_lc(ch_rel='10_11',ch_abs='15-17',outname='iron_band_2',binsize=binsize) #6.045788-7.270824 keV
        gen_lc(ch_rel='12_17',ch_abs='18-29',outname='hard_band',binsize=binsize) #(7.270824- 12.22 keV)
    elif ObsID in ObsList_decay:
        gen_lc(ch_rel='11',ch_abs='15',outname='iron_band_1',binsize=binsize) #6.042514-6.450047 keV
        gen_lc(ch_rel='11_12',ch_abs='15-16',outname='iron_band_2',binsize=binsize) #6.042514-6.85 keV
        gen_lc(ch_rel='14_25',ch_abs='18-29',outname='hard_band',binsize=binsize)  #(7.26-12.21 keV)


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

        #for channel in glob('*lc_bary'):
        #import Misc
        #correct_times(f'{channel}',Misc.doppler_correction.orb_params_v0332)



def make_ccf(second_name='iron_band_1',binsize='0.1'):

    crosscorr=f'''crosscor  cfile1="hard_band.lc_bary" cfile2="{second_name}.lc_bary" window="-" dtnb={binsize} nbint=512 nintfm=INDEF plot=no outfile="{second_name}.ccf" '''
    os.system(crosscorr)

def read_and_plot_ccf(ObsID,second_name='iron_band_1',
                      xlims=[-0.5,10],
                      binsize='0.1'):
    matplotlib.rcParams['figure.figsize'] = 6.6*2, 6.6
    matplotlib.rcParams['figure.subplot.left']=0.15
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=0.85
    matplotlib.rcParams['figure.subplot.top']=0.9

    ccf_fits=fits.open(second_name+'.ccf')

    delay=ccf_fits[1].data['DELAY']
    crosscorr=ccf_fits[1].data['CROSS_CORR']
    error=ccf_fits[1].data['ERROR']
    N=int(len(delay)/2)


    fig,ax=plt.subplots()

    ax.errorbar(delay,crosscorr,error,drawstyle='steps-mid',label=second_name)
    ax.errorbar(delay[N:],crosscorr[0:N+1:][::-1],error[0:N+1:][::-1],alpha=0.5,color='r',drawstyle='steps-mid',ls=':')


    ax.set_xlim(xlims[0],xlims[1])
    ax.set_xlabel('Delay, s')
    ax.set_ylabel('CCF')
    ax.set_title(f'{second_name}: {ObsID}, binsize={binsize} s')
    ax.grid()
    fig.tight_layout()

    if ObsID in ObsList_top:
        if second_name=='iron_band_1':
            lgnd='7-12 keV vs 6.05-6.45 keV'
        elif second_name=='iron_band_2':
            lgnd='7-12 keV vs 6.05-7.27 keV'

    if ObsID in ObsList_decay:
        if second_name=='iron_band_1':
            lgnd='7-12 keV vs 6.04-6.45 keV'
        elif second_name=='iron_band_2':
            lgnd='7-12 keV vs 6.04-6.85 keV'

    ax.legend([lgnd])

    lc_fits=fits.open('hard_band.lc')
    mjdstart=lc_fits[1].header['MJDREFI']+lc_fits[1].header['MJDREFF']+(lc_fits[1].header['TSTARTI'])/86400
    ccf_fits.close()
    lc_fits.close()

    plt.savefig(f'Day{mjdstart}_ccf_{second_name}_{binsize}_{ObsID}.png')
    plt.close()


def int_plot_ccf(name):
    matplotlib.rcParams['figure.figsize'] = 6.6*2, 6.6
    matplotlib.rcParams['figure.subplot.left']=0.15
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=0.85
    matplotlib.rcParams['figure.subplot.top']=0.9

    ccf_fits=fits.open(name+'.ccf')

    delay=ccf_fits[1].data['DELAY']
    crosscorr=ccf_fits[1].data['CROSS_CORR']
    error=ccf_fits[1].data['ERROR']
    N=int(len(delay)/2)


    fig,ax=plt.subplots()

    ax.errorbar(delay,crosscorr,error,drawstyle='steps-mid')
    ax.errorbar(delay[N:],crosscorr[0:N+1:][::-1],error[0:N+1:][::-1],alpha=0.5,color='r',drawstyle='steps-mid',ls=':')



    ax.set_xlabel('Delay, s')
    ax.set_ylabel('CCF')
    ax.set_title(f'{name}')
    ax.grid()
    fig.tight_layout()
    plt.show()



STOP
#%% run loop


binsize='0.02'

err=[]
msg=[]
for k,ObsID in enumerate(ObsList):
    print(' =============== Obs {0} out of {1} ================'.format(str(k+1),str(len(ObsList))))
    try:
        os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')

        #make_light_curves(ObsID,binsize=binsize)
        os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/sa_data_lc_{binsize}sec')
        #make_ccf(second_name='iron_band_1',binsize=binsize)
        #make_ccf(second_name='iron_band_2',binsize=binsize)

        read_and_plot_ccf(ObsID,second_name='iron_band_1',binsize=binsize)
        read_and_plot_ccf(ObsID,second_name='iron_band_2',binsize=binsize)



    except Exception as e:
        print(e)
        print('ERROR OCCURED WITH', ObsID)
        err.append(ObsID)
        msg.append(e)
for e,m in zip(err,msg):
    print(e,m)

err_make_fb=err
msg_make_fb=msg


