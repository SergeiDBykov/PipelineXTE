#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:06:23 2020

@author: s.bykov
"""


from pipeline_core import *

ObsID='90089-11-05-00G'
xte_obs=ObservationXTE(ObsID)
os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')
open_dir_in_term()

binsize=16


#%% make light curves

create_dir(f'products/pcu2top_{binsize}s_spectra')

#ch   21      29  for 7-12 keV
extract=f'saextrct infile=@std2.xdf gtiorfile=APPLY gtiandfile=maketime_file_pcu2.fits outroot=products/pcu2top_{binsize}s_spectra/pcu2top_lc timecol=TIME columns=@pcu2_top.col accumulate=ONE binsz={binsize} printmode=lightcurve lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint=INDEF  chmin = 21 chmax=29 chint=INDEF chbin=INDEF'

run_command(extract,out_path=0,dull=0)




extract=f'saextrct infile=@std2.xdf gtiorfile=APPLY gtiandfile=maketime_file_pcu2.fits outroot=products/pcu2top_{binsize}s_spectra/pcu2top_lc_all_chans timecol=TIME columns=@pcu2_top.col accumulate=ONE binsz={binsize} printmode=lightcurve lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint=INDEF  chmin = INDEF chmax=INDEF chint=INDEF chbin=INDEF'

run_command(extract,out_path=0,dull=0)



os.chdir(f'products/pcu2top_{binsize}s_spectra')
#create_dir('timestamps')
stop
#%% read lc and make chantrans

lcfile=fits.open('pcu2top_lc_all_chans.lc')

tstart=lcfile[1].header['TSTARTI']+lcfile[1].header['TSTARTF']

time_trans=lcfile[1].data['time']

assert tstart==time_trans[0], 'TZERO IS NOT ZERO!'

# for k,time_center in enumerate(time_trans):
#     left_stamp=time_center-5
#     right_stamp=time_center+5

#     cmd=f'timetrans -i pcu2top_lc.lc -f "timestamps/time_sp{k}.txt" -t "{left_stamp}-{right_stamp}"'
#     os.system(cmd)

#%% create a test spectrum, the first one (sp0)


os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')
create_dir(f'products/pcu2top_{binsize}s_spectra/spectra')


for k in range(len(time_trans)):
    print(f'k={k} out of {len(time_trans)}')
    #tmp=[time_trans[k]-binsize/2,time_trans[k]+binsize/2]
    tmp=[time_trans[k],time_trans[k]+binsize]


    extract=f'saextrct infile=@std2.xdf gtiorfile=APPLY gtiandfile=maketime_file_pcu2.fits outroot=products/pcu2top_{binsize}s_spectra/spectra/sp{k} timecol=TIME columns=@pcu2_top.col accumulate=ONE binsz={binsize} printmode=spectrum lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint={tmp[0]}-{tmp[1]}  chmin = INDEF chmax=INDEF chint=INDEF chbin=INDEF'

    run_command(extract,out_path=0,dull=0)

    if k in [0,1,2,4,5,6]:
        ff=fits.open(f'products/pcu2top_{binsize}s_spectra/spectra/sp{k}.pha')
        assert np.abs(ff[1].header['exposure']-binsize)<1, f'Spectra number {k} has wrong expo'
        ff.close()



#%% xspec part

name='std2_spectra_series_ign_ch10_14keV.txt'  #std2_spectra_series_ign_ch10_14keV.txt, std2_spectra_series_spectra_fit.txt

if name=='std2_spectra_series_ign_ch10_14keV.txt':
    outname='std2_spe_ch10_14'
elif name=='std2_spectra_series_spectra_fit.txt':
    outname='std2_spe_fit'
#xspec - /Users/s.bykov/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/std2_spectra....

#%% check rates from xspec and lc8rve

data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2top_{binsize}s_spectra/spectra/{outname}.dat')

time=time_trans
rate=data[:,2]
rate_error=data[:,3]
exposure=data[:,1]

print(exposure)

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

ax.errorbar(lcfile[1].data['time'],lcfile[1].data['rate'],lcfile[1].data['error'],fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Light curve rate',color='c')
ax.set_xlabel('Time, sec')


ax2=ax.twinx()
ax2.errorbar(time,rate,rate_error,fmt='.',color='k',marker='s',ms=2,alpha=0.5)

ax2.set_ylabel('Rate, count/s')

plt.show()




#%% xspec plot, eqw direct evaluation


matplotlib.rcParams['figure.figsize'] = 6.6*2, 6.6*2
matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.85
matplotlib.rcParams['figure.subplot.top']=0.95


import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')

data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2top_{binsize}s_spectra/spectra/{outname}.dat')

time=time_trans-time_trans[0]
rate=data[:,2]
rate_error=data[:,3]

norm_line=data[:,7]
norm_line_low=norm_line-data[:,8]
norm_line_hi=data[:,9]-norm_line
norm_line_err=np.vstack((norm_line_low, norm_line_hi))


eqw=data[:,4]
eqw_low=eqw-data[:,5]
eqw_hi=data[:,6]-eqw
eqw_err=np.vstack((eqw_low, eqw_hi))

exposure=data[:,1]

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

ax.errorbar(time,eqw*1e3,eqw_err*1e3,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line equivalent , eV',color='c')
ax.set_xlabel('Time, sec')


ax2=ax.twinx()
ax2.errorbar(time,rate,rate_error,fmt='.',color='k',marker='s',ms=4,alpha=0.8)

ax2.set_ylabel('Rate, count/s')

fig.tight_layout()
sns.despine(fig,top=1,right=0)

plt.show()



fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

ax.errorbar(time,norm_line*1e3,norm_line_err*1e3,fmt='.',color='r',marker='s',ms=4,alpha=0.8)

ax.set_ylabel('Iron line norm*1e3 ',color='r')
ax.set_xlabel('Time, sec')


ax2=ax.twinx()
ax2.errorbar(time,rate,rate_error,fmt='.',color='k',marker='s',ms=4,alpha=0.8)

ax2.set_ylabel('Rate, count/s')

fig.tight_layout()
sns.despine(fig,top=1,right=0)
plt.show()


fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

CCF=cross_correlation.CrossCorrelation(time,norm_line,rate,circular=0)
lag,ccf=CCF.calc_ccf()
ax.plot(lag,ccf,'r.-',alpha=0.6)
ax.set_xlabel('Delay, sec',fontsize=8)
ax.set_ylabel('Pearson r',fontsize=8)
plt.show()




fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

CCF=cross_correlation.CrossCorrelation(time,eqw,rate,circular=0)
lag,ccf=CCF.calc_ccf()
ax.plot(lag,ccf,'c.-',alpha=0.6)
ax.set_xlabel('Delay, sec',fontsize=8)
ax.set_ylabel('Pearson r',fontsize=8)
plt.show()





#%% xspec plot, eqw estimation


matplotlib.rcParams['figure.figsize'] = 6.6*2, 6.6*2
matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.85
matplotlib.rcParams['figure.subplot.top']=0.95


import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')


oudata=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2top_{binsize}s_spectra/spectra/{outname}.dat')



time=time_trans-time_trans[0]
rate=data[:,2]
rate_error=data[:,3]

eqw=data[:,4]
eqw_err=data[:,5]

chi2=data[:,6]

fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

ax.errorbar(time,eqw,eqw_err,fmt='.',color='c',marker='s',ms=4,alpha=0.8)

ax.plot(time[chi2<=1.5],eqw[chi2<=1.5],marker='s',color='none',mec='r' )

ax.set_ylabel('Iron line equivalent , eV',color='c')
ax.set_xlabel('Time, sec')


ax2=ax.twinx()
ax2.errorbar(time,rate,rate_error,fmt='.',color='k',marker='s',ms=4,alpha=0.8)

ax2.set_ylabel('Rate, count/s')

fig.tight_layout()
sns.despine(fig,top=1,right=0)

plt.show()


