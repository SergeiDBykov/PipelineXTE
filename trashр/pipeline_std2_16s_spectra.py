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

create_dir('products/pcu2top_16s_spectra')

#ch   21      29  for 7-12 keV
extract=f'saextrct infile=@std2.xdf gtiorfile=APPLY gtiandfile=maketime_file_pcu2.fits outroot=products/pcu2top_16s_spectra/pcu2top_lc timecol=TIME columns=@pcu2_top.col accumulate=ONE binsz=16 printmode=lightcurve lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint=INDEF  chmin = 21 chmax=29 chint=INDEF chbin=INDEF'

run_command(extract,out_path=0,dull=0)

os.chdir('products/pcu2top_16s_spectra')
#create_dir('timestamps')
stop
#%% read lc and make chantrans


lcfile=fits.open('pcu2top_lc.lc',skip_header=3)

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
create_dir('products/pcu2top_16s_spectra/spectra')


for k in range(len(time_trans)):
    print(f'k={k} out of {len(time_trans)}')
    tmp=[time_trans[k]-5,time_trans[k]+5]


    extract=f'saextrct infile=@std2.xdf gtiorfile=APPLY gtiandfile=maketime_file_pcu2.fits outroot=products/pcu2top_16s_spectra/spectra/sp{k} timecol=TIME columns=@pcu2_top.col accumulate=ONE binsz=16 printmode=both lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint={tmp[0]}-{tmp[1]}  chmin = INDEF chmax=INDEF chint=INDEF chbin=INDEF'

    run_command(extract,out_path=0,dull=0)



#%% xspec part


#xspec - /Users/s.bykov/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/std2_16s_spectra_series_ign67keV.txt

#std2_16s_spectra_series_ign6_7keV.txt
#std2_16s_spectra_series_ign55_75keV.txt
#std2_16s_spectra_series_ign575_725keV.txt



#%% xspec plot, eqw estimation




matplotlib.rcParams['figure.figsize'] = 6.6*2, 6.6*2
matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.85
matplotlib.rcParams['figure.subplot.top']=0.95


import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')



data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2top_16s_spectra/spectra/16s_spe_5.75_7.25.dat')

#16s_spe_5.75_7.25

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




#%% xspec plot, eqw direct evaluation


matplotlib.rcParams['figure.figsize'] = 6.6*2, 6.6*2
matplotlib.rcParams['figure.subplot.left']=0.15
matplotlib.rcParams['figure.subplot.bottom']=0.15
matplotlib.rcParams['figure.subplot.right']=0.85
matplotlib.rcParams['figure.subplot.top']=0.95


import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')

data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2top_16s_spectra/spectra/16s_spe_fit.dat')

#16s_spe_5.75_7.25

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

# #%% check rate and LC


# matplotlib.rcParams['figure.figsize'] = 6.6*2, 6.6*2
# matplotlib.rcParams['figure.subplot.left']=0.15
# matplotlib.rcParams['figure.subplot.bottom']=0.15
# matplotlib.rcParams['figure.subplot.right']=0.85
# matplotlib.rcParams['figure.subplot.top']=0.95


# import seaborn as sns
# sns.set(style='ticks', palette='deep',context='notebook')

# data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/pcu2top_16s_spectra/spectra/16s_spe.dat')

# time=time_trans
# rate=data[:,2]
# rate_error=data[:,3]


# fig, ax = plt.subplots(1,gridspec_kw={'hspace': 0, 'wspace': 0})

# ax.errorbar(lcfile[1].data['time'],lcfile[1].data['rate'],lcfile[1].data['error'],fmt='.',color='c',marker='s',ms=4,alpha=0.8)

# ax.set_ylabel('Light curve rate',color='c')
# ax.set_xlabel('Time, sec')


# ax2=ax.twinx()
# ax2.errorbar(time,rate,rate_error,fmt='.',color='k',marker='s',ms=4,alpha=0.8)

# ax2.set_ylabel('Rate, count/s')

# fig.tight_layout()
# sns.despine(fig,top=1,right=0)

# plt.show()