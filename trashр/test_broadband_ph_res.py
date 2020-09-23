#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 12:29:25 2020

@author: s.bykov
"""
from PipelineXTE.pipeline_core import *


ObsID='90014-01-03-01'# 90089-11-02-04 ####90089-11-04-04 90089-11-03-01G 90427-01-03-05 90014-01-04-00 90014-01-03-01

datapath=f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cyclabs'

data=np.genfromtxt(datapath+'/ph_res_cyclabs.dat')
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
eqw_err=np.vstack((eqw_low, eqw_hi))

flux712=data[:,7]

flux712=flux712/1e-8

chi2=data[:,2]

ecycle=data[:,9]

po=data[:,8]


norm_line=data[:,10]*1000
norm_line_low=norm_line-data[:,11]*1000
norm_line_hi=data[:,12]*1000-norm_line
norm_line_err=np.vstack((norm_line_low, norm_line_hi))


#%% plot


fig = plt.figure(figsize=(8,6))
rows=5
cols=3
#(rows,cols), (y,x) <- those are coordinates of an axis in subplots
ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=1, colspan=3)
ax_efold=ax_eqw.twinx()
ax_chi=plt.subplot2grid((rows,cols), (1, 0), rowspan=1, colspan=3)
ax_ecycle=plt.subplot2grid((rows,cols), (2, 0), rowspan=1, colspan=3)
ax_po=plt.subplot2grid((rows,cols), (3, 0), rowspan=1, colspan=3)
ax_fe_norm=plt.subplot2grid((rows,cols), (4, 0), rowspan=1, colspan=3)


ax_eqw.set_title(ObsID)
ax_efold.plot(phase,flux712,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
ax_eqw.errorbar(phase,eqw,eqw_err,color='r',drawstyle='steps-mid',alpha=0.8)
ax_eqw.tick_params(axis='y', colors='red')
ax_eqw.spines['left'].set_color('red')

ax_eqw.set_ylabel('Fe Ka Eq. width, eV',color='r')
ax_efold.set_ylabel('Flux 7-12, 1e-8 cgs')
ax_eqw.set_xlabel('Phase')


ax_chi.plot(phase,chi2,drawstyle='steps-mid')
ax_chi.set_ylabel('chi2_red')


ax_ecycle.plot(phase,ecycle,drawstyle='steps-mid')
ax_ecycle.set_ylabel('E_cycle')

ax_po.plot(phase,po,drawstyle='steps-mid')
ax_po.set_ylabel('Gamma')



ax_fe_norm.errorbar(phase,norm_line,norm_line_err,color='r',drawstyle='steps-mid',alpha=0.6)
ax_fe_norm.set_ylabel('iron norm')


plt.show()

fig.savefig(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cyclabs/ph_res_results_{ObsID}.png')



