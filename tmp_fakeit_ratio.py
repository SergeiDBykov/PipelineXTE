#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 10:58:10 2020

@author: s.bykov
"""

from pipeline_core import *



os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90089-11-03-00G_group/products/fasebin/cutoffpl_en_fix_edge_fix')



data1=np.genfromtxt(f'ph9.lda').T

data2=np.genfromtxt(f'ph4.lda').T

fig,ax = plt.subplots(figsize=(8, 6))
label=f"phase 9 lda fakeit / phase 4 lda phakeit"
def ratio_error(a,b,da,db):
    f=a/b
    sigma=np.abs(f)*np.sqrt( (da/a)**2 + (db/b)**2  )
    return f, sigma


rat,delta=ratio_error(data1[2],data2[2],data1[3],data2[3])
energy=data1[0]

ax.errorbar(energy,rat,delta,data1[1],label=label,drawstyle='steps-mid',ls=':',alpha=0.6)
ax.set_xscale('log')

ax.legend(loc='upper left',fontsize=8)
ax.grid('y')
ax.set_ylabel('spectral ratio (lda)',fontsize=8)

plt.show()