#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 15:59:38 2020

@author: s.bykov
"""


from PipelineXTE.pipeline_core import *

os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90427-01-03-11/products/fasebin_timeser_1sec')


def load_spe(num):
    data=np.genfromtxt(f'./spe_data_area/ph_spe_{num}.dat').T
    flux=data[:,2].sum()
    print(f'spectrum {num}, flux {flux}')
    return data,flux
def plot_spe(spe_data):
    plt.loglog(spe_data[:,0],spe_data[:,2])
    plt.errorbar(spe_data[:,0],spe_data[:,2],spe_data[:,3],spe_data[:,1],color='gray',fmt='none',alpha=0.5)
    plt.show()
    plt.grid(True, which="both")


fe_ind=7,8


plt.figure()


data_area=np.genfromtxt('full_spe_area_en.qdp',skip_header=3)
plt.figure()
plt.loglog(data_area[:,0],data_area[:,2],'k+')
plt.xlabel('en, keV')
plt.ylabel('data_area')
plt.show()


data=load_spe(15)
plot_spe(data[0])



data=load_spe(560)
plot_spe(data[0])



# data=load_spe(720)
# plot_spe(data[0])



# data=load_spe(920)
# plot_spe(data[0])









#%% trash

# #old xspec data
# #%% part 1: real effective stuff
# os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90089-11-04-03/products/fasebin_timeser/2sec')

# #get plots: /Users/s.bykov/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/get_xspec_data_area.txt

# #%% spectra #10

# data_area=np.genfromtxt('data_area_en.qdp',skip_header=3)
# data_noarea=np.genfromtxt('data_noarea_en.qdp',skip_header=3)

# data_area_ch=np.genfromtxt('data_area_ch.qdp',skip_header=3)
# data_noarea_ch=np.genfromtxt('data_noarea_ch.qdp',skip_header=3)



# data_noarea=data_noarea[data_noarea[:,2]!=0]
# data_noarea_ch=data_noarea_ch[data_noarea_ch[:,2]!=0]


# plt.figure()
# plt.loglog(data_area[:,0],data_area[:,2])
# test_data=np.genfromtxt('/Users/s.bykov/work/xray_pulsars/rxte/results/out90089-11-04-03/products/fasebin_timeser/2sec/spe_data_area/ph_spe_10.dat')
# test_data=test_data.T

# plt.loglog(test_data[:,0],test_data[:,2],'r.')
# plt.xlabel('en, keV')
# plt.ylabel('data_area')
# plt.show()



# plt.figure()
# plt.loglog(data_noarea[:,0],data_noarea[:,2])
# plt.xlabel('en, keV')
# plt.ylabel('data_noarea')
# plt.show()



# plt.figure()
# plt.plot(data_area[:,0],data_noarea[:,2]/data_area[:,2])
# plt.xlabel('en, keV')
# plt.ylabel('data_noarea/data_area, cm^2')
# plt.show()


# plt.figure()
# plt.plot(data_area_ch[:,0],data_noarea_ch[:,2]/data_area_ch[:,2])
# plt.xlabel('channel')
# plt.ylabel('data_noarea/data_area, cm^2')
# plt.show()



# plt.figure()
# plt.errorbar(data_area[:,0],data_area_ch[:,0],data_area_ch[:,1],data_area[:,1],'b+')
# plt.xlabel('en, keV')
# plt.ylabel('channel')
# plt.show()




# #%% spe 44

# data_noarea_ch_44=np.genfromtxt('data_noarea_ch_spe44.qdp',skip_header=3)
# data_noarea_ch_44=data_noarea_ch_44[data_noarea_ch_44[:,2]!=0]

# data_area_44=np.genfromtxt('data_area_en_spe44.qdp',skip_header=3)
# data_noarea_44=np.genfromtxt('data_noarea_en_spe44.qdp',skip_header=3)
# data_noarea_44=data_noarea_44[data_noarea_44[:,2]!=0]


# plt.figure()
# plt.plot(data_noarea_44[:,0],data_noarea_44[:,2]/data_area_44[:,2])
# plt.plot(data_area[:,0],data_noarea[:,2]/data_area[:,2])

# plt.xlabel('channel')
# plt.ylabel('data_noarea/data_area, cm^2')
# plt.show()


# # plt.figure()
# # plt.loglog(data_area_44[:,0],data_area_44[:,2])
# # plt.xlabel('en, keV')
# # plt.ylabel('data_area')
# # plt.show()

# # plt.loglog(data_area[:,0],data_noarea_ch_44[:,2]/(data_noarea[:,2]/data_area[:,2]))


# #%% mean rates

# os.chdir('/Users/s.bykov/work/xray_pulsars/rxte/results/out90089-11-03-01G/products/ironline_rate_estim')

# class TimeSeries():

#     def __init__(self,lcname):
#         self.fits=fits.open(f'{lcname}.lc_bary')
#         self.time=self.fits[1].data['time']
#         self.rate=self.fits[1].data['rate']
#         self.error=self.fits[1].data['error']
#         self.binsize=np.median(np.diff(self.time))
#         self.fits.close()
#     def divide(self,val):
#         self.rate=self.rate/val
#         self.error=self.error/val

# lc46=TimeSeries('ch1415')
# lc67=TimeSeries('ch1619')
# lc79=TimeSeries('ch2021')
# #lc712=TimeSeries('lc712')
# en=dE=np.array([6.045788+5.232003,7.270824+6.045788,8.09038+7.270824])/2

# dE=np.array([6.045788-5.232003,7.270824-6.045788,8.09038-7.270824])
# en_err=dE/2

# area_fact=np.array([4000,2733,4332])

# for k,lc in enumerate([lc46,lc67,lc79]):
#     lc=lc.divide(dE[k]*area_fact[k])

# rate=np.array([lc.rate.mean() for lc in [lc46,lc67,lc79]])
# rate_error=np.array([lc.rate.std() for lc in [lc46,lc67,lc79]])

# plt.errorbar(en,rate,rate_error,en_err)


# plt.plot(spe[:,0],spe[:,2],'k.',label='data from xspec + 36%')

# plt.xlim(4,9)
# plt.show()
