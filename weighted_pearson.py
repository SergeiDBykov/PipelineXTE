#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 15:37:06 2020

@author: s.bykov
"""


from PipelineXTE.pipeline_core import *
from Misc.TimeSeries import cross_correlation
from scipy import stats

STOP
#%% load data
ObsID='90089-11-03-01G' #90089-11-03-01G  #90427-01-03-02 #90089-11-02-06  #90089-11-01-03 90427-01-04-00 90014-01-02-08
data=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin/cutoffpl/ph_res_cutoffpl.dat')

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
eqw_err=np.vstack((eqw_low, eqw_hi)).max(axis=0)

flux712=data[:,7]
flux712_low=flux712-data[:,8]
flux712_hi=data[:,9]-flux712

flux712=flux712/1e-8
flux712_hi=flux712_hi/1e-8
flux712_low=flux712_low/1e-8

flux712_err=np.vstack((flux712_low, flux712_hi)).max(axis=0)


fig,[ax1,ax2]=plt.subplots(2)
ax1.set_title(ObsID)
ax1.errorbar(phase,eqw,eqw_err,color='c')
ax1.set_ylabel('eqw')
ax1_twin=ax1.twinx()
ax1_twin.errorbar(phase,flux712,flux712_err,color='k')
ax1_twin.set_ylabel('flux')


plt.show()

create_dir('ccf_test_heasoft')

# tmp=np.vstack((phase,flux712,flux712_err)).T
# np.savetxt(f'ccf_test_heasoft/flux712.txt',tmp,delimiter=' ')

# tmp=np.vstack((phase,eqw,eqw_err)).T
# np.savetxt(f'ccf_test_heasoft/eqw.txt',tmp,delimiter=' ')

# tmp=np.vstack((phase,eqw/eqw,eqw_err/eqw_err)).T
# np.savetxt(f'ccf_test_heasoft/const.txt',tmp,delimiter=' ')


#%% ccf of weighted pearson

def m(x, w):
    """Weighted Mean"""
    return np.average(x, weights=w)

def cov(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)

def weighted_pearsonr(x, y, w):
    """Weighted Correlation"""
    return cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))



def cross_correlation_weighted_pearson(x,y1,y2,w1,
                      circular=1,
                      divide_by_mean=1, subtract_mean=1):
    '''
    Computes ccf of two timeseries.
    lags are those of the second array relative to the first. I.E.
    peaking on the negative delays means that the seconds LAGS the first
    peaking on the positive values means that the second precedes the first
    '''
    dx=np.median(np.diff(x))

    if divide_by_mean:
        y1=y1/np.mean(y1)
        y2=y2/np.mean(y2)

    if subtract_mean:
        y1=y1-np.mean(y1)
        y2=y2-np.mean(y2)
    lags_index=np.arange(len(y1))
    if circular:
        #ccf_pos=np.array([stats.pearsonr(y1,np.roll(y2,lag))[0] for lag in lags_index])
        ccf_pos=np.array([weighted_pearsonr(y1,np.roll(y2,lag),w1) for lag in lags_index])
        ccf_neg=np.array([weighted_pearsonr(y1,np.roll(y2,lag),w1) for lag in -lags_index])


        ccf=np.concatenate((ccf_neg,ccf_pos))

        lags=np.concatenate((-lags_index[lags_index],lags_index[lags_index]))*dx

        tmp = lags.argsort()
        return lags[tmp],ccf[tmp]

lags,weighted_ccf=cross_correlation_weighted_pearson(x=phase,
                                                     y1=eqw,w1=1/eqw_err**2,
                                                     y2=flux712)


ax2.plot(lags*4.374,weighted_ccf,color='g',label='weighted')

ax2.grid()




def cross_correlation_pearson(x,y1,y2,
                      circular=1,
                      divide_by_mean=1, subtract_mean=1):
    '''
    Computes ccf of two timeseries.
    lags are those of the second array relative to the first. I.E.
    peaking on the negative delays means that the seconds LAGS the first
    peaking on the positive values means that the second precedes the first
    '''
    dx=np.median(np.diff(x))

    if divide_by_mean:
        y1=y1/np.mean(y1)
        y2=y2/np.mean(y2)

    if subtract_mean:
        y1=y1-np.mean(y1)
        y2=y2-np.mean(y2)
    lags_index=np.arange(len(y1))
    if circular:
        ccf_pos=np.array([stats.pearsonr(y1,np.roll(y2,lag))[0] for lag in lags_index])
        ccf_neg=np.array([stats.pearsonr(y1,np.roll(y2,lag))[0] for lag in -lags_index])
        ccf=np.concatenate((ccf_neg,ccf_pos))

        lags=np.concatenate((-lags_index[lags_index],lags_index[lags_index]))*dx

        tmp = lags.argsort()
        return lags[tmp],ccf[tmp]





lags,weighted_ccf=cross_correlation_pearson(x=phase,
                                                     y1=eqw,
                                                     y2=flux712)


ax2.plot(lags*4.374,weighted_ccf,color='r',label='unweighted')

ax2.grid()
ax2.legend()
plt.show()




def cross_correlation_spearmanr(x,y1,y2,
                      circular=1,
                      divide_by_mean=1, subtract_mean=1):
    '''
    Computes ccf of two timeseries.
    lags are those of the second array relative to the first. I.E.
    peaking on the negative delays means that the seconds LAGS the first
    peaking on the positive values means that the second precedes the first
    '''
    dx=np.median(np.diff(x))

    if divide_by_mean:
        y1=y1/np.mean(y1)
        y2=y2/np.mean(y2)

    if subtract_mean:
        y1=y1-np.mean(y1)
        y2=y2-np.mean(y2)
    lags_index=np.arange(len(y1))
    if circular:
        ccf_pos=np.array([stats.spearmanr(y1,np.roll(y2,lag))[0] for lag in lags_index])
        ccf_neg=np.array([stats.spearmanr(y1,np.roll(y2,lag))[0] for lag in -lags_index])
        ccf=np.concatenate((ccf_neg,ccf_pos))

        lags=np.concatenate((-lags_index[lags_index],lags_index[lags_index]))*dx

        tmp = lags.argsort()
        return lags[tmp],ccf[tmp]



lags,weighted_ccf=cross_correlation_spearmanr(x=phase,
                                                     y1=eqw,
                                                     y2=flux712)


ax2.plot(lags*4.374,weighted_ccf,color='y',label='spearmanr')

ax2.grid()
ax2.legend()
plt.show()




# #test eqw trials
# N=500
# eqw_trials=np.zeros(shape=(N,len(phase)))
# for i in range(N):
#     eqw_trial=eqw+np.random.normal(loc=0,scale=eqw_err)
#     eqw_trials[i]=eqw_trial
#     #plt.plot(phase,eqw_trial,'gray',alpha=0.7)
# plt.errorbar(phase,eqw,eqw_err,color='r',zorder=10,capsize=3)
# plt.errorbar(phase,eqw_trials.mean(axis=0),eqw_trials.std(axis=0),color='c',zorder=15,capsize=3)

# pearsonr_arr=np.zeros(N)

# for i in range(N):
#     pearsonr_arr[i]=pearsonr(eqw_trials[i], flux712)[0]

# N_trials=1000

# ccf_trials=np.zeros(shape=(N_trials,2*len(phase)))

# test_eqw=np.zeros(N_trials)
# test_eqw_ind=12

# test_lag_of_max=np.zeros(N_trials)

# for i in range(N_trials):
#     print(i)
#     eqw_trial=np.random.normal(loc=eqw,scale=eqw_err)
#     test_eqw[i]=eqw_trial[test_eqw_ind]
#     #flux_trial=np.random.normal(loc=flux712,scale=flux712_err)
#     flux_trial=flux712
#     CCF=CrossCorrelation(phase,eqw_trial,flux_trial,circular=1)
#     lag,ccf=CCF.calc_ccf()
#     test_lag_of_max[i]=lag[lag>=0][np.argmax(ccf[lag>=0])]
#     ccf_trials[i]=ccf

# fig,ax=plt.subplots()
# #ax.errorbar(lag,ccf_trials.mean(axis=0),ccf_trials.std(axis=0)*1.645)

# CCF=CrossCorrelation(phase,eqw,flux712,circular=True)
# lag_orig,ccf_orig=CCF.calc_ccf()
# ax.plot(lag_orig,ccf_orig,'b-')
# ax_tmp=ax.twinx()
# ax_tmp.hist(test_lag_of_max,color='m',alpha=0.6,lw=0.8,histtype='step',bins=25)
# ax_tmp.set_ylabel('peak distribution', color='m')

# # plt.figure()
# # plt.plot(lag,ccf-ccf_trials.mean(axis=0),'k:')
# # plt.plot(lag,ccf-np.median(ccf_trials,axis=0),'c-.')

# # fig,ax=plt.subplots()
# # plt.hist(test_eqw,bins=50)
# # plt.axvline(eqw[test_eqw_ind])

# # plt.figure()
# # from scipy import stats
# # import seaborn as sns

# # sns.distplot(ccf_trials[:,11],fit=stats.norm,kde=0)
# # plt.axvline(ccf[11])

# plt.figure()
# from scipy import stats
# import seaborn as sns
# from scipy.stats import norm
# ax.twinx().hist(test_lag_of_max)#,fit=norm, kde=False,bins=25,hist=1)




# if __name__=='main':
# #%%test
# x=np.linspace(0,4,1000)*np.pi
# #y1=np.sin(x)+np.random.normal(x-x,0.5)
# #y2=np.cos(x)+np.random.normal(x-x,0.5)
# def gaussian(x, mu, sig):
#     return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

# y1=gaussian(x, 2, 0.4)
# y2=gaussian(x, 3, 0.4)
# dt=np.median(np.diff(x))
# fig,[ax1,ax2]=plt.subplots(2,1)
# fig.subplots_adjust(hspace=0.5)
# ax1.plot(x,y1,x,y2)


# CCF_obj=CrossCorrelation(x, y1, y2,circular=0)
# CCF_obj.calc_ccf()

# ax2.plot(CCF_obj.lag,CCF_obj.ccf,'g-.',ms=1)

# CCF_obj=CrossCorrelation(x, y1, y2,circular=1)
# CCF_obj.calc_ccf()
# ax2.plot(CCF_obj.lag,CCF_obj.ccf,'b-.',ms=1)
# plt.show()
# max_stuff=CCF_obj.find_max()

# fig,ax=plt.subplots()

# ax.plot(CCF_obj.y1,np.roll(CCF_obj.y2,max_stuff[-1][-1]))

# plt.show()

