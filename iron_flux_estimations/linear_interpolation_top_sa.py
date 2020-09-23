#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:06:23 2020

@author: s.bykov
"""


from pipeline_core import *

ObsID='90089-11-03-01G' # 90089-11-03-00G or 90089-11-03-01G 90089-11-03-02
xte_obs=ObservationXTE(ObsID)
os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')


binsize=0.1



STOP
#%% factors of area and energy width


os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin_spe/conversion_factors')


data_area_en=np.genfromtxt('data_area_en.qdp',skip_header=3)

data_noarea_en=np.genfromtxt('data_noarea_en.qdp',skip_header=3)

data_area_ch=np.genfromtxt('data_area_ch.qdp',skip_header=3)
data_area_ch[:,0]=data_area_ch[:,0]-1
data_noarea_ch=np.genfromtxt('data_noarea_ch.qdp',skip_header=3)
data_noarea_ch[:,0]=data_area_ch[:,0]-1


plt.figure()
plt.plot(data_area_en[:,0],data_noarea_en[:,2]/data_area_en[:,2])
plt.xlabel('en, keV')
plt.ylabel('data_noarea/data_area, cm^2')
plt.show()


plt.figure()
plt.plot(data_area_ch[:,0],data_noarea_ch[:,2]/data_area_ch[:,2])
plt.xlabel('channel')
plt.ylabel('data_noarea/data_area, cm^2')
plt.show()



plt.figure()
plt.errorbar(data_area_en[:,0],data_area_ch[:,0],data_area_ch[:,1],data_area_en[:,1],'b+')
plt.xlabel('en, keV')
plt.ylabel('channel')

plt.show()

#%% print factors
print(f'''
Channel 9:
        energy = {data_noarea_en[6][0]}+-{data_noarea_en[6][1]}
    eff. area = {data_noarea_ch[6][2]/data_area_ch[6][2]}  cm^2
    energy width = {data_noarea_en[6][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[6][1]*2* (data_noarea_ch[6][2]/data_area_ch[6][2]))}


Channel 10:
        energy = {data_noarea_en[7][0]}+-{data_noarea_en[7][1]}
    eff. area = {data_noarea_ch[7][2]/data_area_ch[7][2]}  cm^2
    energy width = {data_noarea_en[7][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[7][1]*2* (data_noarea_ch[7][2]/data_area_ch[7][2]))}


Channel 11:
        energy = {data_noarea_en[8][0]}+-{data_noarea_en[8][1]}
    eff. area = {data_noarea_ch[8][2]/data_area_ch[8][2]}  cm^2
    energy width = {data_noarea_en[8][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[8][1]*2* (data_noarea_ch[8][2]/data_area_ch[8][2]))}

Channel 12:
        energy = {data_noarea_en[9][0]}+-{data_noarea_en[9][1]}
    eff. area = {data_noarea_ch[9][2]/data_area_ch[9][2]}  cm^2
    energy width = {data_noarea_en[9][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[9][1]*2* (data_noarea_ch[9][2]/data_area_ch[9][2]))}


Channel 13:
    energy = {data_noarea_en[10][0]}+-{data_noarea_en[10][1]}
    eff. area = {data_noarea_ch[10][2]/data_area_ch[10][2]}  cm^2
    energy width = {data_noarea_en[10][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[10][1]*2* (data_noarea_ch[10][2]/data_area_ch[10][2]))}
      ''')


#%% make light curves
os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')


'''
channel 7 - 12 (4.825992-5.232003 keV)
channel 8 - 13 (5.232003 - 5.638603 keV)
channel 9 - 14 (5.638603-6.045788 keV)
channel 10 -  15 (6.045788-6.453556 keV)
channel 11 - 16-17 (6.453556-7.270824 keV)
channel 12 - 18-19 (7.270824 - 8.09038 keV)
channel 13 - 20-21  (8.09038-8.912192 keV)
channel 12-17 - 18-29 (7.270824- 12.22 keV)
'''

create_dir(f'products/sa_data_lc_{binsize}sec')



def gen_lc(ch_rel='9',ch_abs='14'):
    cmd=extract=f'saextrct infile=@sa.xdf gtiorfile=APPLY gtiandfile=maketime_file.fits outroot=products/sa_data_lc_{binsize}sec/ch{ch_rel} timecol=TIME columns=GOOD accumulate=ONE binsz={binsize} printmode=lightcurve lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint=INDEF  chmin = INDEF chmax=INDEF chint={ch_abs} chbin=INDEF'
    run_command(extract,out_path=0,dull=0)



gen_lc('9','14')

gen_lc('10','15')

gen_lc('11','16-17')

gen_lc('12','18-19')

gen_lc('13','20-21')

gen_lc('12_17','18-29')


#%% barycorr and orbit correction
os.chdir(f'products/sa_data_lc_{binsize}sec')


orbitfile=glob(xte_obs._orbit_path_+'/*')[0]



for lcname in glob('*lc'):
    faxbary=f'faxbary infile={lcname} outfile={lcname}_bary orbitfiles={orbitfile} ra=5.37500000E+01  dec=5.31730003E+01 barytime=no'
    os.system(faxbary)

for lcname in glob('*lc_bary'):
    import Misc
    correct_times(f'{lcname}',Misc.doppler_correction.orb_params_v0332)


STOP
#%% make new timeseries

'''

Channel 9:
        energy = 5.8421955+-0.20359278
    eff. area = 2646.0869228735896  cm^2
    energy width = 0.40718556
    1/(eff_area*en width) = 0.0009281187047645966


Channel 10:
        energy = 6.249672+-0.20388365
    eff. area = 2816.576150699146  cm^2
    energy width = 0.4077673
    1/(eff_area*en width) = 0.0008706951129767611


Channel 11:
        energy = 6.86219+-0.4086342
    eff. area = 4367.3066357014395  cm^2
    energy width = 0.8172684
    1/(eff_area*en width) = 0.0002801699853467043

Channel 12:
        energy = 7.680602+-0.40977788
    eff. area = 4491.555703678442  cm^2
    energy width = 0.81955576
    1/(eff_area*en width) = 0.00027165937860035115


Channel 13:
    energy = 8.501286+-0.4109063
    eff. area = 4563.682412291262  cm^2
    energy width = 0.8218126
    1/(eff_area*en width) = 0.00026663170529039713

'''

path_to_lc=f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/sa_data_lc_{binsize}sec'

class TimeSeries():

    def __init__(self,lcname):
        self.fits=fits.open(path_to_lc+f'/{lcname}.lc_bary_orb_corr')
        self.time=self.fits[1].data['time']
        self.rate=self.fits[1].data['rate']
        self.error=self.fits[1].data['error']
        self.binsize=np.median(np.diff(self.time))
        self.fits.close()
    def divide(self,val):
        self.rate=self.rate/val
        self.error=self.error/val

#%%read lcurves and plot mean spe
en=np.array([5.8421955,6.249672,6.86219,7.680602,8.501286])
enerr=np.array([0.2034719,0.20388365,0.4086342,0.40977788,0.4109063])

factor=[0.0009281187047645966,0.0008706951129767611,0.0002801699853467043,
        0.00027165937860035115,0.00026663170529039713]

lc9=TimeSeries('ch9')
lc10=TimeSeries('ch10')
lc11=TimeSeries('ch11')
lc12=TimeSeries('ch12')
lc13=TimeSeries('ch13')

lclist=[lc9,lc10,lc11,lc12,lc13]

for k,lc in enumerate(lclist):
    lc=lc.divide(1/factor[k])


rate=np.array([lc.rate.mean() for lc in lclist])
rate_error=np.array([lc.rate.std() for lc in lclist])


plt.figure()

spe=np.genfromtxt(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/fasebin_spe/conversion_factors/data_area_en.qdp',skip_header=3)



plt.errorbar(en,rate,rate_error,enerr,label='spectra_From_lc')
plt.plot(spe[:,0],spe[:,2],'k.',label='spectra_from_xspec (with setpl area, setpl energy)')



a,b=np.polyfit(en[[0,3,4]],rate[[0,3,4]],1)

enaxis=np.linspace(en[0],en[-1],10)
plt.plot(enaxis,a*enaxis+b,label='linear best model ignoring 6-7 kev')

def best_line(x,N):
    return x*a+N

#plt.xscale('log')
plt.show()
plt.ylabel(' counts/s/keV/cm^2')
plt.xlabel('energy')
plt.legend()

plt.show()
plt.savefig(f'lc_spectra.png')
plt.close()




#%% find residual flux

def find_fe_flux(lclist,frac=1):
    N=len(lclist[0].rate)

    diff=[]
    diff_err=[]
    print(f'{frac*100}% OF THE LINEAR APPROX IS USED AS CONTINUUM')


    for i in range(N):
        rate=np.array([x.rate[i] for x in lclist])


        N_opt,_=curve_fit(best_line,en[[0,3,4]],rate[[0,3,4]],
                                  p0=b)
        myline=lambda x: a*x+N_opt


        iron_band_flux=np.average([rate[1],rate[2]],weights=[enerr[1],enerr[2]])
        iron_band_en=np.mean([en[1],en[2]])


        fe_flux=iron_band_flux-myline(iron_band_en)*frac

        diff.append(fe_flux)

    diff=np.asarray(diff)
    diff=np.reshape(diff,N)
    return diff

frac=1
fe_flux=find_fe_flux(lclist,frac=frac)

figure()
plt.plot(lclist[0].time,fe_flux)
plt.xlabel('Time, s')
plt.ylabel('Iron line flux')
plt.show()

plt.savefig(f'fe_flux_lc_frac.png')
plt.close()


figure()
plt.hist(fe_flux,bins=100,label='Fe_flux')
plt.xlabel('Iron line flux ')
plt.show()

plt.savefig(f'fe_flux_hist_{frac}.png')
plt.close()


print(f'{frac*100}% OF THE LINEAR APPROX IS USED AS CONTINUUM')
print(f'Mean Flux: {np.mean(fe_flux)}')
print(f'Std Flux: {np.std(fe_flux)}')
print(f' mean/ std : {np.mean(fe_flux)/np.std(fe_flux)}')



#%% save data

create_dir('fe_line')
os.system(f'cp ch11.lc_bary_orb_corr fe_line/python_lin_approx_fe_line_{frac}.lc_bary_orb_corr')
with fits.open(f'fe_line/python_lin_approx_fe_line_{frac}.lc_bary_orb_corr', mode='update') as hdul:
    hdul[1].data['rate']=fe_flux
    hdul[1].data['error']=fe_flux/10
    hdul.flush()  # changes are written back to original.fits


#%% calc ccfs

crosscorr=f'''crosscor  cfile1="ch12_17.lc_bary_orb_corr" cfile2="fe_line/python_lin_approx_fe_line_{frac}.lc_bary_orb_corr" window="-" dtnb={binsize} nbint=256 nintfm=INDEF plot=no outfile="fe_line/python_lin_approx_fe_line_{frac}.ccf" '''
os.system(crosscorr)



#%% plot ccfs

import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook',rc={"xtick.top" : True,'xtick.direction':'inout','ytick.direction':'inout','xtick.minor.visible':True,'ytick.minor.visible':True})


def read_and_plot_ccf(ObsID,name='python_lin_approx_fe_line_1',
                      xlims=[-0,10],
                      binsize='4.4'):
    matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/2
    matplotlib.rcParams['figure.subplot.left']=0.15
    matplotlib.rcParams['figure.subplot.bottom']=0.15
    matplotlib.rcParams['figure.subplot.right']=0.85
    matplotlib.rcParams['figure.subplot.top']=0.9

    ccf_fits=fits.open('fe_line/'+name+'.ccf')

    delay=ccf_fits[1].data['DELAY']
    crosscorr=ccf_fits[1].data['CROSS_CORR']
    error=ccf_fits[1].data['ERROR']
    N=int(len(delay)/2)


    fig,ax=plt.subplots()

    ax.errorbar(delay,crosscorr,error,drawstyle='steps-mid',label=name)
    ax.errorbar(delay[N:],crosscorr[0:N+1:][::-1],error[0:N+1:][::-1],alpha=0.5,color='r',drawstyle='steps-mid',ls=':')


    ax.set_xlim(xlims[0],xlims[1])
    ax.set_xlabel('Delay, s')
    ax.set_ylabel('CCF')
    #ax.set_title(f'{name}: {ObsID}, binsize={binsize} s')
    ax.grid()
    fig.tight_layout()



    plt.savefig(f'fe_line/python_lin_approx_fe_line_{frac}_{ObsID}_ccf.pdf')
    plt.close()



read_and_plot_ccf(ObsID,name=f'python_lin_approx_fe_line_{frac}')


#%% plot three point

N=lc10.rate.shape[0]

fig,ax=plt.subplots(3)
from random import randint
randi=[randint(0, N) for p in range(0, 3)]
#randi=[100,101,102]
for k,i in enumerate(randi):
    rate=np.array([x.rate[i] for x in lclist])
    rate_error=np.array([x.error[i] for x in lclist])
    ax[k].errorbar(en,rate,rate_error,enerr,label='spectra_From_lc')

    N_opt,N_opt_err=curve_fit(best_line,en[[0,3,4]],rate[[0,3,4]],
                              p0=b,sigma=rate_error[[0,3,4]],absolute_sigma=True)
    myline=lambda x: a*x+N_opt
    ax[k].plot(enaxis,myline(enaxis),'k:',alpha=0.3)



    iron_band_flux=np.average([rate[1],rate[2]],weights=[enerr[1],enerr[2]])
    iron_band_en=np.mean([en[1],en[2]])
    iron_band_flux_err=np.sqrt(rate_error[1]**2+rate_error[2]**2)

    ax[k].errorbar(iron_band_en,iron_band_flux,iron_band_flux_err,[en[2]-en[1]],color='r')

    fe_flux=iron_band_flux-myline(iron_band_en)
    fe_flux_err=iron_band_flux_err

    ax[k].vlines(iron_band_en,iron_band_flux,iron_band_flux-fe_flux,color='c',alpha=0.99)




plt.show()


plt.show()
plt.savefig(f'lc_spectra_few_points.png')
plt.close()
