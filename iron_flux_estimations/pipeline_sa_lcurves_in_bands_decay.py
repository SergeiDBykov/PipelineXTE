#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:06:23 2020

@author: s.bykov
"""


from pipeline_core import *

ObsID='90014-01-03-03' #  90427-01-03-00 90427-01-03-02 90427-01-03-14G 90014-01-02-00 90427-01-03-05 90427-01-03-12 90427-01-03-09 90427-01-03-06 90014-01-03-02(020) 90014-01-03-03
xte_obs=ObservationXTE(ObsID)
os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')
open_dir_in_term()

binsize=5




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
    energy = {data_noarea_en[5][0]}+-{data_noarea_en[5][1]}
    eff. area = {data_noarea_ch[5][2]/data_area_ch[5][2]}  cm^2
    energy width = {data_noarea_en[5][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[5][1]*2* (data_noarea_ch[5][2]/data_area_ch[5][2]))}


Channel 10:
    energy = {data_noarea_en[6][0]}+-{data_noarea_en[6][1]}
    eff. area = {data_noarea_ch[6][2]/data_area_ch[6][2]}  cm^2
    energy width = {data_noarea_en[6][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[6][1]*2* (data_noarea_ch[6][2]/data_area_ch[6][2]))}



Channel 11:
    energy = {data_noarea_en[7][0]}+-{data_noarea_en[7][1]}
    eff. area = {data_noarea_ch[7][2]/data_area_ch[7][2]}  cm^2
    energy width = {data_noarea_en[7][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[7][1]*2* (data_noarea_ch[7][2]/data_area_ch[7][2]))}


Channel 12:
        energy = {data_noarea_en[8][0]}+-{data_noarea_en[8][1]}
    eff. area = {data_noarea_ch[8][2]/data_area_ch[8][2]}  cm^2
    energy width = {data_noarea_en[8][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[8][1]*2* (data_noarea_ch[8][2]/data_area_ch[8][2]))}

Channel 13:
        energy = {data_noarea_en[9][0]}+-{data_noarea_en[9][1]}
    eff. area = {data_noarea_ch[9][2]/data_area_ch[9][2]}  cm^2
    energy width = {data_noarea_en[9][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[9][1]*2* (data_noarea_ch[9][2]/data_area_ch[9][2]))}


Channel 14:
        energy = {data_noarea_en[10][0]}+-{data_noarea_en[10][1]}
    eff. area = {data_noarea_ch[10][2]/data_area_ch[10][2]}  cm^2
    energy width = {data_noarea_en[10][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[10][1]*2* (data_noarea_ch[10][2]/data_area_ch[10][2]))}


Channel 15:
        energy = {data_noarea_en[11][0]}+-{data_noarea_en[11][1]}
    eff. area = {data_noarea_ch[11][2]/data_area_ch[11][2]}  cm^2
    energy width = {data_noarea_en[11][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[11][1]*2* (data_noarea_ch[11][2]/data_area_ch[11][2]))}


Channel 16:
        energy = {data_noarea_en[12][0]}+-{data_noarea_en[12][1]}
    eff. area = {data_noarea_ch[12][2]/data_area_ch[12][2]}  cm^2
    energy width = {data_noarea_en[12][1]*2}
    1/(eff_area*en width) = {1/(data_noarea_en[12][1]*2* (data_noarea_ch[12][2]/data_area_ch[12][2]))}




      ''')

#check:

# print(data_noarea_ch[7][2]*0.00087130476998205-data_area_en[7][2])
# print(data_noarea_ch[8][2]*0.0005685587704151747-data_area_en[8][2])
os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')


#%% make light curves
os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')
create_dir(f'products/sa_data_lc_{binsize}sec')
open_dir_in_term()
'''
channel 11 (6.042514-6.450047 keV)  - 15
channel 12 (6.450047-6.858159 keV) - 16
channel 13 (6.858159-7.266844) - 17
channel 14 (7.266844-7.676101) - 18
channel 15 (7.676101-8.085925) - 19
channel 16 (8.085925-8.496311) -20
channel 17 (8.496311-8.907259) - 21
'''


def gen_lc(ch_rel='10',ch_abs='14'):
    cmd=extract=f'saextrct infile=@sa.xdf gtiorfile=APPLY gtiandfile=maketime_file.fits outroot=products/sa_data_lc_{binsize}sec/ch{ch_rel} timecol=TIME columns=GOOD accumulate=ONE binsz={binsize} printmode=lightcurve lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint=INDEF  chmin = INDEF chmax=INDEF chint={ch_abs} chbin=INDEF'
    run_command(extract,out_path=0,dull=0)

gen_lc('9','13')
gen_lc('10','14')
gen_lc('11','15')
gen_lc('12','16')
gen_lc('13','17')
gen_lc('14','18')
#gen_lc('15','19')
#gen_lc('16','20')
gen_lc('14_25','18-29')


#%% barycorr
os.chdir(f'products/sa_data_lc_{binsize}sec')

orbitfile=glob(xte_obs._orbit_path_+'/*')
if len(orbitfile)==1:
    orbitfile=orbitfile[0]
else:
    raise Exception('more than one orbit files found!!!')


for channel in glob('*lc'):
    bcorr=f"barycorr infile={channel} outfile={channel}_bary orbitfiles={orbitfile} \n"
    #run_command(bcorr,out_path='./',dull=0)
    print(bcorr)

#%% orb corr
for channel in glob('*lc_bary'):
    import Misc
    correct_times(f'{channel}',Misc.doppler_correction.orb_params_v0332)










# #%% combine channels

# lcmath=f'lcmath infile=ch11.lc_bary_orb_corr bgfile=ch12.lc_bary_orb_corr outfile=ch11_and_ch12_tmp.lc_bary_orb_corr multi=0.00087130476998205 multb=0.0005685587704151747 addsubr = yes err_mode=2'

# run_command(lcmath,0,0)

# lcmath=f'lcmath infile=ch11_and_ch12_tmp.lc_bary_orb_corr bgfile=ch13.lc_bary_orb_corr outfile=ch11_and_ch12_and_ch13.lc_bary_orb_corr multi=1 multb=0.0005529993106345093 addsubr = yes err_mode=2'

# run_command(lcmath,0,0)


# lcmath=f'lcmath infile=ch14.lc_bary_orb_corr bgfile=ch15.lc_bary_orb_corr outfile=ch14_and_ch15_tmp.lc_bary_orb_corr multi=0.0005456596986510945 multb=0.0005415316485591997 addsubr = yes err_mode=2'

# run_command(lcmath,0,0)


# lcmath=f'lcmath infile=ch14_and_ch15_tmp.lc_bary_orb_corr bgfile=ch16.lc_bary_orb_corr outfile=ch14_and_ch15_and_ch16.lc_bary_orb_corr multi=1 multb=0.0005357527493898891 addsubr = yes err_mode=2'
# run_command(lcmath,0,0)




# #%% lcmath

# fraction=0.9
# create_dir('fe_line')

# lcmath=f'lcmath infile=ch11_and_ch12_and_ch13.lc_bary_orb_corr bgfile=ch14_and_ch15_and_ch16.lc_bary_orb_corr  outfile=fe_line/fe_line_frac{fraction}.lc_bary_orb_corr multi=1 multb={fraction} addsubr = no'

# run_command(lcmath,0,0)




# fraction_tmp=1
# lcmath=f'lcmath infile=ch11_and_ch12_tmp.lc_bary_orb_corr bgfile=ch14_and_ch15_tmp.lc_bary_orb_corr  outfile=fe_line/fe_line_frac{fraction_tmp}_ch1112_min_14_15.lc_bary_orb_corr multi=1 multb={fraction_tmp} addsubr = no'

# run_command(lcmath,0,0)



# # fraction_try=1.2

# # lcmath=f'lcmath infile=ch11_12_13.lc_bary_orb_corr bgfile=ch14_15.lc_bary_orb_corr  outfile=fe_line/fe_line_frac{fraction_try}_try.lc_bary_orb_corr multi=1 multb={fraction_try} addsubr = no'

# # run_command(lcmath,0,0)




# #%% plot ccfs


# def plot_ccf(filepath,ax):
#     ccf=np.genfromtxt(filepath,skip_header=3)
#     N=int(ccf[:,0].shape[0]/2)
#     norm=np.max(ccf[:,2])
#     ax.errorbar(ccf[:,0],ccf[:,2]/norm,ccf[:,3]/norm,drawstyle='steps-mid',label='data')
#     #ax.errorbar(-ccf[N:,0],ccf[N:,2],ccf[N:,3],alpha=0.5,color='r',drawstyle='steps-mid')
#     ax.errorbar(ccf[N:,0],ccf[0:N+1:,2][::-1]/norm,ccf[0:N+1:,3][::-1]/norm,alpha=0.5,color='m',drawstyle='steps-mid',label='data (neg delay)')
#     ax.set_xlabel('Delay, s')

#     ax.set_ylabel('CCF')
#     fig.tight_layout()
#     sns.despine(fig,top=1,right=0)
#     plt.show()

# fig,ax_ccf=plt.subplots(figsize=(16,6))

# plot_ccf(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}/products/sa_data_lc_{binsize}sec/ccf_{fraction}.qdp',ax_ccf)

# ax_ccf.set_xlabel('iron delay, s')
# ax_ccf.legend()
# ax_ccf.set_title(f'\n\n ccf_{fraction}')

# ax_ccf.set_xlim(-1,15)
# ax_ccf.set_ylim(-0.5,1)

# plt.show()
