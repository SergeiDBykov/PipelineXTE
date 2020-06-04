#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:06:23 2020

@author: s.bykov
"""


from pipeline_core import *

ObsID='90089-11-03-00G' # 90089-11-03-00G or 90089-11-03-01G 90089-11-03-02
xte_obs=ObservationXTE(ObsID)
os.chdir(f'/Users/s.bykov/work/xray_pulsars/rxte/results/out{ObsID}')


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



# #%% combine 10 and 11 channels, and 12 with 13

# lcmath=f'lcmath infile=ch11.lc_bary_orb_corr bgfile=ch10.lc_bary_orb_corr outfile=ch10_and_ch11.lc_bary_orb_corr multi={0.00028} multb={0.00087} addsubr = yes err_mode=2'

# run_command(lcmath,0,0)



# lcmath=f'lcmath infile=ch12.lc_bary_orb_corr bgfile=ch13.lc_bary_orb_corr outfile=ch12_and_ch13.lc_bary_orb_corr multi={0.00027} multb={0.00027} addsubr = yes err_mode=2'

# run_command(lcmath,0,0)


# #%% lcmath

# fraction=0.9
# create_dir('fe_line')

# lcmath=f'lcmath infile=ch10_and_ch11.lc_bary_orb_corr bgfile=ch12_and_ch13.lc_bary_orb_corr outfile=fe_line/fe_line_frac{fraction}.lc_bary_orb_corr multi=1 multb={fraction} addsubr = no'

# run_command(lcmath,0,0)



