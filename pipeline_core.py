#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 16:08:29 2019

@author: s.bykov
"""

#%% imports and definitions
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from glob import glob
import os
from scipy.optimize import curve_fit
import shutil
from Misc.TimeSeries.cross_corr import my_crosscorr
from Misc.plot_spectra import Spectra
from Misc.doppler_correction import correct_times
from subprocess import call

def open_dir():
    call(['open',os.getcwd()])

def open_dir_in_term():
    import appscript
    appscript.app('Terminal').do_script(f'cd {os.getcwd()}')


plt.ioff()
import seaborn as sns
sns.set(style='ticks', palette='deep',context='notebook')



RXTE_path='/Users/s.bykov/work/xray_pulsars/rxte/'
xspec_scripts_path='/Users/s.bykov/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/'

#for beautiful plot in ph_res_results
results_path='/Users/s.bykov/work/xray_pulsars/rxte/plots_results/pandas_data/'
filename='standard_pipeline'
ObsParams=pd.read_pickle(results_path+f'{filename}.pkl')
ObsParams_plot=ObsParams.sort_values('MJD_START')

ObsParams_plot.period_orb_corr= ObsParams.period_orb_corr.replace(to_replace='None',value=np.nan)
ObsParams_plot.period_orb_corr_err= ObsParams.period_orb_corr_err.replace(to_replace='None',value=np.nan)


strange_dets_ID=['90427-01-04-00', '90427-01-04-01', '90427-01-04-04']

def sum_error(a,b,da,db):
    f=a+b
    sigma=np.sqrt(da**2+db**2)
    return f,sigma

def ratio_error(a,b,da,db):
    f=a/b
    sigma=np.abs(f)*np.sqrt( (da/a)**2 + (db/b)**2  )
    return f, sigma

def pulsed_fraction_error(I,I_err):
    arg_max=np.argmax(I)
    arg_min=np.argmin(I)

    Imax=I[arg_max]
    Imax_err=I_err[arg_max]

    Imin=I[arg_min]
    Imin_err=I_err[arg_min]

    tmp_err=sum_error(Imax,Imin,Imax_err,Imin_err)[1]

    numenator=Imax-Imin
    denomenator=Imax+Imin

    numenator_err=tmp_err
    denomenator_err=tmp_err

    return ratio_error(numenator,denomenator,numenator_err,denomenator_err)


def gauss(t,t0,sigma,N):
    return N*np.exp(-(t-t0)**2/(2*sigma**2))/np.sqrt(sigma)

def run_command(cmd,out_path,dull=0):
    print('Running command: ', cmd)
    if dull==0:
        os.system(cmd + '| tee -a {0}/pipeline_log.txt'.format(out_path))
    else:
        print('dull run')
    return 1


def select_pcu(filterfile_name,ObsID,plot=1):
    '''
    this thing choses which pcu has median value of  pcu_on
    equal to 1
    it checks only times then num_pcu_on=the mode of num_pcu_on

    if Obs is in obsid when dets are jumping, return values that i found by eye

    '''
    if ObsID=='90427-01-04-04':
        return np.array([0, 2]),'0,2'
    elif ObsID=='90427-01-04-00':
        return np.array([0, 2, 3]),'0,2,3'
    elif ObsID=='90427-01-04-01':
        return np.array([0, 2]),'0,2'
    else:
        pass

    xfl=fits.open(filterfile_name)

    pcu_on=xfl[1].data['num_pcu_on']
    #set num of pcu as the mode of pcu_on array (most common configuration)
    num_pcu_on=np.bincount(pcu_on[pcu_on>0]).argmax()

    pcu_0_on=xfl[1].data['PCU0_ON']
    pcu_1_on=xfl[1].data['PCU1_ON']
    pcu_2_on=xfl[1].data['PCU2_ON']
    pcu_3_on=xfl[1].data['PCU3_ON']
    pcu_4_on=xfl[1].data['PCU4_ON']

    indeces=pcu_on==num_pcu_on

    working_pcu=[bool(np.mean(x)) for x in [pcu_0_on[indeces],pcu_1_on[indeces],pcu_2_on[indeces],pcu_3_on[indeces],pcu_4_on[indeces]]]
    working_pcu=np.array(working_pcu)
    working_pcu=np.where(working_pcu==True)[0]
    str_working_pcu=str(working_pcu).replace(' ',',').replace('[','').replace(']','')
    if len(working_pcu)!=num_pcu_on:
        print(working_pcu,num_pcu_on)
        raise Exception('number of working pcu is not equal to outpul list length')
    if num_pcu_on==0:
        raise Exception('number of working pcu==0')

    if plot:
        def plot_dets(xfl,num_pcu_on):
            time=xfl[1].data['time']
            fig,axs=plt.subplots(6, 1)
            axs[0].plot(time,xfl[1].data['NUM_PCU_ON'],label='total')
            axs[0].set_ylim(0,6)
            axs[0].axhline(num_pcu_on,color='r',lw=1,ls=':',label='num_pcu_on filter value')
            axs[0].legend()
            for i in range(0,5):
                name='PCU{0}_ON'.format(i)
                axs[i+1].plot(time,xfl[1].data[name],label=name)
                axs[i+1].legend()
                axs[i+1].set_ylim(0,1.5)
            axs[0].set_title(working_pcu)
            fig.savefig(f'dets.png')
        plot_dets(xfl,num_pcu_on)


    return working_pcu,str_working_pcu


def find_max_efsearch_data(efsearcf_fits_file):
    efsearch=fits.open(efsearcf_fits_file)
    period=efsearch[1].data['period']
    chisq=efsearch[1].data['chisqrd1']
    bestper=period[np.argmax(chisq)]
    print(bestper)
    return bestper


def fit_efsearch_data(efsearcf_fits_file,ObsID,ax=None,savefig=1,fit_data=1):
    efsearch=fits.open(efsearcf_fits_file)
    period=efsearch[1].data['period']
    chisq=efsearch[1].data['chisqrd1']
    p0=[period[np.argmax(chisq)],0.001,max(chisq)*np.sqrt(0.001)]
    if ax==None:
        fig,ax=plt.subplots()
    else:
        ax=ax

    if fit_data==False:
        ax.plot(period,chisq)
        ax.set_title(ObsID)
        ax.set_xlabel('Period')
        ax.set_ylabel('\chi^2')
        return None,None
    else:

        try:
            popt,perr=curve_fit(gauss,period,chisq,p0=p0)
            perr=np.sqrt(np.diag(perr))
            ax.plot(period,gauss(period,*p0),'r:',alpha=0.5)
            ax.plot(period,chisq)
            ax.plot(period,gauss(period,*popt),'k:')
            ax.set_title('Period='+str(popt[0])+'  sigma='+str(popt[1])+'\n'+ObsID)
            ax.set_xlabel('Period')
            ax.set_ylabel('\chi^2')
            if savefig:
                fig.savefig(f'efsearch_res_{ObsID}.png')
                plt.close(fig)
            return popt[0],perr[0]
        except:
            ax.plot(period,chisq)
            ax.plot(period,gauss(period,*p0),'r-.')
            ax.set_title('failed fit'+'\n'+ObsID)
            ax.set_xlabel('Period')
            ax.set_ylabel('\chi^2')
            if savefig:
                fig.savefig(f'efsearch_res_{ObsID}.png')
                plt.close(fig)
            return None,None




def create_dir(dir):
    os.system(f'mkdir -p {dir}')


class ObservationXTE():
    def __init__(self,ObsID):
        '''
        we load obs name and write pathes to self
        create if necessary path for results of analysis

        '''
        print('###')
        print(f'Observation {ObsID} loaded')
        self.ObsID=ObsID
        self._data_path_=RXTE_path+'data/AO9/'+ObsID+'/'
        self._pca_path_=self._data_path_+'pca'
        self._orbit_path_=self._data_path_+'orbit'
        os.chdir(RXTE_path+'results/')


        create_dir('out'+self.ObsID)
        os.chdir('out'+ObsID)
        out_path=os.getcwd()
        self.out_path=out_path
        create_dir('products')
        os.chdir('products')
        for folder in ['pcu2_top',
                       'deadtime_corr','std1_lc','fasebin',
                       'bkg_all','fasebin_spe']:
            create_dir(folder)


        os.chdir(self.out_path)
        create_dir('obs_info')
        os.chdir('obs_info')
        self.obs_info_file=open(f"obs_info_{ObsID}.txt","a+")
        self.obs_info_file.close()
        self.obs_info_file=os.path.abspath(f'obs_info_{ObsID}.txt')

        self.spe_info_file=open(f"spe_info_{ObsID}.txt","a+")
        self.spe_info_file.close()
        self.spe_info_file=os.path.abspath(f'spe_info_{ObsID}.txt')

        self.per_info_file=open(f"per_info_{ObsID}.txt","a+")
        self.per_info_file.close()
        self.per_info_file=os.path.abspath(f'per_info_{ObsID}.txt')

        self.fasebin_info_file=open(f"fasebin_info_{ObsID}.txt","a+")
        self.fasebin_info_file.close()
        self.fasebin_info_file=os.path.abspath(f'fasebin_info_{ObsID}.txt')

        self.config_info_file=open(f"config_info_{ObsID}.txt","a+")
        self.config_info_file.close()
        self.config_info_file=os.path.abspath(f'config_info_{ObsID}.txt')



        os.chdir(self.out_path)



    def write_to_obs_info(self,file,name,value):
        '''
        store variables in file
        mode is a+ -append
        save name of variable and its value, for instance
        period 4.37
        '''
        with open(file,'a+') as f:
            f.write(name+'\t'+str(value)+'\n')
            f.close()




    def pandas_series(self,read_obs_info_file=True):
        '''
        lots of those series in a dataframe give an opportunity to
        choose observation for general analysis (period finding, spectral approximation, etc)
        and for phase resolved spectroscopy, and to print datamodes for sa/se data
        '''
        obs_info={}
        if read_obs_info_file:
            for filepath in [self.obs_info_file,self.spe_info_file,
                             self.per_info_file,self.fasebin_info_file,self.config_info_file]:
                if os.stat(filepath).st_size==1: #strange bug: some files are not empty but host 1 carriage return, hence the error. i empty such files by hand.
                    print('1 byte filesize')
                    os.system(f'>{filepath}')
                with open(filepath) as f:
                    for line in f:
                        (key,val) = line.split()            #in case of normal entry as 'chi2 23' or 'se_config E123um'
                        try:
                            obs_info[key] = float(val)
                        except:
                            obs_info[key] = str(val)
                    f.close()

        Series=pd.Series(obs_info,name='obs'+self.ObsID)
        return Series


    def get_configs(self,rewrite=1):
        '''
        it creates a configuration dictionary for NOT standart or housekeeping or not SA\SE modes.
        I am trying to create self.xdf_filelist,  the filelist which is to be analysed with standart pipeline,
        based on some critetion for SE and SA.
        For SA: num of en channels >10, and lower energy is less than ch. 5
        for SE: lower energy<chnnel 5.
        I am checking if all configs of SA or SE data are equal

        '''
        def check_list(x):
            '''
            checks if all elements in an array are equal
            '''
            if x==[]:
                return True
            else:
                return x.count(x[0]) == len(x)
        self.config_dict={}
        data_files=glob(self._pca_path_+'/*.gz')

        self.sa_exists=False
        self.se_exists=False
        self.se_config=[]
        self.sa_config=[]
        self.std1_config='None'
        self.std2_config='None'

        self._xdf_list_std1_=[]
        self._xdf_list_std2_=[]
        self._xdf_list_se_=[]
        self._xdf_list_sa_=[]

        for filepath in data_files:
            with fits.open(filepath) as fits_file:
                data_type=fits_file[1].name
                data_mode=fits_file[1].header['datamode']
                if 'Good_Xenon' in data_mode:
                    raise Exception('good xenon datamode detected!')
                if data_type!='XTE_SA' and data_type!='XTE_SE':
                    #exit if it is, for example, housekeeping files
                    fits_file.close()
                else:
                    tddes2=fits_file[1].header['tddes2']
                    try:
                        #because sa do not have tevtb2 header
                        tevtb2=fits_file[1].header['TEVTB2']
                    except:
                        tevtb2=None
                    #save to self.config
                    self.config_dict[filepath.split('/')[-1]]=[data_mode,data_type,tddes2,tevtb2]

                    #the next section checks whether std1,std2,se,sa formats are presented and available
                    #se and sa (not standart mode but binned or event mode) are subject of abovementionad criteria
                    if data_mode=='Standard1b':
                        self._xdf_list_std1_.append(filepath)
                        self.std1_config=True
                    elif data_mode=='Standard2f':
                        self._xdf_list_std2_.append(filepath)
                        self.std2_config=True
                    #analyse keywords to return if enough energy channels and binning etc
                    if data_type=='XTE_SE':
                        D,E,C=tddes2.replace(' ','').split('&')
                        #channels
                        C=C[2:]
                        C=C.replace(']','')
                        ch_low,ch_hi=C.split('~')
                        ch_low=int(ch_low)
                        ch_hi=int(ch_hi)
                        if ch_low<=5:
                            self._xdf_list_se_.append(filepath)
                            self.se_exists=True
                            self.se_config.append(data_mode)
                    if data_type=='XTE_SA' and data_mode!='Standard1b' and data_mode!='Standard2f':
                        D,E,C,_=tddes2.replace(' ','').split('&')
                        #channels
                        C=C[2:]
                        C=C.replace(']','')
                        ch_num=len(C.split(','))
                        ch_low=int(C[0])
                        ch_hi=int(C[-1])
                        if ch_num>10 and ch_low<=5:
                            self._xdf_list_sa_.append(filepath)
                            self.sa_exists=True
                            self.sa_config.append(data_mode)
                    fits_file.close()
        #checking if sa/se configs are the same. if they are, then assign sa/se config as the first element of an array

        if check_list(self.se_config):
            try:
                self.se_config=self.se_config[0]
            except:
                self.se_config=None
        else:
            print(self.se_config)
            print(self._xdf_list_se_)
            #self._xdf_list_se_=[]
            self.se_config='inhomogeneity'
            print('different configs in SE data')
        if check_list(self.sa_config):
            try:
                self.sa_config=self.sa_config[0]
            except:
                self.sa_config=None
        else:
            print(self.sa_config)
            print(self._xdf_list_sa_)
            #self._xdf_list_sa_=[]
            self.sa_config='inhomogeneity'
            print('different configs in SA data')


        #write configs to a file
        if rewrite:
            fname=self.config_info_file
            os.system(f'> {fname}') #remove content from per_info
        names=['ObsID','SA_config','SE_config','STD1_config','STD2_config']
        values=[self.ObsID,self.sa_config,self.se_config,self.std1_config,self.std2_config]
        for name,val in zip(names,values):
            self.write_to_obs_info(self.config_info_file,name,val)


    #%% ========= filters =========

    def filter_time(self):
        '''
        easy: create filter files for observations
        '''
        appidlist_path='/Users/s.bykov/work/xray_pulsars/rxte/appidlist'
        os.chdir(RXTE_path+'results/')
        out_path=self.out_path
        os.chdir(out_path)

        #path_to_fmi=self._data_path_+'/FMI'
        path_to_fmi=self._data_path_
        xtefilt=' xtefilt -c -a {0} -o {1} -p {2} -t 0.1 -f FILTERFILE'.format(appidlist_path,self.ObsID,path_to_fmi)
        run_command(xtefilt,out_path)


        working_pcu,_=select_pcu(self.out_path+'/FILTERFILE.xfl',self.ObsID)
        num_pcu_on=str(len(working_pcu))
        if self.ObsID=='90427-01-04-04':
            #return np.array([0,  2]),'0,2'
            maketime='maketime FILTERFILE.xfl maketime_file.fits "elv.gt.10.and.offset.lt.0.02.and.(time_since_saa.lt.0.or.time_since_saa.gt.30).and.PCU0_ON.eq.1.and.PCU2_ON.eq.1" compact=no time=TIME'
        elif self.ObsID=='90427-01-04-00':
            #return np.array([0, 2, 3]),'0,2,3'
            maketime='maketime FILTERFILE.xfl maketime_file.fits "elv.gt.10.and.offset.lt.0.02.and.(time_since_saa.lt.0.or.time_since_saa.gt.30).and.PCU0_ON.eq.1.and.PCU2_ON.eq.1.and.PCU3_ON.eq.1" compact=no time=TIME'

        elif self.ObsID=='90427-01-04-01':
            #return np.array([0, 2]),'0,2'
            maketime='maketime FILTERFILE.xfl maketime_file.fits "elv.gt.10.and.offset.lt.0.02.and.(time_since_saa.lt.0.or.time_since_saa.gt.30).and.PCU0_ON.eq.1.and.PCU2_ON.eq.1" compact=no time=TIME'


        else:
            maketime='maketime FILTERFILE.xfl maketime_file.fits "elv.gt.10.and.offset.lt.0.02.and.(time_since_saa.lt.0.or.time_since_saa.gt.30).and.num_pcu_on.eq.{0}" compact=no time=TIME'.format(num_pcu_on)

        run_command(maketime,out_path)


        '''
        easy: create filter files for observations
        uses only PCU2_ON
        RUn only after filter_time(self) created maketime.fits
        '''
        os.chdir(RXTE_path+'results/')
        out_path=self.out_path
        os.chdir(out_path)

        maketime='maketime FILTERFILE.xfl maketime_file_pcu2.fits "elv.gt.10.and.offset.lt.0.02.and.(time_since_saa.lt.0.or.time_since_saa.gt.30).and.PCU2_ON.eq.1" compact=no time=TIME'

        run_command(maketime,out_path)


    def make_file_lists(self):
        '''
        creating filelists as does xdf routine

        '''
        out_path=self.out_path+'/'

        for filename,filelist in zip(['std1.xdf','std2.xdf','se.xdf','sa.xdf'],
                                     [self._xdf_list_std1_,self._xdf_list_std2_,self._xdf_list_se_,self._xdf_list_sa_]):
            if len(filelist)==0:
                pass
            else:
                print('creating a list of files {0}'.format(out_path+filename))
                with open(out_path+filename, 'w') as f:
                    for item in filelist:
                        f.write("{0}\n".format(item[0:len(item)-3])) #deleting .gz part of a string
                    f.close()
        #this part for SE of SA config for fasebin
        xdflist=glob('*.xdf')
        for xdf in ['std1.xdf','std2.xdf']:
            xdflist.remove(xdf)
        if len(xdflist)==0:
            fasebin_cfg='None'
        elif len(xdflist)==1:
            fasebin_cfg=xdflist[0][0:2] #[0:2] because of the extension
        else:
            fasebin_cfg='se' #in case both present
        self.write_to_obs_info(self.obs_info_file,'fasebin_cfg',fasebin_cfg)


    def recreate_se_xdf(self,bitmask='/Users/s.bykov/work/xray_pulsars/rxte/dets_2.bit'):
        '''
        creates new SE files when only pcu2 is left. Returns True if new se.xdf is created and False if no se data present
        '''
        os.chdir(self.out_path)
        se_xdf_file='se.xdf'
        if not os.path.exists(se_xdf_file):
            print('no se xdf file exists')
            return False
        elif os.path.exists('se.xdf_orig'):
            print('already done')
            return False
        try:
            os.system('mkdir se_pcu2')
            for k,filename in enumerate(open(se_xdf_file).readlines(),1):
                filename=filename[:-1]+'.gz'
                print(filename)
                cmd=f'fselect infile={filename} outfile=./se_pcu2/se_{k}.fits expr=@{bitmask}'
                os.system(cmd)
            os.chdir('se_pcu2')
            os.system('rm -f se_pcu2.xdf')
            abspath_se=[os.path.abspath(x) for x in glob("*.fits")]
            with open('se_pcu2.xdf','w+') as f:
                for line in abspath_se:
                    f.write(line+'\n')
            os.chdir('../')
            os.system('mv se.xdf se.xdf_orig')
            os.system('cp ./se_pcu2/se_pcu2.xdf ./se.xdf')
            return True
        except:
            print('se_pcu folder exists, skip recreate_se_xdf function')
            return True





    def filter_and_files(self):
       self.get_configs()
       self.make_file_lists()
       self.recreate_se_xdf()
       self.filter_time()



    #%%=========== spectra =========


    def make_spectrum_and_rsp_std2(self):
        '''
        self-explanatory: spe+response+backgrnd from std2 data.
        mode - 'pcu2': using top layers of pcu2
        '''
        os.chdir(self.out_path)

        datalist='std2'
        binsize='16'

        layers='LR1'
        detectors='2'

        with open(self.out_path+'/pcu2_top.col','w' ) as f:
            f.write(f'X1LSpecPcu2\nX1RSpecPcu2\n')
            f.close()


        extract=f'saextrct infile=@{datalist}.xdf gtiorfile=APPLY gtiandfile=maketime_file_pcu2.fits outroot=products/pcu2_top/{datalist}_{binsize}s timecol=TIME columns=@pcu2_top.col accumulate=ONE binsz={binsize} printmode=SPECTRUM lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint=INDEF  chmin = INDEF chmax=INDEF chint=INDEF chbin=INDEF'
        run_command(extract,out_path=self.out_path)

        os.chdir(self.out_path+'/products/pcu2_top')
        pcarsp=f'pcarsp -f {datalist}_{binsize}s.pha -a ../../FILTERFILE.xfl -l {layers} -j y -p {detectors} -m y -n {datalist}_{binsize}s.rsp'

        run_command(pcarsp,self.out_path)




    #%% ========= background =========
    def calc_bkg_files(self,datalist='std2',binsize='16'):
        '''
        run pcabackest to create baclground filelist
        moronic thing: one can not use arguments for running (run)pcabackest - i need to create
        temporary file with name, say, temp_list.xdf, which may be std2 or std1 or whatnever

        '''
        bkgmodel='/Users/s.bykov/work/xray_pulsars/rxte/pca_bkgd_cmbrightvle_eMv20051128.mdl'
        saahist='/Users/s.bykov/work/xray_pulsars/rxte/pca_saa_history.gz'
        os.chdir(self.out_path)

        try:
            os.system(f'cp FILTERFILE.xfl ./products/bkg_all')   #this is needed because of background files renaming
        except:
            pass
        os.system(f'cp {datalist}.xdf ./products/bkg_all/templist.xdf')

        os.chdir('./products/bkg_all')

        pcabackest=f'runpcabackest infile=@{datalist}.xdf  outsuffix=back outlist="bkg.xdf" modelfile={bkgmodel} filterfile=FILTERFILE.xfl  interval={binsize} layers=Yes gai ncorr=no fullspec=yes saahfile={saahist}'
        run_command(pcabackest,out_path=self.out_path)


        prefix=os.getcwd()+'/'

        f = open("outlist")                  #ultra moronic thing: pcabackest does not write full path to a bkg files
        o = open(f"{datalist}_bkg.xdf","a")
        while 1:
          line = f.readline()
          if not line: break
          line = prefix+line
          o.write(line)
        o.close()


        #os.system(f'cp outlist {datalist}_bkg.xdf')
        os.system(f'cp {datalist}_bkg.xdf ../../')
        os.system('rm outlist')
        os.system('rm templist.xdf')


    def calc_bkg_spe(self):
        binsize='16'
        datalist='std2'
        os.chdir(self.out_path)


        extract=f'saextrct infile=@{datalist}_bkg.xdf gtiorfile=APPLY gtiandfile=maketime_file_pcu2.fits outroot=products/pcu2_top/{datalist}_{binsize}s_bkg timecol=TIME columns=@pcu2_top.col accumulate=ONE binsz={binsize} printmode=SPECTRUM lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint=INDEF  chmin = INDEF chmax=INDEF chint=INDEF chbin=INDEF'
        run_command(extract,out_path=self.out_path)

    #%% =========  deadtime =========
    def calc_deadtime(self):
        '''
        calculate deadtime fraction with formula
        update exposure of source and bkg spectrum
        '''
        os.chdir(self.out_path)
        datalist='std1'
        binsize='1'
        deadtime1_cols='/Users/s.bykov/work/xray_pulsars/rxte/deadtime_1_cols.col'
        extract=f'saextrct infile=@{datalist}.xdf gtiorfile=APPLY gtiandfile=maketime_file.fits outroot=products/deadtime_corr/{datalist}_{binsize}s_deadtime1 timecol=TIME columns=@{deadtime1_cols} accumulate=ONE binsz={binsize} printmode=SPECTRUM lcmode=RATE spmode=RATE timemin=INDEF timemax=INDEF timeint=INDEF  chmin = INDEF chmax=INDEF chint=INDEF chbin=INDEF'
        run_command(extract,out_path=self.out_path,dull=0)

        deadtime2_cols='/Users/s.bykov/work/xray_pulsars/rxte/deadtime_2_cols.col'
        extract=f'saextrct infile=@{datalist}.xdf gtiorfile=APPLY gtiandfile=maketime_file.fits outroot=products/deadtime_corr/{datalist}_{binsize}s_deadtime2 timecol=TIME columns=@{deadtime2_cols} accumulate=ONE binsz={binsize} printmode=SPECTRUM lcmode=RATE spmode=RATE timemin=INDEF timemax=INDEF timeint=INDEF  chmin = INDEF chmax=INDEF chint=INDEF chbin=INDEF'
        run_command(extract,out_path=self.out_path,dull=0)

        C1=fits.open('./products/deadtime_corr/std1_1s_deadtime1.pha')[1].data['rate'][0]
        C2=fits.open('./products/deadtime_corr/std1_1s_deadtime2.pha')[1].data['rate'][0]
        working_pcu,_=select_pcu(self.out_path+'/FILTERFILE.xfl',self.ObsID)
        num_pcu_on=len(working_pcu)
        N=num_pcu_on

        DTF=C1*1e-5/N+C2*1.5*1e-4/N
        DCOR=1/(1-DTF)
        self.DCOR=DCOR
        self.write_to_obs_info(self.obs_info_file,'dead_time_corr',DCOR)

    def update_exposure_for_deadtime(self,path='./products/pcu2_top/std2_16s_sys.pi'):
        os.chdir(self.out_path)
        #DCOR=self.DCOR
        #DCOR=self.get_from_obs_info(self.obs_info_file,'dead_time_corr')
        DCOR=self.pandas_series()['dead_time_corr']
        with fits.open(path,mode='update') as f:
            if 'i changed exposure from' in f[1].header.comments['exposure']:
                raise Exception('correction has alrady been applied')
            else:
                exp=f[1].header['exposure']
                f[1].header['exposure']=exp/DCOR
                f[1].header.comments['exposure']=f'i changed exposure from {exp} to {exp}/{DCOR} (deadtime correction)'
                f.close()



    #%% ========= spectral fitting =============

    def fit_std2_spe(self,model='cutoffpl',spe_path='/products/pcu2_top',
                     chmin=0,chmax=8,error=0.01):

        def apply_systematic_error(self,spe_file='./products/pcu2_top/std2_16s.pha'):
            os.chdir(self.out_path)
            out_file=spe_file.replace('.pha','_sys.pi')
            grp=f'''grppha infile="{spe_file}" outfile="{out_file}"  clobber=yes comm="systematics {chmin}-{chmax} {error} & exit" '''
            run_command(grp,out_path=self.out_path)
            self.update_exposure_for_deadtime()


        apply_systematic_error(self)

        os.chdir(self.out_path+spe_path)
        ObsID=self.ObsID

        os.chdir(self.out_path+'/obs_info')
        #os.system('> {0}'.format(self.spe_info_file)) #remove content from spe_info
        os.system(f'''sed  '/^{model}_/ d' spe_info_{ObsID}.txt > stripped.txt ''')
        #remove lines that start with model name, for instance po_chi, po_tot_err etc
        os.system(f'rm spe_info_{ObsID}.txt')
        os.system(f'mv stripped.txt ./spe_info_{ObsID}.txt')

        os.chdir(self.out_path+spe_path)
        for filePath in ['cutoffpl.xcm','pohi.xcm','mean_spe.dat',f'mean_sp_{ObsID}.ps',f'mean_sp.ps','mean_sp_pohi.ps']:
            if os.path.exists(os.getcwd()+'/'+filePath):
                os.remove(os.getcwd()+'/'+filePath)

        if os.path.exists(model):
            shutil.rmtree(model)
            os.system(f'mkdir -p {model}')
        else:
            os.system(f'mkdir -p {model}')

        os.system(f'xspec - {xspec_scripts_path}{model}.txt')

        os.chdir(model)
        #write spectral data to obs_info
        sp_data=np.genfromtxt('mean_spe.dat')
        sp_data_pars=['chi2','dof','eqw','eqw_lo','eqw_hi',
                      'tot_flux','tot_flux_lo','tot_flux_hi',
                      'fe_flux','fe_flux_lo','fe_flux_hi',
                      'po','efold','ecut','eline','norm_line']
        sp_data_pars=[model+'_'+par for par in sp_data_pars]
        for name,val in zip(sp_data_pars,sp_data):
            self.write_to_obs_info(self.spe_info_file,name,val)

        os.system(f'cp mean_sp_{model}.ps mean_sp_{model}_{ObsID}.ps')
        os.remove(f'mean_sp_{model}.ps')
        os.chdir(self.out_path+spe_path)
        ff=fits.open(f'std2_16s.pha')
        expo=ff[1].header['EXPOSURE']
        self.write_to_obs_info(self.spe_info_file,'EXPOSURE',expo)

        mjd_start=ff[1].header['MJDREFI']+ff[1].header['MJDREFF']+(ff[1].header['TSTART'])/86400
        self.write_to_obs_info(self.obs_info_file,'MJD_START',mjd_start)
        ff.close()


    #%% ========= period finding =============

    def make_std1_lc(self):
        '''
        self-explanatory: lc+ barycentric correction
        '''
        datalist='std1'
        binsize='0.1'

        os.chdir(self.out_path)


        extract=f'saextrct infile=@{datalist}.xdf gtiorfile=APPLY gtiandfile=maketime_file.fits outroot=products/std1_lc/{datalist}_{binsize}s timecol=TIME columns=GOOD accumulate=ONE binsz={binsize} printmode=LIGHCURVE lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint=INDEF  chmin = INDEF chmax=INDEF chint=INDEF chbin=INDEF'


        orbitfile=glob(self._orbit_path_+'/*')
        if len(orbitfile)==1:
            orbitfile=orbitfile[0]
        else:
            raise Exception('more than one orbit files found!!!')

        faxbary=f'faxbary infile=products/std1_lc/{datalist}_{binsize}s.lc outfile=products/std1_lc/{datalist}_{binsize}s_bary.lc orbitfiles={orbitfile} ra=5.37500000E+01  dec=5.31730003E+01 barytime=no'
        run_command(extract,out_path=self.out_path,dull=0)
        run_command(faxbary,out_path=self.out_path,dull=0)


    def std1_lc_orb_corr(self):
        os.chdir(self.out_path)
        os.chdir('products/std1_lc/')
        import Misc
        correct_times('std1_0.1s_bary.lc',Misc.doppler_correction.orb_params_v0332)
        os.system('mv std1_0.1s_bary.lc ./std1_0.1s_bary.lc_original')
        os.system('cp std1_0.1s_bary.lc_orb_corr ./std1_0.1s_bary.lc')

    def make_efsearch(self,lcfile='std1_0.1s_bary',p0='4.3745',rewrite=1,nper='64'):
        '''
        this task runs efsearch with some parameters, saves the file
        my period is between 3.373 and 3.377
        nper is a number of periods to check

        '''
        os.chdir(self.out_path+'/products/std1_lc')
        if rewrite:
            fname=self.per_info_file
            os.system(f'> {fname}') #remove content from per_info
            for filePath in [f'{lcfile}.efs',f'{lcfile}_init.efs']:
                if os.path.exists(filePath):
                    os.remove(filePath)

        efsearch=f'efsearch cfile1="{lcfile}.lc" dper={p0} nphase=32 dres=0.0001 nper=25 outfile="{lcfile}_init.efs" window="-" sepoch=INDEF nbint=INDEF plot=no'
        run_command(efsearch,out_path=self.out_path)

        p0=find_max_efsearch_data('std1_0.1s_bary_init.efs')

        efsearch=f'efsearch cfile1="{lcfile}.lc" dper={p0} nphase=16 dres=0.00005 nper={nper} outfile="{lcfile}.efs" window="-" sepoch=INDEF nbint=INDEF plot=no'
        run_command(efsearch,out_path=self.out_path)



    def find_period(self,lcfile='std1_0.1s_bary',rewrite=1):
        '''
        efsearch files are fitted with gaussian.
        it also writes period and  period sigma to per info file

        '''
        os.chdir(self.out_path+'/products/std1_lc')
        if rewrite:
            fname=self.per_info_file
            os.system(f'> {fname}') #remove content from per_info
            for filePath in ['efsearch_res.png',f'efsearch_res_{self.ObsID}.png']:
                if os.path.exists(filePath):
                    os.remove(filePath)
        per,err=fit_efsearch_data(f'{lcfile}.efs',self.ObsID)
        self.write_to_obs_info(self.per_info_file,'period_orb_corr',per)
        self.write_to_obs_info(self.per_info_file,'period_orb_corr_err',err)




    #%% ========= phase resolved spectroscopy =============

    def make_spe(self):
        os.chdir(self.out_path)


        binsize='16'

        datalist=self.pandas_series()['fasebin_cfg']#self.get_from_obs_info(self.obs_info_file,'fasebin_cfg')
        if datalist=='None':
            raise Exception('No SE and SA: fasebin impossible')
        else:
            pass


        if datalist=='sa':
            detectors_bool,detectors=select_pcu(self.out_path+'/FILTERFILE.xfl',self.ObsID,plot=0)
            layers='L1,R1,L2,R2,L3,R3'

            extract=f'saextrct infile=@{datalist}.xdf gtiorfile=APPLY gtiandfile=maketime_file.fits outroot=products/fasebin_spe/{datalist} timecol=TIME columns=GOOD accumulate=ONE binsz={binsize} printmode=SPECTRUM lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint=INDEF  chmin = INDEF chmax=INDEF chint=INDEF chbin=INDEF'
            run_command(extract,out_path=self.out_path)

            #in SA config B_16ms_46M_0_49_H (for v0332+53)
            #D[0~4] & E[X1L^X1R^X2L^X2R^X3L^X3R]
            #all layers are added and detectors are added, so we cannot filter on detectors

            os.chdir(self.out_path+'/products/fasebin_spe')
            pcarsp=f'pcarsp -f {datalist}.pha -a ../../FILTERFILE.xfl -l {layers} -j y -p {detectors} -m y -n {datalist}.rsp'

            run_command(pcarsp,self.out_path)

            os.system('cp sa.rsp ../fasebin/response.rsp')



        if datalist=='se':
            os.chdir(self.out_path)

            #detectors_bool,detectors=select_pcu(self.out_path+'/FILTERFILE.xfl',plot=0)
            detectors='2'
            layers='L1,R1,L2,R2,L3,R3'


            extract=f'seextrct infile=@{datalist}.xdf gtiorfile=APPLY gtiandfile=maketime_file_pcu2.fits outroot=products/fasebin_spe/{datalist} timecol=TIME columns=EVENT binsz={binsize} printmode=SPECTRUM lcmode=RATE spmode=SUM timemin=INDEF timemax=INDEF timeint=INDEF  chmin = INDEF chmax=INDEF chint=INDEF chbin=INDEF'
            run_command(extract,out_path=self.out_path)

            #in SE config 'E_125us_64M_0_1s' (for v0332+53)
            #D[0~4] & E[X1L^X1R^X2L^X2R^X3L^X3R] & C[0~249]
            #(M[1]{1},D[0:4]{3},C[0~4,5~6,7,8,9,10,11,12,13,14,15,16~17,1 etc
            #i can chose detectors only, not layers
            #i hence use  all layers of pcu2

            os.chdir(self.out_path+'/products/fasebin_spe')
            pcarsp=f'pcarsp -f {datalist}.pha -a ../../FILTERFILE.xfl -l {layers} -j y -p {detectors} -m y -n {datalist}.rsp'

            run_command(pcarsp,self.out_path)

            os.system('cp se.rsp ../fasebin/response.rsp')

    def make_fasebin(self,nph=16):
        os.chdir(self.out_path)
        ser=self.pandas_series()
        datalist=ser['fasebin_cfg']
        if datalist=='None':
            raise Exception('No SE and SA: fasebin impossible')
        else:
            pass

        datafiles=datalist+'.xdf'

        period=ser['period_orb_corr']
        if np.isnan(period):
            raise Exception('period is not determined. ')


        freq=1/period
        mjd_start=ser['MJD_START']

        orbitfile=glob(self._orbit_path_+'/*')
        if len(orbitfile)==1:
            orbitfile=orbitfile[0]
        else:
            raise Exception('more than one orbit files found!!!')

        os.chdir(self.out_path+'/products/fasebin')

        fname=self.fasebin_info_file
        os.system(f'> {fname}') #remove content from fasebin_info


        fasebin_path=os.getcwd()
        for filename in ['/fasebin.sh','/psrtime.dat','/fasebin.pha','/fasebin.dat','/maketime_file.fits','/'+datafiles,'cutoffpl.xcm','ph_res.dat']:
            if os.path.exists(fasebin_path+filename):
                os.remove(fasebin_path+filename)

        for filename in [datafiles,'maketime_file.fits']:
                os.system(f'cp ../../{filename} ./')


        mjd_start=str(mjd_start)
        freq=str(freq)


        psrtime_file_path='/Users/s.bykov/work/xray_pulsars/rxte/psrtime_msv.dat'   #because Molkov's file has appropriate formatting, i cannot echo  v0332+53 params to standart psrtime.dat I dont know why. I will replace times of observation in that file with mine dates
        os.system(f'cp {psrtime_file_path} ./')

        psrbinfile_path='/Users/s.bykov/work/xray_pulsars/rxte/psrbin.dat'   #because Molkov's file has appropriate formatting, i cannot echo  v0332+53 params to standart psrtime.dat I dont know why. I will replace times of observation in that file with mine dates
        os.system(f'cp {psrbinfile_path} ./')


        #orig_line='0332+53  03 34 59.910  53 10 23.30 45658 57346 53360.000000000   0.2285547100000  0.00000D+00   0.00D+00  2.2 P   0332+53'


        f = open("psrtime_msv.dat")
        o = open(f"psrtime.dat","w")
        for line in f.readlines():
            if line.startswith('0332+53'):
                  print(line)
                  line = line.replace('0.2285547100000',freq)
                  #line = line.replace('45658',mjd_start)
                  #line = line.replace('57346',mjd_stop)
                  line=line.replace('2.2 P','0.0 P')
                  line = line.replace('53360.000000000',mjd_start)
                  #line = line.replace('2.2','0.0')
                  print(line)
            else:
              line=line
            o.write(line)
        o.close()


        with open ('fasebin.sh', 'w+') as rsh:
            rsh.write(f'''\
        #! /bin/bash
        TIMING_DIR={fasebin_path}
        export TIMING_DIR
        fasebin orbitfile={orbitfile} sourcename='0332+53' datafile=@{datafiles} outfile='fasebin_orig.pha' gtifile=maketime_file.fits nph={nph} binary='yes'

        ''')

        os.system('chmod +x fasebin.sh')
        os.system('bash fasebin.sh')





    def fit_ph_res(self,model='cutoffpl',chmin=0,chmax=8,error=0.01):
        '''
        adds sys_err to each spectrum
        '''
        os.chdir(self.out_path+'/products/fasebin')

        if not os.path.exists('fasebin_orig.pha'):
            raise Exception('fasebin file not found')

        if not os.path.exists('response.rsp'):
            os.system('cp ../fasebin_spe/s[e,a].rsp ./response.rsp')

        #for filePath in [f'./{model}/ph_res_{model}.dat',f'./{model}/ph_res_{self.ObsID}_{model}.png']:
        #    if os.path.exists(filePath):
        #        os.remove(filePath)

        os.system(f'rm -rf {model}')

        def apply_systematic_error_fasebin(error,chmin,chmax):
            os.system('cp fasebin_orig.pha fasebin_tmp.pha')
            with fits.open('fasebin_tmp.pha',mode='update') as ff:
                N=len(ff[1].data['counts'][0])
                nph=len(ff[1].data)
                #sys_err=astropy.io.fits.column.Column(name = 'SYS_ERR', format = '46D')
                ff[1].data.columns['BASELINE'].name='SYS_ERR'
                ff[1].data.columns['SYS_ERR'].unit=None
                #format is still the same

                sys=np.zeros(N)
                sys[chmin-1:chmax]=error
                for i in range(nph):
                    ff[1].data['SYS_ERR'][i]=sys
                ff.close()
            os.system('mv fasebin_tmp.pha fasebin_sys.pha')

        apply_systematic_error_fasebin(error=error,chmin=chmin+1,chmax=chmax)


        conf=glob('*.xdf')[0][0:2]
        if conf=='se':
            model=model+'_ign11'
        else:
            model=model
        os.system(f'xspec - /Users/s.bykov/work/xray_pulsars/rxte/python_pca_pipeline/xspec_scripts/ph_res_{model}.txt')


    def ph_res_results(self,model='cutoffpl',rewrite=1):



            os.chdir(self.out_path)
            ser=self.pandas_series()
            datalist=ser['fasebin_cfg']
            period=ser['period_orb_corr']
            mjd=ser['MJD_START']
            if datalist=='None':
                raise Exception('No SE and SA: fasebin impossible')
            else:
                pass


            if datalist=='se':
                datamode='E_125us_64M_0_1s'
            elif datalist=='sa':
                datamode='B_16ms_46M_0_49_H'

            os.chdir(self.out_path+'/products/fasebin')
            matplotlib.rcParams['figure.figsize'] = 6.6*2, 6.6
            matplotlib.rcParams['figure.subplot.left']=0.05
            matplotlib.rcParams['figure.subplot.bottom']=0.1
            matplotlib.rcParams['figure.subplot.right']=0.95
            matplotlib.rcParams['figure.subplot.top']=0.85
            if rewrite:
                fname=self.fasebin_info_file
                os.system(f'> {fname}') #remove content from fasebin_info
            #period=ser['period']
            mjd_obs=ser['MJD_START']
            expo=ser['EXPOSURE']



            os.chdir(model)
            filename=f'ph_res_{model}.dat'

            data=np.genfromtxt(filename)

            N_sp=(data[0,1]-1)/2
            spe_num=data[:,0]

            if len(spe_num)!=N_sp:
                raise Exception('Not all spectra were fitted')

            chi2_red=data[:,2]

            bad_phase_num=int(spe_num[np.argmax(chi2_red)])
            better_phase_num=int(spe_num[np.argsort(chi2_red)[-2]])
            godder_phase_num=int(spe_num[np.argsort(chi2_red)[1]])
            good_phase_num=int(spe_num[np.argmin(chi2_red)])



            data=np.vstack((data,data))
            nph=data[0,1]
            data[:,0]=np.arange(1,nph) #create new spe_num
            spe_num=data[:,0]
            phase=((spe_num-1)/(N_sp))

            chi2_red=data[:,2]

            if model=='cutoffpl_no_gauss':
                fig = plt.figure()
                rows=6
                cols=5

                ax_chi=plt.subplot2grid((rows,cols), (1, 0), rowspan=1, colspan=3)
                ax_flux=plt.subplot2grid((rows,cols), (0, 3), rowspan=2, colspan=2)
                ax_per=ax_flux.twinx()
                #ax_chi_hist=plt.subplot2grid((rows,cols), (2, 3), rowspan=1, colspan=1)
                ax_per_find=plt.subplot2grid((rows,cols), (2, 4), rowspan=1, colspan=1)


                fit_efsearch_data(self.out_path+'/products/std1_lc/std1_0.1s_bary.efs',
                                 self.ObsID,ax=ax_per_find,savefig=0,fit_data=0)
                ax_per_find.axvline(period,color='r')
                ax_flux.plot(ObsParams_plot.MJD_START,ObsParams_plot.cutoffpl_tot_flux/1e-8,'b.')
                ax_per.plot(ObsParams_plot.MJD_START,ObsParams_plot.period_orb_corr,'g.')
                ax_per.set_ylim(4.372,4.378)
                ax_chi.set_title(self.ObsID+f'\n Exp: {int(expo)} s \n model: {model}; config: {datamode}')
                ax_chi.step(phase,chi2_red,'b-',where='mid')
                ax_chi.axhline(1,color='b')
                ax_chi.set_ylabel('xi2',color='b')
                ax_chi.set_xlabel('Phase')
                ax_chi.grid(1,'both')

                ax_del_worst=plt.subplot2grid((rows,cols), (3, 0), rowspan=1, colspan=2)
                ax_del_mid=plt.subplot2grid((rows,cols), (4, 0), rowspan=1, colspan=2)
                ax_del_worst.axhline(0,color='k')
                ax_del_mid.axhline(0,color='k')

                ax_del_worst.axhline(0,color='k')
                ax_del_worst.set_xlabel('E, keV')
                ax_del_worst.set_ylabel('$\chi$')
                ax_del_mid.axhline(0,color='k')
                ax_del_mid.set_xlabel('E, keV')
                ax_del_mid.set_ylabel('$\chi$')

                sp1=Spectra(f'./spe_plots/ph_spe_{bad_phase_num}.dat')
                sp1.plot_del(ax_del_worst,mfc='r',label=f'phase {bad_phase_num}')

                sp2=Spectra(f'./spe_plots/ph_spe_{good_phase_num}.dat')
                sp2.plot_del(ax_del_mid,mfc='b',label=f'phase {good_phase_num}')


                sp3=Spectra(f'./spe_plots/ph_spe_{better_phase_num}.dat')
                sp3.plot_del(ax_del_worst,mfc='m',label=f'phase {better_phase_num}')

                sp4=Spectra(f'./spe_plots/ph_spe_{godder_phase_num}.dat')
                sp4.plot_del(ax_del_mid,mfc='g',label=f'phase {godder_phase_num}')


                for axis in [ax_del_worst,ax_del_mid]:
                    axis.legend(loc='upper right')
                try:
                    ax_flux.axvline(mjd_obs,zorder=-100,color='r',ls='-.',alpha=0.7)
                except:
                    pass

                ObsID=self.ObsID
                fig.savefig(f'Day{mjd}_ph_res_{ObsID}_{model}_orb_corr.png')
                plt.close(fig)


                return None


            eqw=data[:,4]
            eqw_low=eqw-data[:,5]
            eqw_hi=data[:,6]-eqw

            eqw=eqw*1e3
            eqw_low=eqw_low*1e3
            eqw_hi=eqw_hi*1e3
            eqw_err=np.vstack((eqw_low, eqw_hi))

            eqw_pf=pulsed_fraction_error(eqw,np.max(eqw_err,axis=0))

            po=data[:,10]
            efold=data[:,11]
            ecut=data[:,12]
            eline=data[:,13]
            norm_line=data[:,14]

            flux712=data[:,7]
            flux712_low=flux712-data[:,8]
            flux712_hi=data[:,9]-flux712

            flux712=flux712/1e-8
            flux712_hi=flux712_hi/1e-8
            flux712_low=flux712_low/1e-8

            flux712_err=np.vstack((flux712_low, flux712_hi))


            fig = plt.figure()
            rows=6
            cols=5
            #(rows,cols), (y,x) <- those are coordinates of an axis in subplots
            ax_eqw = plt.subplot2grid((rows,cols), (0, 0), rowspan=1, colspan=3)
            ax_efold=ax_eqw.twinx()
            ax_chi=plt.subplot2grid((rows,cols), (1, 0), rowspan=1, colspan=3)
            ax_ccf=plt.subplot2grid((rows,cols), (3, 2), rowspan=2, colspan=3)
            ax_ecutpars=plt.subplot2grid((rows,cols), (2, 0), rowspan=1, colspan=3)
            ax_flux=plt.subplot2grid((rows,cols), (0, 3), rowspan=2, colspan=2)
            ax_per=ax_flux.twinx()
            #ax_chi_hist=plt.subplot2grid((rows,cols), (2, 3), rowspan=1, colspan=1)
            ax_per_find=plt.subplot2grid((rows,cols), (2, 4), rowspan=1, colspan=1)
            ax_po_and_norm=plt.subplot2grid((rows,cols), (5, 0), rowspan=1, colspan=3)



            ax_del_worst=plt.subplot2grid((rows,cols), (3, 0), rowspan=1, colspan=2)
            ax_del_mid=plt.subplot2grid((rows,cols), (4, 0), rowspan=1, colspan=2)
            ax_del_worst.axhline(0,color='k')
            ax_del_mid.axhline(0,color='k')

            ax_del_worst.axhline(0,color='k')
            ax_del_worst.set_xlabel('E, keV')
            ax_del_worst.set_ylabel('$\chi$')
            ax_del_mid.axhline(0,color='k')
            ax_del_mid.set_xlabel('E, keV')
            ax_del_mid.set_ylabel('$\chi$')

            sp1=Spectra(f'./spe_plots/ph_spe_{bad_phase_num}.dat')
            sp1.plot_del(ax_del_worst,mfc='r',label=f'phase {bad_phase_num}')

            sp2=Spectra(f'./spe_plots/ph_spe_{good_phase_num}.dat')
            sp2.plot_del(ax_del_mid,mfc='b',label=f'phase {good_phase_num}')


            sp3=Spectra(f'./spe_plots/ph_spe_{better_phase_num}.dat')
            sp3.plot_del(ax_del_worst,mfc='m',label=f'phase {better_phase_num}')

            sp4=Spectra(f'./spe_plots/ph_spe_{godder_phase_num}.dat')
            sp4.plot_del(ax_del_mid,mfc='g',label=f'phase {godder_phase_num}')


            for axis in [ax_del_worst,ax_del_mid]:
                axis.legend(loc='upper right')


            ax_efold.errorbar(phase,flux712,flux712_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
            ax_eqw.errorbar(phase,eqw,eqw_err,color='r',drawstyle='steps-mid',alpha=0.8)
            ax_eqw.tick_params(axis='y', colors='red')
            ax_eqw.spines['left'].set_color('red')

            ax_eqw.set_ylabel('Fe Ka Eq. width, eV',color='r')
            ax_efold.set_ylabel('Flux 7-12, 1e-8 cgs')
            ax_eqw.set_xlabel('Phase')

#            eqw_max=np.where(eqw==max(eqw))[0]
#            flux_max=np.where(flux712==max(flux712))[0]
#            #dT=(phase[eqw_max[eqw_max>flux_max[0]]]-phase[flux_max[0]])*period
#            ax_eqw.hlines(y=max(eqw),
#                          xmin=phase[flux_max[0]],xmax=phase[eqw_max[eqw_max>flux_max[0]]],
#                          linestyles=':',color='m',lw=2,alpha=0.8)
            ax_eqw.set_title(self.ObsID+f'\n Exp: {int(expo)} s \n model: {model}; config: {datamode}')#+f'max-max lag is {dT} s')

            #ax_chi_hist.hist(chi2_red,color='b',alpha=0.7,bins=10)
            #ax_chi_hist.set_xlabel('chi^2/dof hist')

            fit_efsearch_data(self.out_path+'/products/std1_lc/std1_0.1s_bary.efs',
                              self.ObsID,ax=ax_per_find,savefig=0,fit_data=0)
            ax_per_find.axvline(period,color='r')

            ax_chi.step(phase,chi2_red,'b-',where='mid')
            ax_chi.axhline(1,color='b')
            ax_chi.set_ylabel('xi2',color='b')
            ax_chi.set_xlabel('Phase')
            ax_chi.grid(1,'both')


            deltat,_,_=my_crosscorr(phase*period,eqw,flux712,None,ax_ccf,
                         subtract_mean=1,divide_by_mean=1,only_pos_delays=0,
                         y1label='eqw',y2label='F(7-12)',my_only_pos_delays=0)
            deltat_err=period/N_sp

            ax_ecutpars.plot(phase,efold,color='g',label='efold (left)')
            #ax_ecutpars_twin=ax_ecutpars.twinx()
            #ax_ecutpars_twin.plot(phase,ecut,color='y',label='ecut (right)')
            ax_ecutpars.legend(loc='upper right')
            #ax_ecutpars_twin.legend(loc='upper left')
            ax_ecutpars.plot(phase,eline,color='y',label='eline (left)')
            ax_ecutpars.legend(loc='upper left')
            #ax_ecutpars_twin=ax_ecutpars.twinx()
            #ax_ecutpars_twin.plot(phase,po,color='k',label='po gamma (right)')
            #ax_ecutpars_twin_2=ax_ecutpars.twinx()
            #ax_ecutpars_twin_2.plot(phase,norm_line,color='m')

            #ax_ecutpars_twin.legend(loc='upper right')
            ax_po_and_norm_twin=ax_po_and_norm.twinx()
            ax_po_and_norm_twin.plot(phase,po,color='k',label='po gamma (right)')
            ax_po_and_norm.plot(phase,norm_line,color='m',label='iron line norm (left)')
            ax_po_and_norm_twin.legend(loc='upper right')
            ax_po_and_norm.legend(loc='upper left')

            ax_flux.plot(ObsParams_plot.MJD_START,ObsParams_plot.cutoffpl_tot_flux/1e-8,'b.')
            ax_per.plot(ObsParams_plot.MJD_START,ObsParams_plot.period_orb_corr,'g.')
            ax_per.set_ylim(4.372,4.378)
            try:
                ax_flux.axvline(mjd_obs,zorder=-100,color='r',ls='-.',alpha=0.7)
            except:
                pass
            plt.subplots_adjust(wspace=0.7)
            plt.subplots_adjust(hspace=0.3)



            self.write_to_obs_info(self.fasebin_info_file,'eqw_PF_err',eqw_pf[1])
            self.write_to_obs_info(self.fasebin_info_file,'eqw_PF',eqw_pf[0])
            self.write_to_obs_info(self.fasebin_info_file,'deltat',deltat[0])
            self.write_to_obs_info(self.fasebin_info_file,'wrap_deltat',deltat[1])
            self.write_to_obs_info(self.fasebin_info_file,'deltat_err',deltat_err)

            ObsID=self.ObsID
            fig.savefig(f'Day{mjd}_ph_res_{ObsID}_{model}_orb_corr.png')
            plt.close(fig)


            matplotlib.rcParams['figure.figsize'] = 6.6, 6.6/3
            matplotlib.rcParams['figure.subplot.left']=0.1
            matplotlib.rcParams['figure.subplot.bottom']=0.1
            matplotlib.rcParams['figure.subplot.right']=0.9
            matplotlib.rcParams['figure.subplot.top']=0.85
            fig,ax_eqw=plt.subplots()
            ax_efold=ax_eqw.twinx()

            ax_efold.errorbar(phase,flux712,flux712_err,color='k',label='Flux 7-12',drawstyle='steps-mid',ls=':',alpha=0.6)
            ax_eqw.errorbar(phase,eqw,eqw_err,color='r',drawstyle='steps-mid',alpha=0.6)

            ax_eqw.set_ylabel('Fe Ka Eq. width, eV',color='r')
            ax_efold.set_ylabel('Flux 7-12, 1e-8 cgs')
            ax_eqw.set_xlabel('Phase')
            ax_eqw.set_title(self.ObsID+f' ({datamode})')


            fig.tight_layout()
            sns.despine(fig,top=1,right=0)
            fig.savefig(f'Day{mjd}_ph_res_{ObsID}_{model}_report.png',dpi=500)


            plt.close(fig)
            plt.close('all')




