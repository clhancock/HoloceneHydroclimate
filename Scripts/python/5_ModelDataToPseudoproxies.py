
#Load packages
#import pyleoclim as pyleo               # Packages for analyzing LiPD files 
#import lipd                             # Packages for analyzing LiPD files 
import pickle
import csv

import numpy  as np                     # Package with useful numerical functions
import xarray as xr                     #Package for handling file types
import pandas as pd
#import math
import regionmask
from scipy import stats                 # Packages for calculations  
import pymannkendall as mk              # Package for trend detection
import matplotlib.pyplot as plt         # Packages for making figures
import matplotlib as mpl 
#import matplotlib.gridspec as gridspec
#import matplotlib.colors as pltcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs              # Packages for mapping in python
import cartopy.feature as cfeature
#import cartopy.util as cutil
#
#
from scipy.stats import pearsonr

#data_T =  pd.read_csv(gitHubDir+'DataFiles/proxyT.csv')


#%%
#Load Model Data
#
#Set time variables and resolution of data
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/2021_HoloceneHydroclimate/2021_HoloceneHydroclimate/'
modelData = {}
for model in ['hadcm','trace','cmip6']:
    modelData[model] = {}
    for szn in ['ANN','JJA','DJF']:
        modelData[model][szn] = xr.open_dataset(dataDir+'Data/Model/'+model+'/'+model+'_'+szn+'.nc',decode_times=False)

proxyDF =  pd.read_csv(dataDir+'/Data/Proxy/proxyMetaData_HC.csv')

save=False




#%%
#%%
proxyVar    = 'HC'
modelVar    = 'pre'
season      = 'proxy'
standardize = True

regions = []
siteNum = []
rvalueT = []
rvalueH = []
rvalMax = []
rvalAvg = []
pctP = []
pctAnn = []


for reg in np.unique(proxyDF['ipccReg']):
    if reg not in regionmask.defined_regions.ar6.land.abbrevs: continue
    print(reg)
    #Select proxy sites for region
    regData  = proxyDF.loc[proxyDF['ipccReg']==reg] 
    regions.append(reg)
    siteNum.append(np.shape(regData)[0])
    for model in ['hadcm','trace']:
        #Empty dataframe for pseudoproxy timeseries
        pesudoDF = pd.DataFrame()
        ann_n = 0
        p_n = 0
        count=0
        for i in list(regData.index):
            #climate variable
            if   proxyVar == 'pre':                        var = 'pre'
            elif proxyVar == 'p-e':                        var = 'p-e'
            elif proxyVar == 'HC': 
                if regData['climInterp'][i] == 'P':        var = 'pre'
                else:                                      var = 'p-e'
            elif proxyVar == 'tas':                        var = 'tas'
            #seasonality
            if   season == 'annual':                       szn = 'ANN'
            elif season == 'summer':                       szn = 'JJA'
            elif season == 'winter':                       szn = 'DJF'
            else: 
                if   'ummer' in str(regData['season'][i]): szn = 'JJA'
                elif 'inter' in str(regData['season'][i]): szn = 'DJF'
                else:                                      szn = 'ANN'
            if regData['latitude'][i] < 0:         
                if   szn == 'JJA':                         szn = 'DJF'
                elif szn == 'DJF':                         szn = 'JJA'
            if szn == 'ANN': ann_n += 1
            if var == 'pre': p_n += 1
            count += 1
            #geography
            lat = np.argmin(np.abs(modelData[model][szn]['lat'].data-regData['latitude'][i]))
            lon = np.argmin(np.abs(modelData[model][szn]['lon'].data-regData['longitude'][i]))
            #site timeseries
            modelvalues = modelData[model][szn][var][:,lat,lon]
            if standardize: modelvalues = stats.zscore(modelvalues)
            #mask = [idx for idx, val in enumerate(dataModel[model]['time']) if (val < ror val > site['ageMax'])]
            pesudoDF[regData['tsid'][i]] = modelvalues
        pesudoTS = pesudoDF.mean(axis=1)
        regionTS = pd.read_csv(dataDir+'/Data/Model/regionalTS/regional_'+modelVar+'_ANN_'+model+'_land.csv')[reg]
        if   model == 'hadcm': corrH = round(pearsonr(pesudoTS,regionTS)[0],2)
        elif model == 'trace': corrT = round(pearsonr(pesudoTS,regionTS)[0],2)
    rvalueT.append(corrT)
    rvalueH.append(corrH)
    rvalMax.append(max([corrT,corrH]))
    rvalAvg.append(np.mean([corrT,corrH]))
    pctP.append(p_n/count)
    pctAnn.append(ann_n/count)
    

df = pd.DataFrame([regions,siteNum,rvalueT,rvalueH,pctP,pctAnn]).transpose()
df.to_csv(str(dataDir+'Data/Model/'+'pseudoProxyCorr'+'_'+modelVar+'.csv'))


