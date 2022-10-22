#Load packages
import pyleoclim as pyleo               # Packages for analyzing LiPD files 
#import lipd                             # Packages for analyzing LiPD files 
import numpy  as np                     # Package with useful numerical functions
#import xarray as xr                     #Package for handling file types
import pandas as pd
#import math
from scipy import stats                 # Packages for calculations  
import pymannkendall as mk              # Package for trend detection
import matplotlib.pyplot as plt         # Packages for making figures
#import matplotlib as mpl 
import matplotlib.gridspec as gridspec
#import matplotlib.colors as pltcolors
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs              # Packages for mapping in python
import cartopy.feature as cfeature
#import cartopy.util as cutil
#
dataDir='/Volumes/GoogleDrive/My Drive/zResearch/Data/'
save=False
gitHub ='/Volumes/GoogleDrive/My Drive/zResearch/HoloceneHydroclimate/'
#%%
#timebin = timebin; setTime = setTime
clusterComposite = {}
for climVar in ['HC']:
    clusterComposite[climVar] = {}
    clusterData = pd.read_csv(gitHub+'DataSummary/proxyHC_regions.csv')
    for cluster in clusterData['Acronym'].unique():
        if sum(clusterData['Acronym']==cluster) >= 6:        
            print(cluster)
            clusterDict = {}
            clusterDF = clusterData.loc[(clusterData['Acronym'] == cluster)]
            clusterDict['ProxyN']    = len(clusterDF['Lat'])
            clusterDict['Lat']       = list(clusterDF['Lat'])
            clusterDict['Lon']       = list(clusterDF['Lon'])
            clusterDict['Composite'] = pd.read_csv(gitHub+'DataSummary/RegionalComposites/'+climVar+'/'+cluster+'.csv')
            clusterDict['Median']    = list(clusterDict['Composite'].median(axis=1))
            clusterDict['Time']      = timebin
            clusterDict['dataCalc']  = {}
            for time in setTime: 
                clusterDict['dataCalc'][time+'CompSlope'] = calcTrend(clusterDict['Time'],
                                                                      clusterDict['Median'],
                                                                      setTime[time]['min'],
                                                                      setTime[time]['max'],
                                                                      1) 
            for time in range(0,12000,1000):   
                z = [idx for idx, age in enumerate(clusterDict['Time']) if age >= time and age <= time+1000]
                clusterDict['dataCalc']['bin'+str(int(time/1000))+'ka'] = np.nanmean(clusterDict['Median'][min(z):max(z)])





#%%
        
      #  if climateVariable == 'HC': clusterDict['proxyData2'] = data_hc_df.loc[data_hc_df['TSid'].isin(list(clusterDF['TSid']))]
     #   else: clusterDict['proxyData2'] = data_temp_df.loc[data_temp_df['TSid'].isin(list(clusterDF['TSid']))]
        
        
        for time in setTime: 
            clusterDict[time] = {}
            clusterDict[time]['CompSlope'] = [] 
            for i in list(clusterDict['proxyComposite'].columns):
                variable = calcTrend(clusterDict['Time'],list(clusterDict['proxyComposite'][i]),setTime[time]['min'],setTime[time]['max'],1)
                clusterDict[time]['CompSlope'].append(variable['Direction'])   
            variable = calcTrend(clusterDict['Time'],clusterDict['Median'],setTime[time]['min'],setTime[time]['max'],1)
            clusterDict[time]['CompSlopes%'] = list(clusterDict[time]['CompSlope']).count(variable['Direction'])/len(list(clusterDict[time]['CompSlope']))  
        pseudoproxy = {'trace':{},'hadcm':{}}
        for model in pseudoproxy:
            pseudoproxy[model]['array']=[]
            clusterDict[model] = {}
            for tsid in clusterDict['proxyData']['TSid']: 
                try:
                    mask = [idx for idx, val in enumerate(modelData[model]['pseudoproxyTS'][climateVariable]['TSid']) if val == tsid][0] 
                    pseudoproxy[model]['array'].append(modelData[model]['pseudoproxyTS'][climateVariable]['modelTS'][mask])
                except: print('problem with finding TSid: '+tsid)
            pseudoproxy[model]['array']  = np.column_stack(pseudoproxy[model]['array'])
            pseudoproxy[model]['array'] = pseudoproxy[model]['array'].transpose()
            pseudoproxy[model]['median'] = np.nanmedian(pseudoproxy[model]['array'],axis=0)
            clusterDict[model]['array'] = pseudoproxy[model]['array']
            clusterDict[model]['median'] = pseudoproxy[model]['median'] 
        for time in setTime:   
            clusterDict[time]['proxySlope'] = calcTrend(clusterDict['Time'],clusterDict['Median'],setTime[time]['min'],setTime[time]['max'],1)
            for model in pseudoproxy:
                clusterDict[time][model+'Slope'] = calcTrend(clusterDict['Time'],list(pseudoproxy[model]['median']),setTime[time]['min'],setTime[time]['max'],1)
                clusterDict[time][model+'Slopelist'] = []
                for site in pseudoproxy[model]['array']:
                    clusterDict[time][model+'Slopelist'].append(calcTrend(clusterDict['Time'],list(site),setTime[time]['min'],setTime[time]['max'],1))
                clusterDict[time][model+'Slopelist'] = pd.DataFrame(clusterDict[time][model+'Slopelist'])
                calcTrend(clusterDict['Time'],list(pseudoproxy[model]['median']),setTime[time]['min'],setTime[time]['max'],1)
                z = list(clusterDict['proxyData2'][time+'Direction'])
            if clusterDict[time]['proxySlope']['Direction'] == np.NaN:
                clusterDict[time]['proxySlopes%'] = np.NaN
            elif  clusterDict[time]['proxySlope']['Direction'] == 'no trend':
                try: clusterDict[time]['proxySlopes%'] = z.count(clusterDict[time]['proxySlope']['Direction'])/(len(z)-z.count(np.NaN))
                except: clusterDict[time]['proxySlopes%'] = np.NaN
                try: clusterDict[time]['proxySlopes%%'] = clusterDict[time]['proxySlopes%'] 
                except: clusterDict[time]['proxySlopes%%'] = np.NaN
            else:
                try: clusterDict[time]['proxySlopes%'] = z.count(clusterDict[time]['proxySlope']['Direction'])/(len(z)-z.count(np.NaN))
                except: clusterDict[time]['proxySlopes%'] = np.NaN
                try: clusterDict[time]['proxySlopes%%'] = z.count(clusterDict[time]['proxySlope']['Direction'])/(len(z)-z.count(np.NaN)-z.count("no trend"))
                except: clusterDict[time]['proxySlopes%%'] = np.NaN
            clusterDict[time]['agreement'] = 0
            for model in ['trace','hadcm']:
                if np.sign(clusterDict[time]['proxySlope']['linSlope']) == np.sign(clusterDict[time][model+'Slope']['linSlope']): clusterDict[time]['agreement'] += 1
            if clusterDict[time]['proxySlope']['Direction'] == 'decreasing':clusterDict[time]['agreement'] *= -1
            clusterDict[time]['agreement+'] = 0
            for model in ['trace','hadcm']:
                if clusterDict[time]['proxySlope']['Direction'] == clusterDict[time][model+'Slope']['Direction']: clusterDict[time]['agreement+'] += 1
            if clusterDict[time]['proxySlope']['Direction'] == 'decreasing':clusterDict[time]['agreement+'] *= -1
        clusterComposite[climateVariable][clusterName] = clusterDict