
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
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate/'
modelData = {}
for model in ['hadcm','trace','cmip6']:
    modelData[model] = {}
    for szn in ['ANN','JJA','DJF']:
        modelData[model][szn] = xr.open_dataset(dataDir+'Data/Model/'+model+'_'+szn+'.nc',decode_times=False)

proxyDF =  pd.read_csv(dataDir+'/Data/proxyMetaData_HC.csv')

save=False




#%%
#%%
variable    = 'p-e'
season      = 'annual'
standardize = False

regions = []
siteNum = []
rvalueT = []
rvalueH = []
rvalMax = []
rvalAvg = []

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
        for i in list(regData.index):
            #climate variable
            if   variable == 'pre':                        var = 'pre'
            elif variable == 'p-e':                        var = 'p-e'
            elif variable == 'HC': 
                if regData['climInterp'][i] == 'P':        var = 'pre'
                else:                                      var = 'p-e'
            elif variable == 'tas':                        var = 'tas'
            #seasonality
            if   season == 'annual':                       szn = 'ANN'
            elif season == 'summer':                       szn = 'JJA'
            elif season == 'winter':                       szn = 'DJF'
            else: 
                if   'ummer' in str(data_HC['season'][i]): szn = 'JJA'
                elif 'inter' in str(data_HC['season'][i]): szn = 'DJF'
                else:                                      szn = 'ANN'
            if regData['latitude'][i] < 0:         
                if   szn == 'JJA':                         szn = 'DJF'
                elif szn == 'DJF':                         szn = 'JJA'
            #geography
            lat = np.argmin(np.abs(modelData[model][szn]['lat'].data-regData['latitude'][i]))
            lon = np.argmin(np.abs(modelData[model][szn]['lon'].data-regData['longitude'][i]))
            #site timeseries
            modelvalues = modelData[model]['ANN']['p-e'][0:121,lat,lon]
            if standardize: modelvalues = stats.zscore(modelvalues)
            #mask = [idx for idx, val in enumerate(dataModel[model]['time']) if (val < ror val > site['ageMax'])]
            pesudoDF[regData['tsid'][i]] = modelvalues
        pesudoTS = pesudoDF.mean(axis=1)
        regionTS = pd.read_csv(dataDir+'/Data/Model/regionTS/regional_'+var+'_ANN_'+model+'.csv')[reg][0:121]
        if   model == 'hadcm': corrH = round(pearsonr(pesudoTS,regionTS)[0],2)
        elif model == 'trace': corrT = round(pearsonr(pesudoTS,regionTS)[0],2)
    rvalueT.append(corrT)
    rvalueH.append(corrH)
    rvalMax.append(max([corrT,corrH]))
    rvalAvg.append(np.mean([corrT,corrH]))

df = pd.DataFrame([regions,siteNum,rvalueT,rvalueH,rvalMax,rvalAvg]).transpose()
df.to_csv(str(dataDir+'Data/Model/'+'pseudoProxyCorr'+'_'+var+'.csv'))






#%%
def calculatePseudoProxy(proxyDF=data_HC,
                         dataModelDict=dataModel,
                         variable   ='P',#Or T # Change to 'P' or 'P-E' to test assumptions about proxies
                         season     ='annual', # Change to 'ANN' or 'Summer' to test assumptions about proxies
                         standardize=False, #Will calculate z scores for records relative to Holocene timeseries
                         lm=landmask):
    output={}
    for model in dataModelDict.keys():
        output[model] = dict()
        PseudoProxy = {'medianTS':pd.DataFrame()}
        modelLandTS = {'medianTS':pd.DataFrame()}
        for region in np.unique(proxyDF['ipccReg']):
            print(region)
            regData = proxyDF.loc[proxyDF['ipccReg']==region]
            regDF   = pd.DataFrame()
            
            PseudoProxy[region] = regDF
            PseudoProxy['medianTS'][region] = regDF.mean(axis=1)
            #Full region model estimate
            refReg = regionmask.defined_regions.ar6.all
            refReg = refReg.mask_3D(dataModel[model]['lon'],dataModel[model]['lat'])
            refReg = refReg.isel(region=(refReg.abbrevs == region))
            refReg = np.array( refReg.values,dtype='bool')#*lm[model]*1
            values = dataModel[model][var][sea][refReg[0,:,:]]
            weight = np.cos(np.deg2rad(dataModel[model]['lat']))  
            weight = np.concatenate(refReg*weight[:, np.newaxis])
            if sum(sum(weight)) == 0: continue
            weight = [i for i in list(np.concatenate(weight))if i != 0]
            wavg   = list(np.average(values,axis=0,weights=weight))
            modelLandTS['medianTS'][region] = wavg
        output[model]['PseudoProxy'] = PseudoProxy
        output[model]['modelLandTS'] = modelLandTS
    return(output)


pseudoProxy = dict()
for v in ['P-E']:#,'P-E','HC']:
    for s in ['annual']:#,'summer','proxy']:
        name = str(v+'_'+s)
        print(name)
        pseudoProxy[name] =  calculatePseudoProxy(variable=v,season=s)
        for model in  pseudoProxy[name].keys():
            df = pseudoProxy[name][model]['modelLandTS']['medianTS']
            df.to_csv(str(dataDir+'HoloceneHydroclimate/Data/Model/PseudoProxy_FullReg/'+name+'_'+model+'.csv'))

from scipy.stats import pearsonr

#%%
reg = []
num = []
corT = []
corH = []
corMax = []
corAvg = []

for region in np.unique(data_HC['ipccReg']):
    n = np.shape(pseudoProxy[name][model]['PseudoProxy'][region])[1]
    if n >= 6:
        reg.append(region)
        model='trace'
        x1 = list(pseudoProxy[name][model]['PseudoProxy']['medianTS'][region])
        y1 = list(pseudoProxy[name][model]['modelLandTS']['medianTS'][region])
        corT.append(pearsonr(x1,y1)[0])
        model='hadcm'
        x2 = list(pseudoProxy[name][model]['PseudoProxy']['medianTS'][region])
        y2 = list(pseudoProxy[name][model]['modelLandTS']['medianTS'][region])
        corH.append(pearsonr(x2,y2)[0])
        corMax.append(max(pearsonr(x2,y2)[0],pearsonr(x1,y1)[0]))
        corAvg.append(np.mean([pearsonr(x2,y2)[0],pearsonr(x1,y1)[0]]))
        num.append(n)

df = pd.DataFrame([reg,num,corT,corH,corMax,corAvg]).transpose()

#%%
#pseudoProxy
#modelLandTS


z = refReg.values
zz =   landmask['trace']     
np.shape(z*zz)                   
zzz = z*zz
            model='hadcm'
            refReg = regionmask.defined_regions.ar6.all
            refReg = refReg.mask_3D(dataModel[model]['lon'],dataModel[model]['lat'])
            refReg = refReg.isel(region=(refReg.abbrevs == "EAN"))
            refReg = np.array( refReg.values*landmask[model],dtype='bool')
            values = dataModel[model]['Precip']['ANN'][refReg[0,:,:]]
            weight = np.cos(np.deg2rad(dataModel[model]['lat']))  
            weight = np.concatenate(refReg*weight[:, np.newaxis])
            weight = [i for i in list(np.concatenate(weight))if i != 0]
            wavg   = list(np.average(values,axis=0,weights=weight))


for v in ['P','P-E','HC']:
    for s in ['annual','summer','proxy']:
        name = str(v+'_'+s)
        print(name)
        pseudoProxy['name'] =  calculatePseudoProxy(variable=v,season=s)



pseudoProxy['P']['annual']
#traceMask <- nc_open(file.path(githubDir,'Data','Model_transient','netcdf',
 #                              'trace','trace.01-36.22000BP.cam2.LANDFRAC.22000BP_decavg_400BCE.nc'))
#traceMask0 <- ncvar_get(traceMask,'LANDFRAC')[,,2204]
#traceMask0[which(traceMask0>=0.5)] <- 1
#traceMask0[which(traceMask0< 0.5)] <- NA

#hadcmMask <- nc_open(file.path(githubDir,'Data','Model_transient','netcdf',
 #                              'hadcm','deglh.vn1_0.ht_mm_srf.monthly.ANN.100yr.nc'))
#hadcmMaskMask0 <- ncvar_get(hadcmMask,'ht_mm_srf')[,,250]
#hadcmMaskMask0[which(hadcmMaskMask0> 0)] <- 1
#hadcmMaskMask0[which(hadcmMaskMask0<=0)] <- NA




tas_regional = data.weighted(regionMask * weights).mean(dim=("lat", "lon"))
regionVals = z[regionMask,]
tas_CNA = tas.where(r1)

#%%

    for model in dataModelDict.keys():
            pseudoproxy={}
            names_key    = ['TSid','Lat','Lon','Category','CategorySpecific','Interp']
            variable_key = ['Direction','Slope','SlopeSig']
            pseudoproxy['modelTS'] = []
            #Set up dictionary lists to add to 
            for time in setTime:  
                for var in variable_key: pseudoproxy[time+var] = []
            for name in names_key: pseudoproxy[name] = [] 
            for age in range(ageMin,ageMax,binRes): 
                pseudoproxy['bin'+str(int(age/binRes))+'ka'] = []
            #Populate data from model values
            for index, site in proxyDF.iterrows():
                #Add metadata
                for name in names_key: pseudoproxy[name].append(site[name])
                #ID seasonality from proxy metadata
                if   seasonality == 'Ann':      season = 'ANN' #assume everying annual 
                elif seasonality == 'Summer': #summer or annual are summer
                    if site['Season'] == 'winterOnly':
                        if site['Lat'] > 0:     season = 'DJF'
                        else:                   season = 'JJA'
                    else: 
                        if site['Lat'] > 0:     season = 'JJA'
                        else:                   season = 'DJF'
                elif seasonality == 'Proxy': #Use proxy metadata
                    if site['Season'] == 'summerOnly': 
                        if site['Lat'] > 0:     season = 'JJA'
                        else:                   season = 'DJF'
                    if site['Season'] == 'winterOnly':
                        if site['Lat'] > 0:     season = 'DJF'
                        else:                   season = 'JJA'
                    else:                       season = 'ANN'
                #ID appropriate climate variable for pseudoproxy
                if   climVar == 'T':          variable = 'Temp'
                elif climVar == 'HC':
                    if site['Interp'] == 'P': variable = 'Precip'
                    else:                     variable = 'EffM'
                elif climVar == 'P':          variable = 'Precip'
                elif climVar == 'P-E':        variable = 'EffM'
                #ID Model grid Location for proxy site
                lat=np.argmin(np.abs(dataModel[model]['lat']-site['Lat']))
                lon=np.argmin(np.abs(dataModel[model]['lon']-site['Lon']))
                #Standardize and get age range
                modelvalues = dataModel[model][variable][season][lat,lon,:]
                if standardize: modelvalues = stats.zscore(modelvalues)
                mask = [idx for idx, val in enumerate(dataModel[model]['time']) if (val < site['ageMin'] or val > site['ageMax'])]
                modelvalues[mask] = np.NaN
                #if len(pseudoproxy['modelTS']) < 1:
                pseudoproxy['modelTS'].append(list(modelvalues))
                #Calculate trends
                for time in setTime:   
                    out = calcTrend(dataModel[model]['time'],list(modelvalues),setTime[time]['min'],setTime[time]['max'],1)               
                    for var in variable_key: 
                        pseudoproxy[time+var].append(out[var])
                mask = [idx for idx, timeval in enumerate(dataModel[model]['time']) if timeval >= 0 and timeval <= 1000]
                lipdTS_std = np.nanmean(modelvalues[mask])
                #Calcaulte Bin values
                for age in range(ageMin,ageMax,binRes):
                    mask = [idx for idx, val in enumerate(dataModel[model]['time']) if val >= age and val < age+ageRes]
                    value_ts = np.nanmean(modelvalues[mask])
                    pseudoproxy['bin'+str(int(age/binRes))+'ka'].append((value_ts-lipdTS_std))
            output[model] = pd.DataFrame.from_dict(pseudoproxy)
    return(output)



#%%
def calcTrend(timeValues,proxyValues,ageMin,ageMax,direction):    
    divisions = 2; minset = 3; sigLevel = 0.05; CI = 0.95; method = 'lin'
    if len(proxyValues) < 12: minset = 2
    ageRange = (ageMax - ageMin)/divisions
    valuedata = pd.DataFrame([timeValues,proxyValues]).transpose()
    valuedata.columns = ['time','values']
    valuedata = valuedata.loc[(valuedata['time'] <= ageMax) & (valuedata['time'] >= ageMin)]
    valuedata = valuedata.dropna() 
    check = []    
    for timeperiod in range(divisions):
        checkValue = len(valuedata.loc[(valuedata['time'] <= ageMin+(timeperiod+1)*ageRange) 
                                     & (valuedata['time'] >= ageMin+(timeperiod)*ageRange)]['time'])
        if checkValue >= minset: check.append(True)
        else: check.append(False)   
    if check.count(True) == divisions:
        valuedata['time']   = valuedata['time']*-1 #To get correct forward direction
        valuedata['values'] = valuedata['values']*direction
        trendDirection = mk.original_test(valuedata['values']*-1,alpha=sigLevel)
        if   method == 'lin': Slope = stats.linregress( valuedata['time'],valuedata['values'])
        elif method == 'sen': Slope = stats.theilslopes(valuedata['values'],valuedata['time'],CI)
        output = {'Direction':trendDirection[0],'Slope':Slope[0]*1000}
        if output['Direction'] == 'no trend': output['SlopeSig'] = 0
        else:                                 output['SlopeSig'] = output['Slope']
    else: output = {'Direction':np.NaN, 'Slope':np.NaN, 'SlopeSig':np.NaN}
    return(output)

#%%
#Calc model pesuodporxy timeseries based on location, season, length, variable
#
setTime = {'Total':{'min':ageMin,                  'max':ageMax},
           'Early':{'min':ageMin+(ageMax-ageMin)/2,'max':ageMax},
           'Late' :{'min':ageMin,                  'max':ageMin+(ageMax-ageMin)/2}}
#ageMin,ageMax,ageRes
def calculatePseudoProxy(proxyDF,
                         dataModelDict=dataModel,
                         climVar='HC',#Or T # Change to 'P' or 'P-E' to test assumptions about proxies
                         seasonality='Proxy', # Change to 'ANN' or 'Summer' to test assumptions about proxies
                         standardize=True, #Will calculate z scores for records relative to Holocene timeseries
                         binRes=1000):
    output={}
    for model in dataModelDict.keys():
            pseudoproxy={}
            names_key    = ['TSid','Lat','Lon','Category','CategorySpecific','Interp']
            variable_key = ['Direction','Slope','SlopeSig']
            pseudoproxy['modelTS'] = []
            #Set up dictionary lists to add to 
            for time in setTime:  
                for var in variable_key: pseudoproxy[time+var] = []
            for name in names_key: pseudoproxy[name] = [] 
            for age in range(ageMin,ageMax,binRes): 
                pseudoproxy['bin'+str(int(age/binRes))+'ka'] = []
            #Populate data from model values
            for index, site in proxyDF.iterrows():
                #Add metadata
                for name in names_key: pseudoproxy[name].append(site[name])
                #ID seasonality from proxy metadata
                if   seasonality == 'Ann':      season = 'ANN' #assume everying annual 
                elif seasonality == 'Summer': #summer or annual are summer
                    if site['Season'] == 'winterOnly':
                        if site['Lat'] > 0:     season = 'DJF'
                        else:                   season = 'JJA'
                    else: 
                        if site['Lat'] > 0:     season = 'JJA'
                        else:                   season = 'DJF'
                elif seasonality == 'Proxy': #Use proxy metadata
                    if site['Season'] == 'summerOnly': 
                        if site['Lat'] > 0:     season = 'JJA'
                        else:                   season = 'DJF'
                    if site['Season'] == 'winterOnly':
                        if site['Lat'] > 0:     season = 'DJF'
                        else:                   season = 'JJA'
                    else:                       season = 'ANN'
                #ID appropriate climate variable for pseudoproxy
                if   climVar == 'T':          variable = 'Temp'
                elif climVar == 'HC':
                    if site['Interp'] == 'P': variable = 'Precip'
                    else:                     variable = 'EffM'
                elif climVar == 'P':          variable = 'Precip'
                elif climVar == 'P-E':        variable = 'EffM'
                #ID Model grid Location for proxy site
                lat=np.argmin(np.abs(dataModel[model]['lat']-site['Lat']))
                lon=np.argmin(np.abs(dataModel[model]['lon']-site['Lon']))
                #Standardize and get age range
                modelvalues = dataModel[model][variable][season][lat,lon,:]
                if standardize: modelvalues = stats.zscore(modelvalues)
                mask = [idx for idx, val in enumerate(dataModel[model]['time']) if (val < site['ageMin'] or val > site['ageMax'])]
                modelvalues[mask] = np.NaN
                #if len(pseudoproxy['modelTS']) < 1:
                pseudoproxy['modelTS'].append(list(modelvalues))
                #Calculate trends
                for time in setTime:   
                    out = calcTrend(dataModel[model]['time'],list(modelvalues),setTime[time]['min'],setTime[time]['max'],1)               
                    for var in variable_key: 
                        pseudoproxy[time+var].append(out[var])
                mask = [idx for idx, timeval in enumerate(dataModel[model]['time']) if timeval >= 0 and timeval <= 1000]
                lipdTS_std = np.nanmean(modelvalues[mask])
                #Calcaulte Bin values
                for age in range(ageMin,ageMax,binRes):
                    mask = [idx for idx, val in enumerate(dataModel[model]['time']) if val >= age and val < age+ageRes]
                    value_ts = np.nanmean(modelvalues[mask])
                    pseudoproxy['bin'+str(int(age/binRes))+'ka'].append((value_ts-lipdTS_std))
            output[model] = pd.DataFrame.from_dict(pseudoproxy)
    return(output)

#modelHC    = calculatePseudoProxy(data_HC)
#%%
#save=True
modelProxies = {}
for Season in ['Proxy','Ann','Summer']:
    for climateInterp in ['HC','P','P-E']:
        name = climateInterp+'_'+Season
        modelProxies[name] = calculatePseudoProxy(data_HC,
                                                  climVar=climateInterp,
                                                  seasonality=Season)
        if save:
            for model in modelProxies[name].keys():
                modelProxies[name][model].to_csv(gitHub+'DataSummary/'+model+'_Pseudoproxy_'+name+'.csv')


#%%


trace.01-36.22000BP.cam2.PRECT.22000BP_decavgANN_400BCE.nc calculatePseudoProxy(data_HC,climVar='HC',seasonality='Proxy')



#%%
#%%
#Compare models pseudoproxies with proxies
#Set up data to compare
names_key = ['Category', 'CategorySpecific','TSid', 'Lat', 'Lon']
variable_key = []
for time in setTime:
    for variable in ['Slope','SlopeSig']: 
        variable_key.append(time+variable)
for time in range(1,12): 
        variable_key.append('bin'+str(time)+'ka')

#Compare
def dataAgreeSign (proxyDF,modelDF,metric):
    output={}
    for name in ['Category', 'CategorySpecific','TSid', 'Lat', 'Lon']:
        output[name] = proxyDF[name]
    #Make sure comparing same sites
    if sum(proxyDF['Lat']-modelDF['trace']['Lat']) != 0: print('ERROR!!! '); output[name]=[]
    #
    #Get values and sign
    proxyVal  = proxyDF[metric];          proxySign = np.sign(proxyVal)
    traceVal  = modelDF['trace'][metric]; traceSign = np.sign(traceVal)
    hadcmVal  = modelDF['hadcm'][metric]; hadcmSign = np.sign(hadcmVal)
    #Compare agreement
    output['proxy_trace'] = proxySign*traceSign
    output['proxy_hadcm'] = proxySign*hadcmSign
    output['model_model'] = hadcmSign*traceSign
    output['agree_count'] = ((output['proxy_trace']+1)+(output['proxy_hadcm']+1))/2 #0=both models disagree, 1 = 1 model agrees, 2 = both models agree
    output['agree_countSign'] = output['agree_count']*proxySign
    output = pd.DataFrame.from_dict(output)
    return(output)

def agreementByCat(df,var,sortBy='Category',agreement='general'):
    output = {}
    df_dropna = df.dropna(subset=[var])
    if agreement == 'general': 
        output['all'] = round((list(df_dropna[var]).count(1)/len(df_dropna[var]))*100,)    
    else: 
        try: output['all'] = round((list(df_dropna[var]).count(1)/(list(df_dropna[var]).count(1)+list(df_dropna[var]).count(-1)))*100,)    
        except: output['all'] = np.NaN
    #
    for category in df[sortBy].unique(): 
        df_select = df.loc[(df[sortBy] == category)]
        df_dropna = df_select.dropna(subset=[var])
        if agreement == 'general': 
            try: output[category] = round((list(df_dropna[var]).count(1)/len(df_dropna[var]))*100,)    
            except: output[category] = np.NaN    
        else: 
            try: output[category] = round((list(df_dropna[var]).count(1)/(list(df_dropna[var]).count(1)+list(df_dropna[var]).count(-1)))*100,)    
            except: output[category] = np.NaN    
    return(output)


agreeOut={}
for model in ['proxy_hadcm','proxy_trace','model_model']:
    agreeOut[model]={}
    for variable in variable_key: 
        agreeOut[model][variable]={}
        for method in modelProxies.keys():
            modelProxyAgreement = dataAgreeSign(data_HC,modelProxies[method],variable)
            modelProxyAgreeCat  = agreementByCat(modelProxyAgreement,model,sortBy='CategorySpecific',agreement='general')
            agreeOut[model][variable][method] = modelProxyAgreeCat
        agreeOut[model][variable] = pd.DataFrame(agreeOut[model][variable])

#%%
var1 = 'LateSlope'
var2 = 'EarlySlope'



agreeMatrix= np.mean([agreeOut['proxy_hadcm'][var1],agreeOut['proxy_trace'][var1],
                                    agreeOut['proxy_hadcm'][var2],agreeOut['proxy_trace'][var2]],axis=0)

my_norm = mpl.colors.Normalize(30,70) 
my_cmap = plt.cm.get_cmap('RdBu')
color_vals = my_cmap(my_norm(agreeMatrix))

#plt.title('Agreement by proxy type & category with Psuedoproxies (avg trace/hadcm early/late trends)',loc='left')
plt.table(cellText=agreeMatrix,
          colLabels=agreeOut['proxy_hadcm'][var1].columns,
          rowLabels=agreeOut['proxy_hadcm'][var1].index,
          cellColours=color_vals,
          loc='center').scale(1.5, 2)
plt.axis('off')
plt.show()


agreeDiff = agreeMatrix.transpose()
agreeDiff = agreeDiff-agreeDiff[0]
agreeDiff = agreeDiff.transpose()
my_norm = mpl.colors.Normalize(-10,10) 
my_cmap = plt.cm.get_cmap('RdBu')
color_vals = my_cmap(my_norm(agreeDiff))

plt.table(cellText=agreeDiff,
          colLabels=agreeOut['proxy_hadcm'][var1].columns,
          rowLabels=agreeOut['proxy_hadcm'][var1].index,
          cellColours=color_vals,
          loc='center').scale(1.5, 2)
#plt.title('Agreement by proxy type & category with Psuedoproxies (avg trace/hadcm early/late trends)')

plt.axis('off')
plt.show()
#%%
for variable in ['EarlySlopeSig','LateSlopeSig','bin6ka']: #plot proxy values to check calculations
    modelProxyAgreement = dataAgreeSign(data_HC,modelProxies['HC_Proxy'],variable)    
    plot_df = modelProxyAgreement; plot_df = plot_df[plot_df['agree_countSign'].notna()]
    plt.style.use('ggplot')
    plt.figure(figsize=(20,10)); plt.rcParams['axes.facecolor'] ='white'
    plt.rcParams['axes.linewidth'] = 1; plt.rcParams['axes.edgecolor'] = 'k'
    ax1 = plt.subplot(projection=ccrs.Robinson()) 
    ax1.spines['geo'].set_edgecolor('black')
    ax1.set_global(); ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
    ax1.coastlines(); ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
    proxy_scatter = ax1.scatter(plot_df.Lon,plot_df.Lat,c=plot_df['agree_countSign'],
            marker='o',s=130,edgecolor='k',lw=1,alpha=0.8,transform=ccrs.PlateCarree(),
            cmap='BrBG',vmin=-2,vmax=2)
    plt.title(variable,fontsize=30)
    plt.show()

#%%
FigKeyCategory={'Speleothem':   {'s':1,  'marker':'^', 'c':'firebrick'}, 
        'Speleothem (d18O)':    {'s':1,  'marker':'^', 'c':'firebrick'},
        'Speleothem (d13C)':    {'s':1,  'marker':'v', 'c':'indianred'},
        'Speleothem (other)':   {'s':1,  'marker':'>', 'c':'r'},
        'Lake Deposits':        {'s':1.5,'marker':'.', 'c':'cornflowerblue'},
        'Glacier Ice':          {'s':1,  'marker':'D', 'c':'powderblue'},
        'Lake Sediment (d18O)': {'s':1,  'marker':'o', 'c':'darkblue'},
        'Leaf Wax (dD)':        {'s':1,  'marker':'p', 'c':'darkorchid'}, 
        'Pollen':               {'s':1,  'marker':'P', 'c':'forestgreen'},  
        'Pollen (calibrated)':  {'s':1,  'marker':'P', 'c':'seagreen'},    
        'Pollen (uncalibrated)':{'s':1,  'marker':'X', 'c':'olive'}, 
        'Other':                {'s':1,  'marker':'s', 'c':'grey'}, 
        'Other (calibrated)' :  {'s':1,  'marker':'s', 'c':'dimgrey'},  
        'Other (uncalibrated)': {'s':1,  'marker':'s', 'c':'lightgrey'}
        }  

FigKeyExtent = {'Global':{'pltExtent':[0,8,0,10],'GeoExtent':[-180,180,-90,90],'Name':'Global'},
                 'NA':    {'pltExtent':[8,12,2,6],'GeoExtent':[-130,-50,23,55], 'Name':'North America'},
                 'Eur':   {'pltExtent':[8,12,6,9],'GeoExtent':[-12,40,35,70],    'Name':'Europe'}
                 }

for variable in ['EarlySlope','LateSlope','bin6ka']:
    modelProxyAgreement = dataAgreeSign(data_HC,modelProxies['HC_Proxy'],variable)    
    plt.style.use('ggplot')
    plt.figure(figsize=(20,10))
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['axes.facecolor'] ='white'; plt.rcParams['axes.linewidth'] = 1; plt.rcParams['axes.edgecolor'] = 'k'
    gs = gridspec.GridSpec(12,12,wspace=0,hspace=0.2)
    for ax in FigKeyExtent: 
        extent = FigKeyExtent[ax]['pltExtent']
        ax1 = plt.subplot(gs[extent[0]:extent[1],extent[2]:extent[3]],projection=ccrs.PlateCarree())
        ax1.set_extent(FigKeyExtent[ax]['GeoExtent'])
        ax1.coastlines()
        ax1.spines['geo'].set_edgecolor('black')
        ax1.add_feature(cfeature.LAND, facecolor='whitesmoke',edgecolor='k')
        ax1.add_feature(cfeature.LAKES,facecolor='none',      edgecolor='k')
        name = 'CategorySpecific'
        if ax == 'Global': 
              ax1.add_feature(cfeature.BORDERS,facecolor='none',edgecolor='grey',lw=0.2)
              plt.title(variable)
        elif ax == 'NA':
              ax1.add_feature(cfeature.BORDERS,facecolor='none',edgecolor='k',lw=0.3)        
              ax1.add_feature(cfeature.STATES,facecolor='none',edgecolor='grey',lw=0.2)        
        else: ax1.add_feature(cfeature.BORDERS,facecolor='none',edgecolor='k',lw=0.3)            
        for category in modelProxyAgreement[name].unique(): 
            plot_df = modelProxyAgreement.loc[(modelProxyAgreement[name] == category)]
            proxy_scatter = ax1.scatter(plot_df.Lon,plot_df.Lat,c=plot_df['agree_countSign'],
                edgecolor='k',lw=0.5,alpha=0.8,transform=ccrs.PlateCarree(),cmap='BrBG',vmin=-2,vmax=2,
                marker=FigKeyCategory[category]['marker'],s=FigKeyCategory[category]['s']*100)
    if save: plt.savefig(gitHub+'Figures/HC12k_'+name+'.png', dpi=400,format='png')
    else: plt.show()
#%%