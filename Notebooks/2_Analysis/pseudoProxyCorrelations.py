
import numpy  as np                     # Package with useful numerical functions
import xarray as xr                     #Package for handling file types
import pandas as pd
import regionmask
from scipy import stats                 # Packages for calculations  
from scipy.stats import pearsonr
import random
import matplotlib.pyplot   as plt


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

save=True



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
            lon = regData['longitude'][i]
            if lon < 0: lon = 360+lon
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



#%%
proxyVar    = 'HC'
modelVar    = 'pre'
season      = 'proxy'
standardize = True


h_min=[]
h_med=[]
h_max=[]
h_mean=[]
t_min=[]
t_med=[]
t_max=[]
t_mean=[]
for n in range(1,26):
    print(n)
    regions = []
    siteNum = []
    rvalueT = []
    rvalueH = []
    rvalMax = []
    rvalAvg = []
    pctP = []
    pctAnn = []
    for repeat in range(50):
        for reg in np.unique(proxyDF['ipccReg']):
            if reg not in regionmask.defined_regions.ar6.land.abbrevs: continue
            if reg == 'WAN': continue
            if reg == 'EAN': continue
            #Select proxy sites for region
            regData  = proxyDF.loc[proxyDF['ipccReg']==reg]
            if len(regData) >= n:
                regData=regData.sample(n)
            else: continue
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
            #rvalMax.append(max([corrT,corrH]))
            #rvalAvg.append(np.mean([corrT,corrH]))
            #pctP.append(p_n/count)
            #pctAnn.append(ann_n/count)
    #df = pd.DataFrame([regions,siteNum,rvalueT,rvalueH,pctP,pctAnn]).transpose()
    h_min.append(float(np.nanquantile(rvalueH,[0.25])))
    h_med.append(float(np.nanquantile(rvalueH,[0.5])))
    h_max.append(float(np.nanquantile(rvalueH,[0.75])))
    h_mean.append(np.nanmean(rvalueH))
    t_min.append(float(np.nanquantile(rvalueT,[0.25])))
    t_med.append(float(np.nanquantile(rvalueT,[0.5])))
    t_max.append(float(np.nanquantile(rvalueT,[0.75])))
    t_mean.append(np.nanmean(rvalueT))
    #df.to_csv(str(dataDir+'Data/Model/'+'pseudoProxyCorr'+'_'+modelVar+'.csv'))
#%%
df2 = pd.DataFrame([range(1,26),t_mean,h_mean,t_med,h_med]).transpose()
df2.to_csv(str(dataDir+'Data/Model/'+'pseudoProxyCorr'+'_'+modelVar+'_byCount.csv'))
#plt.plot(range(1,26),h_min,color='orange')
#plt.plot(range(1,26),h_med,color='orange')
#plt.plot(range(1,26),h_max,color='orange')
plt.plot(range(1,26),h_mean,color='orange',label='HadCM')
#plt.plot(range(1,26),t_min,color='red')
#plt.plot(range(1,26),t_med,color='red')
#plt.plot(range(1,26),t_max,color='red')
plt.plot(range(1,26),t_mean,color='red',label='TraCE')
plt.legend()
plt.axvline(x=6)
plt.ylabel('correlation coefficient')
plt.xlabel('number of proxy sites')
#plt.ylim([0,1])
plt.xlim([0,25])
plt.show()

#%%
szn='ANN'
model= 'trace'
var='p-e'

loc = 'N_Africa'
if loc == 'N_Africa':
    lat = 20
    lon = 0
elif loc == "C_Asia":
   lat = 45
   lon = 60
elif loc == "E_Asia1":
    lat = 30
    lon = 110
elif loc == "E_Asia2":
    lat = 30
    lon = 90
elif loc == "WNA":
   lat = 40
   lon = -115
   
if lon < 180: lon = 360+lon

lat = np.argmin(np.abs(modelData[model][szn]['lat'].data-lat))
lon = np.argmin(np.abs(modelData[model][szn]['lon'].data-lon))
modelvalues = modelData[model][szn][var][:,lat,lon]

#print(modelvalues)

data = np.flip(stats.zscore(modelvalues).data)
ages = np.flip(stats.zscore(modelvalues).age.data)

data2=np.cumsum(data)


vals=np.mean(modelData[model][szn][var][55:66,:,:],axis=0)-np.mean(modelData[model][szn][var][0:12,:,:],axis=0)

plt.figure(figsize=(6,6))
gs = gridspec.GridSpec(12,8)
ax = plt.subplot(gs[7:12,0:8],projection=ccrs.Robinson()) 
data_cyclic,lon_cyclic = cutil.add_cyclic_point(vals,coord=vals.lon)
model_contour = plt.contourf(lon_cyclic,vals.lat, data_cyclic,transform=ccrs.PlateCarree(),
                             vmin=-0.5,vmax=0.5,cmap='BrBG',levels=20)  
ax.scatter(modelvalues.lon.data,modelvalues.lat.data,c='red',transform=ccrs.PlateCarree(),s=100,ec='k',lw=2)
ax.set_global()
ax.coastlines()
ax1 = plt.subplot(gs[0:6,0:8]) 
ax1.plot(ages,stats.zscore(data),label='P-E')
ax1.plot(ages,stats.zscore(data2),label='cumsum(P-E)')
ax1.legend()
ax1.invert_xaxis()
plt.title("Z-Scores [Trace-21ka 100yr Average Annual P-E]")
plt.show()


data2 = []
for i in range(len(data)):
    if i == 0:
        data2.append(data[i])
    else: data2.append(data2[i-1]+(data[i]-data[i-1]))