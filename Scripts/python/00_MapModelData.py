#%% 1 Load Packages
import numpy      as np
import pandas     as pd
import regionmask as rm
import xarray     as xr
from scipy        import stats
from scipy.stats  import pearsonr
#For Figures:
import matplotlib.pyplot as plt         # Packages for making figures
import matplotlib.gridspec as gridspec
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs              # Packages for mapping in python
import cartopy.feature as cfeature
import cartopy.util as cutil

#%% 2 Load Data
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate/'
modelData = {}
for model in ['hadcm','trace','cmip6']:
    modelData[model] = {}
    for szn in ['ANN','JJA','DJF']:
        modelData[model][szn] = xr.open_dataset(dataDir+'Data/Model/'+model+'_'+szn+'.nc',decode_times=False)

#%% 3 MH anoms
#Settings
save = False
var = 'p-e'
ka = [0.5,6]
models  = ['trace','hadcm','cmip6']
seasons = ['ANN','JJA','DJF'] 
modelAnom = {}
s = 3
mlevels = np.array([i /100 for i in list(range(-55,56,10))])
#Calculate anoms
for model in ['hadcm','trace']:
    for szn in seasons:
        data = modelData[model][szn][var]
        data0 = data.groupby_bins('time',[ka[0]*1000-500,ka[0]*1000+500]).mean(dim='time')
        data1 = data.groupby_bins('time',[ka[1]*1000-500,ka[1]*1000+500]).mean(dim='time')
        modelAnom[model+'_'+szn] = data0
        modelAnom[model+'_'+szn].data = data1.data-data0.data
        modelAnom[model+'_'+szn] = modelAnom[model+'_'+szn].squeeze('time_bins')
for model in ['cmip6']:
    for szn in seasons:
        data = modelData[model][szn][var+'_regrid']
        data = data[0,:,:,:] - data[1,:,:,:] 
        data = data.rename({'lat_regrid':'lat','lon_regrid':'lon'})
        modelAnom[model+'_'+szn] = np.mean(data,axis=0)
#Plot Figure
plt.style.use('ggplot')
plt.figure(figsize=(6.5,5),dpi=400); 
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 0.5; 
plt.rcParams['axes.edgecolor'] = 'k'
gs = gridspec.GridSpec(len(seasons)*s+1,len(models)*s)
for model in models:
    for szn in seasons:
        data = modelAnom[model+'_'+szn]
        i = [x for x, val in enumerate(seasons) if val == szn][0]
        j = [x for x, val in enumerate(models) if val == model][0]
        ax1 = plt.subplot(gs[i*s:i*s+s,j*s:j*s+s],projection=ccrs.Robinson()) 
        ax1.spines['geo'].set_edgecolor('black')
        ax1.set_global()
        ax1.coastlines()
        ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
        ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
        data_cyclic,lon_cyclic = cutil.add_cyclic_point(data,coord=data.lon.data)
        model_contour = plt.contourf(lon_cyclic, data.lat.data, 
                                     data_cyclic,transform=ccrs.PlateCarree(),
                                     levels=mlevels,extend='both',cmap='BrBG')
        plt.title(model+' '+szn,fontsize=8)
ax1 = plt.subplot(gs[len(seasons)*s:len(seasons)*s+1,0:len(models)*s])
ax1.axis('off')
plt.colorbar(model_contour,cax=inset_axes(ax1,width='90%',height="20%",loc="center"),
            orientation="horizontal").set_label('mm/day',fontsize=8,c='black')
plt.tick_params(labelsize=8)
plt.title(var+' Mid-Holocone (6ka) Anomalies Relative to PI (0.5ka)',fontsize=10)
#Save or show
if save: plt.savefig(dataDir+'Figures/midHoloceneModeledAnoms_'+var+'.png',
                     dpi=400,format='png',bbox_inches='tight')       
else: plt.show()

#%% 4 Plot Agreement
save = False
var = var
szn='ANN'
#Calculate sign of each gridcell and flatten into 2d
data = modelData['cmip6'][szn][var+'_regrid']
data = data[0,:,:,:] - data[1,:,:,:] 
data = data.rename({'lat_regrid':'lat','lon_regrid':'lon'})
data = np.sum(np.sign(data),axis=0)
for model in ['hadcm','trace']:
    add = modelAnom[model+'_'+szn]
    for i in range(len(data.lat)):
        for j in range(len(data.lon)): 
            lat=np.argmin(np.abs(add.lat.data-data.lat[i].data))
            lon=np.argmin(np.abs(add.lon.data-data.lon[i].data))
            data[i,j] += np.sign(add[lat,lon]) 
#Calculate %
dataPct = (((data+14)/2)/14)*100
#Load Proxy Data and get geographic position of regionmask labels
refReg = regionmask.defined_regions.ar6.land
proxy = pd.read_csv(dataDir+'Data/midHCpct.csv')[['V1','V2','V3']]
plats = []
plons = []
for i in proxy['V1']: 
    loc = refReg.centroids[refReg.abbrevs.index(i)]
    plats.append(loc[1])
    plons.append(loc[0])
#Plot
plt.style.use('ggplot')
plt.figure(figsize=(6,4)); 
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams['font.family'] = 'Arial'
ax1 = refReg.plot(projection=ccrs.Robinson(),label="abbrev", add_ocean=True, 
                  text_kws=dict(bbox=dict(color="none"),fontsize=0),
                  line_kws=dict(linewidth=0.8))
ax1.spines['geo'].set_edgecolor('black')
ax1.set_global()
ax1.coastlines()
ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
model_contour = plt.pcolormesh(dataPct.lon, dataPct.lat, 
                    dataPct,transform=ccrs.PlateCarree(),
                    cmap='BrBG',alpha=0.8)
proxy_scatter = ax1.scatter(plons,plats,c=list((proxy['V3']/proxy['V2'])*100),
                            s=list(proxy['V2']**0.5*20),edgecolor='k',
                            lw=0.8,alpha=1,cmap=plt.cm.get_cmap('BrBG',5),
                            transform=ccrs.PlateCarree(),vmin=0,vmax=100)
plt.title('Model and Proxy agreement for sign of 6ka-0.5ka difference',fontsize=8)
plt.colorbar(model_contour,cax=inset_axes(ax1,width='70%',height="4%",loc="lower center"),
             orientation="horizontal").set_label('% positive (wet '+var+') mid-Holocene anomaly',fontsize=8,c='black')
plt.tick_params(labelsize=8)
#Save or show
if save: plt.savefig(dataDir+'Figures/midHoloceneAgreement_'+var+'.png',
                     dpi=400,format='png',bbox_inches='tight')       
else: plt.show()



#%% 5 Transient Trends
#Settings
save = False
var = 'p-e'
ka = [0,6,12]
models  = ['trace','hadcm']
szn = 'ANN'
modelTrends = {}
s = 3
mlevels = np.array([i /100 for i in list(range(-11,12,1))])
#Calculate model trends
for model in ['hadcm','trace']:
    for t in range(len(ka)-1):
        data = modelData[model][szn][var]
        i = [x for x, val in enumerate(data['time']) if 
             val >= 1000*ka[t] and val <= 1000*ka[t+1]] 
        data = data[i,:,:]
        dataSl = np.full(np.shape(data)[1:3],np.NaN)
        for lat in range(len(data.lat)):
            for lon in range(len(data.lon)):
                modelReg = stats.linregress(data.time.data,data.data[:,lat,lon])
                dataSl[lat,lon] = modelReg[0]*-1000
        modelTrends[model+'_'+str(ka[t+1])+'to'+str(ka[t])] = dataSl
#Plot FIugres
plt.style.use('ggplot')
plt.figure(figsize=(6.5,5),dpi=400); 
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 0.5; 
plt.rcParams['axes.edgecolor'] = 'k'
gs = gridspec.GridSpec((len(ka)-1)*s+1,len(models)*s)
for model in models:
    lats = modelData[model][szn][var].lat.data
    lons = modelData[model][szn][var].lon.data
    for t in range(len(ka)-1):
        print(model+'_'+str(ka[t+1])+'to'+str(ka[t]))
        data = modelTrends[model+'_'+str(ka[t+1])+'to'+str(ka[t])] 
        i = [x for x, val in enumerate(ka) if val == ka[t]][0]
        j = [x for x, val in enumerate(models) if val == model][0]
        ax1 = plt.subplot(gs[i*s:i*s+s,j*s:j*s+s],projection=ccrs.Robinson()) 
        ax1.spines['geo'].set_edgecolor('black')
        ax1.set_global()
        ax1.coastlines()
        ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
        ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
        data_cyclic,lon_cyclic = cutil.add_cyclic_point(data,coord=lons)
        model_contour = plt.contourf(lon_cyclic, lats, 
                                     data_cyclic,transform=ccrs.PlateCarree(),
                                     levels=mlevels,extend='both',cmap='BrBG')
        plt.title(model+'_'+str(ka[t+1])+'to'+str(ka[t]),fontsize=8)
ax1 = plt.subplot(gs[6:7,0:len(models)*s])
ax1.axis('off')
plt.colorbar(model_contour,cax=inset_axes(ax1,width='90%',height="20%",loc="center"),
            orientation="horizontal").set_label('mm/day/ka',fontsize=8,c='black')
plt.tick_params(labelsize=8)
#plt.set_ticklabels(mlevels)
plt.title(var+' '+szn+' Holocene Trends',fontsize=10)
if save: plt.savefig(dataDir+'Figures/HoloceneTrends_'+var+'_'+szn+'.png',
                     dpi=400,format='png',bbox_inches='tight')       
else: plt.show()

#%% 6 Plot correlation between 2 transient models? 
times = [6,12]
var = 'p-e'
szn = 'ANN'
#Load data
dataHa = modelData['hadcm'][szn][var]
dataTr = modelData['trace'][szn][var]
modelR = np.full(np.shape(dataTr)[1:3],np.NaN)
land = rm.defined_regions.natural_earth.land_110.mask_3D(data.lon,data.lat).squeeze('region').data
t0 = np.argmin(np.abs(times[0]*1000 - dataHa.time.data))   
t1 = np.argmin(np.abs(times[1]*1000 - dataHa.time.data))  
#reGrid 
for lat in range(len(dataTr.lat)):
    for lon in range(len(dataTr.lon)):
        #if land.data[lat,lon] == True:
        i = np.argmin(np.abs(dataTr.lat.data[lat] - dataHa.lat.data))                  
        j = np.argmin(np.abs(dataTr.lon.data[lon] - dataHa.lon.data))
        modelR[lat,lon] = pearsonr(dataTr[t0:t1,lat,lon].data,dataHa[t0:t1,i,j].data)[0]
        
mlevels = np.array([i /10 for i in list(range(-10,11,2))])

plt.style.use('ggplot')
plt.figure(figsize=(6,4))
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.edgecolor'] = 'k'
gs = gridspec.GridSpec(6,2)
ax1 = plt.subplot(gs[0:6,0:2],projection=ccrs.Robinson()) 
ax1.spines['geo'].set_edgecolor('black')
ax1.set_global()
ax1.coastlines()
ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k') #dataSlP= dataSl
data_cyclic,lon_cyclic = cutil.add_cyclic_point(modelR,coord=dataTr.lon.data)
model_contour = plt.contourf(lon_cyclic, dataTr.lat.data, 
                             data_cyclic,transform=ccrs.PlateCarree(),
                             levels=mlevels,cmap='RdBu_r') 
plt.title('Correlation between hadcm/trace ('+var+' '+szn+') ('+str(times[0])+'-'+str(times[1])+'ka)',fontsize=12)
plt.colorbar(model_contour,cax=inset_axes(ax1,width='60%',height="4%",loc=8),
             orientation="horizontal").set_label('r-value (pearsons correlation)',fontsize=12,c='black')
plt.show()
#%% 
mlevels = np.array([-1,-2/3,-1/3,-0.1,0.1,1/3,2/3,1])

for climVar in ['p-e_ANN']:#,'pre_JJA','p-e_ANN','p-e_JJA']:#,'pre_DJF','pre_JJA']:
    var  = climVar[:3]
    szn  = climVar[len(climVar)-3:]
    dataTr = modelData['trace'][szn][var].groupby_bins('time',range(0,6001,100)).mean(dim="time")
    dataHa = modelData['hadcm'][szn][var].groupby_bins('time',range(0,6001,100)).mean(dim="time")
    hadaRg = np.full(np.shape(dataTr)[1:3],np.NaN)
    land   = regionmask.defined_regions.natural_earth.land_110
    land   = land.mask_3D(dataTr.lon,dataTr.lat)  
    for i in range(np.shape(dataTr)[1]):
        for j in range(np.shape(dataTr)[2]):
            if land.data[0,i,j] == True:
                lat = np.argmin(np.abs(dataTr.lat.data[i] - dataHa.lat.data))                  
                lon = np.argmin(np.abs(dataTr.lon.data[j] - dataHa.lon.data))
                hadaRg[i,j] = pearsonr(dataTr[:,i,j].data,dataHa[:,lat,lon].data)[0]
    plt.style.use('ggplot')
    plt.figure(figsize=(12,20))
    plt.rcParams['axes.facecolor'] ='white'
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['axes.edgecolor'] = 'k'
    plt.title(climVar,fontsize=12)
    gs = gridspec.GridSpec(6,2)
    ax1 = plt.subplot(gs[n-1:n,0:1],projection=ccrs.Robinson()) 
    ax1.spines['geo'].set_edgecolor('black')
    ax1.set_global()
    ax1.coastlines()
    ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
    ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
    #dataSlP= dataSl
    data_cyclic,lon_cyclic = cutil.add_cyclic_point(hadaRg,
                                                            coord=dataTr.lon.data)
    model_contour = plt.contourf(lon_cyclic, dataTr.lat.data, 
                                         data_cyclic,transform=ccrs.PlateCarree(),
                                         levels=mlevels,extend='both',cmap='BrBG') 
    plt.title('agreement between hadcm/trace '+climVar,fontsize=12)
    plt.colorbar(model_contour,cax=inset_axes(ax1,width='60%',height="4%",loc=8),ticks=range(0,13),
                 orientation="horizontal").set_label('ka',fontsize=12,c='black')
    plt.show()
print(np.nanmean(hadaRg))

#%% 7 When is the wettest century
save = True
models = ['trace','hadcm']
var = 'pre'
szn = 'ANN'
#Locate wettest century for each models
wetPeak = {}
for model in models:
    data = modelData[model][szn][var]
    wetPeak[model+'_all']  = np.full(np.shape(data)[1:3],np.NaN)
    wetPeak[model+'_land'] = np.full(np.shape(data)[1:3],np.NaN)
    land = rm.defined_regions.natural_earth.land_110.mask_3D(data.lon,data.lat).squeeze('region').data
    for lat in range(len(data.lat)):
        for lon in range(len(data.lon)):
            dataTS = data[0:121,lat,lon].data
            n = [i for i, j in enumerate(dataTS) if j == np.nanmax(dataTS)]
            wetPeak[model+'_all'][lat,lon] = data[:,lat,lon].time[n].data[0]
            if land[lat,lon]: 
                wetPeak[model+'_land'][lat,lon] = wetPeak[model+'_all'][lat,lon]
plt.style.use('ggplot')
plt.figure(figsize=(3.5,3.5),dpi=400)
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 0.5; 
plt.rcParams['axes.edgecolor'] = 'k'
plt.title('Distribution of wettest century in model data')
gs = gridspec.GridSpec(len(models),1)
for model in models:
    n = [i for i, j in enumerate(models) if j == model][0]
    ax = plt.subplot(gs[n:n+1,0:1])
    ax.hist(wetPeak[model+'_all'].flatten(),  bins=12, alpha=0.6, density=False,
            label='All Gridcells', color = 'royalblue', edgecolor='k',lw= 1)
    ax.hist(wetPeak[model+'_land'].flatten(), bins=12, alpha=0.6, density=False,
            label='Land Only', color = 'olivedrab', edgecolor='k', lw= 1)
    if n == 0: 
        ax.legend(loc='upper right',fontsize=8)
        ax.set_title('Distribution of wettest century in model data ('+var+' '+szn+')', 
                     y=1.1, x = 0.42, fontsize = 10)
        ax.tick_params(labelbottom=False) 
    else:  plt.xlabel('Age (yr BP)',fontsize=8)
    plt.ylabel(model,fontsize=8)
    ax.tick_params(labelsize=8)
    ax.set_xlim([12000,0])
    ax.set_ylim([0,2000])

if save: plt.savefig(dataDir+'Figures/HistWettestCentury_'+var+'_'+szn+'.png',
                     dpi=400,format='png',bbox_inches='tight')       
else: plt.show()

#%%
