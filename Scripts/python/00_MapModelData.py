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

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 0.5; 
plt.rcParams['axes.edgecolor'] = 'k'
plt.tick_params(labelsize=8)

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
var = 'tas'
if var == 'tas': 
    cramp, units = 'RdBu_r','degC'
    mlevels = np.array([i /100 for i in list(range(-220,221,40))])
else: 
    cramp, units = 'BrBG', 'mm/day'
    mlevels = np.array([i /100 for i in list(range(-55,56,10))])
    
ka = [0.5,6]
models  = ['cmip6']#['trace','hadcm','cmip6']
seasons = ['ANN','JJA','DJF'] 
modelAnom = {}
s = 3

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
plt.figure(figsize=(2.15*len(models),5),dpi=400); 

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
                                     levels=mlevels,extend='both',cmap=cramp)
        plt.title(model+' '+szn,fontsize=8)
ax1 = plt.subplot(gs[len(seasons)*s:len(seasons)*s+1,0:len(models)*s])
ax1.axis('off')
plt.colorbar(model_contour,cax=inset_axes(ax1,width='90%',height="20%",loc="center"),
            orientation="horizontal").set_label(units,fontsize=8,c='black')
plt.title(var+' (Mid-Holocone - PI)',fontsize=10)
#Save or show
if save: plt.savefig(dataDir+'Figures/Model/midHoloceneModeledAnoms_'+var+'.png',
                     dpi=400,format='png',bbox_inches='tight')       
plt.show()



#%% 5 Transient Trends
#Settings
save = True
var = 'tas'
if var == 'tas': 
    cramp = 'RdBu_r'
    units = 'degC'
    mlevels = np.array([i /100 for i in list(range(-110,111,20))])
else: 
    cramp = 'BrBG'
    units = 'mm/day'
    mlevels = np.array([i /100 for i in list(range(-11,12,1))])
ka = [0,6,12]
models  = ['trace','hadcm']
szn = 'ANN'
modelTrends = {}
s = 3

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
                                     levels=mlevels,extend='both',cmap=cramp)
        plt.title(model+'_'+str(ka[t+1])+'to'+str(ka[t]),fontsize=8)
ax1 = plt.subplot(gs[6:7,0:len(models)*s])
ax1.axis('off')
plt.colorbar(model_contour,cax=inset_axes(ax1,width='90%',height="20%",loc="center"),
            orientation="horizontal").set_label(units+'/ka',fontsize=8,c='black')
plt.tick_params(labelsize=8)
#plt.set_ticklabels(mlevels)
plt.title(var+' '+szn+' Holocene Trends',fontsize=10)
if save: plt.savefig(dataDir+'Figures/Model/HoloceneTrends_'+var+'_'+szn+'.png',
                     dpi=400,format='png',bbox_inches='tight')       
else: plt.show()

#%% 6 Plot correlation between 2 transient models? 


#%% 7 When is the wettest century
save = False
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
