#%% 1 Load Packages
import numpy      as np
import pandas     as pd
import regionmask as rm
import xarray     as xr 
from scipy        import stats
#from scipy.stats  import pearsonr
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
#plt.tick_params(labelsize=8)

#%% 2 Load Data
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate/'
modelData = {}
for model in ['hadcm','trace','cmip6']:
    modelData[model] = {}
    for szn in ['ANN','JJA','DJF']:
        modelData[model][szn] = xr.open_dataset(dataDir+'Data/Model/'+model+'_'+szn+'.nc',decode_times=False)


#%%

data = modelData['cmip6']['ANN']
land = rm.defined_regions.natural_earth.land_110.mask_3D(data.lon_regrid,data.lat_regrid)
land = land.squeeze('region').data
pd.DataFrame(land).to_csv(dataDir+'Data/Model/cmip6LandMask.csv')

#%% 3 MH anoms
#Settings
var = 'p-e'
if var == 'tas': 
    cramp, units = 'RdBu_r','degC'
    mlevels = np.array([i /100 for i in list(range(-220,221,40))])
else: 
    cramp, units = 'BrBG', 'mm/day'
    mlevels = np.array([i /100 for i in list(range(-75,76,10))])
    
seasons = ['ANN','JJA','DJF'] 
s = 3
#Plot Figure
save = True
plt.style.use('ggplot')
plt.figure(figsize=(3.25,6),dpi=600); 
h=len(seasons)
gs = gridspec.GridSpec(h*s+1,s)
for szn in seasons:
    data = modelData['cmip6'][szn][var+'_regrid']
    data = data[0,:,:,:] - data[1,:,:,:] 
    data = data.rename({'lat_regrid':'lat','lon_regrid':'lon'})
    data = np.mean(data,axis=0)
    i = [x for x, val in enumerate(seasons) if val == szn][0]
    ax = plt.subplot(gs[(i)*s:(i)*s+s,0:s],projection=ccrs.Robinson()) 
    ax.spines['geo'].set_edgecolor('black')
    ax.set_global()
    ax.coastlines()
    ax.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
    data_cyclic,lon_cyclic = cutil.add_cyclic_point(data,coord=data.lon.data)
    model_contour = plt.contourf(lon_cyclic, data.lat.data, 
                                 data_cyclic,transform=ccrs.PlateCarree(),
                                 levels=mlevels,extend='both',cmap=cramp)
    ax.annotate(szn, xy=(0, 0), xycoords='data', xytext=(-0.05, 0.5), 
                textcoords='axes fraction', fontsize=8, fontfamily = 'Arial',
                rotation=90, horizontalalignment='center', verticalalignment='center')  
    if szn == 'ANN':
        ax.annotate('CMIP Ensemble Mean ('+var+')',
                    xy=(0, 0), xycoords='data',fontsize=8, xytext=(0.5, 1.1), 
                    textcoords='axes fraction',horizontalalignment='center', verticalalignment='center')  
ax = plt.subplot(gs[(h*s):(h*s+1),0:s])
ax.axis('off')
cbar = plt.colorbar(model_contour,orientation="horizontal",
                    cax=inset_axes(ax,width='90%',height="30%",loc="upper center"),
            ticks=[-0.6,-0.3,0,0.3,0.6]).set_label('Mid-Holocene - Pre-Industrial ('+units+')',
                                                   fontsize=8,c='black')
plt.tick_params(labelsize=8)
#Save or show
if save: plt.savefig(dataDir+'Figures/Model/Anomalies/CMIP_MH-PI_bySeason_'+var+'.png',
                     dpi=600,format='png',bbox_inches='tight')       
else: plt.show()


#%%transient anoms
ka = [0,6]
models  = ['hadcm','trace']#,'cmip6']
seasons = ['ANN','JJA','DJF'] 
modelAnom = {}
s = 3

#Plot Figure
save = True
plt.style.use('ggplot')
plt.figure(figsize=(3.25*len(models),6),dpi=400); 
h=len(seasons)*(len(ka)-1)
gs = gridspec.GridSpec(h*s+1,len(models)*s)
for t in range(0,len(ka)-1):
    for model in models:
        for szn in seasons:
            if model in ['trace','hadcm']:
                data = modelData[model][szn][var]
                data0 = data.groupby_bins('time',[ka[t]*1000-500,ka[t]*1000+500]).mean(dim='time')
                data1 = data.groupby_bins('time',[ka[t+1]*1000-500,ka[t+1]*1000+500]).mean(dim='time')
                data = data0
                data.data = (data1.data-data0.data)
                data = data.squeeze('time_bins')
            else:
                data = modelData[model][szn][var+'_regrid']
                data = data[0,:,:,:] - data[1,:,:,:] 
                data = data.rename({'lat_regrid':'lat','lon_regrid':'lon'})
                data = np.mean(data,axis=0)
            i = [x for x, val in enumerate(seasons) if val == szn][0]
            j = [x for x, val in enumerate(models) if val == model][0]
            ax = plt.subplot(gs[(i+t)*s:(i+t)*s+s,j*s:j*s+s],projection=ccrs.Robinson()) 
            ax.spines['geo'].set_edgecolor('black')
            ax.set_global()
            ax.coastlines()
            plt.ylabel('\n'+'TraCE\n(% of grid cells)',fontsize=8)
            ax.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
            ax.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
            data_cyclic,lon_cyclic = cutil.add_cyclic_point(data,coord=data.lon.data)
            model_contour = plt.contourf(lon_cyclic, data.lat.data, 
                                         data_cyclic,transform=ccrs.PlateCarree(),
                                         levels=mlevels,extend='both',cmap=cramp)
            ax.annotate(szn,
                xy=(0, 0), xycoords='data',
                xytext=(-0.05, 0.5), textcoords='axes fraction',fontsize=8,fontfamily = 'Arial',
                rotation=90,horizontalalignment='center', verticalalignment='center')  
            if szn == 'ANN':
                ax.annotate('CMIP Ensemble Mean ('+var+')',
                    xy=(0, 0), xycoords='data',fontsize=8,
                    xytext=(0.5, 1.1), textcoords='axes fraction',
                    horizontalalignment='center', verticalalignment='center')  
ax = plt.subplot(gs[(h*s):(h*s+1),0:len(models)*s])
ax.axis('off')
#cbar = plt.colorbar(model_contour,cax=inset_axes(ax,width='80%',height="30%",loc="upper center"),
 #           orientation="horizontal",ticks=[-0.55,-0.25,0,0.25,0.55]).set_label(var+' Difference ('+units+')',fontsize=8,c='black')
cbar = plt.colorbar(model_contour,cax=inset_axes(ax,width='90%',height="30%",loc="upper center"),
            orientation="horizontal",ticks=[-0.6,-0.3,0,0.3,0.6]).set_label('Mid-Holocene - Pre-Industrial ('+units+')',fontsize=8,c='black')
plt.tick_params(labelsize=8)
#Save or show
if save: plt.savefig(dataDir+'Figures/Model/Anomalies/CMIP_MH-PI_bySeason_'+var+'.png',
                     dpi=600,format='png',bbox_inches='tight')       
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




#%%
