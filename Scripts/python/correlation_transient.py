#%% 1 Load Packages
import numpy      as np
#import pandas     as pd
import regionmask as rm
import xarray     as xr
#from scipy        import stats
from scipy.stats  import pearsonr
#For Figures:
import matplotlib.pyplot as plt         # Packages for making figures
import matplotlib.gridspec as gridspec
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs              # Packages for mapping in python
import cartopy.util as cutil
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate/'
#%% 
def calcCorrelation (inData,model1,model2,var1,var2,season,t0=0,t1=12,mask=False):
    if model1 == model2:
        modelData1 = inData[model1][season][var1]
        modelData2 = inData[model2][season][var2]
        lats,lons = modelData1.lat,modelData1.lon
    else: 
        modelData1 = inData[model1][season][var1+'_regrid']
        modelData2 = inData[model2][season][var2+'_regrid']
        lats,lons = modelData1.lat_regrid,modelData1.lon_regrid
    t_0 = np.argmin(np.abs(t0*1000 - modelData1.time.data))   
    t_1 = np.argmin(np.abs(t1*1000 - modelData1.time.data))
    out = np.full(np.shape(modelData1)[1:3],np.NaN)
    for lat in range(len(lats)):
        for lon in range(len(lons)):
            out[lat,lon] = pearsonr(modelData1[t_0:t_1,lat,lon].data,modelData2[t_0:t_1,lat,lon].data)[0]

    return(out)

#%% 
#Settings
save=False
times = [0,12]
szn = 'ANN'
v_1 = 'p-e'
v_2 = 'tas'
m_1 = 'trace'
m_2 = 'trace'
#import seaborn as sns
#sns.palplot(sns.color_palette("BrBG", 7))
from matplotlib.colors import LinearSegmentedColormap
cramp = plt.cm.get_cmap('RdGy_r',10)
cramp = LinearSegmentedColormap.from_list('cramp',['#7f3b08','white','#2d004b'],N=10)

#Load data
modelData = {}
for model in ['hadcm','trace']:
    modelData[model] = {}
    if m_1 == m_2:
        modelData[model][szn]   = xr.open_dataset(dataDir+'Data/Model/'+model+'_'+szn+'.nc',decode_times=False)
    else: modelData[model][szn] = xr.open_dataset(dataDir+'Data/Model/'+model+'_'+szn+'_regrid.nc',decode_times=False)

#Calculate using function
rVals = calcCorrelation(modelData,m_1,m_2,v_1,v_2,szn,times[0],times[1])

#%% #Plot 
save = True
mlevels = np.array([i /10 for i in list(range(-10,11,2))])
if m_1 == m_2: 
    lats,lons = modelData[m_1][szn]['lat'],modelData[m_1][szn]['lon']
    name = m_1+'_'+v_1+'_'+v_2+'_'+szn
else: 
    lats,lons = modelData[m_1][szn]['lat_regrid'],modelData[m_1][szn]['lon_regrid']
    name = 'multi'+'_'+v_1+'_'+v_2+'_'+szn

proxy=i#NEED TO CHANGE
refReg = rm.defined_regions.ar6.all
    #Calculate Proxy Percents for regions
pRegs, pVals, plats, plons, = [],[],[],[]
for reg in np.unique(proxy['ipccReg']): 
    #pVals.append(100*(sum(regData>0)/(sum(np.isnan(regData)==False)+sum(regData==0))))#NEED TO CHANGE
    plats.append(refReg.centroids[refReg.abbrevs.index(reg)][1])
    plons.append(refReg.centroids[refReg.abbrevs.index(reg)][0])
    pRegs.append(reg)
    #
plt.style.use('default')
gs = gridspec.GridSpec(12,8)
ax = plt.subplot(gs[0:12,0:6],projection=ccrs.Robinson()) 
plt.title('r-values: '+m_1+' ('+szn+') ('+v_1+':'+v_2+')')
ax.set_global()
ax.coastlines()
data_cyclic,lon_cyclic = cutil.add_cyclic_point(rVals,coord=lons)
model_contour = plt.contourf(lon_cyclic,lats, data_cyclic,transform=ccrs.PlateCarree(),
                             levels=mlevels,cmap=cramp) 
cm = plt.cm.get_cmap(cramp,10)
plt.colorbar(model_contour,cax=inset_axes(ax,width='70%',height="6%",loc=8),orientation="horizontal")


mask = rm.defined_regions.natural_earth.land_110.mask_3D(lons,lats).squeeze('region').data

ax1 = plt.subplot(gs[3:9,6:8]) 
rValsLand = np.full(np.shape(rVals),np.NaN) 
for lat in range(len(lats)):
    for lon in range(len(lons)):
        if mask[lat,lon]: rValsLand[lat,lon] = rVals[lat,lon]
              
ax1.axvline(x=0,color='grey',lw=0.5)
ax1.axhline(y=0,color='grey',lw=0.5)
ax1.axis([mlevels[0],mlevels[-1],min(lats),max(lats)])

ax1.plot(np.nanmean(rVals, axis=1) ,lats,c='k',lw=1.5)

x = np.nanmean(rValsLand, axis=1) 
x[np.isnan(x)] = 0
ax1.plot(x,lats,c='darkgreen',lw=1.5)

ax1.set_xticks([])
ax1.set_yticks(range(-90,91,30))
ax1.set_yticklabels([])
if save: plt.savefig(dataDir+'Figures/Model/TransientCorrelations/'+name+'.png',
                     dpi=400,format='png',bbox_inches='tight')       
plt.show()

#%%
for land in ['all','land']:
    if land == 'all':
        ax1 = plt.subplot(gs[3:9,6:7]) 
        x = np.nanmean(rVals, axis=1) 
    else:
        ax1 = plt.subplot(gs[3:9,7:8]) 
        rValsLand = np.full(np.shape(rVals),np.NaN) 
        for lat in range(len(lats)):
            for lon in range(len(lons)):
                if mask[lat,lon]: rValsLand[lat,lon] = rVals[lat,lon]
        x = np.nanmean(rValsLand, axis=1) 
        x[np.isnan(x)] = 0
    ax1.axvline(x=0,color='grey',lw=0.5)
    ax1.axhline(y=0,color='grey',lw=0.5)
    ax1.axis([mlevels[0],mlevels[-1],min(lats),max(lats)])
    lc = LineCollection([np.column_stack([x[i:i+2],lats[i:i+2]]) for i in range(len(x)-1)], 
                        array=x,cmap=cramp,lw=2)
    ax1.add_collection(lc)
    ax1.plot(x,lats,c='k',lw=2)
    ax1.set_xticks([])
    ax1.set_yticks(range(-90,91,30))
    ax1.set_yticklabels([])
plt.show()
#%%

#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

#plt.colorbar(model_contour,cax=inset_axes(ax,width='70%',height="4%",loc=8),orientation="horizontal")
if save: plt.savefig(dataDir+'Figures/Model/TransientCorrelations/'+name+'.png',
                     dpi=400,format='png',bbox_inches='tight')       
plt.show()