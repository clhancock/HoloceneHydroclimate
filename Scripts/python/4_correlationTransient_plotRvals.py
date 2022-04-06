#Script to calculate and plot correlations between climate fields and model simulations
#Input: Transient Model netcdf files
#Output: Plot of correlation coefficients
#%% 1 Load Packages
import numpy       as np
import pandas      as pd
import regionmask  as rm
import xarray      as xr
from   scipy.stats import pearsonr 

#For Figures:
import matplotlib.pyplot   as plt         # Packages for making figures
import matplotlib.gridspec as gridspec
from   matplotlib.colors   import LinearSegmentedColormap
import cartopy.crs         as ccrs        # Packages for mapping in python
import cartopy.util        as cutil

dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate/'
#%% 2 Define Function for calculating grid cell correlation values
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

#%% Calculate correlation coefficients based on model settings
#Settings
times = [0,12]
szn   = 'ANN'
v_1   = 'p-e'
v_2   = 'tas'
m_1   = 'trace'
m_2   = 'trace'

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
save = False
cramp = plt.cm.get_cmap('RdGy_r',10)
cramp = LinearSegmentedColormap.from_list('cramp',['#7f3b08','white','#2d004b'],N=10)
mlevels = np.array([i /10 for i in list(range(-10,11,2))])

if m_1 == m_2: 
    lats,lons = modelData[m_1][szn]['lat'],modelData[m_1][szn]['lon']
    name = m_1+'_'+v_1+'_'+v_2+'_'+szn
else: 
    lats,lons = modelData[m_1][szn]['lat_regrid'],modelData[m_1][szn]['lon_regrid']
    name = 'multi'+'_'+v_1+'_'+v_2+'_'+szn
import cmasher as cmr

cramp = LinearSegmentedColormap.from_list('cramp',['#40004b','white','#00441b'],N=30)
cramp = cmr.get_sub_cmap(cramp,0.1,0.9,N=20)
 
proxy= pd.read_csv(dataDir+'Data/'+'HC_T_RegionalProxyEnsCorrelations.csv')
refReg = rm.defined_regions.ar6.all
    #Calculate Proxy Percents for regions
pltRegs, pltVals, pltlats, pltlons, = [],[],[],[]
for reg in proxy.columns.values.tolist()[1:]: 
    pltVals.append(np.nanmean(proxy[reg]))
    pltlats.append(refReg.centroids[refReg.abbrevs.index(reg)][1])
    pltlons.append(refReg.centroids[refReg.abbrevs.index(reg)][0])
    pltRegs.append(reg) 
    #
    
plt.figure(figsize=(6,3.01))
gs = gridspec.GridSpec(12,8)
ax = plt.subplot(gs[0:12,0:6],projection=ccrs.Robinson()) 
refRegLand = rm.defined_regions.ar6.land
refRegLand.plot_regions(ax=ax,add_label=False,line_kws=dict(linewidth=0.7))
data_cyclic,lon_cyclic = cutil.add_cyclic_point(rVals,coord=lons)
model_contour = plt.contourf(lon_cyclic,lats, data_cyclic,transform=ccrs.PlateCarree(),
                             vmin=-1,vmax=1,cmap=cramp,levels=20)  
ax.scatter(pltlons,pltlats,c=pltVals,transform=ccrs.PlateCarree(),
           cmap=cramp,vmin=-1,vmax=1,s=40,ec='k',lw=2)
ax.set_global()
ax.annotate('(a)',xy=(0, 0), xycoords='data', xytext=(0.05, 0.95), 
            textcoords='axes fraction', fontsize=8, fontfamily = 'Arial')
cbar = plt.colorbar(model_contour,orientation="horizontal",ticks=[-1,-0.6,-0.3,0,0.3,0.6,1],
                        fraction=0.04, pad=0.04,aspect=30)
cbar.set_label('Pearson Correlation Coefficient')
cbar.ax.set_xticklabels([-1,-0.6,-0.3,0,0.3,0.6,1],fontsize=8)


mask = rm.defined_regions.natural_earth.land_110.mask_3D(lons,lats).squeeze('region').data

ax1 = plt.subplot(gs[2:11,6:8]) 
rValsLand = np.full(np.shape(rVals),np.NaN) 
for lat in range(len(lats)):
    for lon in range(len(lons)):
        if mask[lat,lon]: rValsLand[lat,lon] = rVals[lat,lon]
ax1.annotate('(b)',xy=(0, 0), xycoords='data', xytext=(-0.2, 0.95), 
            textcoords='axes fraction', fontsize=8, fontfamily = 'Arial')
ax1.axvline(x=0,color='black',lw=0.5)
for i in range(-90,91,30):
    ax1.axhline(y=i,xmin=0,xmax=0.05,color='grey',lw=0.3)
    ax1.axhline(y=i,xmin=0.95,xmax=1,color='grey',lw=0.3)

ax1.axis([mlevels[0],mlevels[-1],lats[-2],lats[1]])

x = np.nanmean(rValsLand, axis=1) 
x[np.isnan(x)] = 0
ax1.plot(np.nanmean(rVals,axis=1) ,lats,c='darkslategrey',lw=1.9)
ax1.plot(x,lats,c='darkgoldenrod',lw=1.9)
ax1.scatter(pltVals,pltlats,c='k',s=7)#plt.cm.get_cmap('BrBG',5),
ax1.spines['right'].set_visible(False);ax1.spines['left'].set_visible(False)
ax1.set_xticks([-1,0,1])
ax1.set_yticklabels([-1,0,1],fontsize=8)
ax1.set_yticks([])
ax1.set_yticklabels([])

#Save or show
if save: plt.savefig(dataDir+'Figures/Model/TransientCorrelations/'+name+'.png',
                     dpi=400,format='png',bbox_inches='tight')       
plt.show()

