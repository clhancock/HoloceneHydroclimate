#%%Script to determine the wettest and dryest century of model and proxy data
#Input: csv of proxy metadata and model netcdf files
#Output: Figure of Hustograms

import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec
import numpy               as np
import pandas              as pd
import regionmask          as rm
import seaborn             as sns
import xarray              as xr
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = '8'
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 0.5; 
plt.rcParams['axes.edgecolor'] = 'k'
#plt.tick_params(labelsize=8)

Dir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/2021_HoloceneHydroclimate/2021_HoloceneHydroclimate/'


binSize = 1000
binN = int(12000/binSize)
#Load Proxy Data
proxy = pd.read_csv(Dir+'Data/Proxy/proxyMetaData_HC.csv')
proxyData = proxy#.loc[(proxy['recordRange'] > 0)] #& (proxy['archive'] !='LakeDeposits')]

proxyDataWet = np.histogram(proxyData['maxValAge'],bins=binN,range=[0,12000])[0]
proxyDataDry = np.histogram(proxyData['minValAge'],bins=binN,range=[0,12000])[0]

#Scale proxy data relative to number of available records 
for i in range(0,binN):
    scale = proxyData.loc[(proxy['minAge'] < i*binSize+binSize) & (proxy['maxAge'] > i*binSize)] 
    #print(len(subset))
    proxyDataWet[i] = 100*proxyDataWet[i]/len(scale)
    proxyDataDry[i] = 100*proxyDataDry[i]/len(scale)
    



#%% Calculate Wettest and Dryest Century of transient model
models = ['trace','hadcm']
szn = 'ANN'
var = 'pre'
vals = {}
for model in models:
    data = xr.open_dataset(Dir+'Data/Model/'+model+'/'+model+'_'+szn+'.nc',decode_times=False)
    land = rm.defined_regions.natural_earth.land_110.mask_3D(data)
    vals[model+'_Dry_all'] = np.argmin(data[var].data,axis=0)*abs(data.age[0]-data.age[1]).data
    vals[model+'_Wet_all'] = np.argmax(data[var].data,axis=0)*abs(data.age[0]-data.age[1]).data
    for i in ['Wet','Dry']:
        name  = model+'_'+i
        scale = len(vals[name+'_all'].flatten())
        vals[name+'_land'] = vals[name+'_all']*land.where(land == True).data[0]
        vals[name+'_all']  = np.histogram(vals[name+'_all'].flatten(), 
                                          bins=binN,range=[0,12000])[0]*100/scale
        vals[name+'_land'] = np.histogram(vals[name+'_land'].flatten(),
                                          bins=binN,range=[0,12000])[0]*100/scale

    



#%%Plot
n=40
d=4
palette = sns.color_palette("BrBG_r", n).as_hex()

save=True
plt.rcParams["axes.labelpad"] = 8.0
plt.figure(figsize=(6.5,4)) ###
plt.title('Distribution of wettest century in model data')
gs = gridspec.GridSpec(3,2,hspace=0.1,wspace=0.1)
for i in ['Wet','Dry']:
    if i == 'Wet': col,c,name=0,palette[0+d],  'maxValAge'
    else:          col,c,name=1,palette[n-d-1],'minValAge'
    for model in models:
        row = [i for i, j in enumerate(models) if j == model][0]+1
        total = len(vals[model+'_'+i+'_all'].flatten())
        ax = plt.subplot(gs[row:row+1,col:col+1])
        ax.bar(np.linspace(0.5, 11.5,binN),vals[model+'_'+i+'_all'],
               width=binSize/1000,label='All Grid Cells', color = 'lightgrey', edgecolor='k',lw= 1)
        ax.bar(np.linspace(0.5, 11.5,binN),vals[model+'_'+i+'_land'],
               width=binSize/1000,label='Land Only', color = c, alpha=0.8,edgecolor='k',lw= 1)
        if row == 1: ax.tick_params(labelbottom=False) 
        else: 
            plt.xlabel('Age (ka BP)',fontsize=8)
            ax.legend(loc='upper center',fontsize=8)
        if col == 1: 
            ax.tick_params(labelleft=False) 
            plt.ylabel('',fontsize=0)
        else: 
            if model == 'hadcm': plt.ylabel('\n'+'HadCM\n(% of grid cells)',fontsize=8)
            else: plt.ylabel('\n'+'TraCE\n(% of grid cells)',fontsize=8)
        ax.tick_params(labelsize=8)
        ax.set_xlim([12,0])
        ax.set_ylim([0,50])
        ax.spines['right'].set_visible(True);ax.spines['top'].set_visible(True)
    ax = plt.subplot(gs[0:1,col:col+1])
    if i == 'Wet': proxydata = proxyDataWet
    else:           proxydata = proxyDataDry
    ax.bar(np.linspace(0.5, 11.5,binN),proxydata,
               width=binSize/1000,label='All Grid Cells', color = c,alpha=0.8, edgecolor='k',lw= 1)
    plt.ylabel('Proxy\n(% of records)',fontsize=8)
    if col == 1: plt.ylabel('',fontsize=0)
    ax.tick_params(labelsize=8)
    ax.set_xlim([12,0])
    ax.set_ylim([0,50])
    #
    ax.tick_params(labelbottom=False) 
    ax.set_title('Largest '+i+' '+var.upper()+' Anomaly',fontsize=8)  #y=1.1, x = 0.42,loc='left', fontsize = 10)
    ax.spines['right'].set_visible(True);ax.spines['top'].set_visible(True)
    if col == 1: ax.tick_params(labelleft=False) 

    
if save: plt.savefig(Dir+'Figures/Model/HistWettestCentury_'+var+'_'+szn+'.png',
                         dpi=400,format='png',bbox_inches='tight')       
else: plt.show()

#%% Map
from   matplotlib.colors   import LinearSegmentedColormap
import cartopy.crs         as ccrs        # Packages for mapping in python
import cartopy.util        as cutil

save = False
cramp = plt.cm.get_cmap('RdGy_r',10)
cramp = LinearSegmentedColormap.from_list('cramp',['#7f3b08','white','#2d004b'],N=10)
mlevels = np.array([i /10 for i in list(range(-10,11,2))])

models = ['trace','hadcm']
szn = 'ANN'
var = 'pre'
vals = {}
data = xr.open_dataset(Dir+'Data/Model/'+model+'/'+model+'_'+szn+'.nc',decode_times=False)
    #land = rm.defined_regions.natural_earth.land_110.mask_3D(data)
    #vals[model+'_Dry_all'] = np.argmin(data[var].data,axis=0)*abs(data.age[0]-data.age[1]).data
vals = np.argmax(data[var].data,axis=0)*abs(data.age[0]-data.age[1]).data
    

cramp = plt.cm.get_cmap("PuOr",8)
ax = plt.axes(projection=ccrs.Robinson())
model_contour=plt.pcolormesh(data.lon, data.lat, vals,transform=ccrs.PlateCarree(),
                              vmin=0,vmax=12000,cmap=cramp)

ax.coastlines()
ax.set_global()
cbar = plt.colorbar(model_contour,orientation="horizontal",#ticks=[-1,-0.6,-0.3,0,0.3,0.6,1],
                        fraction=0.04, pad=0.04,aspect=30)
plt.show()





#cbar.set_label('Pearson Correlation Coefficient',fontsize=8, fontfamily = 'Times New Roman')
#cbar.ax.set_xticklabels([-1,-0.6,-0.3,0,0.3,0.6,1],fontsize=8)
plt.show()



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
            textcoords='axes fraction', fontsize=8, fontfamily = 'Times New Roman')



  