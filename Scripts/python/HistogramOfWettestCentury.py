import numpy      as np
import pandas     as pd
import regionmask as rm
import xarray     as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate/'
models = ['trace','hadcm']
szn = 'ANN'
var = 'pre'
modelData = {}
for model in models:
    modelData[model] = {}
    name = 'Data/Model/'+model+'_'+szn+'.nc'
    modelData[model] = xr.open_dataset(dataDir+name,decode_times=False)
proxy = pd.read_csv(dataDir+'Data/proxyMetaData_HC.csv')
proxyData = proxy.loc[proxy['recordRange'] > 9000]


#%%
vals = {}
for model in models:
    data = modelData[model][var]
    land = rm.defined_regions.natural_earth.land_110.mask_3D(data.lon,data.lat)
    land = land.squeeze('region').data
    for i in ['_dry','_wet']:
        vals[model+i+'_all']  = np.full(np.shape(data)[1:3],np.NaN) 
        vals[model+i+'_land'] = np.full(np.shape(data)[1:3],np.NaN) 
        for lat in range(len(data.lat)):
            for lon in range(len(data.lon)):
                ts = data[0:121,lat,lon].data
                if i == '_dry': 
                      n = [i for i, j in enumerate(ts) if j == np.nanmax(ts)]
                else: n = [i for i, j in enumerate(ts) if j == np.nanmin(ts)]
                vals[model+i+'_all'][lat,lon] = data.time[n].data[0]
                if land[lat,lon]: 
                    vals[model+i+'_land'][lat,lon] = vals[model+i+'_all'][lat,
                                                                          lon]

                
 #%%
  
save=True

plt.figure(figsize=(5,4)) ###
plt.title('Distribution of wettest century in model data')
gs = gridspec.GridSpec(3,2)
for i in ['wet','dry']:
    if i == 'wet': col,c,name=0,'olivedrab','maxValAge'
    else: col,c,name=1,'darkgoldenrod','minValAge'
    for model in models:
        row = [i for i, j in enumerate(models) if j == model][0]+1
        ax = plt.subplot(gs[row:row+1,col:col+1])
        ax.hist(vals[model+'_'+i+'_all'].flatten(),  bins=12, alpha=1,density=False,
                label='All Gridcells', color = 'gray', edgecolor='k',lw= 1)
        ax.hist(vals[model+'_'+i+'_land'].flatten(), bins=12, alpha=1,density=False,
                label='Land Only', color = c, edgecolor='k', lw= 1)
        if row == 1:     ax.tick_params(labelbottom=False) 
        else: 
            plt.xlabel('Age (yr BP)',fontsize=8)
            ax.legend(loc='upper right',fontsize=8)
        if col == 1: ax.tick_params(labelleft=False) 
        plt.ylabel(model,fontsize=8)
        ax.tick_params(labelsize=8)
        ax.set_xlim([12000,0])
        ax.set_ylim([0,2000+1000*row])
        ax.spines['right'].set_visible(True);ax.spines['top'].set_visible(True)
    ax = plt.subplot(gs[0:1,col:col+1])
    ax.hist(proxyData[name],  bins=12, alpha=1,density=False,
                label='All Gridcells', color = c, edgecolor='k',lw= 1)
    plt.ylabel('proxy',fontsize=8)
    ax.tick_params(labelsize=8)
    ax.set_xlim([12000,0])
    ax.set_ylim([0,80])
    #
    ax.tick_params(labelbottom=False) 
    ax.set_title('Largest '+i+' anomoly',fontsize=8)  #y=1.1, x = 0.42,loc='left', fontsize = 10)
    ax.spines['right'].set_visible(True);ax.spines['top'].set_visible(True)
    if col == 1: ax.tick_params(labelleft=False) 

    
if save: plt.savefig(dataDir+'Figures/HistWettestCentury_'+var+'_'+szn+'.png',
                         dpi=400,format='png',bbox_inches='tight')       
else: plt.show()