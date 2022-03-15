import numpy      as np
import pandas     as pd
import regionmask as rm
import xarray     as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate/'
models = ['trace','hadcm']
szn = 'ANN'
var = 'p-e'
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
    for i in ['_Dry','_Wet']:
        vals[model+i+'_all']  = np.full(np.shape(data)[1:3],np.NaN) 
        vals[model+i+'_land'] = np.full(np.shape(data)[1:3],np.NaN) 
        for lat in range(len(data.lat)):
            for lon in range(len(data.lon)):
                ts = data[0:121,lat,lon].data
                if i == '_Dry': 
                      n = [i for i, j in enumerate(ts) if j == np.nanmax(ts)]
                else: n = [i for i, j in enumerate(ts) if j == np.nanmin(ts)]
                vals[model+i+'_all'][lat,lon] = data.time[n].data[0]
                if land[lat,lon]: 
                    vals[model+i+'_land'][lat,lon] = vals[model+i+'_all'][lat,
                                                                          lon]

                
 #%%
import seaborn as sns
n=40
d=4
palette = sns.color_palette("BrBG_r", n).as_hex()

save=True
plt.rcParams["axes.labelpad"] = 8.0
plt.figure(figsize=(6.5,4)) ###
plt.title('Distribution of wettest century in model data')
gs = gridspec.GridSpec(3,2,hspace=0.1,wspace=0.1)
for i in ['Wet','Dry']:
    if i == 'Wet': col,c,name=0,palette[0+d],'maxValAge'
    else: col,c,name=1,palette[n-d-1],'minValAge'
    for model in models:
        row = [i for i, j in enumerate(models) if j == model][0]+1
        total = len(vals[model+'_'+i+'_all'].flatten())
        ax = plt.subplot(gs[row:row+1,col:col+1])
        ax.bar(np.linspace(0.5, 11.5,12),
               100*np.histogram(vals[model+'_'+i+'_all'].flatten(),bins=12,range=[0,12000])[0]/total,  
               width=1,label='All Gridcells', color = 'lightgrey', edgecolor='k',lw= 1)
        ax.bar(np.linspace(0.5, 11.5,12),
               100*np.histogram(vals[model+'_'+i+'_land'].flatten(),bins=12,range=[0,12000])[0]/total,  
               width=1,label='Land Only', color = c, alpha=0.8,edgecolor='k',lw= 1)
        if row == 1:     ax.tick_params(labelbottom=False) 
        else: 
            plt.xlabel('Age (ka BP)',fontsize=8)
            ax.legend(loc='upper right',fontsize=8)
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
    ax.bar(np.linspace(0.5, 11.5,12),
               100*np.histogram(proxyData[name],bins=12,range=[0,12000])[0]/len(proxyData[name]),  
               width=1,label='All Grid Cells', color = c,alpha=0.8, edgecolor='k',lw= 1)
    plt.ylabel('Proxy\n(% of records)',fontsize=8)
    if col == 1: plt.ylabel('',fontsize=0)
    ax.tick_params(labelsize=8)
    ax.set_xlim([12,0])
    ax.set_ylim([0,25])
    #
    ax.tick_params(labelbottom=False) 
    ax.set_title('Largest '+i+' Anomaly',fontsize=8)  #y=1.1, x = 0.42,loc='left', fontsize = 10)
    ax.spines['right'].set_visible(True);ax.spines['top'].set_visible(True)
    if col == 1: ax.tick_params(labelleft=False) 

    
if save: plt.savefig(dataDir+'Figures/HistWettestCentury_'+var+'_'+szn+'.png',
                         dpi=400,format='png',bbox_inches='tight')       
else: plt.show()