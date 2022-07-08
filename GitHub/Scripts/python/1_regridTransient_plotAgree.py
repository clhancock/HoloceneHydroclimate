# This scrips regrids transient model data to CMIP resolution using the xESMF package
# It also creates maps for pre/p-e/tas showing mid-Holocene - picontrol agreement
# adpated from code provided by Michael Erb

import cmasher           as cmr
import cartopy.crs       as ccrs
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd 
import regionmask        as rm
import xarray            as xr

dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate/'



#%%Load Regridded Model Data For Plotting
modelData = {}
for model in ['hadcm','trace','cmip6']:
    modelData[model] = {}
    for szn in ['ANN']:
        modelData[model][szn] = xr.open_dataset(dataDir+'Data/Model/'+model+'_'+szn+'_regrid.nc',decode_times=False)


#%% 4. 
#%% 5. Calculate regional timeseries and cmip mh anomolies

for model in ['hadcm','trace']:
    for szn in ['ANN','JJA','DJF']:
        data = xr.open_dataset(saveDir+'Data/Model/'+model+'_'+szn+'.nc',decode_times=False)
        for var in ['pre','p-e','tas']:
            df = maskmodeldata(data[var],data[var].lat,data[var].lon)
            df.to_csv(saveDir+'Data/Model/RegionTS/'+'regional_'+var+'_'+szn+'_'+model+'.csv')
            


#%%Plot agreement
plt.rcParams['font.family'   ] = 'Arial'
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 1; 
plt.rcParams['axes.edgecolor'] = 'k'

save = False
for var in ['pre','p-e','tas']:
    refReg = rm.defined_regions.ar6.all
    refRegLand = rm.defined_regions.ar6.land
    if var == 'tas': 
        cramp, units, n = 'RdBu_r','degC', ['Cooler at 6ka','Warmer at 6ka']
        proxy = pd.read_csv(dataDir+'Data/proxyMetaData_T.csv')
    else:            
        cramp, units, n = 'BrBG', 'mm/day', ['Drier at 6ka','Wetter at 6ka',]
        proxy = pd.read_csv(dataDir+'Data/proxyMetaData_HC.csv')
    #
    #Calculate Proxy Percents for regions
    pRegs, pPcts, plats, plons, = [],[],[],[]
    for reg in np.unique(proxy['ipccReg']): 
        regData = proxy.loc[proxy['ipccReg'] == reg]
        if np.shape(regData)[0] <= 5: continue
        direction = 2*((regData['direction'] != 'negative')-0.5)
        regData = (regData['ka_6']-regData['ka_0.5'])*direction
        pPcts.append(100*(sum(regData>0)/(sum(np.isnan(regData)==False)+sum(regData==0))))
        plats.append(refReg.centroids[refReg.abbrevs.index(reg)][1])
        plons.append(refReg.centroids[refReg.abbrevs.index(reg)][0])
        pRegs.append(reg)
    #
    #Calculate sign of each gridcell and flatten into 2d
    modelVals = modelData['cmip6'][szn][var+'_regrid']
    modelVals = modelVals[0,:,:,:] - modelVals[1,:,:,:]
    modelN    = len(modelVals.model)
    modelAnom = np.sum(np.sign(modelVals),axis=0) 
    for model in ['hadcm','trace']:       
        modelVals = modelData[model][szn][var+'_regrid']
        modelVals = modelVals.groupby_bins('time',[0,1000,5500,6500]).mean(dim="time")
        modelAnom += np.sign(modelVals[2,:,:]-modelVals[1,:,:])
        modelN +=1
    #
    #Plot
    cramp = cmr.get_sub_cmap(cramp,0.1,0.9,N=modelN+1)
    plt.figure(figsize=(5,3))
    ax = plt.subplot(projection=ccrs.Robinson())
    refRegLand.plot_regions(ax=ax,add_label=False,line_kws=dict(linewidth=0.7))
    refReg[pRegs].plot_regions(ax=ax,add_label=False,line_kws=dict(linewidth=1.2))
    model_agree = plt.pcolormesh(modelVals.lon_regrid, modelVals.lat_regrid, 
                                 (modelAnom/modelN)*100,transform=ccrs.PlateCarree(),
                                 cmap=cramp,vmin=-100,vmax=100)
    ax.scatter(plons, plats, c=pPcts, transform=ccrs.PlateCarree(),
               cmap=cramp, vmin=0, vmax=100 ,s=40 ,ec='k', lw=2)
    plt.title('Agreement for sign of MH anomaly within model and proxy data \n(Annual '+var+')',
              fontsize=8)
    cbar = plt.colorbar(model_agree,orientation="horizontal",ticks=range(-100,101,50),
                        fraction=0.04, pad=0.02,aspect=30)
    cbar.ax.set_xticklabels(['100%\n'+n[0],'75%','50%\nEven Split','75%','100%\n'+n[1]],
                            fontsize=8)
    ax.set_global()
    if save: plt.savefig(dataDir+'Figures/Model/midHoloceneAgreement_'+var+'.png',
                         dpi=400,format='png',bbox_inches='tight')       
    plt.show()
    


