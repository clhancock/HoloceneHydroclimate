#=============================================================================
# This regrids transient model data to CMIP resolution using the xESMF package
# It also creates maps for pre/p-e/tas showing mid-Holocene - picontrol agreement
# adpated from code provided by Michael Erb

import numpy             as np
import pandas            as pd 
import regionmask        as rm
import xarray            as xr
#import xesmf             as xe
import cartopy.crs       as ccrs
import matplotlib.pyplot as plt
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes

dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate/'

#%%Regrid and save tansient data to cmip resolution
make_figure = True
for szn in ['ANN','JJA','DJF']:
    dataGrid = xr.open_dataset(dataDir+'Data/Model/cmip6_'+szn+'.nc') #Regrid
    for model in ['hadcm','trace']:
        dataOrig = xr.open_dataset(dataDir+'Data/Model/'+model+'_'+szn+'.nc')
        dataReGr = xr.Dataset({
             'lat': (['lat'],dataGrid['lat_regrid'].data,{'units':'degrees_north'}),
             'lon': (['lon'],dataGrid['lon_regrid'].data,{'units':'degrees_east'}),
            })
        #Weird impact of having values directly at -90 and 90 deg lat
        if model == 'hadcm': dataOrig = dataOrig.isel(lat=list(range(1,len(dataOrig.lat)-1)))
        #Regrid
        regridder = xe.Regridder(dataOrig,dataReGr,'conservative')
        #regridder = xe.Regridder(dataOrig,dataReGr,'conservative',periodic=True,ignore_degenerate=True)
        #regridder = xe.Regridder(dataOrig,dataReGr,'bilinear')         
        dataReGr = regridder(dataOrig,keep_attrs=True)
        #Plot to make sure it looks right
        print(np.max(abs(dataReGr['tas'].data)))
        print(np.max(abs(dataOrig['tas'].data)))
        if make_figure:
            ax1 = plt.subplot2grid((2,1),(0,0),projection=ccrs.PlateCarree())
            ax2 = plt.subplot2grid((2,1),(1,0),projection=ccrs.PlateCarree())
            dataOrig['tas'].isel(time=0).plot.pcolormesh(ax=ax1)
            dataReGr['tas'].isel(time=0).plot.pcolormesh(ax=ax2)
            ax1.coastlines(); ax2.coastlines()
            plt.plot()
            plt.show()
        #Rename to indicate regrid
        for var in dataReGr.keys(): dataReGr = dataReGr.rename({var:var+'_regrid'})
        dataReGr = dataReGr.rename({'lat':'lat_regrid','lon':'lon_regrid'})
        dataReGr.to_netcdf(dataDir+'Data/Model/'+model+'_'+szn+'_regrid'+'.nc')



#%%Load Model Data For Plotting
modelData = {}
for model in ['hadcm','trace','cmip6']:
    modelData[model] = {}
    for szn in ['ANN']:
        modelData[model][szn] = xr.open_dataset(dataDir+'Data/Model/'+model+'_'+szn+'_regrid.nc',decode_times=False)


#%%Plot agreement
#Settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 1; 
plt.rcParams['axes.edgecolor'] = 'k'
plt.tick_params(labelsize=8)

save = True
import cmasher as cmr
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
    ax.scatter(plons,plats,c=pPcts,transform=ccrs.PlateCarree(),cmap=cramp,
               vmin=0,vmax=100,s=40,ec='k',lw=2)
    plt.title('Agreement for sign of MH anomaly within model and proxy data \n(Annual '+var+')',
              fontsize=8)
    cbar = plt.colorbar(model_agree,orientation="horizontal",ticks=range(-100,101,50),
                        fraction=0.04, pad=0.02,aspect=30)
    cbar.ax.set_xticklabels(['100%\n'+n[0],'75%','50%\nEven Split','75%','100%\n'+n[1]],
                            fontsize=8)
    ax.set_global()
    #cbar.ax.xaxis.set_ticks_position("top")
    #cbar.set_label(''+var, labelpad=-35, x=-0.1)#  '% positive (annual '+var+') mid-Holocene anomaly',
    if save: plt.savefig(dataDir+'Figures/Model/midHoloceneAgreement_'+var+'.png',
                         dpi=400,format='png',bbox_inches='tight')       
    plt.show()
    


