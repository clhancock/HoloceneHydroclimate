#==============================================================================
# This regrids model data using the xESMF package.
#    author: Michael P. Erb
#    date  : 2/3/2022
#==============================================================================

import numpy             as np
import regionmask        as rm
import xarray            as xr
import xesmf             as xe
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
save = False
var = 'tas'
if var == 'tas': cramp, units = 'RdBu_r','degC'
else:            cramp, units = 'BrBG', 'mm/day'

#Calculate sign of each gridcell and flatten into 2d
modelVals = modelData['cmip6'][szn][var+'_regrid']
modelVals = modelVals[0,:,:,:] - modelVals[1,:,:,:]
modelN    = len(modelVals.model)
modelAnom = np.sum(np.sign(modelVals),axis=0)
dataPct   = (((modelAnom+modelN)/2)/modelN)*100

#Plot
plt.figure(figsize=(6,4))
ax = rm.defined_regions.ar6.land.plot(projection=ccrs.Robinson(),
                                      add_label=False,line_kws=dict(linewidth=1))
model_agree = plt.pcolormesh(modelVals.lon_regrid, modelVals.lat_regrid, 
                             dataPct,transform=ccrs.PlateCarree(),
                             cmap=cramp,alpha=1,vmin=0,vmax=100)
ax.coastlines()
ax.set_global()
plt.title('Model & Proxy Agreement for sign of MH-PI difference',fontsize=10)
plt.colorbar(model_agree,
             cax=inset_axes(ax,width='70%',height="4%",loc="lower center"),
             orientation="horizontal").set_label(
                 '% positive ('+var+') mid-Holocene anomaly',fontsize=10,c='k')
if save: plt.savefig(dataDir+'Figures/Model/midHoloceneAgreement_'+var+'.png',
                     dpi=400,format='png',bbox_inches='tight')       
plt.show()




