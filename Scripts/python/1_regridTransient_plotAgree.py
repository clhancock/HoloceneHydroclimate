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

Dir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/2021_HoloceneHydroclimate/2021_HoloceneHydroclimate/'



#%%Load Regridded Model Data For Plotting

modelData = {}
for model in ['hadcm','trace','cmip6']:
    modelData[model] = {}
    for szn in ['ANN','JJA','DJF']:
        if model == 'cmip6': end = '.nc'
        else: end =  '.nc'
        modelData[model][szn] = xr.open_dataset(Dir+'Data/Model/'+model+'/'+model+'_'+szn+end,decode_times=False)


#%% 4. #Function for masking model data based on IPCC regions and land

def maskmodeldata(values,lats,lons,geo='all'):
    #Lat weights
    wghts  = np.cos(np.deg2rad(lats))
    #IPCC ref regions
    refReg = rm.defined_regions.ar6.all
    refReg = refReg.mask_3D(lons,lats)  
    if geo == 'land':
        #Land mask
        land   = rm.defined_regions.natural_earth.land_110
        land   = land.mask_3D(lons,lats)   
        land   = np.array([land.squeeze('region').data]*np.shape(refReg)[0])
        #3d mask with land and refRegions
        mask = refReg*land
    else: mask = refReg*1
    #Average value by region
    out = values.weighted(mask * wghts).mean(dim=("lat", "lon")).data
    if np.shape(out)[1] == len(refReg):#46: 
        out = pd.DataFrame(out) 
        try: out.index  = list(values.age.data)
        except:  out.index  = list(values.model.data)
        out.columns= list(mask.abbrevs.data)
        return out
    else: print("Error: out df does not have correct dims")
    
    
    
#%% 5. Calculate regional timeseries and cmip mh anomolies

for model in ['hadcm','trace','cmip6']:
    for szn in ['ANN','JJA','DJF']:
        data =  modelData[model][szn]
        for var in ['pre','p-e','tas']:
            for geovar in ['all','land']:
                df = maskmodeldata(data[var],data[var].lat,data[var].lon,geo=geovar)
                df.to_csv(Dir+'Data/Model/RegionalTS/'+'regional'+'_'+var+'_'+szn+'_'+model+'_'+geovar+'.csv')


#%%Plot agreement
modelData = {}
for model in ['hadcm','trace','cmip6']:
    modelData[model] = {}
    for szn in ['ANN','JJA','DJF']:
        if model == 'cmip6': end = '.nc'
        else: end =  '_regrid.nc'
        modelData[model][szn] = xr.open_dataset(Dir+'Data/Model/'+model+'/'+model+'_'+szn+end,decode_times=False)

plt.rcParams['font.family'   ] = 'Arial'
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 1; 
plt.rcParams['axes.edgecolor'] = 'k'

szn = 'ANN'
save = False
for var in ['pre','p-e','tas']:
    refReg = rm.defined_regions.ar6.all
    refRegLand = rm.defined_regions.ar6.land
    if var == 'tas': 
        cramp, units, n = 'RdBu_r','degC', ['Cooler at 6ka','Warmer at 6ka']
        proxy = pd.read_csv(Dir+'Data/proxy/proxyMetaData_T.csv')
    else:            
        cramp, units, n = 'BrBG', 'mm/day', ['Drier at 6ka','Wetter at 6ka',]
        proxy = pd.read_csv(Dir+'Data/proxy/proxyMetaData_HC.csv')
    #
    #Calculate Proxy Percents for regions
    #Calculate Proxy Percents for regions
    pRegs, pPcts, plats, plons, = [],[],[],[]
    for reg in np.unique(proxy['ipccReg']): 
        regData = proxy.loc[proxy['ipccReg'] == reg]
        if np.shape(regData)[0] < 4: continue
        direction = 2*((regData['direction'] != 'negative')-0.5)
        regData = (regData['ka_6']-regData['ka_0.5'])*direction
        if sum(~np.isnan(regData)) < 4: continue
        pPcts.append(100*(sum(regData>0)/(sum(np.isnan(regData)==False)+sum(regData==0))))
        plats.append(refReg.centroids[refReg.abbrevs.index(reg)][1])
        plons.append(refReg.centroids[refReg.abbrevs.index(reg)][0])
        pRegs.append(reg)
    #Calculate sign of each gridcell and flatten into 2d
    modelVals = modelData['cmip6'][szn][var]
    modelN    = len(modelVals.model)
    modelAnom = np.sum(np.sign(modelVals),axis=0) 
    for model in ['hadcm','trace']:       
        modelVals = modelData[model][szn]
        modelVals = modelVals.rename({var+'_regrid':var,'lat_regrid':'lat','lon_regrid':'lon'})
        modelVals = modelVals[var].groupby_bins('time',[0,1000,5500,6500]).mean(dim="time")
        modelAnom += np.sign(modelVals[2,:,:]-modelVals[0,:,:])
        modelN +=1
    #
    #Plot
    cramp = cmr.get_sub_cmap(cramp,0.1,0.9,N=modelN+1)
    #cramp = cmr.get_sub_cmap(cramp,0.1,0.9,N=int(modelN/2)) 
    plt.figure(figsize=(5,3))
    ax = plt.subplot(projection=ccrs.Robinson())
    refRegLand.plot_regions(ax=ax,add_label=False,line_kws=dict(linewidth=1))
    refReg[pRegs].plot_regions(ax=ax,add_label=False,line_kws=dict(linewidth=1.2))
    model_agree = plt.pcolormesh(modelVals.lon, modelVals.lat, 
                                 (modelAnom/modelN)*100,transform=ccrs.PlateCarree(),
                                 cmap=cramp,vmin=-100,vmax=100)
    ax.scatter(plons, plats, c=pPcts, transform=ccrs.PlateCarree(),
               cmap=cramp, vmin=0, vmax=100 ,s=40 ,ec='k', lw=2)
    cbar = plt.colorbar(model_agree,orientation="horizontal",ticks=range(-100,101,50),
                        fraction=0.04, pad=0.145,aspect=30)
    cbar.ax.set_xticklabels(['100%\n'+n[0],'75%','50%\nEven Split','75%','100%\n'+n[1]],
                            fontsize=8,
                            verticalalignment='baseline')
    cbar.ax.xaxis.set_ticks_position("top")    
    ax.set_global()
    cbar.ax.text(0,-1.6,'Agreement for sign of MH anomaly within model and proxy data \n(Annual '+var.upper()+')',
                 horizontalalignment='center',
        verticalalignment='center',fontsize=8)
    plt.tight_layout()
    if save: plt.savefig(Dir+'Figures/Model/Agreement/midHoloceneAgreement_'+var+'.png',
                         dpi=400,format='png',bbox_inches='tight')       
    plt.show()
    #















#%%

model='trace'
szn='ANN'
data =  modelData[model][szn]
values = data[var]
lats = values.lat
lons = values.lon
geo='all'
wghts  = np.cos(np.deg2rad(lats))                                                                                          
#IPCC ref regions
refReg = rm.defined_regions.ar6.all.mask_3D(lons,lats)  



#%%
import regionmask
data =  modelData[model][szn]
values = data
airtemps=values
mask_3D = regionmask.defined_regions.srex.mask_3D(airtemps)
weights = np.cos(np.deg2rad(airtemps.lat))
ts_airtemps_regional = airtemps.weighted(mask_3D * weights).mean(dim=("lat", "lon"))
ts_airtemps_regional.pre.plot(col="region", col_wrap=6);

#%%
var='tas'

filename = 'trace.01-36.22000BP.cam2.'+key[model][var]['varName']+'.22000BP_decavg'+szn+'_400BCE'
data = xr.open_dataset(dataDir+model+'/'+filename+'.nc',decode_times=False)
#Change variable names
data = data.rename({key[model]['age']['varName']: 'age'})
data = data.rename({key[model]['lat']:            'lat'})
data = data.rename({key[model]['lon']:            'lon'})
data = data.rename({key[model][var]['varName']:    var})
data = data[var]
#conversion   = key[model]['age']['conv']
#data = data.assign_coords(age = (data['age'] * conversion[0] + conversion[1]) )  

airtemps = data# xr.tutorial.load_dataset("air_temperature")

proj = ccrs.Robinson()

ax = plt.subplot(111, projection=proj)

(airtemps.isel(age=1)).plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree())

ax.coastlines();
#%%

mask_3D   = rm.defined_regions.ar6.all.mask_3D(airtemps)
mask_full = rm.defined_regions.ar6.all.mask_3D(airtemps, drop=False)
r2 = mask_3D.isel(region=(mask_3D.abbrevs == "EAS"))
airtemps_EAS = airtemps.where(r2)
proj = ccrs.Robinson()
ax = plt.subplot(111, projection=proj)
(airtemps_EAS.isel(age=5)-airtemps_EAS.isel(age=4))[:,:,0].plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree())
ax.coastlines()
#%%
airtemps_EAS = airtemps.where(r2)
weights = np.cos(np.deg2rad(airtemps_EAS.lat))
ts_airtemps_regional = airtemps_EAS.weighted(mask_3D * weights).mean(dim=("lat", "lon"))
ts_airtemps_regional.plot()

#%%

mask = regionmask.defined_regions.srex.mask_3D(lon, lat)
mask
r2 = refReg.isel(region=(refReg.abbrevs == "EAS"))
airtemps_cna = values.where(r2)

mask = refReg*1
z = values.weighted(mask * wghts).mean(dim=("lat", "lon")).data
z = pd.DataFrame(z) 
try: z.index  = list(values.age.data)
except:  z.index  = list(values.model.data)
z.columns= list(mask.abbrevs.data)
