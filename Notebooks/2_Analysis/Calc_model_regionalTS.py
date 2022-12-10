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

#Plot Settings
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 0.5; 
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams.update({'font.size': 10})

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








