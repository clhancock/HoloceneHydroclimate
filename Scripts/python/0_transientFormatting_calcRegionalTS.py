#This script formats model data to compare with each other and model simulations 
#Input:  HadCM & TraCE netcdf files with non-standard nameing/unit conventions
#Output: Each model and season a netcdf with standard naming conventions is produced
#        With both original and new (regrid) spatial resolution
#1 Load Packages
import numpy             as np
#import xesmf             as xe
import xarray            as xr
import matplotlib.pyplot as plt
import cartopy.crs       as ccrs

dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/'
saveDir = dataDir + '2021_HoloceneHydroclimate/2021_HoloceneHydroclimate/'
dataDir = dataDir + 'Model/'

#%% 2 Transient Models 
#Model variable names and conversion to common time/variable units (generally mm/day)
key = {'trace':{'lat':'lat',
                'lon':'lon',
                'age':{'varName':'time', 'conv':[-1000,0]}, #To 0-12ka Holocene
                'pre':{'varName':'PRECT', 'conv':[(60*60*24*1000),0]}, #converts m/s to mm/day
                'evp':{'varName':'QFLX',  'conv':[((1/1000)*(60*60*24*1000)),0]}, #converts kg/m2/s to m/s to mm/day
                'tas':{'varName':'TREFHT','conv':[1,-273.15]}}, #converts K to degC
       'hadcm':{'lat':'latitude',
                'lon':'longitude',
                'age':{'varName':'t',              'conv':[-1,0]}, #To 0-12ka Holocene
                'pre':{'varName':'precip_mm_srf',   'conv':[((1/1000)*(60*60*24*1000)),0]}, #converts kg/m2/s to m/s to mm/day
                'evp':{'varName':'totalEvap_mm_srf','conv':[1,0]}, #already in mm/day
                'tas':{'varName':'temp_mm_srf',     'conv':[1,-273.15]}}} #converts K to degC

binsize = 100
binvec = list(range(0-int(binsize/2),int(12000+1+binsize/2),binsize))
binyrs = list(range(0,12000+1,binsize))

for model in ['trace','hadcm']:
    for szn in ['ANN','DJF','JJA']:
         print(szn)
         newData = {}
         for var in ['pre','evp','tas']:
             #Load Data
             varName = key[model][var]['varName']
             if   model== 'hadcm': filename = 'deglh.vn1_0.'+varName+'.monthly.'+szn+'.010yr'
             elif model== 'trace': filename = 'trace.01-36.22000BP.cam2.'+varName+'.22000BP_decavg'+szn+'_400BCE'
             data = xr.open_dataset(dataDir+model+'/'+filename+'.nc',decode_times=False)
             #Change variable names
             data = data.rename({key[model]['age']['varName']: 'age'})
             data = data.rename({key[model]['lat']:            'lat'})
             data = data.rename({key[model]['lon']:            'lon'})
             data = data.rename({key[model][var]['varName']:    var})
             data = data[var]
             print(data.units)
             #Make uniform shape    
             if model == 'hadcm': data = data.squeeze('surface')
             #Uniform Time Conversion
             conversion   = key[model]['age']['conv']
             data = data.assign_coords(age = (data['age'] * conversion[0] + conversion[1]) )  
             #Bin to 100 yr resolution
             data = data.groupby_bins('age',binvec).mean(dim="age")
             data = data.rename({'age_bins': 'age'})
             data['age'] = binyrs  
                 #minAge = np.argmin(abs(data.age.data-12000))
                 #maxAge = np.argmin(abs(data.age.data-0))
                 #data = data[var][minAge:(maxAge+1),::]
             #Uniform Value Conversion 
             conversion = key[model][var]['conv']
             data  = data * conversion[0] + conversion[1]   
             #Store Data
             newData[var] = data
         #Merge into a single xarray
         newData = xr.merge([newData['pre'],newData['evp'],newData['tas']])  
         newData['p-e'] = newData.pre - newData.evp
         #Add metadata
         newData.tas.attrs['unit']                                  = szn+' mean degC'
         for var in ['pre','p-e','evp']: newData[var].attrs['unit'] = szn+' mean mm/day'
         for var in ['pre','evp','tas']: newData[var].attrs['original_name'] = key[model][var]['varName']
         newData.age.attrs['unit']                                  = 'yr BP'
         #Save Data
         newData.to_netcdf(saveDir+'Data/Model/'+model+'/'+model+'_'+szn+'.nc')

#%%Regrid and save tansient data with cmip spatial resolution
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
        #regridder = xe.Regridder(dataOrig,dataReGr,'conservative_norm',periodic=True,ignore_degenerate=True)
        #regridder = xe.Regridder(dataOrig,dataReGr,'bilinear')         
        dataReGr = regridder(dataOrig,keep_attrs=True)
        #Plot to make sure it looks right
        print(np.max(abs(dataReGr['tas'].data)))
        print(np.max(abs(dataOrig['tas'].data)))
        if make_figure:
            ax1 = plt.subplot2grid((2,1),(0,0),projection=ccrs.PlateCarree())
            ax2 = plt.subplot2grid((2,1),(1,0),projection=ccrs.PlateCarree())
            dataOrig['tas'].isel(age=0).plot.pcolormesh(ax=ax1)
            dataReGr['tas'].isel(age=0).plot.pcolormesh(ax=ax2)
            ax1.coastlines(); ax2.coastlines()
            plt.plot()
            plt.show()
        #Rename to indicate regrid
        for var in dataReGr.keys(): dataReGr = dataReGr.rename({var:var+'_regrid'})
        dataReGr = dataReGr.rename({'lat':'lat_regrid','lon':'lon_regrid'})
        dataReGr.to_netcdf(dataDir+'Data/Model/'+model+'_'+szn+'_regrid'+'.nc')


    
