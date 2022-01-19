#This script uses hadcm & trace netcdf files with non-standard nameing/unit conventions
#For model & variable, a csv is produced with rows=age, cols=regions, & vals=mm/a or degC
#Data are not on the same grid

#%% Load Packages
import numpy  as np
import pandas as pd
import regionmask 
import xarray as xr

#%% Settings and names
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Data/Model/'
saveDir = dataDir

#Set time variables and resolution of data
ageMin=0; ageMax=12000; ageRes=100
binvec = list(range(ageMin-50,ageMax+51,ageRes))
binyrs = list(range(ageMin,ageMax+1,ageRes))

#Model variable names and conversion to common time/variable units (generally mm/day)
key = {'trace':{'lat':'lat',
                'lon':'lon',
                'time':{'varName':'time', 'conv':[-1000,0]}, #To 0-12ka Holocene
                'pre':{'varName':'PRECT', 'conv':[(60*60*24*1000),0]}, #converts m/s to mm/day
                'evp':{'varName':'QFLX',  'conv':[((1/1000)*(60*60*24*1000)),0]}, #converts kg/m2/s to m/s to mm/day
                'tas':{'varName':'TREFHT','conv':[1,-273.15]}}, #converts K to degC
       'hadcm':{'lat':'latitude',
                'lon':'longitude',
                'time':{'varName':'t',              'conv':[-1,1950]}, #To 0-12ka Holocene
                'pre':{'varName':'precip_mm_srf',   'conv':[((1/1000)*(60*60*24*1000)),0]}, #converts kg/m2/s to m/s to mm/day
                'evp':{'varName':'totalEvap_mm_srf','conv':[1,0]}, #already in mm/day
                'tas':{'varName':'temp_mm_1_5m',    'conv':[1,-273.15]}}} #converts K to degC

#%% Load Michael's netcdf files and convert to uniform units/nameing
modelData = {}
for model in key.keys():
    print(model)
    modelData[model] = {}
    for szn in ['ANN','DJF','JJA']:
        print(szn)
        sznData = {}
        for var in ['pre','evp','tas']:
            #Load Data
            varName = key[model][var]['varName']
            if   model== 'hadcm': filename = 'deglh.vn1_0.'+varName+'.monthly.'+szn
            elif model== 'trace': filename = 'trace.01-36.22000BP.cam2.'+varName+'.22000BP_decavg'+szn+'_400BCE'
            data = xr.open_dataset(dataDir+model+'/'+filename+'.nc',decode_times=False)
            #Change variable names
            data = data.rename({key[model]['time']['varName']: 'time',
                                key[model]['lat']:             'lat', 
                                key[model]['lon']:             'lon', 
                                key[model][var]['varName']   : var})
            #Make uniform size    
            if model == 'hadcm': 
                if var == 'tas': data = data.squeeze('ht')
                else:            data = data.squeeze('surface')
            #Value Conversion
            conversion = key[model][var]['conv']
            data[var]  = data[var] * conversion[0] + conversion[1]        
            #Time Conversion
            conversion   = key[model]['time']['conv']
            data['time'] = data['time'] * conversion[0] + conversion[1]    
            #Bin to 100 yr resolution
            databin = data.groupby_bins('time',binvec).mean(dim="time")
            databin = databin.rename({'time_bins': 'time'})
            databin['time'] = binyrs
            #Save
            sznData[var+'_'+szn] = databin
        #Merge into a single xarray
        sznData = xr.merge([sznData['pre_'+szn],sznData['evp_'+szn],sznData['tas_'+szn]])  
        sznData['p-e'] = sznData.pre - sznData.evp
        modelData[model][szn] = sznData

#np.save(str(dataDir+'HoloceneHydroclimate/'+'Data/'+'Model/modelPy.npy'),modelData)

#%% Save regional timeseries for each climate variable as an individual csv
szn = 'JJA'
for model in key.keys():
    for var in ['pre','p-e','tas']:
        vals   = modelData[model][szn][var]
        #Lat weights
        wghts  = np.cos(np.deg2rad(vals.lat))
        #IPCC ref regions
        refReg = regionmask.defined_regions.ar6.land
        refReg = refReg.mask_3D(vals.lon,vals.lat)   
        #Land mask
        land   = regionmask.defined_regions.natural_earth.land_110
        land   = land.mask_3D(vals.lon,vals.lat)   
        land   = np.array([land.squeeze('region').data]*np.shape(refReg)[0])
        #3d mask with land and refRegions
        refRegLand = refReg*land
        #Average value by region
        df = vals.weighted(refRegLand * wghts).mean(dim=("lat", "lon")).data
        df = pd.DataFrame(df, columns = list(refRegLand.abbrevs.data))
        df.index = binyrs
        df.to_csv(saveDir+var+'_'+szn+'_'+model+'.csv')



