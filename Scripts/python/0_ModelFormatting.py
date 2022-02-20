#This script uses cmip, hadcm & trace netcdf files with non-standard nameing/unit conventions
#For each model and season a netcdf with standard naming conventions is produced
#Data are not on the same grid (other than CMIP)

#Load Packages
import numpy      as np
import pandas     as pd
import regionmask as rm
import xarray     as xr

#%% 2 Settings and names
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/'
saveDir = dataDir + 'HoloceneHydroclimate/'
dataDir = dataDir + 'Model/'

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
                'time':{'varName':'t',              'conv':[-1,0]}, #To 0-12ka Holocene
                'pre':{'varName':'precip_mm_srf',   'conv':[((1/1000)*(60*60*24*1000)),0]}, #converts kg/m2/s to m/s to mm/day
                'evp':{'varName':'totalEvap_mm_srf','conv':[1,0]}, #already in mm/day
                'tas':{'varName':'temp_mm_srf',    'conv':[1,-273.15]}}} #converts K to degC
  
cmip6 = {}  
cmip6['pi'] = {
    'AWI-ESM-1-1-LR': "/AWI-ESM-1-1-LR_piControl_r1i1p1f1_gn_185501-195412_clim_tas_pr_evspsbl.nc",
    'CESM2':          "/CESM2_piControl_r1i1p1f1_gn_000101-120012_clim_tas_pr_evspsbl.nc",
    'EC-Earth3-LR':   "/EC-Earth3-LR_piControl_r1i1p1f1_gr_221901-241912_clim_tas_pr_evspsbl.nc",
    'FGOALS-f3-L':    "/FGOALS-f3-L_piControl_r1i1p1f1_gr_060001-116012_clim_tas_pr_evspsbl.nc",
    'FGOALS-g3':      "/FGOALS-g3_piControl_r1i1p1f1_gn_020001-089912_clim_tas_pr_evspsbl.nc",
    'GISS-E2-1-G':    "/GISS-E2-1-G_piControl_r1i1p1f1_gn_415001-500012_clim_tas_pr_evspsbl.nc",
    #'HadGEM3-GC31-LL':"/HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_185001-234912_clim_tas_pr_evspsbl.nc", #no evap
    'INM-CM4-8':      "/INM-CM4-8_piControl_r1i1p1f1_gr1_185001-238012_clim_tas_pr_evspsbl.nc",
    'IPSL-CM6A-LR':   "/IPSL-CM6A-LR_piControl_r1i1p1f1_gr_185001-209912_clim_tas_pr_evspsbl.nc",
    'MIROC-ES2L':     "/MIROC-ES2L_piControl_r1i1p1f2_gn_185001-234912_clim_tas_pr_evspsbl.nc",
    'MPI-ESM1-2-LR':  "/MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_185001-284912_clim_tas_pr_evspsbl.nc",
    'MRI-ESM2-0':     "/MRI-ESM2-0_piControl_r1i1p1f1_gn_185001-255012_clim_tas_pr_evspsbl.nc",
    #'NESM3':          "/NESM3_piControl_r1i1p1f1_gn_050001-099912_clim_tas_pr_evspsbl.nc", #no evap
    'NorESM1-F':      "/NorESM1-F_piControl_r1i1p1f1_gn_150101-170012_clim_tas_pr_evspsbl.nc",
    #'NorESM2-LM':     "/NorESM2-LM_piControl_r1i1p1f1_gn_160001-210012_clim_tas_pr_evspsbl.nc" #anom
    }
cmip6['mh'] = {
    'AWI-ESM-1-1-LR': "/AWI-ESM-1-1-LR_midHolocene_r1i1p1f1_gn_310601-320512_clim_tas_pr_evspsbl.nc",
    'CESM2':          "/CESM2_midHolocene_r1i1p1f1_gn_000101-070012_clim_tas_pr_evspsbl.nc",
    'EC-Earth3-LR':   "/EC-Earth3-LR_midHolocene_r1i1p1f1_gr_224501-244712_clim_tas_pr_evspsbl.nc",
    'FGOALS-f3-L':    "/FGOALS-f3-L_midHolocene_r1i1p1f1_gr_072001-121912_clim_tas_pr_evspsbl.nc",
    'FGOALS-g3':      "/FGOALS-g3_midHolocene_r1i1p1f1_gn_062701-112612_clim_tas_pr_evspsbl.nc",
    'GISS-E2-1-G':    "/GISS-E2-1-G_midHolocene_r1i1p1f1_gn_290001-299912_clim_tas_pr_evspsbl.nc",
    #'HadGEM3-GC31-LL':"/HadGEM3-GC31-LL_midHolocene_r1i1p1f1_gn_225001-234912_clim_tas_pr.nc",
    'INM-CM4-8':      "/INM-CM4-8_midHolocene_r1i1p1f1_gr1_188001-207912_clim_tas_pr_evspsbl.nc",
    'IPSL-CM6A-LR':   "/IPSL-CM6A-LR_midHolocene_r1i1p1f1_gr_185001-239912_clim_tas_pr_evspsbl.nc",
    'MIROC-ES2L':     "/MIROC-ES2L_midHolocene_r1i1p1f2_gn_800001-809912_clim_tas_pr_evspsbl.nc",
    'MPI-ESM1-2-LR':  "/MPI-ESM1-2-LR_midHolocene_r1i1p1f1_gn_100101-150012_clim_tas_pr_evspsbl.nc",
    'MRI-ESM2-0':     "/MRI-ESM2-0_midHolocene_r1i1p1f1_gn_195101-215012_clim_tas_pr_evspsbl.nc",
    #'NESM3':          "/NESM3_midHolocene_r1i1p1f1_gn_179801-189712_clim_tas_pr.nc",
    'NorESM1-F':      "/NorESM1-F_midHolocene_r1i1p1f1_gn_150101-170012_clim_tas_pr_evspsbl.nc",
    #'NorESM2-LM':     "/NorESM2-LM_midHolocene_r1i1p1f1_gn_210101-220012_clim_tas_pr_evspsbl.nc"
    }

seasons = {'ANN':[0,1,2,3,4,5,6,7,8,9,10,11],'DJF':[0,1,11],'JJA':[5,6,7]}

#%% 3 Load CMIP netcdf files and convert to uniform units/nameing
cmipEns = {'ANN':[],'JJA':[],'DJF':[]}
for time in ['mh','pi']:
    ens = {'ANN':[],'JJA':[],'DJF':[]}
    for model in cmip6[time].keys():
        print(model)
        data = xr.open_dataset((dataDir+'cmip6/'+time+cmip6[time][model]))
        data = data.rename({'pr_regrid':'pre_regrid','pr':'pre',
                            'evspsbl_regrid':'evp_regrid','evspsbl':'evp',
                            'tas_regrid':'tas_regrid','tas':'tas'})
        data = data.drop('climatology_bnds')
        for szn in ens.keys():
            idx = seasons[szn] #index for season
            vals = data*1 #data management
            for var in ['pre','evp','tas']:
                #Id model variable name and conversion to degC or mm/day
                if   var == 'pre': conversion = [((1/1000)*(60*60*24*1000)),0] #converts kg/m2/s to m/s to mm/day
                elif var == 'evp': conversion = [((1/1000)*(60*60*24*1000)),0] #converts kg/m2/s to m/s to mm/day
                elif var == 'tas': conversion = [1,-273.15] #converts K to degC
                for scale in ['','_regrid']:
                    vals[var+scale] = vals[var+scale] * conversion[0] + conversion[1] #convert
                    vals[var+scale] = vals[var+scale][idx,:,:].weighted(vals['days_per_month'][idx]).mean(dim=("month")) #weighted avg for season
                    vals['p-e'+scale] = vals['pre'+scale] - vals['evp'+scale] #calculate p-e
            vals = vals.drop('days_per_month')#data management
            ens[szn].append(vals) #save
    for szn in ens.keys(): # combine models into a single xarray for each season
        ens[szn] = xr.concat(ens[szn],dim='model')
        ens[szn] = ens[szn].assign_coords(model=list(cmip6[time].keys()))
        cmipEns[szn].append(ens[szn])
#%% 
for szn in ens.keys(): #combine pi and mh into single xarray with time dim
    cmipEns[szn] = xr.concat(cmipEns[szn],dim='time')
    cmipEns[szn] = cmipEns[szn].assign_coords(time=['mh','pi'])

#np.save(saveDir+'Data/Model/cmip6_modelPy.npy',cmipEns)
#%% 
for szn in ens.keys(): #combine pi and mh into single xarray with time dim
    cmipEns[szn] = cmipEns[szn].drop(['pre','evp','tas','lat','lon'])
    cmipEns[szn].to_netcdf(saveDir+'Data/Model/cmip6'+'_'+szn+'.nc')


#%% 4. Load Michael's netcdf files and convert to uniform units/nameing
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
            if   model== 'hadcm': filename = 'deglh.vn1_0.'+varName+'.monthly.'+szn+'.010yr'
            elif model== 'trace': filename = 'trace.01-36.22000BP.cam2.'+varName+'.22000BP_decavg'+szn+'_400BCE'
            data = xr.open_dataset(dataDir+model+'/'+filename+'.nc',decode_times=False)
            #Change variable names
            data = data.rename({key[model]['time']['varName']: 'time',
                                key[model]['lat']:             'lat', 
                                key[model]['lon']:             'lon', 
                                key[model][var]['varName']   : var})
            #Make uniform shape    
            if model == 'hadcm': data = data.squeeze('surface')
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
        sznData.to_netcdf(saveDir+'Data/Model/'+model+'_'+szn+'.nc')
        modelData[model][szn] = sznData

#%% 5. Calculate regional timeseries and cmip mh anomolies

def maskmodeldata(values,lats,lons):
    #Lat weights
    wghts  = np.cos(np.deg2rad(lats))
    #IPCC ref regions
    refReg = rm.defined_regions.ar6.land
    refReg = refReg.mask_3D(lons,lats)   
    #Land mask
    land   = rm.defined_regions.natural_earth.land_110
    land   = land.mask_3D(lons,lats)   
    land   = np.array([land.squeeze('region').data]*np.shape(refReg)[0])
    #3d mask with land and refRegions
    mask = refReg*land
    #Average value by region
    out = values.weighted(mask * wghts).mean(dim=("lat", "lon")).data
    if sum(np.shape(out)) == 46: 
        out = pd.DataFrame(out) 
        out.index = list(mask.abbrevs.data)
    else: out = pd.DataFrame(out, columns = list(mask.abbrevs.data))
    return(out)

for model in ['hadcm','trace','cmip6']:
    for szn in ['ANN','JJA','DJF']:
        data = xr.open_dataset(saveDir+'Data/Model/'+model+'_'+szn+'.nc',decode_times=False)
        for var in ['pre','p-e','tas']:
            if model == 'cmip6':
                df = pd.DataFrame()
                for n in range(len(data.model)):
                    vals = data[var][0,n,:,:]-data[var][1,n,:,:]
                    x = [i for i, w in enumerate(np.nanmean(vals,axis=0)) if np.isnan(w) == False]
                    y = [i for i, w in enumerate(np.nanmean(vals,axis=1)) if np.isnan(w) == False]
                    vals = vals[y,x]
                    df[str(data.model[n].data)] = maskmodeldata(vals,vals.lat,vals.lon)[0]
                df=df.transpose()
            else: df = maskmodeldata(data[var],data[var].lat,data[var].lon)
            df.to_csv(saveDir+'Data/Model/RegionTS/'+'regional_'+var+'_'+szn+'_'+model+'.csv')
            
