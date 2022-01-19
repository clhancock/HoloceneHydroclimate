#This script uses CMIP6 netcdf files with standard nameing/unit conventions
#For variable, a csv is produced with rows=models, cols=regions, & vals=6ka-0k
#Data are not on the same grid, but can change the variable names to _regrid

#%% Load Packages
import numpy  as np
import pandas as pd
import regionmask 
import xarray as xr

#%% Settings and names
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Data/Model/cmip6/'
saveDir = dataDir

CMIP6 = {}
CMIP6['pi'] = {
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
CMIP6['mh'] = {
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

#%% Calculate midHolocene - pi differences for each climate variable
mhAnom = {}
for model in CMIP6['mh'].keys():
    print(model)
    mhAnom[model] = {}
    mh = xr.open_dataset((dataDir+'mh'+CMIP6['mh'][model]))
    pi = xr.open_dataset((dataDir+'pi'+CMIP6['pi'][model]))
    for szn in seasons.keys():
        for var in ['pre','evp','tas']:
            #Id model variable name and conversion to degC or mm/day
            if   var == 'pre': 
                varName    = 'pr'
                conversion = [((1/1000)*(60*60*24*1000)),0] #converts kg/m2/s to m/s to mm/day
            elif var == 'evp': 
                varName    = 'evspsbl'
                conversion = [((1/1000)*(60*60*24*1000)),0] #converts kg/m2/s to m/s to mm/day
            elif var == 'tas': 
                varName    = 'tas'
                conversion = [1,-273.15] #converts K to degC
            #Load values for mh - pi using best units 
            mhVals = mh[varName][seasons[szn]] * conversion[0] + conversion[1] 
            piVals = pi[varName][seasons[szn]] * conversion[0] + conversion[1] 
            mhVals = mhVals.weighted(mh['days_per_month'][seasons[szn]]).mean(dim=("month"))
            piVals = piVals.weighted(pi['days_per_month'][seasons[szn]]).mean(dim=("month"))
            #Save difference between mh - pi
            mhAnom[model][var+'_'+szn] = (mhVals-piVals) / (piVals+mhVals)/2
        mhAnom[model]['p-e_'+szn] = mhAnom[model]['pre_'+szn] - mhAnom[model]['evp_'+szn]
    mhAnom[model]['lats'] = mh.lat
    mhAnom[model]['lons'] = mh['lon']
    
#%%Calculate weighted average by region
for climVar in ['pre_ANN','p-e_ANN','tas_ANN','pre_JJA','p-e_JJA','tas_JJA']:
    df = pd.DataFrame()
    for model in CMIP6['mh'].keys():
        vals = mhAnom[model][climVar]
        #Lat weights
        wghts = np.cos(np.deg2rad(vals.lat))
        #IPCC ref regions
        refReg = regionmask.defined_regions.ar6.land
        refReg = refReg.mask_3D(mhAnom[model]['lons'],mhAnom[model]['lats'])   
        #Land mask
        land   = regionmask.defined_regions.natural_earth.land_110
        land   = land.mask_3D(mhAnom[model]['lons'],mhAnom[model]['lats']) 
        land   = np.array([land.squeeze('region').data]*np.shape(refReg)[0])
        #3d mask with land and refRegions
        refRegLand = refReg*land
        #Average value by region
        df[model] = list(vals.weighted(refRegLand * wghts).mean(dim=("lat", "lon")).data)
    df.index      = list(refRegLand.abbrevs.data)
    df = df.transpose()
    df.to_csv(saveDir+climVar+'_cmip.csv')


