#This script formats model data to compare with each other and model simulations 
#Input:  CMIP6, netcdf files with non-standard nameing/unit conventions
#Output: Each model and season a netcdf with standard naming conventions is produced
#Notes:  Data are not on the same grid (other than CMIP6)

#1 Load Packages
import numpy             as np
import xarray            as xr
import matplotlib.pyplot as plt
import cartopy.crs       as ccrs

#%% 2 Settings and names
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/'
saveDir = dataDir + 'HoloceneHydroclimate/'
dataDir = dataDir + 'Model/'

#Dictionary of CMIP file names
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

cmip6['diff'] = {}
seasons = {'ANN':[0,1,2,3,4,5,6,7,8,9,10,11],'DJF':[0,1,11],'JJA':[5,6,7]}
#%%
monthlyMean = []
cmipEns     = {'ANN':[],'JJA':[],'DJF':[]}

for model in cmip6['mh'].keys():
    print(model)
    for time in ['mh','pi']:
        #Load data and save regrid data
        data = xr.open_dataset((dataDir+'cmip6/'+time+cmip6[time][model]))
        data = data.drop(['pr','evspsbl','tas','lat','lon','climatology_bnds'])
        data = data.rename({'pr_regrid':'pre','evspsbl_regrid':'evp','tas_regrid':'tas'})
        data = data.rename({'lat_regrid':'lat','lon_regrid':'lon'})
        for var in ['pre','evp','tas']:
            #Id model variable name and conversion to degC or mm/day
            if var == 'tas': conversion = [1,-273.15] #converts K to degC 
            else:            conversion = [((1/1000)*(60*60*24*1000)),0] #converts kg/m2/s to m/s to mm/day
            data[var] = data[var] * conversion[0] + conversion[1]
        data['p-e'] = data.pre-data.evp
        cmip6[time][model] = data
    #Calculate midHolocene Anomaly
    cmip6['diff'][model] = cmip6['mh'][model] - cmip6['pi'][model]
    cmip6['diff'][model].days_per_month.data = cmip6['mh'][model].days_per_month.data
    #Save for monthly mean
    monthlyMean.append(cmip6['diff'][model])
    #Save Seasonal values
    for szn in seasons.keys():
        vals = cmip6['diff'][model].isel(month=seasons[szn])
        vals = vals.weighted(vals.days_per_month).mean(dim=("month"))
        cmipEns[szn].append(vals.drop('days_per_month'))
        
monthlyMean = xr.concat(monthlyMean,dim='model').mean(dim=("model"))
monthlyMean.to_netcdf(saveDir+'Data/Model/midHolocene/cmip6'+'_MonthlyMean.nc')
for szn in seasons.keys(): 
    cmipEns[szn] = xr.concat(cmipEns[szn],dim='model')
    cmipEns[szn] = cmipEns[szn].assign_coords(model=list(cmip6['diff'].keys()))
    cmipEns[szn].to_netcdf(saveDir+'Data/Model/midHolocene/cmip6'+'_'+szn+'.nc')
    


#%% Plot to confirm using annual temps
pltvals = cmipEns['ANN']['pre'].mean(dim=("model"))
ax = plt.axes(projection=ccrs.PlateCarree())
plt.contourf(pltvals.lon,pltvals.lat,pltvals,transform=ccrs.PlateCarree(),
             extend='both',cmap='RdBu_r',
             levels= np.array([i /10 for i in list(range(-9,10,2))]))
ax.coastlines()
ax.set_global()
plt.show()
