#This script formats model data to compare with each other and model simulations 
#Input:  CMIP6, netcdf files with non-standard naming/unit conventions (provided by MP Erb)
#Output: Each model and season a netcdf with standard naming conventions is produced (Available on clhancock github)
#Notes:  Transient data are not on the same grid (other than CMIP6)
         #Updated on 4.24 to use new files provided by MP Erb

#1 Load Packages
import numpy             as np
import xarray            as xr
import matplotlib.pyplot as plt
import cartopy.crs       as ccrs

#%% 2 Settings and names
dataDir = '/Users/chrishancock/Library/CloudStorage/OneDrive-NorthernArizonaUniversity/Research/Manuscript/'
saveDir = dataDir + 'HoloceneHydroclimate/HoloceneHydroclimate/'
dataDir = dataDir + 'Model/'


#Dictionary of CMIP file names
cmip6 = {}  
cmip6['pi'] = {
    'AWI-ESM-1-1-LR': "/AWI-ESM-1-1-LR_piControl_r1i1p1f1_gn_185501-195412_",
    'CESM2':          "/CESM2_piControl_r1i1p1f1_gn_000101-120012_",
    'EC-Earth3-LR':   "/EC-Earth3-LR_piControl_r1i1p1f1_gr_221901-241912_",
    'FGOALS-f3-L':    "/FGOALS-f3-L_piControl_r1i1p1f1_gr_060001-116012_",
    'FGOALS-g3':      "/FGOALS-g3_piControl_r1i1p1f1_gn_020001-089912_",
    'GISS-E2-1-G':    "/GISS-E2-1-G_piControl_r1i1p1f1_gn_415001-500012_",
    #'HadGEM3-GC31-LL':"/HadGEM3-GC31-LL_piControl_r1i1p1f1_gn_185001-234912_", #no evap
    'INM-CM4-8':      "/INM-CM4-8_piControl_r1i1p1f1_gr1_185001-238012_",
    'IPSL-CM6A-LR':   "/IPSL-CM6A-LR_piControl_r1i1p1f1_gr_185001-384912_",
    'MIROC-ES2L':     "/MIROC-ES2L_piControl_r1i1p1f2_gn_185001-234912_",
    'MPI-ESM1-2-LR':  "/MPI-ESM1-2-LR_piControl_r1i1p1f1_gn_185001-284912_",
    'MRI-ESM2-0':     "/MRI-ESM2-0_piControl_r1i1p1f1_gn_185001-255012_",
    #'NESM3':          "/NESM3_piControl_r1i1p1f1_gn_050001-099912_", #no evap
    'NorESM1-F':      "/NorESM1-F_piControl_r1i1p1f1_gn_150101-170012_",
    #'NorESM2-LM':     "/NorESM2-LM_piControl_r1i1p1f1_gn_160001-210012_" #anom
    }
cmip6['mh'] = {
    'AWI-ESM-1-1-LR': "/AWI-ESM-1-1-LR_midHolocene_r1i1p1f1_gn_310601-320512_",
    'CESM2':          "/CESM2_midHolocene_r1i1p1f1_gn_000101-070012_",
    'EC-Earth3-LR':   "/EC-Earth3-LR_midHolocene_r1i1p1f1_gr_224501-244712_",
    'FGOALS-f3-L':    "/FGOALS-f3-L_midHolocene_r1i1p1f1_gr_072001-121912_",
    'FGOALS-g3':      "/FGOALS-g3_midHolocene_r1i1p1f1_gn_062701-112612_",
    'GISS-E2-1-G':    "/GISS-E2-1-G_midHolocene_r1i1p1f1_gn_290001-299912_",
    #'HadGEM3-GC31-LL':"/HadGEM3-GC31-LL_midHolocene_r1i1p1f1_gn_225001-234912_clim_tas_pr.nc",
    'INM-CM4-8':      "/INM-CM4-8_midHolocene_r1i1p1f1_gr1_188001-207912_",
    'IPSL-CM6A-LR':   "/IPSL-CM6A-LR_midHolocene_r1i1p1f1_gr_185001-239912_",
    'MIROC-ES2L':     "/MIROC-ES2L_midHolocene_r1i1p1f2_gn_800001-809912_",
    'MPI-ESM1-2-LR':  "/MPI-ESM1-2-LR_midHolocene_r1i1p1f1_gn_100101-150012_",
    'MRI-ESM2-0':     "/MRI-ESM2-0_midHolocene_r1i1p1f1_gn_195101-215012_",
    #'NESM3':          "/NESM3_midHolocene_r1i1p1f1_gn_179801-189712_clim_tas_pr.nc",
    'NorESM1-F':      "/NorESM1-F_midHolocene_r1i1p1f1_gn_150101-170012_",
    #'NorESM2-LM':     "/NorESM2-LM_midHolocene_r1i1p1f1_gn_210101-220012_"
    }

#
# Load, standardize, and save data for comparison with transient
monthlyMean = []
cmip6['diff'] = {}
seasons = {'ANN':[0,1,2,3,4,5,6,7,8,9,10,11],'DJF':[0,1,11],'JJA':[5,6,7]}
cmipEns     = {'ANN':[],'JJA':[],'DJF':[]}

for model in cmip6['mh'].keys():
    print(model)
    for time in ['mh','pi']:
        #Load data
        data = xr.open_dataset(dataDir+'cmip6/'+time+cmip6[time][model]+'regridclim.nc')
        print(time)
        #
        print(np.shape(data.pr))
        try: data=data.drop('height')
        except: print(time)
        #Standard Names (for comparison with transient)
        data = data.rename({'pr':'pre','evspsbl':'evp','tas':'tas'})
        #Standard Units (for comparison with transient)
        for var in ['pre','evp','tas']:
            #Id model variable name and conversion to degC or mm/day
            if var == 'tas': conversion = [1,-273.15]                    #converts K to degC 
            else:            conversion = [((1/1000)*(60*60*24*1000)),0] #converts kg/m2/s to m/s to mm/day
            data[var] = data[var] * conversion[0] + conversion[1]
        #Add P-E
        data['p-e'] = data.pre-data.evp
        #Add Units
        data.tas.attrs['unit']                                  = 'mean degC'
        for var in ['pre','p-e','evp']: data[var].attrs['unit'] = 'mean mm/day'
        #Compile data
        cmip6[time][model] = data
        data.close()
    #Calculate midHolocene Anomaly
    cmip6['diff'][model] = cmip6['mh'][model] - cmip6['pi'][model]
    cmip6['diff'][model].days_in_months.data = cmip6['mh'][model].days_in_months.data
    #Save for monthly mean
    monthlyMean.append(cmip6['diff'][model])
    #Save Seasonal means
    for szn in seasons.keys():
        vals = cmip6['diff'][model].isel(month=seasons[szn])
        vals = vals.weighted(vals.days_in_months).mean(dim=("month"))
        cmipEns[szn].append(vals.drop('days_in_months'))

#%%Save to_netcdf
monthlyMean = xr.concat(monthlyMean,dim='model').mean(dim=("model"))
monthlyMean.to_netcdf(saveDir+'Data/Model/cmip6/cmip6'+'_MonthlyMean.nc')
for szn in seasons.keys(): 
    cmipEns[szn] = xr.concat(cmipEns[szn],dim='model')
    cmipEns[szn] = cmipEns[szn].assign_coords(model=list(cmip6['diff'].keys()))
    cmipEns[szn].to_netcdf(saveDir+'Data/Model/cmip6/cmip6'+'_'+szn+'.nc')
    


#%% Plot to confirm using annual temps
pltvals = cmipEns['ANN']['tas'].mean(dim=("model"))
ax = plt.axes(projection=ccrs.PlateCarree())
plt.contourf(pltvals.lon,pltvals.lat,pltvals,transform=ccrs.PlateCarree(),
             extend='both',cmap='RdBu_r',
             levels= np.array([i /10 for i in list(range(-9,10,2))]))
plt.title('MH - PI')
ax.coastlines()
ax.set_global()
plt.show()
