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
scale = '_regrid'

for model in CMIP6['mh'].keys():
    print(model)
    mhAnom[model] = {}
    mh = xr.open_dataset((dataDir+'mh'+CMIP6['mh'][model]))
    pi = xr.open_dataset((dataDir+'pi'+CMIP6['pi'][model]))
    for szn in seasons.keys():
        for var in ['pre','evp','tas']:
            #Id model variable name and conversion to degC or mm/day
            if   var == 'pre': 
                varName    = 'pr'+scale
                conversion = [((1/1000)*(60*60*24*1000)),0] #converts kg/m2/s to m/s to mm/day
            elif var == 'evp': 
                varName    = 'evspsbl'+scale
                conversion = [((1/1000)*(60*60*24*1000)),0] #converts kg/m2/s to m/s to mm/day
            elif var == 'tas': 
                varName    = 'tas'+scale
                conversion = [1,-273.15] #converts K to degC
            #Load values for mh - pi using best units 
            mhVals = mh[varName][seasons[szn]] * conversion[0] + conversion[1] 
            piVals = pi[varName][seasons[szn]] * conversion[0] + conversion[1] 
            mhVals = mhVals.weighted(mh['days_per_month'][seasons[szn]]).mean(dim=("month"))
            piVals = piVals.weighted(pi['days_per_month'][seasons[szn]]).mean(dim=("month"))
            #Save difference between mh - pi
            mhAnom[model][var+'_'+szn] = (mhVals-piVals) 
        mhAnom[model]['p-e_'+szn] = mhAnom[model]['pre_'+szn] - mhAnom[model]['evp_'+szn]
    mhAnom[model]['lats'] = mh['lat'+scale]
    mhAnom[model]['lons'] = mh['lon'+scale]

#%%Calculate % agreement within CMIP


z = pd.read_csv(dataDir+'midHCpct.csv')[['V1','V2','V3']]
text_kws = dict(
    bbox=dict(color="none"),
    #path_effects=[pe.withStroke(linewidth=0, foreground="w")],
    color="#67000d",
    fontsize=0,
)


lats = mhAnom[model]['lats']
lons = mhAnom[model]['lons']
for climVar in ['p-e_DJF']:
    vals = np.zeros((len(lats),len(lons)))
    for model in CMIP6['mh'].keys():
        vals += mhAnom[model][climVar].data
        #vals += np.sign(mhAnom[model][climVar].data)

#land   = regionmask.defined_regions.natural_earth.land_110
#land   = land.mask_3D(lons,lats).data.squeeze()
#vals=(vals/2+6)/12
vals=vals/12

#vals[vals==0]=['nan']
refReg =     regionmask.defined_regions.ar6.land
z = pd.read_csv(dataDir+'midHCpct.csv')[['V1','V2','V3']]
plats = []
plons = []
for i in z['V1']: 
    loc = refReg.centroids[refReg.abbrevs.index(i)]
    plats.append(loc[1])
    plons.append(loc[0])

plt.style.use('ggplot')
plt.figure(figsize=(20,10)); plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 1; plt.rcParams['axes.edgecolor'] = 'k'
ax1 = refReg.plot(projection=ccrs.Robinson(), label="abbrev", add_ocean=True, text_kws=text_kws)
ax1.spines['geo'].set_edgecolor('black')
ax1.set_global()
ax1.coastlines()
ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
#model_contour = plt.pcolormesh(lons, lats,vals,
        #            transform=ccrs.PlateCarree(),vmin=0,vmax=1,
         #           cmap=plt.cm.get_cmap('BrBG',5),alpha=0.8)
#proxy_scatter = ax1.scatter(plons,plats,c=list((z['V3']/z['V2'])*100),
 #       s=list(z['V2']**0.5*80),edgecolor='k',lw=3,alpha=1,cmap=plt.cm.get_cmap('BrBG',5),
  #      transform=ccrs.PlateCarree(),vmin=0,vmax=100)
mlevels = np.array([-1,-2/3,-1/3,-0.1,0.1,1/3,2/3,1])*0.5
data_cyclic,lon_cyclic = cutil.add_cyclic_point(vals,
                        coord=lons)
model_contour = plt.contourf(lon_cyclic, lats, 
                  data_cyclic,transform=ccrs.PlateCarree(),
                  levels=mlevels,extend='both',cmap='BrBG')

plt.title('CMIP p-e DJF'+str(ka)+'ka-0.5ka',fontsize=20)
plt.colorbar(model_contour,cax=inset_axes(ax1,width='60%',height="4%",loc=8),
            orientation="horizontal").set_label('mm/day',fontsize=12,c='black')
           
#ax1.add_feature(cfeature.OCEAN,facecolor='white',edgecolor='k')
plt.colorbar(model_contour,cax=inset_axes(ax1,width='60%',height="4%",loc=8),
                         orientation="horizontal").set_label('% wet',fontsize=12,c='black')
           


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


