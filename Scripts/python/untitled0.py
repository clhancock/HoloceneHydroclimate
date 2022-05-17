#%% 1 Load Packages
import numpy      as np
import pandas     as pd
import regionmask as rm
import xarray     as xr 
from scipy        import stats
#from scipy.stats  import pearsonr
#For Figures:
import matplotlib.pyplot   as plt         # Packages for making figures
import matplotlib.gridspec as gridspec
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs         as ccrs        # Packages for mapping in python
import cartopy.feature     as cfeature
import cartopy.util        as cutil
  
#Plot Settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = '8'
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 0.5; 
plt.rcParams['axes.edgecolor'] = 'k'
#plt.tick_params(labelsize=8)
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


#%% 2 Load Data


  
for region in ['EAS','SAH','SAS','WAF']:
    monthly = {}
    spatial = {}  
    for t in ['mh','pi']:
        n=0
        for model in cmip6[t].keys():
            print(model)
            data = xr.open_dataset((dataDir+'cmip6/'+t+cmip6[t][model]))
            #Calculate Monthly Values
            tas  = data.tas_regrid.rename({'lat_regrid':'lat','lon_regrid':'lon'})
            pre  = data.pr_regrid.rename({'lat_regrid':'lat','lon_regrid':'lon'})
            if n == 0: 
                spatial[t+'_tas'] = tas
                spatial[t+'_pre'] = pre
                #months = data.days_per_month
            else: 
                spatial[t+'_tas'] += tas
                spatial[t+'_pre'] += pre
            n+=1
        #
        lats = data.lat_regrid
        lons = data.lon_regrid
        wghts  = np.cos(np.deg2rad(lats))
        #IPCC ref regions
        refReg = rm.defined_regions.ar6.all.mask_3D(lons,lats)   
        land   = rm.defined_regions.natural_earth.land_110.mask_3D(lons,lats)   
        mask = refReg * land.squeeze(drop=True)
        for var in ['tas','pre']:
            spatial[t+'_'+var] /= n
            spatial[t+'_'+var] = spatial[t+'_'+var] * mask
            spatial[t+'_'+var] = spatial[t+'_'+var].where(mask.isel(region=(mask.abbrevs == region)))
            monthly[t+'_'+var] = spatial[t+'_'+var].mean(dim=("lat","lon","region"))
            spatial[t+'_'+var+'1'] = spatial[t+'_'+var][[0,1,11],:,:,:].mean(dim=("month","region"))
            spatial[t+'_'+var+'7'] = spatial[t+'_'+var][[5,6,7],:,:,:].mean(dim=("month","region"))
            spatial[t+'_'+var] = spatial[t+'_'+var].mean(dim=("month","region"))
    tas = (monthly['mh_tas']-monthly['pi_tas'])
    pre = (monthly['mh_pre']-monthly['pi_pre'])*((1/1000)*(60*60*24*1000))
    var = 'pre'
    val1 = (spatial['mh_'+var+'1']-spatial['pi_'+var+'1'])*((1/1000)*(60*60*24*1000))
    val2 = (spatial['mh_'+var+'7']-spatial['pi_'+var+'7'])*((1/1000)*(60*60*24*1000))
    months = [1,2,3,4,5,6,7,8,9,10,11,12]
    if var == 'tas': c = 'RdBu_r' 
    else: c = 'BrBG' 
    plt.style.use('seaborn-dark')
    plt.figure(figsize=(10,5))
    gs = gridspec.GridSpec(2,10)
    ax1 = plt.subplot(gs[0:2,0:4]) 
    ax3 = plt.subplot(gs[0:1,4:11],projection=ccrs.PlateCarree()) 
    ax4 = plt.subplot(gs[1:2,4:11],projection=ccrs.PlateCarree()) 
    ax3.coastlines(); ax4.coastlines()
    levels=[-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8]
    val1.plot.pcolormesh(ax=ax3,transform=ccrs.PlateCarree(),cmap=c,#add_colorbar=False,
                         levels=levels,extend='both')
    val2.plot.pcolormesh(ax=ax4,transform=ccrs.PlateCarree(),cmap=c,#add_colorbar=False,
                         levels=levels,extend='both')
    
    ax3.set_title('DJF Precip')
    ax4.set_title('JJA Precip')
    
    ax2 = ax1.twinx()
    ax2.bar(months, pre, color ='royalblue',width = 1,alpha=0.8)
    ax1.plot(months, tas, '-o',color='Red')
    ax1.set_ylim(-2.5,2.5)
    ax2.set_ylim(-2.5,2.5)
    ax1.set_xlim(0.5,12.5)
    ax2.set_xlim(0.5,12.5)
    ax1.set_ylabel('Temperature (degC)')
    ax2.set_ylabel('Precipitation (mm/day)')
    ax1.set_xlabel('Month')
    if   region == 'EAS': 
        ax3.set_extent((98, 150, 16, 48), crs=ccrs.PlateCarree())
        ax4.set_extent((98, 150, 16, 48), crs=ccrs.PlateCarree())
    elif region == 'SAS': 
        ax3.set_extent((58, 102, 0, 35), crs=ccrs.PlateCarree())
        ax4.set_extent((58, 102, 0, 35), crs=ccrs.PlateCarree())
    elif region == 'SAH': 
        ax3.set_extent((-22, 40, 3, 37), crs=ccrs.PlateCarree())
        ax4.set_extent((-22, 40, 3, 37), crs=ccrs.PlateCarree())
    elif region == 'WAF': 
        ax3.set_extent((-22, 17, -5, 17), crs=ccrs.PlateCarree())
        ax4.set_extent((-22, 17, -5, 17), crs=ccrs.PlateCarree())
    else: 
        ax3.set_global()
        ax4.set_global()
    refReg = rm.defined_regions.ar6.all
    refReg[region,region].plot_regions(ax=ax3,add_label=False,line_kws=dict(linewidth=1.5))
    refReg[region,region].plot_regions(ax=ax4,add_label=False,line_kws=dict(linewidth=1.5))
    ax1.plot([0,12.5],[0,0],'k',lw=2)
    plt.title('(MidHolocene minus PIcontrol)',size=14,loc='left',fontweight="bold")
    plt.suptitle(region,size=20,x=0.16,fontweight="bold")
    plt.show()
