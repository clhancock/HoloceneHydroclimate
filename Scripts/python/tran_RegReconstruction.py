#This script uses hadcm & trace netcdf files with non-standard nameing/unit conventions
#For model & variable, a csv is produced with rows=age, cols=regions, & vals=mm/a or degC
#Data are not on the same grid

#%% 1 Load Packages
import numpy  as np
import pandas as pd
import regionmask 
import xarray as xr
import matplotlib.pyplot as plt         # Packages for making figures
import matplotlib.gridspec as gridspec
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs              # Packages for mapping in python
import cartopy.feature as cfeature
import cartopy.util as cutil
from scipy import stats

#%% 2 Settings and names
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Data/Model/'
saveDir = dataDir

#Set time variables and resolution of data
ageMin=0; ageMax=18000; ageRes=100
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


#%% 3 Load Michael's netcdf files and convert to uniform units/nameing
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
        modelData[model][szn] = sznData

np.save(str(dataDir+'HoloceneHydroclimate_modelPy.npy'),modelData)

#%% 4 Save regional timeseries for each climate variable as an individual csv
modelData = np.load(str(dataDir+'HoloceneHydroclimate_modelPy.npy'),allow_pickle=True).item()
szn = 'ANN'
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


#%% 5 MAP DATA BY 1KA BINS
for ka in [*range(6,12,6)]:
    for model in ['hadcm','trace']:
        for climVar in ['pre_ANN','pre_JJA','pre_DJF']:
            name = model+' '+climVar+' '+str(ka)+'ka-0.5ka'
            var  = climVar[:3]
            szn  = climVar[len(climVar)-3:]
            data = modelData[model][szn][var]
            dataLM = data.groupby_bins('time',[0,1000]).mean(dim='time')
            dataKa = data.groupby_bins('time',[ka*1000-500,ka*1000+500]).mean(dim='time')
            dataAn = (dataKa.data-dataLM.data)#/((dataKa.data+dataLM.data))
            plt.style.use('ggplot')
            plt.figure(figsize=(20,10)); plt.rcParams['axes.facecolor'] ='white'
            plt.rcParams['axes.linewidth'] = 1; 
            plt.rcParams['axes.edgecolor'] = 'k'
            ax1 = plt.subplot(projection=ccrs.Robinson()) 
            ax1.spines['geo'].set_edgecolor('black')
            ax1.set_global()
            ax1.coastlines()
            ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
            ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
            data_cyclic,lon_cyclic = cutil.add_cyclic_point(dataAn[0,:,:],
                                                            coord=data.lon.data)
            mlevels = np.array([-1,-2/3,-1/3,-0.1,0.1,1/3,2/3,1])*0.5
            model_contour = plt.contourf(lon_cyclic, data.lat.data, 
                                         data_cyclic,transform=ccrs.PlateCarree(),
                                         levels=mlevels,extend='both',cmap='BrBG')
            plt.title(name,fontsize=30)
            plt.colorbar(model_contour,cax=inset_axes(ax1,width='60%',height="4%",loc=8),
                         orientation="horizontal").set_label('mm/day',fontsize=12,c='black')
            plt.colorbar(model_contour,cax=inset_axes(ax1,width='60%',height="4%",loc=8),
                         orientation="horizontal").set_label('mm/day',fontsize=12,c='black')
            ###########################
            plt.show()
            #plt.savefig(dataDir+'figs/'+name[:-6]+'.png',dpi=400,format='png',bbox_inches='tight')       

#%% 6 Plot Slopes
mlevels = np.array([-1,-2/3,-1/3,-0.1,-0.06,0.06,0.1,1/3,2/3,1])*0.1

for climVar in ['p-e_JJA']:#,'pre_JJA','p-e_ANN','p-e_JJA']:#,'pre_DJF','pre_JJA']:
    var  = climVar[:3]
    szn  = climVar[len(climVar)-3:]
    plt.style.use('ggplot')
    plt.figure(figsize=(10,22))
    plt.rcParams['axes.facecolor'] ='white'
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['axes.edgecolor'] = 'k'
    plt.title(climVar,fontsize=8)
    gs = gridspec.GridSpec(6,2)
    n = 0
    for ka in [*range(6,18,6)]:
        #name = model+' '+climVar+' '+str(ka)+'ka-0.5ka'
        n += 1
        for model in ['trace','hadcm']:
            if   model == 'hadcm':
                ax1 = plt.subplot(gs[n-1:n,1:2],projection=ccrs.Robinson()) 
            elif model == 'trace': 
                ax1 = plt.subplot(gs[n-1:n,0:1],projection=ccrs.Robinson()) 
            ax1.spines['geo'].set_edgecolor('black')
            ax1.set_global()
            ax1.coastlines()
            ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
            ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
            #ax1 = refReg.plot(projection=ccrs.Robinson(), label="abbrev", add_ocean=True, text_kws=text_kws)
            data = modelData[model][szn][var]
            dataSl = np.full(np.shape(data)[1:3],np.NaN)
            for i in range(np.shape(data)[1]):
                for j in range(np.shape(data)[2]):
                    modelReg = stats.linregress(data[(n-1)*60:(n*60)+1,i,j].time.data,data.data[(n-1)*60:(n*60)+1,i,j])
                    dataSl[i,j] = modelReg[0]*-1000
            data_cyclic,lon_cyclic = cutil.add_cyclic_point(dataSl,
                                                            coord=data.lon.data)
            model_contour = plt.contourf(lon_cyclic, data.lat.data, 
                                         data_cyclic,transform=ccrs.PlateCarree(),
                                         levels=mlevels,extend='both',cmap='BrBG')                            
            plt.title(model+' '+str(ka)+'ka-'+str(ka-6)+'ka '+climVar,fontsize=8)
            #plt.colorbar(model_contour,cax=inset_axes(ax1,width='60%',height="8%",loc=8),
                 #        orientation="horizontal").set_label('mm/day/1ka',fontsize=8,c='black')
    plt.show()
    #plt.savefig(dataDir+'figs/'+climVar+'.png',dpi=400,format='png',bbox_inches='tight')       

#%% 7 Plot correlation between 2 transient models? 
mlevels = np.array([-1,-2/3,-1/3,-0.1,0.1,1/3,2/3,1])

from scipy.stats import pearsonr
for climVar in ['p-e_ANN']:#,'pre_JJA','p-e_ANN','p-e_JJA']:#,'pre_DJF','pre_JJA']:
    var  = climVar[:3]
    szn  = climVar[len(climVar)-3:]
    dataTr = modelData['trace'][szn][var].groupby_bins('time',range(0,6001,100)).mean(dim="time")
    dataHa = modelData['hadcm'][szn][var].groupby_bins('time',range(0,6001,100)).mean(dim="time")
    hadaRg = np.full(np.shape(dataTr)[1:3],np.NaN)
    land   = regionmask.defined_regions.natural_earth.land_110
    land   = land.mask_3D(dataTr.lon,dataTr.lat)  
    for i in range(np.shape(dataTr)[1]):
        for j in range(np.shape(dataTr)[2]):
            if land.data[0,i,j] == True:
                lat = np.argmin(np.abs(dataTr.lat.data[i] - dataHa.lat.data))                  
                lon = np.argmin(np.abs(dataTr.lon.data[j] - dataHa.lon.data))
                hadaRg[i,j] = pearsonr(dataTr[:,i,j].data,dataHa[:,lat,lon].data)[0]
    plt.style.use('ggplot')
    plt.figure(figsize=(12,20))
    plt.rcParams['axes.facecolor'] ='white'
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['axes.edgecolor'] = 'k'
    plt.title(climVar,fontsize=12)
    gs = gridspec.GridSpec(6,2)
    ax1 = plt.subplot(gs[n-1:n,0:1],projection=ccrs.Robinson()) 
    ax1.spines['geo'].set_edgecolor('black')
    ax1.set_global()
    ax1.coastlines()
    ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
    ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
    #dataSlP= dataSl
    data_cyclic,lon_cyclic = cutil.add_cyclic_point(hadaRg,
                                                            coord=dataTr.lon.data)
    model_contour = plt.contourf(lon_cyclic, dataTr.lat.data, 
                                         data_cyclic,transform=ccrs.PlateCarree(),
                                         levels=mlevels,extend='both',cmap='BrBG') 
    plt.title('agreement between hadcm/trace '+climVar,fontsize=12)
    plt.colorbar(model_contour,cax=inset_axes(ax1,width='60%',height="4%",loc=8),ticks=range(0,13),
                 orientation="horizontal").set_label('ka',fontsize=12,c='black')
    plt.show()
print(np.nanmean(hadaRg))
#%% When is the wettest century

for climVar in ['p-e_ANN']:#,'pre_JJA','p-e_ANN','p-e_JJA']:#,'pre_DJF','pre_JJA']:
        var  = climVar[:3]
        szn  = climVar[len(climVar)-3:]
        plt.style.use('ggplot')
        plt.figure(figsize=(12,20))
        plt.rcParams['axes.facecolor'] ='white'
        plt.rcParams['axes.linewidth'] = 1
        plt.rcParams['axes.edgecolor'] = 'k'
        plt.title(climVar,fontsize=12)
        gs = gridspec.GridSpec(6,2)
        n = 0
            #name = model+' '+climVar+' '+str(ka)+'ka-0.5ka'
        n += 1
        nn = 0
        nnn = 0
        for model in ['trace','hadcm']:
            if   model == 'trace':
                ax1 = plt.subplot(gs[n-1:n,1:2],projection=ccrs.Robinson()) 
            elif model == 'hadcm': 
                ax1 = plt.subplot(gs[n-1:n,0:1],projection=ccrs.Robinson()) 
            ax1.spines['geo'].set_edgecolor('black')
            ax1.set_global()
            ax1.coastlines()
            ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
            ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
            data = modelData[model][szn][var][0:121,:,:]
            dataPk = np.full(np.shape(data)[1:3],np.NaN)
            #dataSlP= dataSl
            for i in range(np.shape(data)[1]):
                for j in range(np.shape(data)[2]):
                    dataTS = data[0:121,i,j].data
                    dataPk[i,j] = np.mean([c for c, d in enumerate(dataTS) if d == np.nanmax(dataTS)])                    
                    nnn+=1
                    if dataPk[i,j] < 10: nn+=1
                    if dataPk[i,j] > 110: nn+=1
                    #dataSlP[i,j] = modelReg[3]
            data_cyclic,lon_cyclic = cutil.add_cyclic_point(dataPk/10,
                                                            coord=data.lon.data)
            model_contour = plt.contourf(lon_cyclic, data.lat.data, 
                                         data_cyclic,transform=ccrs.PlateCarree(),
                                         cmap=plt.cm.get_cmap('PuOr',12))
            plt.title(model+' timing of wettest century '+climVar,fontsize=12)
            plt.colorbar(model_contour,cax=inset_axes(ax1,width='60%',height="4%",loc=8),ticks=range(0,13),
                         orientation="horizontal").set_label('ka',fontsize=12,c='black')
        plt.show()
        print(nn/nnn)

    #plt.savefig(da
  #%% Agreement between sign of change 
ka=6
hadcm = modelData['hadcm']['ANN']['pre']
hadcm0k = hadcm.groupby_bins('time',[ka*1000-500,ka*1000+500]).mean(dim='time')
hadcm6k = hadcm.groupby_bins('time',[(ka+6)*1000-500,(ka+6)*1000+500]).mean(dim='time')
hadcm  = hadcm6k.data-hadcm0k.data

trace = modelData['trace']['ANN']['pre']
trace0k = trace.groupby_bins('time',[ka*1000-500,ka*1000+500]).mean(dim='time')
trace6k = trace.groupby_bins('time',[(ka+6)*1000-500,(ka+6)*1000+500]).mean(dim='time')
trace   = trace6k.data-trace0k.data

hadcm.resize(np.shape(trace))

agree = (np.sign(trace) + np.sign(hadcm))

plt.style.use('ggplot')
plt.figure(figsize=(20,10)); plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 1; plt.rcParams['axes.edgecolor'] = 'k'
ax1 = plt.subplot(projection=ccrs.Robinson()) 
ax1.spines['geo'].set_edgecolor('black')
ax1.set_global()
ax1.coastlines()
ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
mlevels = np.array([-0.6,0.6])
model_contour = plt.pcolormesh(trace0k.lon.data, trace0k.lat.data, 
                    agree[0,:,:],transform=ccrs.PlateCarree(),
                    cmap='BrBG',alpha=0.8)

plt.title('ANN pre '+str(ka)+'ka-0.5ka',fontsize=20)
#plt.title(name,fontsize=30)


#%%

