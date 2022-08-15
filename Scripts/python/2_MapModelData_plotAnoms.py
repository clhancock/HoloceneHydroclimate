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
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 0.5; 
plt.rcParams['axes.edgecolor'] = 'k'
#plt.tick_params(labelsize=8)

#%% 2 Load Data
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/2021_HoloceneHydroclimate/2021_HoloceneHydroclimate/'
modelData = {}
for model in ['hadcm','trace','cmip6']:
    modelData[model] = {}
    for szn in ['ANN','JJA','DJF']:
        modelData[model][szn] = xr.open_dataset(dataDir+'Data/Model/'+model+'/'+model+'_'+szn+'.nc',decode_times=False)


#%%Land Mask
data = modelData['cmip6']['ANN']
land = rm.defined_regions.natural_earth.land_110.mask_3D(data.lon,data.lat)
land = land.squeeze('region').data
pd.DataFrame(land).to_csv(dataDir+'Data/Model/cmip6LandMask.csv')

#%% 3 MH anoms
#Settings
for var in ['pre','p-e','tas']:
    if var == 'tas': 
        cramp, units = 'RdBu_r','degC'
        mlevels = np.array([i /100 for i in list(range(-250,251,20))])
        ticklabels = [-2.0,-1.0,0,1.0,2.0]
    else: 
        cramp, units = 'BrBG', 'mm/day'
        mlevels = np.array([i /100 for i in list(range(-75,76,6))])
        ticklabels = [-0.6,-0.3,0,0.3,0.6]
    labels=['(a)','(b)','(c)']
    label=0
    seasons = ['ANN','JJA','DJF'] 
    s = 3
    #Plot Figure
    save = True
    plt.style.use('ggplot')
    plt.figure(figsize=(3.25,6),dpi=600); 
    h=len(seasons)
    gs = gridspec.GridSpec(h*s+1,s)
    for szn in seasons:
        data = modelData['cmip6'][szn][var]
        #data = data.rename({'lat_regrid':'lat','lon_regrid':'lon'})
        data = np.mean(data,axis=0)
        i = [x for x, val in enumerate(seasons) if val == szn][0]
        ax = plt.subplot(gs[(i)*s:(i)*s+s,0:s],projection=ccrs.Robinson()) 
        ax.spines['geo'].set_edgecolor('black')
        ax.set_global()
        ax.coastlines()
        ax.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
        data_cyclic,lon_cyclic = cutil.add_cyclic_point(data,coord=data.lon.data)
        model_contour = plt.contourf(lon_cyclic, data.lat.data, 
                                     data_cyclic,transform=ccrs.PlateCarree(),
                                     levels=mlevels,extend='both',cmap=cramp) 
        if szn == 'ANN':
            ax.annotate('CMIP Ensemble Mean (Mid-Holocene - Pre-Industrial)  ',
                        xy=(0, 0), xycoords='data',fontsize=8, xytext=(0.5, 1.1), 
                        textcoords='axes fraction',horizontalalignment='center', verticalalignment='center')  
            ax.annotate('Annual', xy=(0, 0), xycoords='data', xytext=(-0.05, 0.5), 
                    textcoords='axes fraction', fontsize=8, fontfamily = 'Arial',
                    rotation=90, horizontalalignment='center', verticalalignment='center') 
        else: 
            ax.annotate(szn, xy=(0, 0), xycoords='data', xytext=(-0.05, 0.5), 
                    textcoords='axes fraction', fontsize=8, fontfamily = 'Arial',
                    rotation=90, horizontalalignment='center', verticalalignment='center') 
        ax.annotate(labels[label], 
                        xy=(0, 0), xycoords='data', xytext=(0.05, 0.95), 
                        textcoords='axes fraction', fontsize=8, fontfamily = 'Arial')
        label+=1
    ax = plt.subplot(gs[(h*s):(h*s+1),0:s])
    ax.axis('off')
    cbar = plt.colorbar(model_contour,orientation="horizontal",
                        cax=inset_axes(ax,width='90%',height="30%",loc="upper center"),
                ticks=ticklabels).set_label(var.upper()+' ('+units+')',
                                                       fontsize=8,c='black')
    plt.tick_params(labelsize=8)
    #Save or show
    if save: plt.savefig(dataDir+'Figures/Model/Anomalies/CMIP_MH-PI_bySeason_'+var+'.png',
                         dpi=600,format='png',bbox_inches='tight')       
    else: plt.show()


#%% 4 transient anoms
ka = [0.5,6,12]
models  = ['trace','hadcm']#,'cmip6']
seasons = ['ANN'] 
modelAnom = {}
s = 3
for var in ['pre','p-e','tas']:
    if var == 'tas': 
        cramp, units = 'RdBu_r','degC'
        mlevels = np.array([i /100 for i in list(range(-250,251,20))])
        ticklabels = [-2.0,-1.0,0,1.0,2.0]
    else: 
        cramp, units = 'BrBG', 'mm/day'
        mlevels = np.array([i /100 for i in list(range(-75,76,6))])
        ticklabels = [-0.6,-0.3,0,0.3,0.6]
    labels=['(a)','(b)','(c)','(d)']
    label=0
    #Plot Figure
    save = True
    plt.style.use('ggplot')
    plt.figure(figsize=(3.25*len(models),4),dpi=400); 
    h=len(seasons)*(len(ka)-1)
    gs = gridspec.GridSpec(h*s+1,len(models)*s)
    for t in range(0,len(ka)-1):
        for model in models:
            for szn in seasons:
                data = modelData[model][szn][var]
                data0 = data.groupby_bins('age',[ka[t]*1000-100,ka[t]*1000+100]).mean(dim='age')
                data1 = data.groupby_bins('age',[ka[t+1]*1000-100,ka[t+1]*1000+100]).mean(dim='age')
                data = data0
                data.data = (data1.data-data0.data)
                data = data.squeeze('age_bins')
                i = [x for x, val in enumerate(seasons) if val == szn][0]
                j = [x for x, val in enumerate(models) if val == model][0]
                ax = plt.subplot(gs[(i+t)*s:(i+t)*s+s,j*s:j*s+s],projection=ccrs.Robinson()) 
                ax.spines['geo'].set_edgecolor('black')
                ax.set_global()
                ax.coastlines()
                ax.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
                ax.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
                data_cyclic,lon_cyclic = cutil.add_cyclic_point(data,coord=data.lon.data)
                model_contour = plt.contourf(lon_cyclic, data.lat.data, 
                                             data_cyclic,transform=ccrs.PlateCarree(),
                                             levels=mlevels,extend='both',cmap=cramp)
                if j == 0:
                    ax.annotate(str(ka[t+1])+'ka - '+str(ka[t])+'ka', 
                        xy=(0, 0), xycoords='data', xytext=(-0.05, 0.5), 
                        textcoords='axes fraction', fontsize=8, fontfamily = 'Arial',
                        rotation=90, horizontalalignment='center', verticalalignment='center')  
                if t == 0:
                    if   model == 'hadcm': modelName = 'HadCM'
                    elif model == 'trace': modelName = 'TraCE'
                    ax.annotate(modelName, xy=(0, 0), xycoords='data',fontsize=8,
                        xytext=(0.5, 1.1), textcoords='axes fraction',
                        horizontalalignment='center', verticalalignment='center')  
                ax.annotate(labels[label], 
                        xy=(0, 0), xycoords='data', xytext=(0.05, 0.95), 
                        textcoords='axes fraction', fontsize=8, fontfamily = 'Arial')
                label+=1
    ax = plt.subplot(gs[(h*s):(h*s+1),0:len(models)*s])
    ax.axis('off')
    #cbar = plt.colorbar(model_contour,cax=inset_axes(ax,width='80%',height="30%",loc="upper center"),
     #           orientation="horizontal",ticks=[-0.55,-0.25,0,0.25,0.55]).set_label(var+' Difference ('+units+')',fontsize=8,c='black')
    cbar = plt.colorbar(model_contour,cax=inset_axes(ax,width='90%',height="30%",loc="upper center"),
                orientation="horizontal",ticks=ticklabels).set_label('Annual '+var.upper()+' ('+units+')',fontsize=8,c='black')
    plt.tick_params(labelsize=8)
    #Save or show
    if save: plt.savefig(dataDir+'Figures/Model/Anomalies/Trans_Anoms_byAge_byModel_'+var+'.png',
                         dpi=600,format='png',bbox_inches='tight')       
    plt.show()
    
    
    


 #%% 5 Transient Trends
#Settings
save = True
for var in ['pre','p-e','tas']:
    if var == 'tas': 
        cramp = 'RdBu_r'
        units = 'degC'
        mlevels = np.array([i /100 for i in list(range(-110,111,20))])
    else: 
        cramp = 'BrBG'
        units = 'mm/day'
        mlevels = np.array([i /100 for i in list(range(-11,12,1))])
    ka = [0,6,12]
    models  = ['trace','hadcm']
    szn = 'ANN'
    modelTrends = {}
    s = 3
    
    #Calculate model trends
    for model in ['hadcm','trace']:
        for t in range(len(ka)-1):
            data = modelData[model][szn][var]
            i = [x for x, val in enumerate(data['age']) if 
                 val >= 1000*ka[t] and val <= 1000*ka[t+1]] 
            data = data[i,:,:]
            dataSl = np.full(np.shape(data)[1:3],np.NaN)
            for lat in range(len(data.lat)):
                for lon in range(len(data.lon)):
                    modelReg = stats.linregress(data.age.data,data.data[:,lat,lon])
                    dataSl[lat,lon] = modelReg[0]*-1000
            modelTrends[model+'_'+str(ka[t+1])+'to'+str(ka[t])] = dataSl
    #Plot FIugres
    plt.style.use('ggplot')
    plt.figure(figsize=(6.5,5),dpi=400); 
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['axes.facecolor'] ='white'
    plt.rcParams['axes.linewidth'] = 0.5; 
    plt.rcParams['axes.edgecolor'] = 'k'
    gs = gridspec.GridSpec((len(ka)-1)*s+1,len(models)*s)
    for model in models:
        lats = modelData[model][szn][var].lat.data
        lons = modelData[model][szn][var].lon.data
        for t in range(len(ka)-1):
            print(model+'_'+str(ka[t+1])+'to'+str(ka[t]))
            data = modelTrends[model+'_'+str(ka[t+1])+'to'+str(ka[t])] 
            i = [x for x, val in enumerate(ka) if val == ka[t]][0]
            j = [x for x, val in enumerate(models) if val == model][0]
            ax1 = plt.subplot(gs[i*s:i*s+s,j*s:j*s+s],projection=ccrs.Robinson()) 
            ax1.spines['geo'].set_edgecolor('black')
            ax1.set_global()
            ax1.coastlines()
            ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
            ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
            data_cyclic,lon_cyclic = cutil.add_cyclic_point(data,coord=lons)
            model_contour = plt.contourf(lon_cyclic, lats, 
                                         data_cyclic,transform=ccrs.PlateCarree(),
                                         levels=mlevels,extend='both',cmap=cramp)
            plt.title(model+'_'+str(ka[t+1])+'to'+str(ka[t]),fontsize=8)
    ax1 = plt.subplot(gs[6:7,0:len(models)*s])
    ax1.axis('off')
    plt.colorbar(model_contour,cax=inset_axes(ax1,width='90%',height="20%",loc="center"),
                orientation="horizontal").set_label(units+'/ka',fontsize=8,c='black')
    plt.tick_params(labelsize=8)
    #plt.set_ticklabels(mlevels)
    plt.title(var+' '+szn+' Holocene Trends',fontsize=10)
    if save: plt.savefig(dataDir+'Figures/Model/TransientTrends/HoloceneTrends_'+var+'_'+szn+'.png',
                         dpi=400,format='png',bbox_inches='tight')       
    else: plt.show()


