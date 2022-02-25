#%% 1 Load Packages
import numpy      as np
#import pandas     as pd
#import regionmask as rm
import xarray     as xr
#from scipy        import stats
from scipy.stats  import pearsonr
#For Figures:
import matplotlib.pyplot as plt         # Packages for making figures
import matplotlib.gridspec as gridspec
#from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs              # Packages for mapping in python
import cartopy.util as cutil
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/HoloceneHydroclimate/HoloceneHydroclimate/'
#%% 
def calcCorrelation (inData,model1,model2,var1,var2,season,t0=0,t1=12,mask=False):
    if model1 == model2:
        modelData1 = inData[model1][season][var1]
        modelData2 = inData[model2][season][var2]
        lats,lons = modelData1.lat,modelData1.lon
    else: 
        modelData1 = inData[model1][season][var1+'_regrid']
        modelData2 = inData[model2][season][var2+'_regrid']
        lats,lons = modelData1.lat_regrid,modelData1.lon_regrid
    t_0 = np.argmin(np.abs(t0*1000 - modelData1.time.data))   
    t_1 = np.argmin(np.abs(t1*1000 - modelData1.time.data))
    out = np.full(np.shape(modelData1)[1:3],np.NaN)
    for lat in range(len(lats)):
        for lon in range(len(lons)):
            out[lat,lon] = pearsonr(modelData1[t_0:t_1,lat,lon].data,modelData2[t_0:t_1,lat,lon].data)[0]

    return(out)

#%% 
#Settings
times = [0,12]
szn = 'ANN'
v_1 = 'p-e'
v_2 = 'tas'
m_1 = 'hadcm'
m_2 = 'hadcm'
#import seaborn as sns
#sns.palplot(sns.color_palette("BrBG", 7))
from matplotlib.colors import LinearSegmentedColormap
cramp = plt.cm.get_cmap('RdGy_r',10)
cramp = LinearSegmentedColormap.from_list('cramp',['Black','white','#4A148C'],N=10)

#Load data
modelData = {}
for model in ['hadcm','trace']:
    modelData[model] = {}
    if m_1 == m_2:
        modelData[model][szn]   = xr.open_dataset(dataDir+'Data/Model/'+model+'_'+szn+'.nc',decode_times=False)
    else: modelData[model][szn] = xr.open_dataset(dataDir+'Data/Model/'+model+'_'+szn+'_regrid.nc',decode_times=False)

#Calculate using function
rVals = calcCorrelation(modelData,m_1,m_2,v_1,v_2,szn,times[0],times[1])

#Plot 
mlevels = np.array([i /10 for i in list(range(-10,11,2))])
if m_1 == m_2: 
    lats,lons = modelData[m_1][szn]['lat'],modelData[m_1][szn]['lon']
    name = m_1+'_'+v_1+'_'+v_2+'_'+szn
else: 
    lats,lons = modelData[m_1][szn]['lat_regrid'],modelData[m_1][szn]['lon_regrid']
    name = 'multi'+'_'+v_1+'_'+v_2+'_'+szn
gs = gridspec.GridSpec(12,8)
ax = plt.subplot(gs[0:10,0:8],projection=ccrs.Robinson()) 
plt.title('Correlation between '+m_1+'('+szn+'_'+v_1+') & '+m_2+'('+szn+'_'+v_2+')')
ax.set_global()
ax.coastlines()
data_cyclic,lon_cyclic = cutil.add_cyclic_point(rVals,coord=lons)
model_contour = plt.contourf(lon_cyclic,lats, data_cyclic,transform=ccrs.PlateCarree(),
                             levels=mlevels,cmap=cramp) 
cm = plt.cm.get_cmap(cramp,10)
ax1 = plt.subplot(gs[10:12,1:7]) 
n, bins, patches = ax1.hist(rVals.flatten(), range = [-1,1],
                           bins=10, alpha=1, density=True,
                           label='All Gridcells',edgecolor='k',lw= 1)
for i, p in enumerate(patches): plt.setp(p, 'facecolor', cm(i/10)) 
ax1.set_xlim([-1,1])

ax1.set_xticks(mlevels)
ax1.axvline(x=np.mean(rVals),c='k')
ax1.spines['right'].set_visible(True);ax1.spines['top'].set_visible(True)
#plt.colorbar(model_contour,cax=inset_axes(ax,width='70%',height="4%",loc=8),orientation="horizontal")

if save: plt.savefig(dataDir+'Figures/Model/TransientCorrelations/'+name+'.png',
                     dpi=400,format='png',bbox_inches='tight')       
plt.show()