import geopandas as gpd
from   geopandas.tools  import sjoin
from   shapely.geometry import Point
import geoplot
import matplotlib as mpl
import matplotlib.pyplot as plt         # Packages for making figures
import cartopy.crs as ccrs              # Packages for mapping in python
import cartopy.feature as cfeature
import pandas as pd
import numpy  as np
#
gitHubDir  = '/Volumes/GoogleDrive/My Drive/zResearch/HoloceneHydroclimate/'
save = True # for figures

# Outputs
# csv of proxy table with added region columns based on spatial join 
# map of proxy site within IPCC regions 
# maps of regions colored by proxy trends

#%%
#
#Set variables used in region ID 
climVar='T'   #HC or T
threshold = 3  #Minimum number of records for a region to be included in analysis
Oceans = False #Option to only get Land or Land-Ocean regions
#
#Load proxy data as a point geoDataframe
proxyData    = pd.read_csv(gitHubDir+'DataFiles/proxy'+climVar+'.csv')
proxyData['geometry'] = proxyData.apply(lambda x: Point((float(x.Lon), float(x.Lat))), axis=1)
proxyData    = gpd.GeoDataFrame(proxyData, geometry='geometry')
proxyData    = proxyData.set_crs(epsg=4326)   #WGS84
proxyDataPrj = proxyData.to_crs("EPSG:32662") #Project to  #WGS_1984_Plate_Carre
#
#Load ippc regions shapefile as a polygon geoDataframe
ipccRegions = gpd.read_file(gitHubDir+'DataFiles/IPCC-WGI-reference-regions-v4_shapefile/IPCC-WGI-reference-regions-v4.shp')
if Oceans == False:
    ipccRegions = ipccRegions.loc[(ipccRegions['Type'] != 'Ocean')] 
ipccRegionsPrj = ipccRegions.to_crs("EPSG:32662") #Project to  #WGS_1984_Plate_Carre
#
#Spatial join to ID which polygon each point is within
proxyDataByRegionPrj = sjoin(proxyDataPrj, ipccRegionsPrj, how="left", op="within")
proxyDataByRegion    = proxyDataByRegionPrj.to_crs("EPSG:4326") #for mapping
if save: proxyDataByRegion.to_csv(gitHubDir+'DataFiles/proxy'+climVar+'_regions.csv')
#
#Map regions and proxy data 
crs = ccrs.PlateCarree()
fig, ax = plt.subplots(subplot_kw={'projection': crs})
plt.style.use('ggplot')
plt.figure(figsize=(8,15))
#Figure Settings and basemap items
plt.rcParams['axes.facecolor'] ='white'
plt.rcParams['axes.linewidth'] = 0.2
ax.spines['geo'].set_edgecolor('black')
ax.set_global(); 
ax.coastlines()
ax.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k',lw=0.2)
#Plot proxy sites and ippc regions
ax.add_geometries(ipccRegions['geometry'], crs=crs, 
                  facecolor='none', edgecolor='orange')
proxy_scatter = ax.scatter(proxyDataByRegion.Lon,proxyDataByRegion.Lat,
                  c=proxyDataByRegion['bin6ka'],
                  marker='o',s=5,edgecolor='k',lw=0.2,alpha=0.8,
                  cmap='BrBG',vmin=-3,vmax=3)
plt.show()

#%%
#
#
varList = ['EarlySlope','LateSlope','bin6ka']
proxyDataByRegion=proxyDataByRegion
regionsProxyData = {'Name':[],'Count':[],'geometry':[]}
for var in varList: regionsProxyData[var] = []  
             
for name in (proxyDataByRegion['Name'].unique()):
    if type(name) != str: continue
    dataSubset = proxyDataByRegion.loc[(proxyDataByRegion['Name'] == name)]
    regionsProxyData['Name'].append(name)
    regionsProxyData['Count'].append(np.shape(dataSubset)[0])
    for var in varList: regionsProxyData[var].append(np.nanmedian(dataSubset[var]))
    regionsProxyData['geometry'].append(ipccRegions['geometry'][[index for index, item in enumerate(list(ipccRegions['Name'])) if item == name][0]])
    
regionsProxyData = pd.DataFrame.from_dict(regionsProxyData)
regionsProxyData = regionsProxyData.loc[(regionsProxyData['Count'] >= threshold)]
regionsProxyData = regionsProxyData.loc[(regionsProxyData['Name'] != 'N.Atlantic-Ocean')]

#Save file to folder
if save: regionsProxyData.to_csv(gitHubDir+'DataFiles/proxy'+climVar+'_IPCCregions.csv')

#%%
#Plot values for each region to investigate regional variations
#
regionsProxyData = gpd.GeoDataFrame(regionsProxyData, geometry='geometry')
norm=0.5
#Plot climate variables 
for var in varList:
    regionsMap = regionsProxyData.loc[regionsProxyData[var].notna()] 
    geoplot.choropleth(regionsMap,figsize=(8, 4),edgecolor='k',
                       hue=regionsMap[var],cmap='BrBG',
                       norm=mpl.colors.Normalize(norm*-1,norm),
                       legend=True)
    plt.title(var+" median value for regions with at least "+str(threshold)+" records")
    if save: plt.savefig(gitHubDir+'Figures/Test_IPCCregionMedianValues/'+var+'.png', dpi=400,format='png')
    else: plt.show()

#Plot Number of records in each regions
var='Count'
regionsMap = regionsProxyData.loc[regionsProxyData[var].notna()] 
geoplot.choropleth(regionsMap,figsize=(8,4),edgecolor='k',
                       hue=regionsMap[var],cmap='Greens',
                       norm=mpl.colors.Normalize(0,50),
                       legend=True)
plt.title(var+" of records within each region")
if save: plt.savefig(gitHubDir+'Figures/Test_IPCCregionMedianValues/'+var+'.png', dpi=400,format='png')
else: plt.show()
