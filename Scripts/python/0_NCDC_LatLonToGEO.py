#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ### Takes Longitude and Latitude, Finds NASA GCMD Location Keyword
# #### Files needed:  
# 1) "lon-lat.txt": two columns [lon, lat] with header  
# 2) "ne_10m_admin_0_map_units" country shapefile from https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-details/  
# 3) "tl_2018_us_state" US states shapefile from https://www2.census.gov/geo/tiger/TIGER2018/STATE/  
# 4) "gpr_000a11a_e" Canadian province shapefile from https://open.canada.ca/data/en/dataset/35ee219c-a3b0-448b-a952-3e195cb40b70  
# 5) "World_Seas_IHO_v1" shapefile from http://www.marineregions.org/downloads.php  
# 6) "GCMD-country-lookup.txt" provides shapefile IDs and corresponding GCMD location keywords for country  
# 7) "GCMD-state-lookup.txt" provides shapefile IDs and corresponding GCMD location keywords for US state  
# 8) "GCMD-province-lookup.txt" provides shapefile IDs and corresponding GCMD location keywords for Canadian province  
# 9) "GCMD-ocean-lookup.txt" provides shapefile IDs and corresponding GCMD location keywords for ocean basin  
# 
# Notes: If site lat/lon cannot be placed on shapefile, a "Not found" warning will be returned. These sites will need to be handled manually. In addition, the shapefile does not demarcate large lakes (e.g., Tanganyika, Malawi), which do have GCMD location keywords. These I also did by hand.



import pandas as pd
import fiona
import shapely
from shapely import geometry
import numpy as np



# READ IN SITE LONS AND LATS
#sites = pd.read_csv("lon-lat.txt", delim_whitespace=True)
sites = pd.read_csv("/Users/chrishancock/Desktop/latlon.csv",sep=',')
#sites = location_data[['lon','lat']]
sites.head(5)


#data_dir = '/projects/pd_lab/data/modern_datasets/location_files/'
data_dir = '/Users/chrishancock/Desktop/location_files/'

# READ IN GCMD LOCATION KEYWORD LOOK-UP TABLES
GCMD_state = {}
with open(data_dir+"GCMD-state-lookup.txt") as f:
    for line in f:
       (key, val) = line.strip().split(None,1)
       GCMD_state[key] = val

GCMD_province = {}
with open(data_dir+"GCMD-province-lookup.txt") as f:
    for line in f:
       (key, val) = line.strip().split(None,1)
       GCMD_province[key] = val

GCMD_country = {}
with open(data_dir+"GCMD-country-lookup.txt") as f:
    for line in f:
       (key, val) = line.strip().split(None,1)
       GCMD_country[key] = val

GCMD_ocean = {}
with open(data_dir+"GCMD-ocean-lookup.txt") as f:
    for line in f:
       (key, val) = line.strip().split(None,1)
       GCMD_ocean[key] = val



location_txt = []

# GEOLOCATE AND OUTPUT NASA GCMD LOCATION KEYWORD
for index, row in sites.iterrows():  # loop through sites
    
    found = False
    point = shapely.geometry.Point(row['lon'],row['lat'])  # grab site lon and lat
    
    with fiona.open(data_dir+"ne_10m_admin_0_map_units/ne_10m_admin_0_map_units.shp") as countries:  # open country shapefile
        for country in countries:  # loop through countries
            shape = shapely.geometry.asShape(country['geometry'])
            if shape.contains(point):
                if (country['properties']['NAME']=='United States of America'):  # If US, find state and output
                    with fiona.open(data_dir+"tl_2018_us_state/tl_2018_us_state.shp") as states:  # open state shapefile
                        for state in states:   # loop through states
                            shape = shapely.geometry.asShape(state['geometry'])
                            if shape.contains(point):
                                print(row['dataSetName'],GCMD_state.get(state['properties']['GEOID']))
                                location_txt.append(GCMD_state.get(state['properties']['GEOID']))
                                found = True
                                break
                elif (country['properties']['NAME']=='Canada'):   # If Canada, find province and output
                    with fiona.open(data_dir+"gpr_000a11a_e/gpr_000a11a_e.shp") as provinces:  # open province shapefile
                        for province in provinces:  # loop through provinces
                            shape = shapely.geometry.asShape(province['geometry'])
                            if shape.contains(point):
                                print(row['dataSetName'],GCMD_province.get(province['properties']['PRUID']))
                                location_txt.append(GCMD_province.get(province['properties']['PRUID']))
                                found = True
                                break
                else:   # not US or Canada, output country name
                    print(row['dataSetName'],GCMD_country.get(country['id']))
                    location_txt.append(GCMD_country.get(country['id']))
                    found = True
                    break
                break
                
    if (found==False):   # point not found on land, try oceans
        with fiona.open(data_dir+"World_Seas_IHO_v1/World_Seas.shp") as basins:  # open ocean basin shapefile
            for basin in basins:   # loop through ocean basins
                shape = shapely.geometry.asShape(basin['geometry'])
                if shape.contains(point):
                    print(row['dataSetName'],GCMD_ocean.get(basin['id']))
                    location_txt.append(GCMD_ocean.get(basin['id']))
                    found = True
                    break
                    
    if (found==False):    # point not found on land or ocean due to some vagary of shapefile
        location_txt.append('')
        print(row['dataSetName'],"Not found")


### OUTPUT

# Write the output to a new file.
output_array = np.concatenate((np.array(sites['dataSetName'])[:,None],np.array(location_txt)[:,None]),axis=1)
np.savetxt(data_dir+'location_output.txt',output_array,delimiter=',',fmt='%s')

