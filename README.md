# HoloceneHydroclimate
Repository for data, code, and figures used by Hancock et al. (2023).

The proxy records in the most recent Holocene Hydroclimate dataset are provided here:  
https://lipdverse.org/HoloceneHydroclimate/current_version/

A table listing the proxy records is provided here: 
https://raw.githack.com/clhancock/HoloceneHydroclimate/main/Figures/Proxy/TableS1/TableS1.html 
and described by: 
https://raw.githack.com/clhancock/HoloceneHydroclimate/main/Figures/Proxy/TableS1/TableS1_Key.pdf

A complete list of references is provided here:
https://lipdverse.org/HoloceneHydroclimate/current_version/references.html

Cite: Hancock, C. L., McKay, N. P., Erb, M. P., Kaufman, D. S., Routson, C. R., Ivanovic, R. F., et al. (2023). Global synthesis of regional Holocene hydroclimate variability using proxy and model data. Paleoceanography and Paleoclimatology, 38, e2022PA004597. https://doi.org/10.1029/2022PA004597

[ArticleCoverPage.pdf](https://github.com/clhancock/HoloceneHydroclimate/files/12729999/Cover.pdf)


## Data ================================

### Models
```
Trace/Hadcm:

-NetCDF files (Annual/JJA/DJF) binned to 100-year resolution.
-Each dataset includes temperature, precipitation, evaporation, and P-E. 
-Examples of the original file names that these were created from are "trace.01-36.22000BP.cam2.PRECT.22000BP_decavgANN_400BCE.nc" and "deglh.vn1_0.precip_mm_srf.monthly.ANN.010yr_s.nc" respectively. 
-Data for each model are provided in the native resolution and a regridded spatial resolution. 
-The original files were modified to standardize units/names between different models according to the /Notebooks/1_DataPrep/transientFormatting.py script. 
(Dimensions:  (lon: 135, lat: 90, age: 121))

These results are shown in Fig. 5
```

```
CMIP6:

-NetCDF files of mid-Holocene minus preindustrial anomalies. 
-Each dataset includes temperature, precipitation, evaporation, and P-E
-Data for each model are provided in the native resolution and a regridded spatial resolution. 
-The original files were modified to standardize units/names between different models according to the /Notebooks/1_DataPrep/ cmip6Formatting.py script.
(Dimensions:  (lon: 135, lat: 90, model: 12))

These results are shown in Fig. 4
```

```
RegionalTS:

-csv files containing the mean of data binned by IPCC region for pre, tas, and p-e.  
-For hadcm/trace, each region is a column and each row an age.
-For cmip6, each row is a model
-Land/all files distinguish if the regional mean includes an ocean mask. 
-Data calculated using /Notebooks/1_DataPrep/CalculateModelValuesByRegion.ipynb.

These results are shown in Fig. 6
```

```
pseudoProxyCorr_pre.csv & pseudoProxyCorr_pre_byCount.csv:

-csv files of correlations between pseudoproxies and regional mean (RegionalTS data)
-Data calculated using /Notebooks/2_Analysis/pseudoProxyCorrelations.py

These results are shown in the SM
```



### Proxies

```
lipdData.rds

-rds file of LiPD extracted to ts objects (list of proxy records and their data) 
-Contains hydroclimate and temperature files
-Includes additional standardization / metadata calculation not included in the lipdverse files (such as IPCC region of each site)
-Data downloaded from https://lipdverse.org/ and modified by /Notebooks/1_DataPrep/1_StandardizeLiPDs.Rmd.

These results are shown in Fig. 1
```

```
Proxy_MetaData

-csv files containing key metadata which are used for plotting
-Different files for hydroclimate and temp12k
-Unlike lipdData.rds, this file can be opened in R and python. 
-Created by /Notebooks/1_DataPrep/1_StandardizeLiPDs.Rmd.
```


### Regional Composites
```
-csv files for hydroclimate and temperature proxy composites for each region
-Hydroclimate composites are standardized anomalies (z-scores). Temperatures are degC anomalies. Both relative to the Holocene mean. 
-Each region has a csv file which includes the entire composite ensemble
-'MedianTS_byRegion.csv' list the ensemble median for each region in a single file
-Results of /Notebooks/2_Analysis/Composite.R
-Results of /Notebooks/2_Analysis/CompositeCorrelations.R (zip file)

These results are shown in Fig. 2,3,6,7
```


## Notebooks ================================

```
1_DataPrep

1_StandardizeLiPDs.Rmd:      for loading/standardizing data LiPD data (proxies)
Convert2calAge.R:            for converting radiocarbon years to calendar years
cmip6Formatting.py:          for standardizing CMIP6 data
transientFormatting.py:      for standardizing TraCE/HadCM data
```

```
2_Analysis

Composite.R:                 for compositing the data
CompositeCorrelations.R:     for calculating the correlation between HC and T composite ensembles
pseudoProxyCorrelations.py:  for testing the correlation between pseudoproxies and the regional mean.
```

```
3_Figures

For creating figures in Hancock et al. (2023). 
Each script named to indicate which figure it was used to create (Fig1_...)
```


## Figures ================================

```
Figure used in publication are identified in '/Notebooks/3_Figures/README.md'
Figures created by /Notebooks/3_Figures/...
```

