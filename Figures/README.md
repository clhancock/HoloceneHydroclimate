# HoloceneHydroclimate Figures
Repository for data, code, and figures used by Hancock et al. (2023).

The proxy records in the Holocene Hydroclimate dataset are provided here:  
https://lipdverse.org/HoloceneHydroclimate/current_version/

A table listing the proxy records is provided here: 
https://raw.githack.com/clhancock/HoloceneHydroclimate/main/Figures/Proxy/TableS1/TableS1.html
and described by:
https://raw.githack.com/clhancock/HoloceneHydroclimate/main/Figures/Proxy/TableS1/TableS1_Key.pdf

Hancock, C., McKay, N. P., Erb, M. P., Kaufman, D. K., Routson, C., Ivanovic, R. F., Gregoire, L. J., and Valdes, P.: Global synthesis of regional Holocene hydroclimate variability using proxy and model data, Paleoceanography and Paleoclimatology


## Figure 1 - Proxy Map
Geographic distribution of proxy data (top). Sites symbolized by proxy category (Table 1). IPCC reference regions (Iturbide et al., 2020) shown with solid black lines. Temporal coverage (bottom) shows the number of available records through the Holocene. Values are based on the range between the minimum and maximum age of each record. Colors symbolize proxy categories and match the map legend.
```
/Figures/Proxy/Figure1_HC_ProxyMapwithTime.png
    created by /Notebooks/3_Figures/Fig1_ProxyDataSummary.Rmd
```

## Figure 2 - Composite Grid
Proxy composite timeseries. Y-axes (right side) show z-score anomalies. Gray shading signifies portions of the timeseries when fewer than 50% of proxy records contribute to the timeseries. The dark and light blue and gray shaded envelopes indicate the 50% and 95% quantiles of the ensemble range, respectively. Colors and symbols for the maps and temporal coverage plots match Figure 1. Bars under the maps represent the proportion of annual (black), summer (red), and winter (blue) records. The regions are arranged according to geographic proximity and divided between tropical (a-p) and extratropical (q-af) locations. Region domains and names are according to Iturbide et al. (2021).
```
/Figures/RegionComposites/compositeSummaryGrid_HC_1.png
/Figures/RegionComposites/compositeSummaryGrid_HC_2.png
    created by /Notebooks/3_Figures/Fig2_CompositeGrid.Rmd
```

## Figure 3 - Monsoon Composite
Proxy record composites from monsoon regions. Summer insolation anomalies at 30° north and south are plotted for comparison (Berger & Loutre, 1991). Composites are median ensemble members plotted relative to their last millennium average (0-1 ka). Gray shading shows the 10-90% quantiles of the aggregated ensembles for each hemisphere.
```
/Figures/RegionComposites/Monsoon_HC.png 
    created by /Notebooks/3_Figures/Fig3_MonsoonRegions.Rmd
```

## Figure 4 - mid-Holocene CMIP6 Map
Mid-Holocene precipitation anomalies for the CMIP6 ensemble mean. Green colors signify greater precipitation during the mid-Holocene relative to preindustrial and brown colors indicate drier conditions. Values represent the mean during (a) annual, (b) boreal summer, and (c) boreal winter.
```
/Figures/Model/Anomalies/MH_Anoms_pre.png
    created by /Notebooks/3_Figures/Fig4&5_modelData.ipynb
```

## Figure 5 - Holocene Maps (TraCE/HadCM)
Annual precipitation anomalies for two time slices from two transient simulations. Maps show the mean difference between 1,000-year bins centered on (a,c) 12 and 6 ka (i.e., 12 ± 0.5 ka minus 6 ± 0.5 ka), and (b,d) 6 and 0.5 ka. Green and brown colors indicate wetter and drier conditions, respectively, during the earlier period. The 6-0.5 ka values (right column) are analogous to the CMIP6 values shown in Figure 4a. The same figure, but showing P-E anomalies, is in Figure S3 for comparison.
```
/Figures/Model/Anomalies/Trans_Anoms_pre.png
    created by /Notebooks/3_Figures/Fig4&5_modelData.ipynb
```

## Figure 6 - Composite Map
Comparison between proxy composite and simulated regional annual precipitation timeseries. TraCE (red) and HadCM (yellow) timeseries of mean annual precipitation (mm/year). Boxplot at 6 ka is based on the 12 model CMIP6 ensemble. Model timeseries are calculated as the latitude-weighted mean of all grid cells over land within each region. The variance of each proxy (blue) ensemble is scaled to match the average centennial-scale variance of the two transient timeseries. The y axes for each plot are adjusted to fit the full range of variability within that region over the Holocene. In all plots, each y-axis tick indicates a change of 50 mm/year. All data are shown as the difference relative to the last millennium average, which is shown by the black horizontal line, and positive values (up) indicate wetter conditions. Letters in the corner of each plot indicate the region as marked in Figure 2. Similar to Figure S5, which shows temperature. 
```
/Figures/RegionComposites/GlobalTSmap/global_pre_ANN_compBandPlt.png
    created by /Notebooks/3_Figures/Fig6_CompositeMap.Rmd
```

## Figure 7 - Holocene Correlation Map
Correlation between precipitation and temperature in Holocene climate simulations and proxy records. (a) Map colors show the mean correlation between the two climate variables from the HadCM and TraCE simulations. Circle colors indicate the correlation from the proxy composites. More saturated colors represent stronger correlations. (b) Zonal average of correlations from simulations (lines) and from regional proxy composites (dots) shown in panel (a). All correlations calculated using data binned to 100 years.
```
/Figures/Model/TransientCorrelations/multi_pre_tas_ANN_0to12ka.png
created by /Notebooks/Model/3_Figures/Fig7_correlations.ipynb
```

## Figure 8 - mid-Holocene Agreement Map
Agreement of the sign of mid-Holocene annual (a) precipitation and (b) surface air temperature anomalies relative to preindustrial. Model-to-model agreement is indicated by the color of the base map, and proxy-to-proxy agreement is indicated by colors of circles within regions, where more saturated colors on either end of the scale represent stronger agreement. Proxy–model agreement is indicated by the similarity between the map and circle colors. Models include TraCE, HadCM, and the 12 CMIP6 simulations, all regridded to a common spatial resolution (2.66° longitude by 2° latitude). Each color step from left to right on the color bar represents one additional model among the 14 with a positive (wetter or warmer) mid-Holocene anomaly. For the two transient models, signs represent the difference in mean values between 6.5-5.5 ka and 1-0 ka. The agreement among proxy records (colored circle within each region) is represented as the proportion of available proxy records with positive anomalies during the mid-Holocene (6.5-5.5 ka average) relative to the last millennium (1-0 ka average). Available records are defined as those with at least one data point within both the mid-Holocene and last millennium age ranges. Circles are only shown for regions with at least four available proxy records.
```
/Figures//Model/Fig8_midHoloceneAgreement_pre_tas_ANN.png
created by /Notebooks/3_Figures/Fig8&9_ProxyModelAgreement.ipynb
```

## Figure 9 - mid-Holocene Agreement Scatter
Scatterplot showing the extent of proxy–model agreement from Figure 8. Percent of proxy records within each region with a (a) wet or (b) warm mid-Holocene anomaly, each plotted against the percent of pseudoproxies with a corresponding anomaly. Perfect agreement is represented by a 1:1 line. Pseudoproxies represent model values calculated using the geographic locations and seasonality metadata of the proxy records. The regional percent is calculated as the mean of the percent of positive anomalies for each proxy site relative to the 14-model ensemble. (c,d) Same as (a) and (b), but showing the mean magnitude of the simulated anomaly for each proxy site. Regional abbreviations are in Figure S6 and Figure 2. 
```
/Figures/Model/Figure9_AgreeScatter.png.png
created by /Notebooks/Model/Fig8&9_ProxyModelAgreement.ipynb
```





