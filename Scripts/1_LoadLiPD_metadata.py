#Load packages
import pyleoclim as pyleo               # Packages for analyzing LiPD files 
import lipd                             # Packages for analyzing LiPD files 
import numpy  as np                     # Package with useful numerical functions
#import xarray as xr                     #Package for handling file types
import pandas as pd
#import math
from scipy import stats                 # Packages for calculations  
import pymannkendall as mk              # Package for trend detection
import matplotlib.pyplot as plt         # Packages for making figures
#import matplotlib as mpl 
import matplotlib.gridspec as gridspec
#import matplotlib.colors as pltcolors
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs              # Packages for mapping in python
import cartopy.feature as cfeature
#import cartopy.util as cutil
#
dataDir = '/Volumes/GoogleDrive/My Drive/zResearch/Data/'
gitHubDir = '/Volumes/GoogleDrive/My Drive/zResearch/HoloceneHydroclimate/'
save    = False # for figures
# Outputs
# csv with table of proxy Tsid, metadata, and caclualted values (slope/binvalues)
# map of proxy site symbolized based on catergory or season 
# maps of calculated variables
# pie chart of data soruces

#%%
#Read all available LiPD files & extrat to timeseries objects
#
LiPDdatabase   = pyleo.Lipd(dataDir+'LiPD/database')
LiPDdatabaseTS = LiPDdatabase.to_tso()
LiPDnew        = pyleo.Lipd(dataDir+'LiPD/NewHoloceneHydroclimate')
LiPDnewTS      = LiPDnew.to_tso()

#%%
# filter timeseries relevant to temp12k and hc12k
#
lipdHC,lipdT,dataSetList = [],[],['Pass Age Control','Dune.Finney.2012','Dunee.Finney.2012']

#Add files from newdata folder
for ts in LiPDnewTS: #These files don't have proper metadata
    dataSetList.append(ts['dataSetName'])
    try: #Not all ts have an interpretation to search for 
        if ts['paleoData_interpretation'][0]['variable'] in ['P','P-E','P_E','M']: 
            lipdHC.append(ts)
    except: pass
#Add files from general database folder without duplicates
for ts in LiPDdatabaseTS:
    if ts['dataSetName'] in dataSetList: continue #To prevent duplicate files
    try: #Not all ts have a compilation data to search for 
        for compilation in ts['paleoData_inCompilationBeta']:
            #Check if in temp12k
            if compilation['compilationName'] == "Temp12k":
                if any(v == '1_0_2' for v in compilation['compilationVersion']): 
                    lipdT.append(ts) ### add record
            #Check if in HC12k
            elif compilation['compilationName'] == "HoloceneHydroclimate":
                ts['paleoData_values'] = [np.nan if x == 'nan' else x for x in ts['paleoData_values']] #Convert Lake Deposit records to integers
                if ts['archiveType'] == 'LakeDeposits':
                    ts['paleoData_values'] = [x-2 for x in ts['paleoData_values']] #Change 1,2,3 values to -1 (high),0,1 (low)
                    if np.count_nonzero(~np.isnan(ts['paleoData_values'][0:11])) >= 9: #Set minimal resolution
                        if np.nanmin(ts['paleoData_values'][0:11]) != np.nanmax(ts['paleoData_values'][0:11]): #has some variability over Holocene
                            lipdHC.append(ts) ### add record
                else: lipdHC.append(ts) ### add record
    except: pass

print(' - Found '+str(len(lipdT))+' temperature proxies -')
print(' - Found '+str(len(lipdHC))+' hydroclimate proxies -')

#%%
fTS = lipd.filterTs(lipdHC,'paleoData_TSid == LPDfdd2a0239')
ts = pyleo.LipdSeries(fTS[0]).clean(True)
plt.plot(ts.time,ts.value,linewidth=1)
plt.show()

#%%
#Create summary table with metadata imformation and statistics for LiPD records 
#
ageMin=0; ageMax=12000; ageRes=1000; ageMinRef = 0; ageMaxRef = 1000
setTime = {'Total':{'min':ageMin,                  'max':ageMax},
           'Early':{'min':ageMin+(ageMax-ageMin)/2,'max':ageMax},
           'Late' :{'min':ageMin,                  'max':ageMin+(ageMax-ageMin)/2}}
#Function to assign proxy catergories based on archive/proxy/and unit type        
def assignProxyCategory(archive,proxy,unit):
    if archive == 'Speleothem': 
        if  proxy == 'd18O':        Category = 'Speleothem';    CategorySpec = Category+' ('+proxy+')'
        elif proxy == 'd13C':       Category = 'Speleothem';    CategorySpec = Category+' ('+proxy+')'
        else:                       Category = 'Speleothem';    CategorySpec = 'Speleothem (other)'
    elif archive == 'LakeDeposits': Category = 'Lake Deposits'; CategorySpec = 'Lake Deposits'
    elif archive == 'GlacierIce':   Category = 'Glacier Ice';   CategorySpec = 'Glacier Ice'
    elif proxy   == 'dDwax':        Category = 'Leaf Wax (dD)'; CategorySpec = 'Leaf Wax (dD)'
    elif archive == 'LakeSediment' and proxy == 'd18O': Category = 'Lake Sediment (d18O)'; CategorySpec = 'Lake Sediment (d18O)'
    elif proxy ==   'pollen': 
        Category = 'Pollen'
        if unit == 'mm' or unit == 'mm/a': CategorySpec = Category+' (calibrated)'
        else: CategorySpec = Category+' (uncalibrated)'
    else:
        Category = 'Other'
        if unit == 'mm' or unit == 'mm/a': CategorySpec = Category+' (calibrated)'
        else: CategorySpec = Category+' (uncalibrated)'
    return([Category,CategorySpec])

#Function to calculate trend of timeseries data
def calcTrend(timeValues,proxyValues,ageMin,ageMax,direction):    
    divisions = 2; minset = 3; sigLevel = 0.05; CI = 0.95; method = 'lin'
    if len(proxyValues) < 12: minset = 2
    ageRange = (ageMax - ageMin)/divisions
    valuedata = pd.DataFrame([timeValues,proxyValues]).transpose()
    valuedata.columns = ['time','values']
    valuedata = valuedata.loc[(valuedata['time'] <= ageMax) & (valuedata['time'] >= ageMin)]
    valuedata = valuedata.dropna() 
    check = []    
    for timeperiod in range(divisions):
        checkValue = len(valuedata.loc[(valuedata['time'] <= ageMin+(timeperiod+1)*ageRange) 
                                     & (valuedata['time'] >= ageMin+(timeperiod)*ageRange)]['time'])
        if checkValue >= minset: check.append(True)
        else: check.append(False)   
    if check.count(True) == divisions:
        valuedata['time']   = valuedata['time']*-1 #To get correct forward direction
        valuedata['values'] = valuedata['values']*direction
        trendDirection = mk.original_test(valuedata['values']*-1,alpha=sigLevel)
        if   method == 'lin': Slope = stats.linregress( valuedata['time'],valuedata['values'])
        elif method == 'sen': Slope = stats.theilslopes(valuedata['values'],valuedata['time'],CI)
        output = {'Direction':trendDirection[0],'Slope':Slope[0]*1000}
        if output['Direction'] == 'no trend': output['SlopeSig'] = 0
        else:                                 output['SlopeSig'] = output['Slope']
    else: output = {'Direction':np.NaN, 'Slope':np.NaN, 'SlopeSig':np.NaN}
    return(output)
#FCE4EC F8BBD0 F48GB1
#Function to create dataframe of relevent data / quick qc sheet to sort data by
def LiPDmetadata (lipd_list,standardize):
    problem_records = []; RecordNo=0
    data_matrix = {'Category':[],'CategorySpecific':[],
                   'ageMin':[],'ageMax':[],'ageRange':[],'ageRes':[]}
    names_metadata = {'TSid':'paleoData_TSid','Lat':'geo_meanLat','Lon':'geo_meanLon','Name':'dataSetName',
                      'Archive':'archiveType','Proxy':'paleoData_proxy','Unit':'paleoData_units'}
    names_interp = {'Interp':'variable','Season':'seasonalityGeneral','Direction':'direction'}
    for key in names_metadata: data_matrix[key] = []
    for key in names_interp:   data_matrix[key] = []
    for age in range(ageMin,ageMax,ageRes): data_matrix['bin'+str(int(age/ageRes))+'ka'] = []
    for time in setTime: 
        for value in ['Direction','Slope','SlopeSig']: data_matrix[time+value] = []
    for ts in lipd_list:
        print(RecordNo); RecordNo += 1
        try:
            lipdTS = pyleo.LipdSeries(ts).slice([ageMin,ageMax])
            data_matrix['ageMin'].append( np.nanmin(lipdTS.time))
            data_matrix['ageMax'].append( np.nanmax(lipdTS.time))
            data_matrix['ageRange'].append(data_matrix['ageMax'][-1]-data_matrix['ageMin'][-1])
            data_matrix['ageRes'].append(data_matrix['ageRange'][-1]/len(lipdTS.value))
            for key in names_metadata: 
                try:    data_matrix[key].append(ts[names_metadata[key]])
                except: data_matrix[key].append('')
            for key in names_interp:
                try:    data_matrix[key].append(ts['paleoData_interpretation'][0][names_interp[key]])
                except: data_matrix[key].append('')
            category = assignProxyCategory(data_matrix['Archive'][-1],data_matrix['Proxy'][-1],data_matrix['Unit'][-1])
            data_matrix['Category'].append(category[0])
            data_matrix['CategorySpecific'].append(category[1])
            if data_matrix['Direction'][-1] == 'negative': direction = -1
            else: direction = 1
            if standardize: lipdTS = lipdTS.standardize() 
            lipdTS_ref = np.nanmean(lipdTS.slice([ageMinRef,ageMaxRef]).value)
            for age in range(ageMin,ageMax,ageRes):
                value_ts = np.nanmean(lipdTS.slice([age,age+ageRes]).value)
                data_matrix['bin'+str(int(age/ageRes))+'ka'].append((value_ts-lipdTS_ref)*direction)
            for time in setTime:
                trendvalues = calcTrend(lipdTS.time,lipdTS.value,setTime[time]['min'],setTime[time]['max'],direction)
                for val in ['Direction','Slope','SlopeSig']: data_matrix[time+val].append(trendvalues[val])  
        except: problem_records.append(ts['dataSetName'])
    print(' --- Problematic records (n='+str(len(problem_records))+') ---')
    for record in problem_records: print(record)
    return(data_matrix)
    
dataHC = LiPDmetadata(lipdHC,True)
dataT  = LiPDmetadata(lipdT,False)

#%%
#Convert to dataFrame and upload to gitHub folder. Next step id regions and model pseudoproxies
#
cal = True
if cal: tsidList = ['GHe74bbbb1','LPDe897b7c7','WEB8783f810','WEB85c68a3a','WEBc0fa285c','WEB9660e756','WEB820982f9']
else:   tsidList = ['GHe74bbbb1','LPDe897b7c7','WEBa0656c03','WEBc9694a07','WEBd1b743fa','WEB801cc768','WEB720c22d6']
data_HC = pd.DataFrame.from_dict(dataHC); 
print("No. of hydroclimate records");                    print(len(data_HC))
data_HC = data_HC.loc[(data_HC['Direction'] != '')];     print(len(data_HC))
data_HC = data_HC.loc[(data_HC['Season'] != 'summer+') & 
                      (data_HC['Season'] != 'winter+')]; print(len(data_HC))
data_HC = data_HC.iloc[[i for i,val in enumerate(data_HC.TSid) if val not in tsidList]]; print(len(data_HC))


data_T = pd.DataFrame.from_dict(dataT); 
print("No. of temperature records");                     print(len(data_T))
data_T = data_T.loc[(data_T['Unit'] == 'degC')];         print(len(data_T))
data_T = data_T.loc[(data_T['Season'] != 'summer+') & 
                    (data_T['Season'] != 'winter+')];    print(len(data_T))

data_HC.to_csv(gitHubDir+'DataFiles/proxyHC.csv')
data_T.to_csv( gitHubDir+'DataFiles/proxyT.csv')






#%%
#Plot figures based on metadata (category or season) with isnet maps
#
FigKeyCategory={'Speleothem':   {'s':1,  'marker':'^', 'c':'firebrick'}, 
        'Speleothem (d18O)':    {'s':1,  'marker':'^', 'c':'firebrick'},
        'Speleothem (d13C)':    {'s':1,  'marker':'v', 'c':'indianred'},
        'Speleothem (other)':   {'s':1,  'marker':'>', 'c':'r'},
        'Lake Deposits':        {'s':1.5,'marker':'.', 'c':'cornflowerblue'},
        'Glacier Ice':          {'s':1,  'marker':'D', 'c':'powderblue'},
        'Lake Sediment (d18O)': {'s':1,  'marker':'o', 'c':'darkblue'},
        'Leaf Wax (dD)':        {'s':1,  'marker':'p', 'c':'darkorchid'}, 
        'Pollen':               {'s':1,  'marker':'P', 'c':'forestgreen'},  
        'Pollen (calibrated)':  {'s':1,  'marker':'P', 'c':'seagreen'},    
        'Pollen (uncalibrated)':{'s':1,  'marker':'X', 'c':'olive'}, 
        'Other':                {'s':1,  'marker':'s', 'c':'grey'}, 
        'Other (calibrated)' :  {'s':1,  'marker':'s', 'c':'dimgrey'},  
        'Other (uncalibrated)': {'s':1,  'marker':'s', 'c':'lightgrey'}
        }  
FigKeySeason = {'summerOnly':   {'s':1,  'marker':'^', 'c':'r'},
                'winterOnly':   {'s':1,  'marker':'v', 'c':'b'},
                'Annual':       {'s':0.5,  'marker':'o', 'c':'dimgrey'},
                'annual':       {'s':0.5,  'marker':'o', 'c':'dimgrey'},
                'not specified':{'s':0.5,  'marker':'o', 'c':'dimgrey'},
                '':             {'s':0.5,  'marker':'o', 'c':'dimgrey'}
        }
FigKeyExtent = {'Global':{'pltExtent':[0,8,0,10],'GeoExtent':[-180,180,-90,90],'Name':'Global'},
                 'NA':    {'pltExtent':[8,12,2,6],'GeoExtent':[-130,-50,23,55], 'Name':'North America'},
                 'Eur':   {'pltExtent':[8,12,6,9],'GeoExtent':[-12,40,35,70],    'Name':'Europe'}
                 }
for FigKey in [FigKeySeason,FigKeyCategory]:
    plt.style.use('ggplot')
    plt.figure(figsize=(20,10))
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['axes.facecolor'] ='white'; plt.rcParams['axes.linewidth'] = 1; plt.rcParams['axes.edgecolor'] = 'k'
    gs = gridspec.GridSpec(12,12,wspace=0,hspace=0.2)
    for ax in FigKeyExtent: 
        extent = FigKeyExtent[ax]['pltExtent']
        ax1 = plt.subplot(gs[extent[0]:extent[1],extent[2]:extent[3]],projection=ccrs.PlateCarree())
        ax1.set_extent(FigKeyExtent[ax]['GeoExtent'])
        ax1.coastlines()
        ax1.spines['geo'].set_edgecolor('black')
        ax1.add_feature(cfeature.LAND, facecolor='whitesmoke',edgecolor='k')
        ax1.add_feature(cfeature.LAKES,facecolor='none',      edgecolor='k')
        if FigKey  == FigKeyCategory: name = 'CategorySpecific'
        elif FigKey == FigKeySeason:  name = 'Season'
        if ax == 'Global': 
              ax1.add_feature(cfeature.BORDERS,facecolor='none',edgecolor='grey',lw=0.2)            
        elif ax == 'NA':
              ax1.add_feature(cfeature.BORDERS,facecolor='none',edgecolor='k',lw=0.3)        
              ax1.add_feature(cfeature.STATES, facecolor='none',edgecolor='grey',lw=0.2)        
        else: ax1.add_feature(cfeature.BORDERS,facecolor='none',edgecolor='k',lw=0.3)            
        for category in sorted(data_HC[name].unique())[::-1]: 
            plot_df = data_HC.loc[(data_HC[name] == category)]
            proxy_scatter = ax1.scatter(plot_df.Lon,plot_df.Lat,c=FigKey[category]['c'],
                edgecolor='k',lw=1,alpha=0.8,transform=ccrs.PlateCarree(),
                marker=FigKey[category]['marker'],s=FigKey[category]['s']*100,
                label=category)
        if ax == 'Global': 
            proxylegend = ax1.legend(loc='center left',bbox_to_anchor=(0,.35))
    proxylegend.get_frame().set_alpha(1)        
    plt.draw()
    if save: plt.savefig(gitHubDir+'Figures/HC12k_'+name+'.png', dpi=400,format='png')
    else: plt.show()
#%%
for variable in ['EarlySlopeSig','LateSlopeSig','bin6ka']: #plot proxy values to check calculations
    plt.style.use('ggplot')
    plt.figure(figsize=(20,10)); plt.rcParams['axes.facecolor'] ='white'
    plt.rcParams['axes.linewidth'] = 1; plt.rcParams['axes.edgecolor'] = 'k'
    ax1 = plt.subplot(projection=ccrs.Robinson()) 
    ax1.spines['geo'].set_edgecolor('black')
    ax1.set_global(); ax1.add_feature(cfeature.LAND,facecolor='whitesmoke',edgecolor='k')
    ax1.coastlines(); ax1.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
    plot_df = data_HC; plot_df = plot_df[plot_df[variable].notna()]
    proxy_scatter = ax1.scatter(plot_df.Lon,plot_df.Lat,c=plot_df[variable],
            marker='o',s=130,edgecolor='k',lw=1,alpha=0.8,transform=ccrs.PlateCarree(),
            cmap='BrBG',vmin=-1,vmax=1)
    plt.title(variable,fontsize=30)
    plt.show()
    
#%%
dataSource = []
for ts in lipdHC:
    if ts['paleoData_TSid'] in dataHC['TSid']:
        test = {}
        try:    test['url'] = ts['originalDataUrl']
        except: test['url'] = ''
        try:    test['createdBy'] = ts['createdBy']
        except: test['createdBy'] = ''
        try:    test['calibration'] = ts['paleoData_calibration_method']
        except: test['calibration'] = ''
        try:    test['doi1'] = ts['pub1_doi']
        except: test['doi1'] = ''
        try:    test['doi2'] = ts['pub2_doi']
        except: test['doi2'] = ''
        if   test['createdBy'] == 'http://github.com/nickmckay/oxfordLakeStatus2Lipd':
                dataSource.append('Oxford Lake Status')
        elif test['createdBy'] == 'sisal2lipd': 
                dataSource.append('SISAL')
        elif test['url'] == 'wNAm': 
                dataSource.append('wNA')
        elif test['calibration'] == 'JM18_MAT':
                dataSource.append('Marsicek')
        elif 'gov/paleo/study/15444' in  test['url'] or '10.5194/cp-10-1605-2014' == test['doi2']: 
                dataSource.append('Arctic')
        elif ts['dataSetName'][0:2] == 'LS':
                dataSource.append('iso2k')
        elif '10.25921/4RY2-G808' in test['url'] or '/paleo/study/27330' in test['url']: 
                dataSource.append('Temp12k')
        else: dataSource.append('Miscellanious')


datasource = {'source':[],'count':[],'sourceNo':[]}
for source in dataSource:
    if source not in datasource['source']:
        datasource['source'].append(source)
        datasource['sourceNo'].append(source+' (n='+str(dataSource.count(source))+')')
        datasource['count'].append(dataSource.count(source))
        
fig1, ax1 = plt.subplots()
ax1.pie(datasource['count'], labels=datasource['sourceNo'],startangle=90,
        colors=['lightgrey','firebrick','k','powderblue','forestgreen','cornflowerblue','rosybrown','darkslategrey'])
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.show()
if save: plt.savefig(gitHubDir+'Figures/HC12k_DataSourceAttempt'+'.png', dpi=400,format='png')
else: plt.show()

