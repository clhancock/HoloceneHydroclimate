#Load packages
import pandas as pd
import pyleoclim as pyleo
#import Lipd as lipd
import numpy as np
import pymannkendall as mk              # Package for trend detection
from scipy import stats    # Packages for calculations  
#import math

dataDir='/Volumes/GoogleDrive/My Drive/zResearch/Data/LiPD/'

#%%
#Load all available LiPD files & 
#
LiPDdatabase   = pyleo.Lipd(dataDir+'database')
LiPDdatabaseTS = LiPDdatabase.to_tso()
LiPDnew        = pyleo.Lipd(dataDir+'NewHoloceneHydroclimate')
LiPDnewTS      = LiPDnew.to_tso()
#%%
# filter timeseries relevant to temp12k and hc12k
#
lipdHC = []; lipdT = []; dataSetList = []
#Add files from newdata folder
for ts in LiPDnewTS: #These files don't have proper metadata
    if ts['paleoData_TSid'] in dataSetList: continue #To prevent duplicate files
    dataSetList.append(ts['paleoData_TSid']); dataSetList.append(ts['dataSetName'])
    try: #Not all ts have an interpretation to search for 
        if ts['paleoData_interpretation'][0]['variable'] in ['P','P-E','P_E','M']: 
            lipdHC.append(ts)
    except: pass
#Add files from general database folder without duplicates
for ts in LiPDdatabaseTS:
    if ts['dataSetName'] in dataSetList or ts['paleoData_TSid'] in dataSetList: continue #To prevent duplicate files
    dataSetList.append(ts['paleoData_TSid'])
    try: #Not all ts have a compilation data to search for 
        for compilation in ts['paleoData_inCompilationBeta']:
            #Check if in temp12k
            if compilation['compilationName'] == "Temp12k":
                if any(v == '1_0_2' for v in compilation['compilationVersion']): 
                    lipdT.append(ts) ###
            #Check if in HC12k
            elif compilation['compilationName'] == "HoloceneHydroclimate":
                ts['paleoData_values'] = [np.nan if x == 'nan' else x for x in ts['paleoData_values']] #Convert Lake Deposit records to integers
                if ts['archiveType'] == 'LakeDeposits':
                    ts['paleoData_values'] = [x-2 for x in ts['paleoData_values']] #Change 1,2,3 values to -1,0,1  
                    if np.count_nonzero(~np.isnan(ts['paleoData_values'][0:11])) >= 9: #Set minimal resolution
                        lipdHC.append(ts) ###
                else: lipdHC.append(ts) ###
    except: pass

print(' - Found '+str(len(lipdT))+' temperature proxies -')
print(' - Found '+str(len(lipdHC))+' hydroclimate proxies -')
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
    divisions = 3; minset = 2; sigLevel = 0.05; CI = 0.95; method = 'lin'
    if len(proxyValues) < 12: minset = 1
    ageRange = (ageMax - ageMin)/divisions
    check = []
    valuedata = pd.DataFrame([timeValues,proxyValues]).transpose()
    valuedata.columns = ['time','values']
    valuedata = valuedata.loc[(valuedata['time'] <= ageMax) & (valuedata['time'] >= ageMin)]
    valuedata = valuedata.dropna() 
    for timeperiod in range(divisions):
        checkValue = len(valuedata.loc[(valuedata['time'] <= ageMin+(timeperiod+1)*ageRange) 
                                     & (valuedata['time'] >= ageMin+(timeperiod)*ageRange)]['time'])
        if checkValue >= minset: check.append(True)
        else: check.append(False)
    if check.count(True) == divisions:
        valuedata['time']   = valuedata['time']*-1 #To get correct forward direction
        valuedata['values'] = valuedata['values']*direction
        trendDirection = mk.original_test(valuedata['values']*-1,alpha=sigLevel)
        if method == 'lin':   Slope = stats.linregress(valuedata['time'],valuedata['values'])
        elif method == 'sen': Slope = stats.theilslopes(valuedata['values'],valuedata['time'],CI)
        output = {'Direction':trendDirection[0],'Slope':Slope[0]*1000}
        if output['Direction'] == 'no trend': output['SlopeSig'] = 0
        else:                                 output['SlopeSig'] = output['slope']
    else: output = {'Direction':np.NaN, 'Slope':np.NaN, 'SlopeSig':np.NaN}
    return(output)

#Function to create dataframe of relevent data / quick qc sheet to sort data by
def df_LiPDmetadata (lipd_list,standardize):
    problem_records = []; RecordNo=0
    data_matrix = {'Category':[],'CategorySpecific':[],
                   'ageMin':[],'ageMax':[],'ageRange':[],'ageRes':[]}
    names_metadata = {'TSid':'paleoData_TSid','Lat':'geo_meanLat','Lon':'geo_meanLon',
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
            lipdTS = ts['lipdSeries'].slice([ageMin,ageMax])
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
            else:                                          direction = 1
            if standardize: lipdTS = lipdTS.standardize() 
            lipdTS_ref = np.nanmean(lipdTS.slice([ageMinRef,ageMaxRef]).value)
            for age in range(ageMin,ageMax,ageRes):
                value_ts = np.nanmean(lipdTS.slice([age,age+ageRes]).value)
                data_matrix['bin'+str(int(age/ageRes))+'ka'].append((value_ts-lipdTS_ref)*direction)
            for time in setTime:
                trendvalues = calcTrend(lipdTS.time,lipdTS.value,setTime[time]['min'],setTime[time]['max'],direction)
                for value in ['Direction','Slope','SlopeSig']: data_matrix[time+value].append(trendvalues[value])  
        except: problem_records.append(ts['dataSetName'])
    print(' --- Problematic records (n='+str(len(problem_records))+') ---')
    for record in problem_records: print(record)
    return(data_matrix)
    
data_hc   = df_LiPDmetadata(lipdHC,True)
data_temp = df_LiPDmetadata(lipdT,False)
#%%
