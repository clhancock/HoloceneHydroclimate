library(maptools)
library(dplyr)
library(lipdR)
library(reactable)
library(reactablefmtr)
library(htmlTable)
library(htmltools)
library(htmlwidgets)
library(webshot)
dir     <- '/Volumes/GoogleDrive/My Drive/zResearch/Manuscript/2021_HoloceneHydroclimate/2021_HoloceneHydroclimate' #
var <- 'HC'

#Load Data
lipdData <- readRDS(file=file.path(dir,'Data','Proxy','LiPD','lipdData.rds'))
lipdTSO <- lipdData[[var]]
#regList <- readShapePoly(file.path(dir,'Data','IPCC_refRegions','IPCC-WGI-reference-regions-v4.shp'))


#Create table with metadata from LiPD files
tbl <- tibble(Region         = pullTsVariable(lipdTSO,'geo_ipccRegion'),
              Dataset        = pullTsVariable(lipdTSO,'dataSetName'),
              TSid           = pullTsVariable(lipdTSO,'paleoData_TSid'),
              Archive        = pullTsVariable(lipdTSO,'archiveType'),
              Category       = pullTsVariable(lipdTSO,'CategorySpecific'),
              Proxy          = pullTsVariable(lipdTSO,'paleoData_proxy'),
              Season         = pullTsVariable(lipdTSO,'climateInterpretation1_seasonalityGeneral'),
              Interp         = pullTsVariable(lipdTSO,'climateInterpretation1_variable'),
              Direction      = pullTsVariable(lipdTSO,'climateInterpretation1_direction'),
              PublicationDOI = pullTsVariable(lipdTSO,'pub1_doi'),
              SourceURL      = pullTsVariable(lipdTSO,'originalDataUrl'),
              AgeRange       = paste(round(pullTsVariable(lipdTSO,'ageMax')/1000,1),
                                     round(pullTsVariable(lipdTSO,'ageMin')/1000,1),sep=' - '),
              Resolution     = round(pullTsVariable(lipdTSO,'ageRes')),
              AgeControlN    = as.character(pullTsVariable(lipdTSO,'chronData_agesN_12k')),
              AgeControlMax  = as.character(pullTsVariable(lipdTSO,'chronData_agesMaxGap_12k')),
              Lat            = round(pullTsVariable(lipdTSO,'geo_latitude'),2),
              Lon            = round(pullTsVariable(lipdTSO,'geo_longitude'),2)) #Source

#Modify url to include https
for (i in 1:nrow(tbl)){
 row <- tbl[i,]
 if(substr(row$PublicationDOI, 1, 3) == '10.'){row$PublicationDOI <- paste('https://doi.org/',row$PublicationDOI,sep='')}
 if(substr(row$PublicationDOI, 1, 3) == '10.'){row$PublicationDOI <- paste('https://doi.org/',row$PublicationDOI,sep='')}
 if(is.na(row$AgeControlN)){
   if(row$Category == "Shoreline (Lake Level)"){
     row$AgeControlN<-"Shoreline Dates"
   } else if (row$Archive %in% c("GlacierIce","Wood")){
     row$AgeControlN<-"Layer Counting"
   } else{row$AgeControlN<-"Not Available"}
 }
 if(is.na(row$AgeControlMax)){row$AgeControlMax<-"Not Available"}
 tbl[i,] <- row
}



#Order data by region and dataset and replace region name from Acronym to full name
tbl2 <- tbl[1,]
#tbl2 <- NA
start <- 1
for (reg in as.character(refregions$Acronym)){
  regData <- tbl %>% filter(tbl$Region==reg) %>% arrange(Dataset)
  regData$Region <- as.character(refregions$Name[which(refregions$Acronym==reg)])
  if (nrow(regData) == 0){next}
  end <- start + nrow(regData) - 1 
  tbl2[start:end,] <- regData
  start <- end + 1
}


#Id Divides between regions
#idx <- c(lengthtbl2$Region))
#idx2 <- c()
#for (i in 1:(length(tbl2$Region)-1)){
 #if(tbl2$Region[i] != tbl2$Region[i+1]){
  #idx <- c(idx,i)
  #idx2 <- c(idx2,i+1)
 #}
#}
filelist <- as.vector(list.files(path = file.path(dir,"Figures","Proxy","Dashboard")))

#Create Table
outTbl <- reactable(tbl2,
 groupBy = "Region",
 defaultPageSize = nrow(tbl2),
 height = 500,#), nrow(tbl2),
 wrap         = FALSE,
 resizable    = TRUE,
 searchable   = FALSE,
 filterable   = TRUE,
 sortable     = FALSE,
 showSortable = FALSE,
 borderless   = TRUE,
 compact      = TRUE,
 theme = reactableTheme(
  borderWidth = 1.2, borderColor='Black',
  headerStyle = list(`border-bottom` = "double",`border-top` = "thin solid", borderColor = "#555",borderWidth = 2,align = 'left')),
 style = list(fontFamily = 'sans',fontSize = 10),
 #rowStyle = function(index) {
  #if (index %in% idx) {list(      `border-bottom` = "thin solid")}
  #else if (index %in% idx2) {list(`border-top   ` = "thin solid")}},
 columns   = list(
  Region          = colDef(width=160,align = 'left'),
  Dataset         = colDef(width=150,align = 'left', name='Dataset Name \n (Site.Author.Year)', cell = function(value, index){
    htmltools::tags$a(href = paste('http://lipdverse.org/HoloceneHydroclimate/current_version/',value,'.html',sep=''), target = "_blank", as.character(value))}),
  TSid            = colDef(width=80, align = 'left', name='TSid \n (unique)', cell = function(value, index){
    htmltools::tags$a(href = paste("https://raw.githack.com/clhancock/HoloceneHydroclimate/main/Figures/Proxy/Dashboard/",
                                   filelist[which(grepl(value,filelist))],sep=""), target = "_blank", as.character(value))}),
  Archive         = colDef(width=80, align = 'left'),
  Category        = colDef(width=80, align = 'left'),
  Proxy           = colDef(width=80, align = 'left'),
  Interp          = colDef(width=45, align = 'center'),
  Season          = colDef(width=60, align = 'center'),
  Direction       = colDef(width=60, align = 'center'),
  PublicationDOI  = colDef(width=150,align = 'left', cell = function(value, index){
    htmltools::tags$a(href = value, target = "_blank", as.character(value))}),
  #SourceURL  = colDef(width=150,align = 'left', cell = function(value, index){
    #htmltools::tags$a(href = value, target = "_blank", as.character(value))}),
  SourceURL = colDef(width=150,align = 'left'),
  AgeRange        = colDef(width=60, align = 'center', name='Age Range \n (ka)',filterable=FALSE),
  Resolution      = colDef(width=60, align = 'center',name='Resolution \n (yrs/sample)'),
  AgeControlN     = colDef(width=80, align = 'center', name='Age Control \n (#)'),
  AgeControlMax   = colDef(width=120, align = 'center', name='Max Age Control Gap \n (yrs)'),
  Lat             = colDef(width=40, align = 'right'),
  Lon             = colDef(width=50, align = 'right')
  ))  %>%
  add_title("Appendix 1. List of proxy records included in the Holocene Hydroclimate dataset.
             Data are grouped by geographical region which are ordered according to Iturbide et al., (2020). 
             Within each region, records are listed alphabetically according to their dataset name.
             Columns can be resized using the column boundaries in the header row. 
             The empty boxes below the heading allow you to search for specific records. A full discription of each row is provided in the appendix file of Hancock et al., 2022",
            font_size=12, font_weight='normal', font_style='italic')

outTbl




html_file <- file.path(dir,'Figures','Proxy','Appendix1_Table','Appendix1.html')
saveWidget(widget = outTbl, file = html_file, selfcontained = TRUE)
