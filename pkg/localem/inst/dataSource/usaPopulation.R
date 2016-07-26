
# install.packages("Pmisc", repos='http://r-forge.r-project.org')
library("Pmisc")


# rates
lungRates = diseasemapping::cancerRates(
		'usa', 2000, 'F', 'Lung'
)

# table of counties and fips codes

countyTable = read.table(
		url('http://www2.census.gov/geo/docs/reference/codes/files/national_county.txt'),
		header=FALSE, stringsAsFactors=FALSE, sep=',', quote="\"",
		col.names=c('ab1','id1','subId2','name2','junk'),
		colClasses=rep('character',5)
)
countyTable$id2 = paste(
		'USA',countyTable$id1,
		'_0', countyTable$subId2,
		sep=''
)		


# Kentucky tract

kentuckyTract = gcensus(
		country = 'USA',
		year = 2000,
		subset = c(name1='Kentucky'),
		level = 4, income=TRUE,
		simplifyTolerance=0.002,
		crs='EPSG:3088'
)


# expected

kentuckyTract = diseasemapping::getSMR(
		kentuckyTract, 
		5*lungRates
)		

kentuckyTract = kentuckyTract[,grep("^(m|f)_", names(kentuckyTract), invert=TRUE)]		

save(kentuckyTract,
		file='/home/patrick/workspace/localem/pkg/localem/data/kentuckyTract.RData',
		compress='xz')		

file.info('/home/patrick/workspace/localem/pkg/localem/data/kentuckyTract.RData')['size']

if(FALSE) {
	library(raster)
	kMapSmall = mapmisc::openmap(kentuckyTract[1,], buffer=10000)
	mapmisc::map.new(kMapSmall)
	plot(kMapSmall, add=TRUE)
	plot(kentuckyTract, add=TRUE, border='#0000FF50',lwd=3)
	
}


# Background map

kMapOrig  = mapmisc::openmap(kentuckyTract, path='stamen-toner')
kMap = mapmisc::tonerToTrans(kMapOrig)

save(kMap,
		file='/home/patrick/workspace/localem/pkg/localem/data/kMap.RData',
		compress='xz')		
file.info('/home/patrick/workspace/localem/pkg/localem/data/kMap.RData')['size']


# kentucky cases


kentuckyCases = diseasemapping:::usCancer(
  site = 'Lung',
  sex = 'F',
  state = 'Kentucky',
  year = c(2003, 2007)
)


kentuckyCases$ab1 = toupper(kentuckyCases$stateCode)
kentuckyCases$name2 = paste(kentuckyCases$County, "County", sep=' ')

kentuckyCases = merge(
		kentuckyCases,
		countyTable[,c('ab1','id2','name2')],
		all.x=TRUE, all.y=FALSE
)

table(is.na(kentuckyCases$id2))	


# Kentucky counties

kentuckyCounty = gcensus(
		country = 'USA',
		year = 2000,
		subset = c(name1='Kentucky'),
		level = 2, income=TRUE,
		simplifyTolerance=0.002,
		crs='EPSG:3088'
)		



kentuckyCounty = sp::merge(
		kentuckyCounty,
		kentuckyCases
)

kentuckyCounty = kentuckyCounty[,
		grep("^(id|name|cases)", names(kentuckyCounty), ignore.case=TRUE)]		


kentuckyCounty = sp::spTransform(
		kentuckyCounty,
		kentuckyTract@proj4string
)			

if(FALSE) {
	library('raster')
	kCol = mapmisc::colourScale(
			kentuckyCounty$Cases,
			breaks=9, dec=-1,style='quantile'
	)
	mapmisc::map.new(kentuckyCounty)
	plot(kentuckyCounty, add=TRUE,
			col=kCol$plot)
	plot(kMap, add=TRUE)
	mapmisc::legendBreaks("topleft", kCol)
	
	
}	

save(kentuckyCounty,
		file='/home/patrick/workspace/localem/pkg/localem/data/kentuckyCounty.RData',
		compress='xz')		
file.info('/home/patrick/workspace/localem/pkg/localem/data/kentuckyCounty.RData')['size']


# California, tract

californiaTract = gcensus(
		country = 'USA',
		year = 2000,
		subset = c(name1='California'),
		level = 4, income=TRUE,
		simplifyTolerance=0.002,
		crs='EPSG:3310'
)


californiaTract = diseasemapping::getSMR(
		californiaTract, 
		lungRates
)		

californiaTract = californiaTract[,grep("^(m|f)_", names(californiaTract), invert=TRUE)]		

save(californiaTract,
		file='/home/patrick/workspace/localem/sandbox/californiaTract.RData',
		compress='xz')		

file.info('/home/patrick/workspace/localem/sandbox/californiaTract.RData')['size']



if(FALSE) {
	library(raster)
	cMapSmall = mapmisc::openmap(californiaTract[1,], buffer=10000)
	mapmisc::map.new(cMapSmall)
	plot(cMapSmall, add=TRUE)
	plot(californiaTract, add=TRUE, border='#0000FF50',lwd=3)
}



# 		california cases


californiaCases = diseasemapping:::usCancer(
  site = 'Lung',
  sex = 'F',
  state = 'California',
  year = 2003
)

californiaCases$ab1 = toupper(californiaCases$stateCode)
californiaCases$name2 = paste(californiaCases$County, "County", sep=' ')

# California counties

californiaCounty = gcensus(
		country = 'USA',
		year = 2000,
		subset = c(name1='California'),
		level = 2, income=TRUE,
		simplifyTolerance=0.002,
		crs='EPSG:3310'
)

calTable = countyTable[countyTable$ab1 == 'CA',]
calTable$name2Combined = calTable$name2

multiCounties = grep("-", californiaCases$County, value=TRUE)
for(D in multiCounties) {
	DnoSp = gsub("[[:space:]]", "", D)
	toGrep = paste("^(", gsub("-", "|", DnoSp), ")County", sep='')
	cHere = grep(toGrep, gsub("[[:space:]]", "", calTable$name2))
	calTable[cHere, 'name2Combined'] = paste(D, 'County')
}


californiaCounty = sp::merge(
		californiaCounty,
		calTable[,c('id2','name2Combined')]
)

californiaMerged = rgeos::gUnaryUnion(
		californiaCounty,
		californiaCounty$name2Combined
)

rownames(californiaCases) = gsub("[[:space:]]", "", californiaCases$name2)
californiaMerged = sp::spChFIDs(californiaMerged, gsub("[[:space:]]", "", names(californiaMerged) ))

californiaCases2 = californiaCases[match(
				gsub("[[:space:]]", "", names(californiaMerged)),
				gsub("[[:space:]]", "", californiaCases$name2)
		), ]				


californiaMerged = sp::SpatialPolygonsDataFrame(
		californiaMerged, californiaCases2[,c('name2','Cases')]
)


californiaMerged = sp::spTransform(
		californiaMerged,
		californiaTract@proj4string
)			

californiaCounty = sp::spTransform(
		californiaCounty,
		californiaTract@proj4string
)			


if(FALSE) {
	cMapOrig = mapmisc::openmap(californiaMerged, path='stamen-toner', fact=3)
	cMap = mapmisc::tonerToTrans(cMapOrig)
	
	library('raster')
	cCol = mapmisc::colourScale(
			californiaMerged$Cases,
			breaks=9, dec=-1,style='quantile'
	)
	mapmisc::map.new(californiaCounty)
	plot(californiaMerged, add=TRUE,
			col=cCol$plot)
	plot(californiaCounty, border='blue',add=TRUE)
	plot(californiaMerged, border='grey',add=TRUE, lwd=2)
	plot(cMap, add=TRUE)
	mapmisc::legendBreaks("topright", cCol)
	
	
}	


save(californiaMerged,
		file='/home/patrick/workspace/localem/sandbox/californiaMerged.RData',
		compress='xz')		
file.info('/home/patrick/workspace/localem/sandbox/californiaMerged.RData')['size']

