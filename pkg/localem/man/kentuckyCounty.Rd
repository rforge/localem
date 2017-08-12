\docType{data}
\name{kentuckyCounty}
\alias{kentuckyCounty}
\alias{kentuckyTract}
\alias{kentuckyPuma}
\alias{kMap}
\title{Cancer case data in Kentucky}
\format{SpatialPolygonsDataFrame of Kentucky census boundaries with cancer cases and population data.}
\usage{
data('kentuckyCounty')
data('kentuckyTract')
data('kentuckyPuma')
data('kMap')
}
\description{
Lung cancer case counts, population, and map of Kentucky.
}
\examples{
# plot case data
data('kentuckyCounty')
data('kMap')

countyCol = mapmisc::colourScale(kentuckyCounty$count, 
							breaks = 7, style = 'quantile', dec = -1)

mapmisc::map.new(kentuckyCounty)
plot(kentuckyCounty, 
	col = countyCol$plot, 
	add = TRUE)
plot(kMap, add = TRUE)
mapmisc::legendBreaks('topleft', countyCol)

# plot population data
data('kentuckyTract')

kCol = mapmisc::colourScale(kentuckyTract$expected, 
						breaks = 9, style = 'equal', dec = 0)

mapmisc::map.new(kentuckyTract)
plot(kentuckyTract, 
	col = kCol$plot, border = '#00000010', 
	add = TRUE)
plot(kMap, add = TRUE)
mapmisc::legendBreaks('topleft', kCol)
}
\keyword{datasets}
