#*********************************************
#*********************************************
#' This package provides code for simulating acoustic-trawl survey estimates from spatial biomass data, developed for the PELFOSS project (Observational system simulator for pelagic fish stocks).
#'
#' @export
#' @rdname pelfoss
#' 
pelfoss <- function(){}

#*********************************************
#*********************************************
#' Read a NORWECOM file.
#'
#' \code{readNORWECOMbiomass} reads a NORWECOM biomass file, with biomass on a irregular grid Long, Latt.
#' \code{readNORWECOMsuperind} reads a NORWECOM super individual file.
#' \code{readNcVarsAndAttrs} extracts variables with long names as attributes.
#' \code{addTimeHoursFrom1950} adds time from the "hours since 1950"-variable of the NORWECOM data.
#' \code{getBiomassXY} converts to Cartesian coordinates for the biomass data.
#' \code{getSuperIndXY} converts to Cartesian coordinates for the superindividual data.
#'
#' @param x	The input NetCDF4 file.
#' 
#' @return
#'
#' @export
#' @rdname readNORWECOMbiomass
#' 
readNORWECOMbiomass <- function(
	ncfile, centroid, 
	vars = c("MACbiom", "Long", "Latt"), 
	rename = list(MACbiom="Biom")){
		
	# Read the NORWECOM superindividual file:
	out <- readNcVarsAndAttrs(ncfile, vars=vars, rename=rename)
	# Add time:
	out <- addTimeHoursFrom1950(out, timedim="T")
	
	# Use only the data object, and discard the nc object:
	out <- out$data
	
	##### NOTE 1: This should be fixed in the data: #####
	#out$Long[168, 206] <- mean(out$Long[ 167:169, 206], na.rm=TRUE)
	
	# Get Cartesian coordinates:
	out <- getBiomassXY(out, centroid=centroid)
	
	# Return:
	out
}
#'
#' @export
#' @rdname readNORWECOMbiomass
#' 
readNORWECOMsuperind <- function(
	ncfile, centroid, 
	vars = c("xpos", "ypos", "Long", "Latt", "female", "age", "inumb", "length", "sweight", "pweight"), 
	rename = list(xpos="gridLongInd", ypos="gridLattInd", Long="gridLong", Latt="gridLatt"), 
	removeIntitialNA = TRUE){
	
	# Read the NORWECOM superindividual file:
	out <- readNcVarsAndAttrs(ncfile, vars=vars, rename=rename)
	# Add time:
	out <- addTimeHoursFrom1950(out, timedim="time")
	
	# Use only the data object, and discard the nc object:
	out <- out$data
	
	# Remove superindividuals which are initially NA:
	if(removeIntitialNA){
		out <- getSuperIndOfDay(out, day=NULL, maxValue=1e10, NAsByFirst=TRUE)
	}
	
	# Get Cartesian coordinates:
	out <- getSuperIndXY(out, centroid=centroid)
	
	# Return:
	out
}
#'
#' @export
#' @importFrom ncdf4 nc_open ncvar_get ncatt_get
#' @keywords internal
#' @rdname readNORWECOMbiomass
#' 
readNcVarsAndAttrs <- function(ncfile, vars, rename=NULL){
	getNcVarAndAttr <- function(y, nc){
		out <- ncvar_get(nc, y)
		att <- ncatt_get(nc, y)$long_name
		attr(out, "long_name") <- att
		out
	}
	# Read the NORWECOM superindividual file:
	nc <- nc_open(ncfile)
	data <- lapply(vars, getNcVarAndAttr, nc=nc)
	
	if(length(rename)){
		rename <- rename[intersect(names(rename), vars)]
		if(length(rename)){
			toBeRenamed <- match(names(rename), vars)
			vars[toBeRenamed] <- unlist(rename)
		}
	}
	
	names(data) <- vars
	
	list(nc=nc, data=data)
}
#'
#' @export
#' @keywords internal
#' @rdname readNORWECOMbiomass
#' 
addTimeHoursFrom1950 <- function(x, timedim="time"){
	# Add time from the dimensions:
	x$data$time <- as.POSIXct("1950-01-01", tz = "UTC") + x$nc$dim[[timedim]]$vals * 60^2
	
	x
}
#'
#' @export
#' @keywords internal
#' @importFrom Rstox geo2xy
#' @rdname readNORWECOMbiomass
#' 
getBiomassXY <- function(biomass, centroid, proj="aeqd", units="kmi", x_0=0, y_0=0, ellps="WGS84", datum="WGS84"){
	# Convert the locations to x,y using the 'aeqd' projection and centered at 0, 68:
	thisLonLat <- data.frame(lon=c(biomass$Long), lat=c(biomass$Latt))
	xy_new <- Rstox::geo2xy(thisLonLat, list(proj=proj, units=units, lon_0=centroid[1], lat_0=centroid[2], x_0=x_0, y_0=y_0, ellps=ellps, datum=datum), data.frame.out=TRUE)

	# Define the grid to interpolate onto:
	biomass$X <- xy_new$x
	biomass$Y <- xy_new$y
	dim(biomass$X) <- dim(biomass$Long)
	dim(biomass$Y) <- dim(biomass$Long)
	
	biomass
}
#'
#' @export
#' @keywords internal
#' @importFrom fields interp.surface
#' @importFrom Rstox geo2xy
#' @rdname readNORWECOMbiomass
#' 
getSuperIndXY <- function(superInd, centroid, proj="aeqd", units="kmi", x_0=0, y_0=0, ellps="WGS84", datum="WGS84"){
	# Get the dimensions of the irregular geographical grid:
	dimxy <- dim(superInd$gridLong)
	# Get the sequence spanning the dimensions of the irregular grid (i.e., indices of the grid, which are the units of the positions xpos and ypos):
	seqx <- seq_len(dimxy[1])
	seqy <- seq_len(dimxy[2])
	# Define two lists of data, where x and y are the sequences of indices defined above in both lists, and the z is the actual longitude and latitude values respectively of the irregular grid:
	dataLon <- list(x=seqx, y=seqy, z=superInd$gridLong)
	dataLat <- list(x=seqx, y=seqy, z=superInd$gridLatt)
	
	# Define the matrix of points given as partial indices in the irregular grid:
	outpos <- cbind(c(superInd$gridLongInd), c(superInd$gridLattInd))
	
	# Interpolate each of longitude and latitude onto the superindividual positions:
	superInd$Long <- fields::interp.surface(dataLon, outpos)
	superInd$Latt <- fields::interp.surface(dataLat, outpos)
	
	# Reset dimensions:
	dimout <- dim(superInd$gridLongInd)
	dim(superInd$Long) <- dimout
	dim(superInd$Latt) <- dimout
	
	# Define the geographical coordinates of the superindividuals:
	#x$superIndLonLat <- data.frame(lon=c(outLon), lat=c(outLat))
	# lonlat <- interp::interpp(x=x$lonlat, z=z, xo=xo, yo=yo, output="points")
	
	# Convert to cartesian coordinates:
	thisLonLat <- data.frame(lon=c(superInd$Long), lat=c(superInd$Latt))
	xy_new <- Rstox::geo2xy(thisLonLat, list(proj=proj, units=units, lon_0=centroid[1], lat_0=centroid[2], x_0=x_0, y_0=y_0, ellps=ellps, datum=datum), data.frame.out=TRUE)
	
	superInd$X <- xy_new[,1]
	dim(superInd$X) <- dimout
	superInd$Y <- xy_new[,2]
	dim(superInd$Y) <- dimout
	
	superInd
}


#*********************************************
#*********************************************
#' Interpolate NORWECOM biomass onto log distances.
#'
#' @param x	The input NetCDF4 file.
#' 
#' @return
#'
#' @export
#' @importFrom interp interp
#' @rdname interpolateBiomassToTransects
#' 
interpolateBiomassToTransects <- function(biomass, transect, centroid, dayvec, margin=0.05, proj="aeqd", units="kmi", x_0=0, y_0=0, ellps="WGS84", datum="WGS84"){
	# Function to add margin to the ranges from which the points in the original locations are selected:
	addMargin <- function(x, margin=0.05){
		rangex <- range(x)
		rangex <- rangex + c(-1, 1) * margin[1] * diff(rangex)
		rangex
	}
	
	# Function for interpolating the biomass values onto the log distances for one specific day 'day', which is looked for in the vector 'dayvec', which is the days for each log distance:
	interpolateOneStratumOneDay <- function(stratum, day, biomass, transect, dayvec, margin){
		# Get the locations of the requested day, and convert NA to 0:
		z <- c(biomass$Biom[,,day])
		z[is.na(z)] <- 0
	
		# Use points in the original data only from the current day and stratum:
		inTransects <- which(dayvec==day & temp$stratum == stratum)
		xo <- transect$x_mid[inTransects]
		yo <- transect$y_mid[inTransects]
		
		# Interpolate only based on the points in proximety to the output points:
		rangex <- addMargin(xo, margin)
		rangey <- addMargin(yo, margin)
		valid <- c(biomass$X) >= rangex[1] & c(biomass$X) <= rangex[2] & c(biomass$Y) >= rangey[1] & c(biomass$Y) <= rangey[2]
		
		# Interpolate onto the output grid:
		#zo <- akima::interp(x=c(coords$x), y=c(coords$y), z=z, xo=xo, yo=yo)
		x <- c(biomass$X)[valid]
		y <- c(biomass$Y)[valid]
		z <- z[valid]
		
		zo <- interp::interp(x=x, y=y, z=z, xo=xo, yo=yo, output="points")
		out <- cbind(inTransects, zo$z)
		
		return(out)
	}
	
	# Function for interpolating the biomass values onto the log distances for one specific day 'day', which is looked for in the vector 'dayvec', which is the days for each log distance:
	interpolateOneDay <- function(day, biomass, transect, dayvec, margin){
		udays <- unique(dayvec)
		bio <- do.call(rbind, lapply(udays, interpolateOneDay, biomass=biomass, transect=transect, dayvec=dayvec, margin=margin))
	
		# Order by the days:
		bio[order(bio[,1]), 2]
	}
	
	# Function for interpolating the biomass values onto the log distances for one specific day 'day', which is looked for in the vector 'dayvec', which is the days for each log distance:
	interpolateOneDayOld <- function(day, biomass, transect, dayvec, margin, thr=10){
		# Get the locations of the requested day, and convert NA to 0:
		z <- c(biomass$Biom[,,day])
		z[is.na(z)] <- 0
		
		# Interpolate only based on the points in proximety to the output points:
		inDay <- which(dayvec==day)
		xo <- transect$x_mid[inDay]
		yo <- transect$y_mid[inDay]
		
		#margin <- 0.1
		#rangex <- range(transect$x_mid[inDay])
		#rangex <- rangex + c(-margin[1], margin[1])# * diff(rangex)
		#rangey <- range(transect$y_mid[inDay])
		#rangey <- rangey + c(-margin[2], margin[2])# * diff(rangey)
		
		valid <- NULL
		
		while(sum(valid)<thr){
			rangex <- addMargin(xo, margin)
			rangey <- addMargin(yo, margin)
			valid <- c(biomass$X) >= rangex[1] & c(biomass$X) <= rangex[2] & c(biomass$Y) >= rangey[1] & c(biomass$Y) <= rangey[2]
			margin <- margin * 4
		}
		
		#rangex <- addMargin(xo, margin)
		#rangey <- addMargin(yo, margin)
		#
		#valid <- c(biomass$X) >= rangex[1] & c(biomass$X) <= rangex[2] & c(biomass$Y) >= rangey[1] & c(biomass$Y) <= rangey[2]
		#if(sum(valid)==0){
		#	warning("")
		#	return(NULL)
		#}
		
	
		
		# Interpolate onto the output grid:
		#zo <- akima::interp(x=c(coords$x), y=c(coords$y), z=z, xo=xo, yo=yo)
		x <- c(biomass$X)[valid]
		y <- c(biomass$Y)[valid]
		z <- z[valid]
		#xo=transect$x_mid[inDay]
		#yo=transect$y_mid[inDay]
		zo <- interp::interp(x=x, y=y, z=z, xo=xo, yo=yo, output="points")
		#out <- cbind(inDay, zo$z)
		out <- cbind(day, zo$z)
		
		return(out)
	}
	
	#margin <- margin * c(diff(range(biomass$X)), diff(range(biomass$Y)))
	
	# Run through the 
	udays <- unique(dayvec)
	bio <- do.call(rbind, lapply(udays, interpolateOneDayOld, biomass=biomass, transect=transect, dayvec=dayvec, margin=margin))
	
	# Order by the days:
	bio <- bio[order(bio[,1]), ]
	
	# Fill in the biomass values at the indices with days present in the interpolation:
	out <- double(length(dayvec))
	out[dayvec %in% bio[, 1]] <- bio[,2]
	out
	#bio[order(bio[,1]), 2]
}


#*********************************************
#*********************************************
#' Extract one day of valid superindividuals.
#'
#' @param x	The input NetCDF4 file.
#' 
#' @return
#'
#' @export
#' @keywords internal
#' @rdname getSuperIndOfDay
#' 
getSuperIndOfDay <- function(superInd, day=NULL, maxValue=1e10, NAsByFirst=FALSE){
	# Get the number of days:
	numDays <- length(superInd$time)
	if(length(day) == 0){
		day <- seq_len(numDays)
	}
	# Get the indices of variables with days at the second dimension:
	withDays <- which(sapply(superInd, function(x) identical(ncol(x), numDays)))
	# Get the indices of superindividuals with valid data of the given day:
	if(NAsByFirst){
		valid <- which(superInd[[withDays[1]]][, 1] < maxValue)
	}
	else{
		valid <- which(superInd[[withDays[1]]][, day] < maxValue)
	}
	
	# Extract the valid values of the given day:
	superInd[withDays] <- lapply(superInd[withDays], "[", valid, day)
	# Return:
	superInd
}
#'
#' @export
#' @keywords internal
#' @rdname getSuperIndOfDay
#' 
getBiomassOfSurvey <- function(x, daysOfSurvey=NULL){
	if(length(x$biomassOfSurvey)==0){
		if(length(x$daysOfSurvey)){
			daysOfSurvey <- x$daysOfSurvey
		}
		# Get the biomass of the days of the survey:
		Biom <- x$biomass$Biom[,,daysOfSurvey]
		dim(Biom) <- c(prod(dim(Biom)[1:2]), length(daysOfSurvey))
		Biom <- rowMeans(Biom)
		x$biomassOfSurvey <- data.frame(Long = c(x$biomass$Long), Latt = c(x$biomass$Latt), Biom = Biom)
	}
	return(x)
}
#'
#' @export
#' @keywords internal
#' @rdname getSuperIndOfDay
#' 
getSuperIndOfSurvey <- function(x, midDayOfSurvey=NULL){
	if(length(x$superIndOfSurvey)==0  && length(x$superInd)){
		if(length(x$midDayOfSurvey)){
			midDayOfSurvey <- x$midDayOfSurvey
		}
		# Get the supeInd of the mid day of the survey:
		x$superIndOfSurvey <- as.data.frame(getSuperIndOfDay(x$superInd, midDayOfSurvey)[c("Long", "Latt")])
	}
	return(x)
}


#*********************************************
#*********************************************
#' Interpolate NORWECOM biomass onto log distances.
#'
#' @param x	The input NetCDF4 file.
#' 
#' @return
#'
#' @export
#' @importFrom sp point.in.polygon
#' @rdname getTotalBiomass
#' 
getTotalBiomass <- function(superInd, stratumPolygons, strataInd=NULL, days=NULL, type=c("biomass", "length"), lonlatFilter=NULL){
	
	sumBiomass <- function(d, inside){
		sum(d$pweight[inside] * d$inumb[inside])
	}
	
	meanLengthByTS <- function(d, inside){
		sqrt(weighted.mean(d$length[inside]^2, d$inumb[inside]))
	}
	
	
	# Function for getting the biomass of all strata of one day:
	getTotalBiomassStratum <- function(day, superInd, stratumPolygons, lonlatFilter){
		# Get the indices in the variable in thisSuperInd which are inside the straum:
		getBiomassInside <- function(stratumPolygon, thisSuperInd){
			if(length(stratumPolygon)){
				inside <- sp::point.in.polygon(
					#point.x = thisSuperInd$X, 
					#point.y = thisSuperInd$Y, 
					point.x = thisSuperInd$Long, 
					point.y = thisSuperInd$Latt, 
					pol.x = stratumPolygon[,1], 
					pol.y = stratumPolygon[,2])
				inside <- which(inside > 0)
			}
			else{
				inside <- seq_along(thisSuperInd$X)
			}
			
			# Apply also the filter given by the funciton lonlatFilter:
			if(is.function(lonlatFilter)){
				insideFilter <- which(lonlatFilter(thisSuperInd))
				inside <- intersect(inside, insideFilter)
			}
			
			# Get the total biomass:
			if(type[1] == "biomass"){
				total <- sumBiomass(thisSuperInd, inside)
			}
			else if(type[1] == "length"){
				total <- meanLengthByTS(thisSuperInd, inside)
			}
			else{
				stop("Invalid type in getTotalBiomass()")
			}
			#total <- sum(thisSuperInd$pweight[inside] * thisSuperInd$inumb[inside])
			total
		}
		
		# Get valid values of the given day:
		thisSuperInd <- getSuperIndOfDay(superInd, day)
	
		# Get the total biomass of the entire stock:
		total <- getBiomassInside(NULL, thisSuperInd=thisSuperInd)
		
		# Get the biomass of each stratum:
		byStratum <- sapply(stratumPolygons, getBiomassInside, thisSuperInd=thisSuperInd)
		
		# Return the biomass in the strata, the sum of all strata, and the total biomass:
		if(length(strataInd)==0){
			strataInd <- seq_along(byStratum)
		}
		c(byStratum, sum(byStratum[strataInd]), sum(byStratum), total)
	}
	
	# Get a sequence of the days:
	if(length(days)==0){
		days <- seq_along(superInd$time)
	}
	
	# Get the biomass for all days and combine to a data frame:
	biomass <- lapply(days, getTotalBiomassStratum, superInd=superInd, stratumPolygons=stratumPolygons, lonlatFilter=lonlatFilter)
	biomass <- do.call("rbind", biomass)
	biomass <- as.data.frame(biomass)
	biomassNames <- c(if(length(names(stratumPolygons))) names(stratumPolygons) else seq_along(stratumPolygons), "NonEmptyStrata", "AllStrata", "Total")
	names(biomass) <- biomassNames
	# Return:
	biomass
}


#*********************************************
#*********************************************
#' Simulate trawl stations given NORWECOM superindividuals.
#'
#' @param x	The input NetCDF4 file.
#' 
#' @return
#'
#' @export
#' @importFrom Rstox getSeedV
#' @rdname simulateTrawl
#'
simulateTrawl <- function(superInd, fishStation, seed=0, radius=10, nn=NULL, N=100, lengthunit=2){
	# Small funciton for calcualting the Eucledian distance given x, y, ...:
	edist <- function(...){
		x <- do.call(cbind, list(...))
		sqrt(rowSums(x^2))
	}
	mround <- function(x, base){
		base * round(x / base)
	}
	getValidLengthunit <- function(lengthunit){
		# From getNMDinfo("Lengdeenhet"):
		# name description
		#    1        1 mm
		#    2        5 mm
		#    3        1 cm
		#    4        3 cm
		#    5        5 cm
		#    6      0.5 mm
		#    7      0.1 mm
		validLengthunit <- c(1e-3, 5e-3, 10e-3, 30e-3, 50e-3, 0.5e-3, 0.1e-3)
		validLengthunit[as.numeric(lengthunit)]
	}
	
	# Function for simulating one trawl sample based on the surrounding superindividuals:
	simulateTrawlOne <- function(ind, superInd, fishStation, seedV, radius=10, nn=NULL, N=100, validLengthunit=2){
		
		# Get the day of the trawl station:
		day <- findInterval(fishStation$starttime[ind], superInd$time)
		#thisSuperInd <- lapply(superInd[c("X", "Y", "inumb", "female", "age", "length", "pweight")], function(x) x[, day])
		thisSuperInd <- getSuperIndOfDay(superInd, day)
		
		# First find the superindividuals within the radius around the station:
		d <- edist(
			thisSuperInd$X - fishStation$x_mid[ind], 
			thisSuperInd$Y - fishStation$y_mid[ind]
		)
		
		# Get the superindividuals inside the radius:
		# Also discard invalid superindividuals:
		if(length(nn)){
			inside <- which(thisSuperInd$pweight < 1e10)
			inside <- order(d[inside])[seq_len(nn)]
		}
		else{
			inside <- which(d <= radius & thisSuperInd$pweight < 1e10)
		}
		superIndCount = length(inside)
		
		# Get the number of fish of the valid superindividuals:
		size <- thisSuperInd$inumb[inside]
		
		# We could consider weighting by the distance in the below probability of selecting a fish from the superindividual:
		prob <- size / sum(size)
		set.seed(seedV[ind])
		if(superIndCount == 0){
			return(data.frame())
		}
		else if(superIndCount == 1){
			s <- rep(inside, N)
		}
		else{
			s <- sample(inside, size=N, prob=prob, replace=TRUE)
		}
		
		# Generate the individual samples:
		individual <- data.frame(
			serialno = fishStation$serialno[ind], 
			superInd = s,
			superIndCount = superIndCount,
			minDist = min(d),
			specimenno = seq_len(N), 
			lengthunit = lengthunit, 
			length = thisSuperInd$length[s],
			sex = 2 - thisSuperInd$female[s], # The definition of biotic xml has male=2 and female=1. In the Norwecom model male=0 and female=1.
			weight.individual = thisSuperInd$pweight[s],
			age..agedetermination.1 = thisSuperInd$age[s]
		)
		
		# Round off length to the lengthunit (given in cm in the NORWECOM model output):
		individual$length <- individual$length * 1e-2
		individual$length <- mround(individual$length, validLengthunit)
		# Round off weight to grams (given in grams in the MORWECOM output):
		individual$weight.individual <- individual$weight.individual * 1e-3
		individual$weight.individual <- round(individual$weight.individual, digits=3)
		
		
		individual
	}
	
	# Get the length unit from the reference data:
	### validLengthunit <- getNMDinfo("Lengdeenhet")
	### validLengthunit <- validLengthunit$description[validLengthunit$name == lengthunit]
	### # Extract the number and unit:
	### validLengthunit <- strsplit(validLengthunit, " ")[[1]]
	### validLengthunit <- as.numeric(validLengthunit[1]) * ifelse(validLengthunit[2] == "cm", 10, 1)
	### # Convert to meters:
	### validLengthunit <- validLengthunit * 1e-3
	validLengthunit <- getValidLengthunit(lengthunit)
	
	# Draw seeds:
	nstations <- nrow(fishStation)
	seedV <- Rstox::getSeedV(seed, nstations)
	
	# Draw trawl samples:
	individual <- lapply(seq_len(nstations), simulateTrawlOne, superInd=superInd, fishStation=fishStation, seedV=seedV, radius=radius, nn=nn, N=N, validLengthunit=validLengthunit)
	superIndCount <- sapply(individual, function(x) if(length(x$superIndCount)) head(x$superIndCount, 1) else 0)
	minDist <- sapply(individual, function(x) if(length(x$minDist)) head(x$minDist, 1) else NA)
	print(summary(superIndCount))
	print(summary(minDist))
	individual <- do.call("rbind", individual)
	
	# Merge the fish stations and individuals
	out <- list(biotic=merge(fishStation, individual), superIndCount=superIndCount)
	out
}
#'
#' @export
#' @rdname simulateTrawl
#'
drawTrawlStationInd <- function(nasc, nstations=1, seed=0, probfact=0){
	# Generate the cummulative probability distribution as the cumsum of the nasc, and draw based on this distribution:
	nasc0 <- nasc
	nasc0[is.na(nasc0)] <- 0
	
	prob <- nasc0 / sum(nasc0)
	# Add a small probability to all positions:
	prob <- prob + probfact / length(prob)
	
	if(nstations > length(prob)){
		warning("The number of trawl stations to draw exceeds the number og log distances. All log distances drawn (without replacement).")
		nstations <- length(prob)
	}
	
	ind <- sample(seq_along(prob), nstations, prob=prob, replace=FALSE)
	# Return:
	ind
}


#*********************************************
#*********************************************
#' Funciton for converting from a biomass density to nautical area scattering coefficient NASC (linear values, not dB).
#'
#' @param Wg	The biomass area density in units of g/m^2 (gram per square meter in the horizontal plane). 
#' @param a,b	Parameters of the length-weight relationship Wg = a * Lcm^b, where Wg is the weight in grams and Lcm is the length in cm. Typical values are e.g., a = 0.01 and b = 3.
#' @param m,TS0	The parameters of the target strength-length relationship TS = m * log10(Lcm) + TS0, typically m = 20 and TS0 = -71.9 (herring, from Foote, K. G. 1987. Fish target strengths for use in echo integrator surveys. Journal of the Acoustical Society of America, 82: 981– 987.)
#' 
#' @return
#'
#' @export
#' @rdname biomass2nasc
#'
biomass2nasc <- function(Wg, LcmOne, a, b, m=20, TS0=-71.9){
	# Function for converting from length in cm to weight in g:
	Lcm2Wg <- function(Lcm, a, b){
		# Standard length - weight relationship:
		# 	w = a * L^b
		a * Lcm^b
	}
	# Target strength of a fish with length given in cm:
	TS <- function(Lcm, m=20, TS0=-71.9){
		TS <- m * log10(Lcm) + TS0
	}
	# Backscattering cross section of a fish with length given in cm:
	sigmabs <- function(Lcm, m=20, TS0=-71.9){
		10^(TS(Lcm=Lcm, m=m, TS0=TS0) / 10)
	}
	
	# Get the weight of one fish:
	WgOne <- Lcm2Wg(Lcm=LcmOne, a=a, b=b)
	
	# Get the fractional number of fish:
	numFish <- Wg / WgOne
	
	# Get the backscattering of one fish:
	sigmabsOne <- sigmabs(LcmOne, m=m, TS0=TS0)
	
	# Convert to NASC:
	nasc <- 4 * pi * 1852^2 * numFish * sigmabsOne
	
	return(nasc)
}


#*********************************************
#*********************************************
#' Plot the output from biomass2tsb() or runPelfoss().
#'
#' @param x					The output from \code{\link{biomass2tsb}} or \code{\link{runPelfoss}}.
#' @param ncolour			The number of colours used to plot the biomass field in layers.
#' @param col				A list of hue, saturation and value ranges, between which a sequence of 'ncolour' is used for the biomass layers, or a color vector such as gray(ncolour, 0.2, 0.8).
#' @param logColscale		Logical: If TRUE, use 10 * log10 for the plotted biomass layers.
#' @param firstcol			The color to use for the first biomass level (the default, NA, omits plotting the first layer, which may include 0-values).
#' @param NASCthr			A value below which NASC values are omitted when plotted as black dots with size proportional to NASC^NASCexp.
#' @param NASCmax_size		Maximum size of the NASC dots.
#' @param biomassAlpha		Transparency of the biomass layers.
#' @param trawlSize			The size to use for the red asteriks representing trawl stations.
#' 
#' @return
#'
#' @export
#' @import ggplot2
#' @importFrom Rstox plotStratum 
#' @rdname plotPELFOSS
#'
plotPELFOSS <- function(
	x, 
	ncolour = 20, 
	col = list(h=c(0.6, 0.98), s=c(0.35, 1), v=c(0.9, 0.3)), 
	logColscale   = TRUE, 
	firstcol = NA, 
	NASCthr = 0.001, 
	NASCexp = 2, 
	NASCmax_size = 2, 
	biomAlpha = 0.1, 
	biomShape = 20, 
	biomSize = 2, 
	biomThr = 1, 
	trawlSize = 2, 
	stratumcol = NULL, 
	plot = c("biom", "nasc", "trawl"), 
	#alldays = TRUE, 
	...){
	
	# Run the original surveyPlanner plot:
	p <- plotStratum(x$transects, ...)
	p0 <- p
	
	# Get the biomass to plot, either as an average of the days of the survey, or at the mid day of the survey:
	x <- getBiomassOfSurvey(x)
	
	maxBiom <- max(x$biomassOfSurvey$Biom, na.rm=TRUE) * biomThr
	if(logColscale){
		BiomSeq <- seq(0, 10*log10(maxBiom), length.out=ncolour)
		BiomSeq <- 10^(BiomSeq/10)
	}
	else{
		BiomSeq <- seq(0, maxBiom, length.out=ncolour)
	}
	
	if(is.list(col) && all(c("h", "s", "v") %in% names(col))){
		col <- hsv(h=seq(col$h[1], col$h[2], l=ncolour), s=seq(col$s[1], col$s[2], l=ncolour), v=seq(col$v[1], col$v[2], l=ncolour))
	}
	
	if(length(firstcol)){
		col[1] <- firstcol
	}
	colorInterval <- findInterval(x$biomassOfSurvey$Biom, BiomSeq)
	
	if("superind" %in%  tolower(plot)){
		x <- getSuperIndOfSurvey(x)
		#if(length(x$superIndOfSurvey)==0  && length(x$superInd)){
		#	x$superIndOfSurvey <- as.data.frame(getSuperIndOfDay(x$superInd, x$midDayOfSurvey)[c("Long", "Latt")])
		#}
		#else if(length(x$superIndOfSurvey)==0){
		#	warning("Super individual data missing, and addSuperInd cannot be TRUE.")
		#}
		p <- p + geom_point(data=x$superIndOfSurvey, aes(x=Long, y=Latt), shape=2, color="black", size=0.2)
	}
	
	# Add biomass to the transect plot:
	if("biom" %in%  tolower(plot)){
		for(i in seq_len(ncolour)){
			p <- p + geom_point(data=x$biomassOfSurvey[colorInterval==i, ], aes(x=Long, y=Latt), colour=col[i], alpha=biomAlpha, shape=biomShape, size=biomSize)
		}
	}
	
	# Add NASC values larger than NASCthr of the fraction of NASC relative to maxNASC:
	maxNasc <- max(x$transects$Transect$nasc, na.rm=TRUE)
	pointSize <- (x$transects$Transect$nasc/maxNasc)^NASCexp
	
	pointSize[pointSize < NASCthr] <- NA
	x$transects$Transect <- cbind(x$transects$Transect, pointSize=pointSize)
	if("nasc" %in%  tolower(plot)){
		p <- p + geom_point(data=x$transects$Transect, aes(x=lon_mid, y=lat_mid, size=pointSize), shape=20)  +  scale_size_area(max_size=NASCmax_size, guide=FALSE)
	}
	
	# Plot the trawl stations
	if("trawl" %in%  tolower(plot)){
		p <- p + geom_point(data=x$transects$Transect[ x$transects$Transect$trawl, ], aes(x=lon_mid, y=lat_mid), shape=42, color="red", size=trawlSize)
	}
	
	if(length(stratumcol)){
		p <- p + scale_fill_manual(values = rep(stratumcol, length.out=nrow(x$transects$Stratum)))
	}
	
	print(p)
	list(p, p0)
}
#'
#' @export
#' @importFrom Rstox getPlottingUnit 
#' @rdname plotPELFOSS
#'
plotTSB <- function(x, unit="mt"){
	scale <- Rstox::getPlottingUnit(unit, var="weight")$scale
	x$totBiom$Total <- x$totBiom$Total / scale
	x$totBiom$NonEmptyStrata <- x$totBiom$NonEmptyStrata / scale
	x$totBiom$AllStrata <- x$totBiom$AllStrata / scale
	TSB <- getTSB(x, unit=unit)$bootstrap$abnd
	
	plot(x$totBiom$Total, ylim=c(0, max(x$totBiom$Total)))
	points(x$totBiom$AllStrata, col=2)
	points(x$totBiom$NonEmptyStrata, col=3)
	abline(h = TSB$Ab.Sum.mean)
	abline(h = TSB$'Ab.Sum.50%', lty=2)
	abline(v = range(x$daysOfSurvey))
}
#'
#' @export
#' @importFrom Rstox getPlottingUnit 
#' @rdname plotPELFOSS
#'
getTSB <- function(x, unit="mt", digits=4){
	
	# Remove the name column:
	x$report$bootstrap$abnd <- x$report$bootstrap$abnd[,-1]
	
	# Scale the Rstox report back to the baseunit, which is grams and the same as used in the input biomass data:
	x$report$bootstrap$abnd[names(x$report$bootstrap$abnd) != "Ab.Sum.cv"] <- x$report$bootstrap$abnd[names(x$report$bootstrap$abnd) != "Ab.Sum.cv"] * x$report$bootstrap$scale
	
	# The scale to the requested unit:
	scale <- getPlottingUnit(unit, var="weight")$scale
	x$report$bootstrap$abnd[names(x$report$bootstrap$abnd) != "Ab.Sum.cv"] <- x$report$bootstrap$abnd[names(x$report$bootstrap$abnd) != "Ab.Sum.cv"] / scale
	#ThSB_inside <- x$totBiom$NonEmptyStrata[x$midDayOfSurvey] / scale
	#ThSB_insideAll <- x$totBiom$AllStrata[x$midDayOfSurvey] / scale
	#ThSB_total <- x$totBiom$Total[x$midDayOfSurvey] / scale
	ThSB_inside <- mean(x$totBiom$NonEmptyStrata[x$daysOfSurvey]) / scale
	ThSB_insideAll <- mean(x$totBiom$AllStrata[x$daysOfSurvey]) / scale
	ThSB_total <- mean(x$totBiom$Total[x$daysOfSurvey]) / scale
	TSB_meanByInside <- x$report$bootstrap$abnd[["Ab.Sum.mean"]] / ThSB_inside
	TSB_meanByInsideAll <- x$report$bootstrap$abnd[["Ab.Sum.mean"]] / ThSB_insideAll
	TSB_meanByTotal <- x$report$bootstrap$abnd[["Ab.Sum.mean"]] / ThSB_total
	TSB_medianByInside <- x$report$bootstrap$abnd[["Ab.Sum.50%"]] / ThSB_inside
	TSB_medianByInsideAll <- x$report$bootstrap$abnd[["Ab.Sum.50%"]] / ThSB_insideAll
	TSB_medianByTotal <- x$report$bootstrap$abnd[["Ab.Sum.50%"]] / ThSB_total
	
	
	x$report$bootstrap$abnd <- cbind(
		data.frame(
			startDayOfSurvey = x$startDayOfSurvey, 
			midDayOfSurvey = x$midDayOfSurvey, 
			endDayOfSurvey = x$endDayOfSurvey, 
			startDateOfSurvey = x$startDateOfSurvey,
			endDateOfSurvey = x$endDateOfSurvey, stringsAsFactors=FALSE
		),
		x$report$bootstrap$abnd[, c("Ab.Sum.5%", "Ab.Sum.50%", "Ab.Sum.95%", "Ab.Sum.mean", "Ab.Sum.sd", "Ab.Sum.cv")], 
		data.frame(
			ThSB_inside = ThSB_inside, 
			ThSB_insideAll = ThSB_insideAll, 
			ThSB_total = ThSB_total, 
			TSB_meanByInside = TSB_meanByInside, 
			TSB_meanByInsideAll = TSB_meanByInsideAll, 
			TSB_meanByTotal = TSB_meanByTotal, 
			TSB_medianByInside = TSB_medianByInside, 
			TSB_medianByInsideAll = TSB_medianByInsideAll, 
			TSB_medianByTotal = TSB_medianByTotal, stringsAsFactors=FALSE
		)
	)
	
	# Set the significant digits:
	areNumeric <- sapply(x$report$bootstrap$abnd, is.numeric)
	x$report$bootstrap$abnd[areNumeric] <- signif(x$report$bootstrap$abnd[areNumeric], digits)
	x$report
}
#'
#' @export
#' @rdname plotPELFOSS
#' 
getAllTSB <- function(dir, outfile=NULL){
	# Get PELFOSS files:
	getPelfossFiles <- function(l, lshort, dir, type="output"){
		# Define keys to grep by:
		keys <- list(
			output = "NORWECOM_Rstox/output/r/report/biomass2tsb.RData", 
			plot = "NORWECOM_Rstox/output/r/report/biomass2tsb.png"
		)
		
		# Get the directory in which to put the files:
		allPelfossDir <- file.path(dir, paste("all", type, sep="_"))
		suppressWarnings(dir.create(allPelfossDir))
		
		# Get the file names as the concatination of the folders:
		atPelfossFiles <- grep(keys[[type]], l)
		pelfossFiles <- l[atPelfossFiles]
		# Split the file names and paste the folders for the survey, resolution, year, direction+timing and seed:
		filetag <- strsplit(lshort[atPelfossFiles], "/")
		filetag <- sapply(filetag, function(x) paste(head(x, 5), collapse="_"))
		# Add the file base names:
		pelfossFiles_new <- paste(filetag, basename(pelfossFiles), sep="_")
		# Expand to the full path:
		pelfossFiles_new <- file.path(allPelfossDir, pelfossFiles_new)
		
		# Copy the files:
		file.copy(pelfossFiles, pelfossFiles_new, overwrite=TRUE)
		
		return(list(orig=pelfossFiles, new=pelfossFiles_new, filetag=filetag))
	}
	
	splitRunName <- function(x){
		getBetween <- function(x, start, end){
			if(is.character(start)){
				start <- regexpr(start, x) + nchar(start)
			}
			if(length(end)==0){
				return(substring(x, start))
			}
			else if(is.character(end)){
				end <- regexpr(end, x) - 1
			}
			return(substr(x, start, end))
		}
		
		survey <- getBetween(x, 1, "_res_")
		res <- getBetween(x, "_res_", "_year_")
		year <- getBetween(x, "_year_", "_direction_")
		dir <- getBetween(x, "_direction_", "_timing_")
		timing <- getBetween(x, "_timing_", "_seed_")
		seed <- getBetween(x, "_seed_", NULL)
		
		
		
		#
		#
		#
		#survey <- substr(x, 1, regexpr("_res_", x) - 1)
		#res <- substr(x, regexpr("_res_", x) + 4, regexpr("_year_", x) - 1)
		#year <- substr(x, regexpr("_year_", x) + 5, regexpr("_direction_", x) - 1)
		#dir <- substr(x, regexpr("_direction_", x) + 10, regexpr("_timing_", x) - 1)
		#timing <- substr(x, regexpr("_timing_", x) + 7, regexpr("_seed_", x) - 1)
		#seed <- substring(x, regexpr("_seed_", x) + 5)
		data.frame(
			Survey = survey, 
			Resolution = res, 
			Year = year, 
			Direction = dir, 
			Timing = timing, 
			Seed = seed, stringsAsFactors=FALSE
		)
	}
	
	# List all projects:
	projectsDir <- file.path(dir, "project")
	l <- list.files(projectsDir, recursive=TRUE, full.names=TRUE)
	lshort <- list.files(projectsDir, recursive=TRUE, full.names=FALSE)
	
	# Get and copy files:
	allOutputFiles <- getPelfossFiles(l, lshort, dir, type="output")
	allPlotFiles <- getPelfossFiles(l, lshort, dir, type="plot")
	
	# Load all output files into a list:
	data <- lapply(allOutputFiles$new, function(x) mget(load(x)))
	names(data) <- allOutputFiles$filetag
	
	# Get the reports:
	reports <- do.call(rbind, lapply(data, function(x) getTSB(x$biomass2tsb.out)$bootstrap$abnd))
	#reports <- do.call(rbind, lapply(seq_along(data), function(ind) {print(ind); print(allOutputFiles$filetag[ind]); getTSB(data[[ind]]$biomass2tsb.out)$bootstrap$abnd}))
	
	reports <- cbind(splitRunName(rownames(reports)), reports, name=rownames(reports), stringsAsFactors=FALSE)
	rownames(reports) <- NULL
	
	
	reports <- reports[order(reports$Survey, reports$Year),]
	#reports <- reports[order(reports$Survey),]
	
	getStratumSysytem <- function(x){
		paste(head(unlist(strsplit(x, "_")), 2), collapse="_")
	}
	stratumSysytem <- sapply(reports$Survey, getStratumSysytem)
	surveyRegion <- paste(stratumSysytem, "res", reports$Resolution, "year", reports$Year, "direction", reports$Direction, "timing", reports$Timing, "seed", reports$Seed, sep="_")
	
	
	ThSB_inside <- c(reports$ThSB_inside[match(surveyRegion, reports$name)])
	reports$Ab.Sum.05_ThSB_inside <- reports$'Ab.Sum.5%' / ThSB_inside
	reports$Ab.Sum.mean_ThSB_inside <- reports$Ab.Sum.mean / ThSB_inside
	reports$Ab.Sum.95_ThSB_inside <- reports$'Ab.Sum.95%' / ThSB_inside
	reports$CaseNr <- seq_len(nrow(reports))
	
	
	
	


	# Write the reports to file:
	if(length(outfile)==0 || !is.character(outfile)){
		outfile <- file.path(dir, "reports", "allReports.txt")
	}
	write.table(reports, file=outfile, sep="\t", dec=".", row.names=FALSE, fileEncoding="UTF-8")
	
	
	return(list(reports=reports, data=data, allOutputFiles=allOutputFiles, allPlotFiles=allPlotFiles))
}


#*********************************************
#*********************************************
#' Simulate survey based on NORWECOM data.
#'
#' @param projectName		The biomass area density in units of g/m^2 (gram per square meter in the horizontal plane). 
#' @param xmlOutputDir		The biomass area density in units of g/m^2 (gram per square meter in the horizontal plane). 
#' @param biomass			The path to a NORWECOM NetCDF4 file with biomass in grams per square meter in a irregular geographical grid. The file should contain the variables "Biom", "Long" and "Latt".
#' @param superInd			The path to a NORWECOM NetCDF4 file with superindividuals. The file should contain the variables "gridLongInd", "gridLattInd", "gridLong", "gridLatt", "female", "age", "inumb", "length", "sweight" and "pweight".
#' @param stratum			The path to a file with the stratum polygons given as a two column matrix with stratum name in the first column and a MULTIPOLYGON wkt string in the second conlumn.
#' @param startdate			The start date of the survey, given as "%d/%m", e.g., 31 January is "31/1".
#' @param dateshift 		An integer number of days to shift the survey by, negative values shifting to earlier in the year.
#' @param centroid			The centroid of the data, given in longitude, latitude.
#' @param seed				A single integer, or a list of seeds used in the funciton, including the following seeds: ('transect') for drawing transects using surveyPlanner(), ('trawl') for drawing trawl stations along the transects with probability as a funciton of the NASC (see \code{probfact}), and ('bootstrap') for getting the final estimate using runBootstrap() (see \code{nboot}).
#' @param xsd				A named list of xsd versions of the acoustic and biotic file format.
#' @param nTrawl			The number of trawls to draw. Implies a penalty on the total transect time by \code{hoursPerTrawl}.
#' @param type,knots,equalEffort,bearing,distsep,margin	See \code{\link{surveyPlanner}}
#' @param tsn				The tsn code of the species.
#' @param m,TS0				The parameters of the target strength-length relationship TS = m * log10(Lcm) + TS0, typically m = 20 and TS0 = -71.9 (herring, from Foote, K. G. 1987. Fish target strengths for use in echo integrator surveys. Journal of the Acoustical Society of America, 82: 981– 987.)
#' @param platform			The platform to use, defaulted to G.O.Sars. Only kept for cosmetic reasons.
#' @param distance,sweepWidth	The trawled distance and the seew width of the simulated trawl.
#' @param probfact			A numeric indicating an addition in the probability of selecting a log distance for trawling. I.e., add probfact / numberOfTransects to each log distance probability NASC / sum(NASC), and normalize to toal probability = 1.
#' @param radius			The radius around the trawl station within which individuals are drawn from the superindividuals.
#' @param N					The number of individuals to draw for each trawl station.
#' @param ...				Arguments passed to \code{\link{surveyPlanner}}. 
#' 
#' @return
#'
#' @export
#' @import Rstox
#' @importFrom maptools unionSpatialPolygons
#' @rdname biomass2tsb
#'
biomass2tsb <- function(
	# For the StoX project and writing XML files:
	projectName, xmlOutputDir, 
	# For the biomass and superindividuals and other global options:
	biomass, superInd, stratum, startdate = "31/1", dateshift = 0, centroid = c(0, 68), seed = 0, nboot=10, xsd = list(acoustic="1", biotic="1.4"),
	# For transects and NASC:
	nTrawl = 50, hours = list(240), type = "RectEnclZZ", knots = 10, equalEffort = FALSE, bearing = "along", distsep = 1, margin = 0.1, tsn = 161722, m = 20, TS0 = -71.9, 
	# For trawl:
	platform = 4174, distance = 5, sweepWidth = 25, probfact = 1, radius=10, 
	N=100, cores=list(biotic=1, acoustic=1, bootstrap=1), unit="mt", lonlatFilter=NULL, 
	...){
		
	# Not used, calculate the pure transect time outside of the function, and do not consider the time used by trawling as a delay of the acoustic logging.
	# hoursPerTrawl = 2, 
	#' @param daysOfSurvey		The number of days reserved for the survey.
	
	# Save the inputs to the function:
	biomass2tsbInputs <- mget(setdiff(names(formals(biomass2tsb)), "..."))
	
	
	if(!is.list(seed) && length(seed) == 1){
		seed <- list(transect=seed, trawl=seed, bootstrap=seed)
	}
		
	# Hard coded:
	pel_ch_thickness = 100
	##### NOTE 4: Find the meaning of acocat in some file from Rolf, but it has no effec unless filtered on: #####
	acocat = 99
	# The TS-length relations require 38 kHz:
	freq = 38000
	# For now we only consider one channel, and that this is a pelagic channel:
	nChannels = 1
	channelType = "P"
	byStratum = FALSE
	keepTransport = FALSE
	missiontype = 4
	missionnumber = 1
	samplenumber = 1
	
	if(!is.list(cores)){
		cores <- as.list(rep(cores, length.out=3))
		names(cores) <- c("biotic", "acoustic", "bootstrap")
	}
	
	
	#########################################################
	##### 1. Read NORWECOM files (convert to Cartesian) #####
	#########################################################

	# Read the biomass:
	##### NOTE 2: Remove the rename parameter once the files have the requested variable names, also in the readNORWECOMsuperind(): #####
	if(is.character(biomass) && file.exists(biomass)){
		cat("Read biomass file...\n")
		biomass <- readNORWECOMbiomass(
			biomass, centroid=centroid, 
			#vars=c("HERbiom", "Long", "Latt"), 
			#rename=list(HERbiom="Biom")
			vars=c("Biom", "Long", "Latt")
		)
	}


	# Read the superindividuals:
	##### NOTE 3: Remove the rename parameter once the files have the requested variable names, also in the readNORWECOMsuperind(): #####
	if(is.character(superInd) && file.exists(superInd)){
		cat("Read superInd file...\n")
		superInd <- readNORWECOMsuperind(
			superInd, centroid=centroid, 
			vars = c("xpos", "ypos", "Long", "Latt", "female", "age", "inumb", "length", "sweight", "pweight"), 
			rename=list(xpos="gridLongInd", ypos="gridLattInd", Long="gridLong", Latt="gridLatt")
		)
	}

	#########################################################
	#########################################################


	############################################################################
	##### 2. Run surveyPlanner(), extract NASC, return acoustic data frame #####
	############################################################################

	cat("Run surveyPlanner()...\n")
	
	##### Parameters: #####
	# For the surveyPlanner():
	#timePenalty <- nTrawl * hoursPerTrawl
	#totalhours <- 24 * totalhours - timePenalty
	#hours = list(24 * daysOfSurvey)

	##### surveyPlanner: #####
	# Get the year, and allow for dates in the biomass crossing over to neighbouring years:
	year = format(as.Date(biomass$time),"%Y")
	year <- median(year)
	# Add the year to the start date:
	startdate <- paste(year, startdate, sep="/")
	startdate <- as.Date(startdate, format="%Y/%d/%m")
	t0 <- as.POSIXct(paste(startdate, "00:00:00"), format="%Y-%m-%d %H:%M:%S", tz="UTC")
	# Add the dateshift:
	t0 <- t0 + dateshift * 24 * 60 * 60

	transects <- Rstox::surveyPlanner(
		projectName = stratum, 
		type = type, 
		bearing = bearing, 
		hours = hours, 
		t0 = t0, 
		knots = knots, 
		seed = seed$transect, 
		equalEffort = equalEffort, 
		distsep = distsep, 
		margin = margin, 
		byStratum = byStratum, 
		keepTransport = keepTransport, 
		centroid = centroid, 
		...
	)
	
	
	##### Create NASC values: #####
	# Interpolate directly onto the log distances:
	dayvec <- findInterval(transects$Transect$time_start, biomass$time)
	daysOfSurvey <- unique(dayvec)
	startDayOfSurvey <- min(daysOfSurvey)
	endDayOfSurvey <- max(daysOfSurvey)
	midDayOfSurvey <- round(mean(c(startDayOfSurvey, endDayOfSurvey)))
	startDateOfSurvey <- min(transects$Transect$time_mid)
	endDateOfSurvey <- max(transects$Transect$time_mid)
	
	
	Bg <- interpolateBiomassToTransects(biomass, transect=transects$Transect, centroid=centroid, dayvec=dayvec)
	# NORWECOM biomass is given in grams
	#Wkg <- Wg * 1e-3

	# Get condition exponent from the superindividuals:
	s <- getSuperIndOfDay(superInd=superInd, day=1)
	Lcm <- s$length
	Wg <- s$sweight
	#pW <- s$pweight * 1e-3
	#plot((Lcm^3), Wg)
	b_true <- 3
	a_true <- Wg / Lcm^3
	summary(a_true)
	a_true <- median(a_true)

	# Get the typical length of the fish, by a weighted average of the superindividuals in the survey region:
	id <- double(length(transects$Input$lonlatSP))
	# Get the union of the SpatialPolygons of the strata:
	stratumUnion <- maptools::unionSpatialPolygons(transects$Input$lonlatSP, id)
	# Select the fist in case there are holes in the union:
	stratumUnionMatrix <- Rstox::getMatrixList(stratumUnion)
	if(is.list(stratumUnionMatrix)){
		stratumUnionMatrix <- stratumUnionMatrix[1]
	}
	else{
		stratumUnionMatrix <- list(stratumUnionMatrix)
	}
	
	## Get the mean length of the fish weighted by the expected backscatter:
	#startDayOfSurvey <- as.numeric(strftime(min(transects$Tra$time_mid), format = "%j"))
	#endDayOfSurvey <- as.numeric(strftime(max(transects$Tra$time_mid), format = "%j"))
	##midDayOfSurvey <- as.numeric(strftime(median(transects$Tra$time_mid), format = "%j"))
	#midDayOfSurvey <- round(mean(c(startDayOfSurvey, endDayOfSurvey)))
	
	
	LcmMean <- getTotalBiomass(superInd=superInd, stratumPolygons=stratumUnionMatrix, day=midDayOfSurvey, type="length")$AllStrata

	
	# Get the NASC values from the biomass horizontal area density:
	#nasc <- biomass2nasc(Wg, a=a_true, b=b_true, m=m, TS0=TS0)
	nasc <- biomass2nasc(Bg, LcmOne=LcmMean, a=a_true, b=b_true, m=m, TS0=TS0)
	if(all(is.na(nasc)) || !isTRUE(any(nasc>0))){
		warning("No biomass along the transects")
		return(NULL)
	}
	
	
	transects$Transect <- cbind(transects$Transect, nasc=nasc)


	# Create the acoustic data frame to write to LUF20 file:
	numt <- nrow(transects$Transect)
	allZeros <- double(numt)
	allOnes <- allZeros + 1
	transceiver <- length(freq)
	
	acoustic <- data.frame(
		#distance_list = NA, 
		#ch_type = NA, 
		cruise = transects$Transect$cruise[1], 
		integrator_dist = transects$Transect$dist_stop - transects$Transect$dist_start, 
		pel_ch_thickness = rep(pel_ch_thickness, length.out=numt),
		log_start = transects$Transect$log_start,
		lon_start = transects$Transect$lon_start,
		lon_stop = transects$Transect$lon_stop,
		lat_start = transects$Transect$lat_start,
		#start_time = transects$Transect$time_start,
		start_time = format(transects$Transect$time_start),
		num_pel_ch = allOnes,
		upper_interpret_depth = allZeros,
		upper_integrator_depth = allZeros, 
		acocat = acocat, 
		freq = freq,
		type = channelType,
		transceiver = transceiver,
		sa..ch.1 = round(nasc, digits=6), 
		nasc = nasc, stringsAsFactors=FALSE
	)
	
	############################################################################
	############################################################################


	############################################
	##### 3. Get total theoretical biomass #####
	############################################

	cat("Get total biomass...\n")
	strataInd <- match(transects$Stratum$stratum, names(transects$Input$lonlat))
	totBiom <- getTotalBiomass(superInd=superInd, stratumPolygons=transects$Input$lonlat, strataInd=strataInd, type="biomass", lonlatFilter=lonlatFilter)
	
	############################################
	############################################


	################################################################
	##### 4. Simulate trawl stations, return biotic data frame #####
	################################################################

	cat("Draw trawl stations...\n")
	
	# Draw trawl stations (get the indices of the log distances with trawl):
	trawlInd <- drawTrawlStationInd(nasc, nTrawl, probfact=probfact)
	# Add to the transect matrix:
	transects$Transect$trawl <- FALSE
	transects$Transect$trawl[trawlInd] <- TRUE
	
	# The following gives the required attributes: readHIXSD("1.4", "biotic")$attrs_required.
	# Here we include the 'mission' info in the fish station data frame:
	fishStation <- data.frame(
		missiontype = missiontype, 
		cruise = transects$Transect$cruise[1], 
		year = format(as.Date(head(transects$Transect$time_start, 1)),"%Y"),
		platform = platform,
		missionnumber = missionnumber, 
		# Ordinary fish station variables and attributes:
		serialno = seq_along(trawlInd), 
		starttime = transects$Transect$time_mid[trawlInd],
		longitudestart = transects$Transect$lon_start[trawlInd],
		latitudestart = transects$Transect$lat_start[trawlInd],
		x_mid = transects$Transect$x_mid[trawlInd],
		y_mid = transects$Transect$y_mid[trawlInd], 
		species = tsn, 
		samplenumber = samplenumber, 
		distance = distance, stringsAsFactors=FALSE
	)

	# Simulate the trawls:
	biotic <- simulateTrawl(superInd, fishStation, seed=seed$trawl, radius=radius, N=N)
	superIndCount <- biotic$superIndCount
	biotic <- biotic$biotic
	
	### Add trawl info to the transects$Transect: ###
	# Add superindividual count for the trawl stations:
	transects$Transect$superIndCount <- NA
	transects$Transect$superIndCount[trawlInd] <- superIndCount
	# Add serialno for the trawl stations:
	transects$Transect$serialno <- NA
	transects$Transect$serialno[trawlInd] <- fishStation$serialno
	
	# Discard strata with no biotic data:
	#sumBiotic <- by(transects$Transect$stratum, transects$Transect$superIndCount, sum, na.rm=TRUE)
	sumBiotic <- by(transects$Transect$superIndCount, transects$Transect$stratum, sum, na.rm=TRUE)
	hasNoBiotic <- sumBiotic == 0
	if(all(hasNoBiotic)){
		warning("All strata had no biotic data (no non-empty trawl stations")
		return(NULL)
	}
	else if(any(hasNoBiotic)){
		badStrata <- names(sumBiotic)[hasNoBiotic]
		warning(paste("The following strata had no biotic data (no non-empty trawl stations), and vere discarded in the acousic.xml file: "), paste(badStrata, collapse=", "))
		transects$Transect <- transects$Transect[! transects$Transect$stratum %in% badStrata, , drop=FALSE]
	}
	
	
	################################################################
	################################################################


	# The plotting should be revised, and possibly put inside the plotStratum(), and maybe put outside of the function:
	
	#########################
	##### 5. Visualize: #####
	#########################

	#p <- plotStratum(transects)
    #
	## Plot with the trals and NASC:
	#p1 <- plotPELFOSS(
	#	p, 
	#	transects,
	#	biomass,
	#	day=222, 
	#	trawlInd=trawlInd, 
	#	ncolour=40, 
	#	maxBiom=700, 
	#	maxNasc=4000, 
	#	#col=list(h=c(0.95, 0.98), s=c(0.35, 1), v=c(0.9, 0.5)), 
	#	NASCthr=0.001, 
	#	NASCexp=2
	#)

	#########################
	#########################


	########################
	##### 6. write xml #####
	########################

	tempfile_biotic <- file.path(xmlOutputDir, "biotic.xml")
	tempfile_acoustic <- file.path(xmlOutputDir, "acoustic.xml")

	cat("Write biotic xml...\nTime used:\n")
	temp <- system.time(Rstox::writeBioticXML(biotic, tempfile_biotic, xsd=xsd$biotic, cores=cores$biotic))
	#cat("Time used:\n")
	#print(temp)
	
	cat("Write acoustic xml...\nTime used:\n")
	temp <- system.time(Rstox::writeAcousticXML(acoustic, tempfile_acoustic, xsd=xsd$acoustic, cores=cores$acoustic))
	#cat("Time used:\n")
	#print(temp)
	
	
	########################
	########################


	#######################################
	##### 7. Generate and run project #####
	#######################################

	cat("Create and run the StoX project...\n")
	
	# Create the NORWECOM project:
	input_full <- file.path(Rstox::getProjectPaths(projectName)$projectPath, "input", c("acoustic", "biotic"))
	input <- file.path("input", c("acoustic", "biotic"))
	names(input) <- basename(input)
	names(input_full) <- basename(input)
	# define the paths to the input files:
	file_biotic <- file.path(input_full["biotic"], "biotic.xml")
	file_acoustic <- file.path(input_full["acoustic"], "acoustic.xml")


	model <- list(
		# Baseline: 
		"ReadProcessData", 
		ReadAcousticXML = list(
			FileName1 = file.path(input["acoustic"], "acoustic.xml")
		), 
		NASC = list(
			AcousticData = "ReadAcousticXML", 
			LayerType = "WaterColumn"
		), 
		ReadBioticXML = list(
			FileName1 = file.path(input["biotic"], "biotic.xml")
		), 
		StationLengthDist = list(
			BioticData = "ReadBioticXML", 
			LengthDistType = "PercentLengthDist"
		), 
		RegroupLengthDist = list(
			LengthDist = "StationLengthDist", 
			LengthInterval = "1.0"
		), 
		StratumArea = list(
			ProcessData = "ReadProcessData", 
			AreaMethod = "Accurate"
		), 
		DefineAcousticPSU = list(
			ProcessData = "ReadProcessData", 
			AcousticData = "ReadAcousticXML", 
			DefinitionMethod = "UseProcessData"
		),
		MeanNASC = list(
			ProcessData = "ReadProcessData", 
			NASC = "NASC", 
			SampleUnitType = "PSU"
		), 
		BioStationAssignment = list(
			ProcessData = "ReadProcessData", 
			BioticData = "ReadBioticXML", 
			AssignmentMethod = "UseProcessData", 
			EstLayers = "1~PELBOT"
		), 
		TotalLengthDist = list(
			ProcessData = "ReadProcessData", 
			LengthDist = "RegroupLengthDist"
		), 
		SweptAreaDensity = list(
			ProcessData = "ReadProcessData", 
			SweptAreaMethod = "LengthDependent", 
			BioticData = "ReadBioticXML", 
			LengthDist = "TotalLengthDist", 
			DistanceMethod = "FullDistance",
			sweepWidthMethod = "Constant", 
			sweepWidth = sweepWidth
		), 
		AcousticDensity = list(
			ProcessData = "ReadProcessData", 
			LengthDist = "TotalLengthDist", 
			NASC = "MeanNASC", 
			m = m, 
			a = TS0
		),  
		MeanDensity = list(
			ProcessData = "ReadProcessData", 
			Density = "AcousticDensity", 
			SampleUnitType = "Stratum"
		), 
		SumDensity = list(
			Density = "MeanDensity"
		), 
		Abundance = list(
			Density = "SumDensity", 
			PolygonArea = "StratumArea"
		), 
		IndividualDataStations = list(
			ProcessData = "ReadProcessData", 
			Abundance = "Abundance"
		), 
		IndividualData = list(
			BioticData = "ReadBioticXML", 
			IndividualDataStations = "IndividualDataStations"
		), 
		SuperIndAbundance = list(
			Abundance = "Abundance", 
			IndividualData = "IndividualData", 
			ProcessData = "ReadProcessData", 
			AbundWeightMethod = "Equal"
		), 
		"WriteProcessData", 
		# Baseline report: 
		#FillMissingData = list(
		#	SuperIndividuals = "SuperIndAbundance", 
		#	FillVariables = "ImputeByAge", 
		#	Seed = "1", 
		#	FillWeight = "Regression"
		#), 
		#EstimateByPopulationCategory = list(
		#	SuperIndividuals = "FillMissingData", 
		#	LengthInterval = "2.0", 
		#	Scale = "1000", 
		#	Dim1 = "LenGrp", 
		#	Dim2 = "age", 
		#	Dim3 = "SpecCat"
		#)
		EstimateByPopulationCategory = list(
			SuperIndividuals = "SuperIndAbundance", 
			LengthInterval = "1.0", 
			Scale = "1000", 
			Dim1 = "LenGrp", 
			Dim2 = "age", 
			Dim3 = "SpecCat"
		)
	)

	# Create the project skeleton:
	Rstox::createProject(projectName, model=model, ow=TRUE)


	# Copy the input files to the project input directory:
	file.copy(tempfile_biotic, file_biotic, overwrite=TRUE)
	file.copy(tempfile_acoustic, file_acoustic, overwrite=TRUE)



	##### Modify and run the StoX project:
	# Get the path to the project.xml file_
	projectXML <- Rstox::getProjectPaths(projectName)$projectXML

	# Write the stratum polygons, the edsupsu, the psustratum, the suassignment, and the bioticassignment to the project.xml file directly:
	insertStratumpolygon(transects$Input$lonlat, file=projectXML)

	insertEdsupsu(transects$Transect, file=projectXML)

	insertPsustratum(transects$Transect, file=projectXML)

	insertSuassignment(transects$Transect, file=projectXML)

	# Write the bioticassignment
	
	#Stations <- data.frame(
	#	stratum = transects$Stratum$stratum, 
	#	stationID = paste(biotic$cruise, biotic$serialno, sep="/"), stringsAsFactors=FALSE
	#)
	
	hasFish <- which(transects$Transect$superIndCount > 0)
	Stations <- data.frame(
		stratum = transects$Transect$stratum[hasFish], 
		stationID = paste(transects$Transect$cruise[hasFish], transects$Transect$serialno[hasFish], sep="/"), stringsAsFactors=FALSE
	)
	insertBioticassignment(Stations, file=projectXML)


	# Reopen and save the project. This is required to order the process data properly:
	Rstox::reopenProject(projectName)
	Rstox::saveProject(projectName)
	
	# Get the baseline:
	#g <- runBaseline(projectName, exportCSV=TRUE)

	# Bootstrap:
	boot <- Rstox::runBootstrap(projectName, nboot=nboot, seed=seed$bootstrap, bootstrapMethod="AcousticTrawl", cores=cores$bootstrap, ...)

	# Get the total weight with precision:
	#reports <- getReports(projectName)
	#reports <- getReports(projectName, grp1=NULL)
	#reports <- getReports(projectName, var="Weight")
	report <- Rstox::getReports(projectName, var="Weight", unit="mt", grp1=NULL)
	
	# Save project data and close the project:
	Rstox::saveProjectData(projectName)
	Rstox::closeProject(projectName)
	
	#   TSB Ab.Sum.5% Ab.Sum.50% Ab.Sum.95% Ab.Sum.mean Ab.Sum.sd Ab.Sum.cv
	# 1 TSB   4124917    6376444    9903694     6864052   2228430 0.3246522
	#######################################
	#######################################
	

	# Output a list of relevant objects:
	out <- list(report=report, totBiom=totBiom, superInd=superInd, transects=transects, biomass=biomass, t0=t0, daysOfSurvey=daysOfSurvey, startDayOfSurvey=startDayOfSurvey, midDayOfSurvey=midDayOfSurvey, endDayOfSurvey=endDayOfSurvey, startDateOfSurvey=startDateOfSurvey, endDateOfSurvey=endDateOfSurvey, projectName=projectName)
	
	# Add total (estimated) stock biomass (TSB) and theoretical stock biomass (ThSB):
	#out$report <- getTSB(out, unit=unit)
	
	# Add function inputs:
	out$input <- biomass2tsbInputs
	
	return(out)
}



#*********************************************
#*********************************************
#' Run the PELFOSS framework.
#'
#' This function runs the entire PELFOSS framework from the NORWECOM model output files (biomass in a grid and superindividuals), through surveyPlanner() and StoX to the estimated total stock biomass compared to the corresponding theoretical values:
#'
#' @param dir			The directory holding the NORWECOM files, the stratum files and the output xml files from the function (temporary stored for insertion into the StoX project).
#' @param survey		The name of the survey to simulate, one of "Herring_IESNS" (International ecosystem survey in the Nordic Seas), "Herring_NASSHS" (Norwegian acoustic spring spawning herring survey) or "Mackerel_IESSNS" (The International Ecosystem Summer Survey in the Nordic Seas)
#' @param year,res		The year and resolution (in km) to run. The NORWECOM files are located in the following folder structure: species / resolution / year, where the resolution folders has names such as "res_4km".
#' @param seed				A single integer, or a list of seeds used in the funciton, including the following seeds: ('transect') for drawing transects using surveyPlanner(), ('trawl') for drawing trawl stations along the transects with probability as a funciton of the NASC (see \code{probfact}), and ('bootstrap') for getting the final estimate using runBootstrap() (see \code{nboot}).
#' @param type			The type of survey. See \code{\link{surveyPlanner}}
#' @param reversed		Logical: If TRUE run the survey in the opposite direciton.
#' @param projectName   The name or full path of the StoX project of the simulated survey.
#' @param nboot			Number of bootstrap replicates.
#' @param xsd			A named list of xsd versions of the acoustic and biotic file format.
#' @param ...			Additional inputs overriding the defaults returned by pelfossDefaults().
#' 
#' @return
#'
#' @export
#' @importFrom Rstox saveProjectData closeProject getProjectPaths
#' @rdname runPelfoss
#' 
runPelfoss <- function(dir, surveyPar, year=2010, res=4, dateshift=0, reversed=FALSE, seed=0, projectName="NORWECOM_Rstox", nboot=10, xsd=list(acoustic="1", biotic="1.4"), format="png", ...){
	
	# Define the path to the temporary XML files (saved here and then copied to the project for safety);
	xmlOutputDir <- file.path(dir, "XMLfiles")
	input <- list(dir=dir, xmlOutputDir=xmlOutputDir, year=year, res=res, dateshift=dateshift, reversed=reversed, seed=seed, projectName=projectName, nboot=nboot, xsd=xsd, format=format)
	
	# Parameters not listed explicitely in the function arguments, but which can be changed through "...":
	defaults <- pelfossDefaults()
	# Replace by parameters given in "...":
	lll <- list(...)
	common <- intersect(names(defaults), names(lll))
	defaults[common] <- lll[common]
	
	files <- getNorwecomPaths(dir, survey=surveyPar$survey, year=year, res=res, reversed=reversed, dateshift=dateshift, seed=seed)
	
	#surveyPar <- setSurveyParameters(survey=survey, reversed=reversed)
	
	# Reverse the survey by reversing the stratum order and orientations:
	if(reversed){
		surveyPar$strata <- rev(surveyPar$strata)
		surveyPar$rev <- rev(!surveyPar$rev)
		if(!is.list(surveyPar$hours)){
			surveyPar$hours <- rev(surveyPar$hours)
		}
	}

	# Run the function taking model output through surveyPlanner() and StoX to the TSB:
	inputs <- c(input, defaults, files, surveyPar)
	#validnames <- names(formals(biomass2tsb))
	#inputs <- inputs[names(inputs) %in% validnames]
	
	# Run the biomass2tsb:
	x <- do.call("biomass2tsb", inputs)
	if(length(x)==0){
		return(NULL)
	}
	
	# Get the paths of the project:
	paths <- Rstox::getProjectPaths(projectName)
	
	# Save the output from biomass2tsb (except from "superInd" and "biomass" but adding "biomassOfSurvey" and "superIndOfSurvey") to the output folder of the project:
	x <- getBiomassOfSurvey(x)
	x <- getSuperIndOfSurvey(x)
	x$files <- c(files, list(year=year, res=res, reversed=reversed, dateshift=dateshift, seed=seed))
	biomass2tsb.out <- x[!names(x) %in% c("superInd", "biomass")]
	biomass2tsb.outfile <- file.path(paths$RReportDir, "biomass2tsb.RData")
	save("biomass2tsb.out", file=biomass2tsb.outfile)
	
	# Save the plot to the output folder of the project:
	# Run the plot:
	p <- plotPELFOSS(x)
	# Get the file names:
	filenamebase <- file.path(paths$RReportDir, c("biomass2tsb", "surveyPlanner"))
	filename <- paste(filenamebase, format, sep=".")
	
	lll <- list(...)
	if(!all(c("width", "height") %in% names(lll))){
		lll$width <- 5000
		lll$height <- 3000
		lll$res <- 500
	}
	
	# Run the ggplots returned by plotPELFOSS:
	for(i in seq_along(p)){
		if(length(format)){
			do.call(format, c(list(filename=filename[i]), Rstox::applyParlist(lll, format)))
			#moveToTrash(filename)
		}
		print(p[[i]])
		if(length(format)){
			dev.off()
		}
	}
	
	# Copy the project to the PELFOSS directory:
	if(!file.exists(files$projectPath)){
		dir.create(files$projectPath, recursive=TRUE)
	}
	cat("Copying project ", paths$projectPath, " to ", files$projectPath, ".\n", sep="")
	file.copy(paths$projectPath, files$projectPath, recursive=TRUE)
	
	list(files$projectPath)
}
#'
#' @export
#' @keywords internal
#' @rdname runPelfoss
#' 
pelfossDefaults <- function(){
	defaults <- list(
		margin = 0.1, # Related to equalEffort
		bearing = "along",
		distsep = 1, # One nmi log distance is standard
		probfact = 1, # This reserves half of the probability of drawing trawls to the NASC and half to chance
		radius = 10, # Increase this to include more superindividuals in the trawl sample
		N = 100, # Draw 100 fish per trawl
		cores = 6, # Use 6 cores for both biotic and acoustic xml and bootstrap
		unit = "mt",
		platform = 4174, # G.O.Sars
		distance = 5, # Trawling distance, equal for all stations, thus ineffective in an acousic-trawl survey
		sweepWidth = 25 # Seewp width, see distance above
	)
	defaults
}
#'
#' @export
#' @keywords internal
#' @rdname runPelfoss
#' 
getNorwecomPaths <- function(dir, survey="Herring_IESNS", year=2010, res=4, reversed=FALSE, dateshift=0, seed=0){
	# Get the NORWECOM file paths:
	
	# Extract species name from the survey:
	species <- strsplit(survey, "_", fixed=TRUE)[[1]][1]
	# List the files and identify the file type (one of biomass and superInd):
	data <- list.files(file.path(dir, "data", species, paste0("res_", res, "km"), paste0("year_", year)), full.names=TRUE)
	type <- sapply(strsplit(basename(data), "_", fixed=TRUE), head, 1)
	
	# Get the biomass and superInd files:
	biomass <- data[tolower(type) == "biomass"]
	superInd <- data[tolower(type) == "superind"]
	# Get stratum file:
	stratum <- list.files(file.path(dir, "stratum", survey), full.names=TRUE)[1]
	
	# Get the directory in which to copy the StoX project containing all files and plots:
	direction <- c("Normal", "Reversed")[reversed + 1]
	projectPath <- file.path(dir, "project", survey, paste0("res_", res, "km"), paste0("year_", year), paste0("direction_", direction, "_timing_", dateshift, "days"), paste0("seed_", seed))
	
	# Get fishery file:
	fishery <- list.files(file.path(dir, "fishery", species), full.names=TRUE)[1]
	
	# Return the paths:
	list(biomass=biomass, superInd=superInd, stratum=stratum, species=species, projectPath=projectPath, survey=survey, fishery=fishery)
}
