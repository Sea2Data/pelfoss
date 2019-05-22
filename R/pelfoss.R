#*********************************************
#*********************************************
#' Run the PELFOSS framework.
#'
#' This function runs the entire PELFOSS framework from the NORWECOM model output files (biomass in a grid and superindividuals), through surveyPlaner() and StoX to the estimated total stock biomass compared to the corresponding theoretical values:
#'
#' @param dir				The directory holding the NORWECOM files, the stratum files and the output xml files from the function (temporary stored for insertion into the StoX project).
#' @param run				The name of the run, which could be used to reflect the input parameters, such as "Run1_parallel" or "Run1_zigzag".
#' @param biomass,superind	The biomass and superind data read from \code{readNcVarsAndAttrs} (used in \code{readNorwecomBiomass} and \code{readNorwecomSuperind}).
#' @param surveyPar			A list of parameters to be used as input to biomass2tsb (see ?biomass2tsb).
#' @param year,res			The year and resolution (in km) to run. The NORWECOM files are located in the following folder structure: species / resolution / year, where the resolution folders has names such as "res_4km".
#' @param dateshift 		An integer number of days to shift the survey by, negative values shifting to earlier in the year.
#' @param reversed			Logical: If TRUE run the survey in the opposite direciton.
#' @param seed				A single integer, or a list of seeds used in the funciton, including the following seeds: ('transect') for drawing transects using surveyPlanner(), ('trawl') for drawing trawl stations along the transects with probability as a funciton of the NASC (see \code{probfact}), and ('bootstrap') for getting the final estimate using runBootstrap() (see \code{nboot}).
#' @param projectName   	The name or full path of the StoX project of the simulated survey.
#' @param nboot				Number of bootstrap replicates.
#' @param xsd				A named list of xsd versions of the acoustic and biotic file format.
#' @param format			The file format of the saved plot, given as a string naming the function to use for saving the plot (such as bmp, jpeg, png, tiff), with \code{filename} as its first argument. Arguments fo the functions are given as \code{...}. Dimensions are defaulted to width=5000, height=3000, resolution to 500 dpi. If \code{format} has length 0, the plot is shown in the graphics window, and not saved to file.
#' @param ...				Additional inputs overriding the defaults returned by pelfossDefaults().
#' 
#' @export
#' @importFrom Rstox saveProjectData closeProject getProjectPaths
#' @importFrom grDevices dev.off
#' @rdname runPelfoss
#' 
runPelfoss <- function(dir, run="Run1", biomass=NULL, superind=NULL, surveyPar=list(), year=2010, res=4, dateshift=0, reversed=FALSE, seed=0, projectName="PelfossProjectRstox", nboot=10, xsd=list(acoustic="1", biotic="1.4"), format="png", ...){

	dirs <- getPelfossSkeleton(dir, run=run)
	createPelfossSkeleton(dir, run=run)

	# Define a list of inputs to the biomass2tsb() function;
	input <- list(dir=dir, xmlOutputDir=dirs$XMLfilesdir, year=year, res=res, dateshift=dateshift, reversed=reversed, seed=seed, projectName=projectName, nboot=nboot, xsd=xsd, format=format)

	# Parameters not listed explicitely in the function arguments, but which can be changed through "...":
	defaults <- pelfossDefaults()
	# Replace by parameters given in "...":
	lll <- list(...)
	common <- intersect(names(defaults), names(lll))
	defaults[common] <- lll[common]

	files <- getPelfossPaths(dirs=dirs, survey=surveyPar$survey, year=year, res=res, reversed=reversed, dateshift=dateshift, seed=seed)

	#surveyPar <- setSurveyParameters(survey=survey, reversed=reversed)

	# Reverse the survey by reversing the stratum order and orientations:
	if(reversed){
		surveyPar$strata <- rev(surveyPar$strata)
		surveyPar$rev <- rev(!surveyPar$rev)
		if(!is.list(surveyPar$hours)){
			surveyPar$hours <- rev(surveyPar$hours)
		}
	}

	# Merge all inputs, starting with 
	inputs <- c(input, defaults, surveyPar, files)
	# Add the biomass and superind data (lists) if given:
	if(length(biomass)){
		inputs$biomass <- biomass
	}
	if(length(superind)){
		inputs$superind <- superind
	}
	
	#validnames <- names(formals(biomass2tsb))
	#inputs <- inputs[names(inputs) %in% validnames]

	# Run the function taking model output through surveyPlanner() and StoX to the TSB:
	x <- do.call("biomass2tsb", inputs)
	if(length(x)==0){
		return(NULL)
	}

	# Get the paths of the project:
	paths <- Rstox::getProjectPaths(projectName)

	# Save the output from biomass2tsb (except from "superind" and "biomass" but adding "biomassOfSurvey" and "superindOfSurvey") to the output folder of the project:
	x <- getBiomassOfSurvey(x)
	#x <- getSuperindOfSurvey(x)
	x$files <- c(files, list(year=year, res=res, reversed=reversed, dateshift=dateshift, seed=seed))
	#biomass2tsb.out <- x[!names(x) %in% c("superind", "biomass")]
	#biomass2tsb.outfile <- file.path(paths$RReportDir, "biomass2tsb.RData")
	#save("biomass2tsb.out", file=biomass2tsb.outfile)
	pelfossOutput <- x[!names(x) %in% c("superind", "biomass")]
	pelfossOutput.file <- file.path(paths$RReportDir, "pelfossOutput.RData")
	save("pelfossOutput", file=pelfossOutput.file)
	
	# Save the plot to the output folder of the project:
	# Run the plot:
	p <- plotPelfoss(x)
	# Get the file names:
	filenamebase <- file.path(paths$RReportDir, c("pelfossOutput", "surveyPlanner"))
	filename <- paste(filenamebase, format, sep=".")

	lll <- list(...)
	if(!all(c("width", "height") %in% names(lll))){
		lll$width <- 5000
		lll$height <- 3000
		lll$res <- 500
	}

	# Run the ggplots returned by plotPelfoss:
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
	file.copy(paths$projectPath, dirname(files$projectPath), recursive=TRUE)

	list(files$projectPath)
}


#*********************************************
#*********************************************
#' Read a NORWECOM file.
#'
#' \code{readNorwecomBiomass} Reads a NORWECOM biomass file, with biomass on a irregular grid Long, Latt. \cr \cr
#' \code{readNorwecomSuperind} Reads a NORWECOM super individual file. \cr \cr
#' \code{readNcVarsAndAttrs} Extracts variables with long names as attributes. \cr \cr
#' \code{getBiomassXY} Converts to Cartesian coordinates for the biomass data. \cr \cr
#' \code{getSuperindXY} Converts to Cartesian coordinates for the superindividual data. \cr \cr
#' \code{interpolateBiomassToTransects} Interpolate NORWECOM biomass onto log distances. \cr \cr
#'
#' @param ncfile							The input NetCDF4 file with a (possibly irregular) grid of biomass in gram per square meters, and geographical positions (see the input \code{vars}).
#' @param vars								A string vector of the variables to read from the \code{ncfile}.
#' @param rename							A list of variableNameInNCfile = newVariableName, renaming the variables read from the \code{ncfile} to the variables required by the function \code{biomass2tsb} ("Biom", "Long", "Latt", for output from \code{readNorwecomBiomass}, and "gridLongInd", "gridLattInd", "gridLong", "gridLatt", "age", "inumb", "length", "sweight", "pweight", for output from \code{readNorwecomSuperind}).
#' @param centroid							The geographical position in decimal degrees (longitude, latitude) used as centoid in the conversion to Cartesian positions.
#' @param na.rm								Logical: If TRUE, remove superindividuals which are initially NA (thus born within the time period of the file). This 
#' @param maxValue							The value above which variables are identified as NA.
#' @param x									The output from \code{readNcVarsAndAttrs}.
#' @param biomass,superind					The biomass and superind data read from \code{readNcVarsAndAttrs} (used in \code{readNorwecomBiomass} and \code{readNorwecomSuperind}).
#' @param proj,units,x_0,y_0,ellps,datum	The proj4 parameters used in the projection from geographical to Cartesian coordinates centered at \code{centroid}.
#' @param transects							The output from \code{\link{surveyPlanner}}.
#' @param logDays							A vector of the days in the biomass data at which the log distances in \code{transects} are timed.
#' @param lonlatmargin						The margin to expand the ranges of longitude and latitude when interpolating each stratum. A large value demands more CPU time and memory.
#' @param superindFilter					An optional function acting on the superind data, such as longitude or latitude, e.g., superindFilter <- function(x) x$longitude < 20.
#' @param m,TS0								The parameters of the target strength-length relationship TS = m * log10(Lcm) + TS0, typically m = 20 and TS0 = -71.9 (herring, from Foote, K. G. 1987. Fish target strengths for use in echo integrator surveys. Journal of the Acoustical Society of America, 82: 981 - 987.)
#'
#' @importFrom ncdf4 nc_open ncvar_get ncatt_get
#' @rdname readBiomass
#' @export
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
#' @rdname readBiomass
#' @export
#' 
readNorwecomBiomass <- function(
	ncfile, centroid, 
	vars = c("Biom", "Long", "Latt"), 
	rename = list(Biom="biomass", Long="longitude", Latt="latitude")
	){
	
	# Read the NORWECOM superindividual file:
	out <- readNcVarsAndAttrs(ncfile, vars=vars, rename=rename)
	# Add time, using the dimension information in the NetCDF4-file:
	out <- addTimeFromHoursSince1950(out, timedim="T")
	
	# Use only the data object, and discard the nc object:
	out <- out$data
	
	# Get Cartesian coordinates:
	out <- projLonLat2XY(out, centroid=centroid)
	
	return(out)
}
#'
#' @rdname readBiomass
#' @export
#' 
readNorwecomSuperind <- function(
	ncfile, centroid, 
	vars = c("xpos", "ypos", "Long", "Latt", "female", "age", "length", "pweight", "inumb", "site", "stage"), 
	rename = list(pweight="weight"), 
	#rename = list(xpos="gridLongInd", ypos="gridLattInd", Long="gridLong", Latt="gridLatt"), 
	na.rm = TRUE, 
	maxValue = 1e10, 
	excludeSite4 = TRUE, 
	stagesToRemove = c(0, 1, 2)
	){
	
	# Function for setting large values to NA:
	setLargeToNA <- function(x, maxValue){
		x[x > maxValue] <- NA
		x
	}
	# Get the variables with days:
	getVarsWithDays <- function(x){
		withDays <- which(sapply(x, function(y) identical(ncol(y), length(x$time))))
		withDays
	}
	# Funciton for removing super individuals with intital NA at the first available variable:
	removeIntitialNA <- function(superind, maxValue){
		# Get the indices of variables with days at the second dimension:
		withDays <- getVarsWithDays(superind)
	
		# Get the indices of superindividuals with valid data of the given day:
		valid <- which(superind[[withDays[1]]][, 1] < maxValue)
		superind[withDays] <- lapply(superind[withDays], "[", valid, )
		return(superind)
	}
	# Function for replacing all data for fish with site==4 by NA:
	insertNAAtSite4 <- function(superind){
		# Get the indices of variables with days at the second dimension:
		withDays <- getVarsWithDays(superind)
	
		# Get the indices of site = 4, and replace by NAs:
		atSite4 <- superind$site == 4
		superind[withDays] <- lapply(superind[withDays], replace, atSite4, NA)
		return(superind)
	}
	# Funciton for removing super individuals with intital NA at the first available variable:
	removeStages <- function(superind, stagesToRemove=c(0, 1, 2)){
		# Get the indices of variables with days at the second dimension:
		withDays <- getVarsWithDays(superind)
		
		# Get the indices of superindividuals with valid data of the given day:
		young <- superind$stage %in% stagesToRemove
		superind[withDays] <- lapply(superind[withDays], replace, young, NA)
		return(superind)
	}
	
	
	# Read the NORWECOM superindividual file:
	out <- readNcVarsAndAttrs(ncfile, vars=vars, rename=rename)
	# Add time, using the dimension information in the NetCDF4-file:
	out <- addTimeFromHoursSince1950(out, timedim="time")
	
	# Use only the data object, and discard the nc object:
	out <- out$data
	
	# Convert too large values to NA (as intended):
	out <- lapply(out, setLargeToNA, maxValue=maxValue)
	
	##### Filters: #####
	# Remove all superindividuals which are not valid at the start of the year (i.e., discard eggs and larva, which may have different definitions of age and other variables than the juveniles and adults):
	if(na.rm){
		out <- removeIntitialNA(out, maxValue=maxValue)
	}
	
	# Remove eggs, larvae and juveniles:
	if(length(stagesToRemove)){
		out <- removeStages(out, stagesToRemove=stagesToRemove)
	}
	
	
	# Change added on 2018-11-27 after comment from Erik Askov Moussing, who informed that super individuals with site=4 are in limbo, thus no longer existing in the model, so these should be removed:
	if(excludeSite4){
		out <- insertNAAtSite4(out)
	}
	
	# Convert from the NORWECOM specified locations to actual longitude-latitude values:
	out <- getNorwecomSuperindLonLat(out)
	
	# Get Cartesian coordinates:
	out <- projLonLat2XY(out, centroid=centroid)
	
	# Return:
	out
}
#'
#' @importFrom fields interp.surface
#' @rdname readBiomass
#' @export
#' 
getNorwecomSuperindLonLat <- function(superind){
	# The NORWECOM superind data are given with location as partial values in an irregular grid of longitude and latitude. Here we convert to actual longitude and latitude values:
	# Get the dimensions of the irregular geographical grid:
	dimxy <- dim(superind$Long)
	# Get the sequence spanning the dimensions of the irregular grid (i.e., indices of the grid, which are the units of the positions xpos and ypos):
	seqx <- seq_len(dimxy[1])
	seqy <- seq_len(dimxy[2])
	# Define two lists of data, where x and y are the sequences of indices defined above in both lists, and the z is the actual longitude and latitude values respectively of the irregular grid:
	dataLon <- list(x=seqx, y=seqy, z=superind$Long)
	dataLat <- list(x=seqx, y=seqy, z=superind$Latt)

	# Define the matrix of points given as partial indices in the irregular grid:
	superindpos <- cbind(c(superind$xpos), c(superind$ypos))

	# Interpolate each of longitude and latitude onto the superindividual positions:
	superind$longitude <- fields::interp.surface(dataLon, superindpos)
	superind$latitude <- fields::interp.surface(dataLat, superindpos)

	# Reset dimensions:
	dimsuperind <- dim(superind$xpos)
	dim(superind$longitude) <- dimsuperind
	dim(superind$latitude) <- dimsuperind
	
	# Remove the obsolete variables:
	superind$Long <- NULL
	superind$Latt <- NULL
	superind$xpos <- NULL
	superind$ypos <- NULL
	
	return(superind)
}
#'
#' @importFrom Rstox geo2xy
#' @rdname readBiomass
#' @export
#' 
projLonLat2XY <- function(x, centroid, proj="aeqd", units="kmi", x_0=0, y_0=0, ellps="WGS84", datum="WGS84"){
	# Convert the locations to x,y using the 'aeqd' projection and centered at 0, 68:
	thisLonLat <- data.frame(lon=c(x$lon), lat=c(x$lat))
	xy_new <- Rstox::geo2xy(thisLonLat, list(proj=proj, units=units, lon_0=centroid[1], lat_0=centroid[2], x_0=x_0, y_0=y_0, ellps=ellps, datum=datum), data.frame.out=TRUE)

	# Define the grid to interpolate onto:
	x$x <- xy_new$x
	x$y <- xy_new$y
	dim(x$x) <- dim(x$lon)
	dim(x$y) <- dim(x$lat)
	
	return(x)
}
#'
#' @importFrom interp interp
#' @rdname readBiomass
#' @export
#' 
interpolateBiomassToTransects <- function(biomass, superind, transects, centroid, logDays, lonlatmargin=0.05, superindFilter=NULL, m=20, TS0=-71.9, insideRange=NULL, condPar=NULL, LcmMean=NULL){
	# Function to add margin to the ranges from which the points in the original locations are selected:
	addMargin <- function(x, lonlatmargin=0.05){
		rangex <- range(x)
		rangex <- rangex + c(-1, 1) * lonlatmargin[1] * diff(rangex)
		rangex
	}
	
	# Function for interpolating the biomass values onto the log distances for one specific day 'superindDay', which is looked for in the vector 'logDays', which is the days for each log distance:
	interpolateOneDayOld <- function(superindDay, biomass, superind, transects, logDays, lonlatmargin, stratumMatrix, superindFilter=NULL, m=20, TS0=-71.9, thr=10, union=FALSE, insideRange=NULL){
		
		# Get the locations of the requested day, and convert NA to 0:
		z <- c(biomass$biomass[,,superindDay])
		z[is.na(z)] <- 0
		
		# Interpolate only based on the points in proximety to the output points:
		inDay <- which(logDays==superindDay)
		xo <- transects$Transect$x_mid[inDay]
		yo <- transects$Transect$y_mid[inDay]
		
		valid <- NULL
		while(length(valid) < thr){
			rangex <- addMargin(xo, lonlatmargin)
			rangey <- addMargin(yo, lonlatmargin)
			valid <- which(c(biomass$x) >= rangex[1] & c(biomass$x) <= rangex[2] & c(biomass$y) >= rangey[1] & c(biomass$y) <= rangey[2])
			lonlatmargin <- lonlatmargin * 4
		}
		
		#rangex <- addMargin(xo, lonlatmargin)
		#rangey <- addMargin(yo, lonlatmargin)
		#
		#valid <- c(biomass$X) >= rangex[1] & c(biomass$X) <= rangex[2] & c(biomass$Y) >= rangey[1] & c(biomass$Y) <= rangey[2]
		#if(sum(valid)==0){
		#	warning("")
		#	return(NULL)
		#}
		
	
		
		# Interpolate onto the output grid:
		#zo <- akima::interp(x=c(coords$x), y=c(coords$y), z=z, xo=xo, yo=yo)
		x <- c(biomass$x)[valid]
		y <- c(biomass$y)[valid]
		z <- z[valid]
		zo <- interp::interp(x=x, y=y, z=z, xo=xo, yo=yo, output="points")
		
		##### Add the NASC, by first fitting the fish length-weight relationship to the superindividual data, and then using this to convert from biomass to NASC: #####
		#  Get first the parameters a and b in the length-weight relationship:
		transectsXyOfDay <- cbind(xo, yo)
		if(length(condPar) == 0){
			#condPar <- fitLengthWeightSuperind(superindDay, superind, superindFilter=superindFilter, b=NULL, transectsXyOfDay=transectsXyOfDay, insideRange=insideRange)$fit
			condPar <- fitLengthWeightSuperind(superindDay, superind, superindFilter=superindFilter, stratumPolygons=stratumMatrix, b=NULL, transectsXyOfDay=transectsXyOfDay, insideRange=insideRange)$fit
		}
		
		# Get the mean of the square length weighted by the number of indivudyals per superindividual:
		temp <- getTotalBiomass(superind=superind, stratumPolygons=stratumMatrix, days=superindDay, type="length", union=union, transectsXyOfDay=transectsXyOfDay, insideRange=insideRange)
		if(length(LcmMean) == 0){
			LcmMean <- temp$AllStrata
		}
		
		# Convert to NASC
		temp <- biomass2nasc(Wg=zo$z, LcmOne=LcmMean, a=condPar[1], b=condPar[2], m=m, TS0=TS0)
	
		if(all(is.na(temp$NASC)) || !isTRUE(any(temp$NASC>0))){
			warning("No biomass along the transects")
			# return(NULL)
		}
		
		# Combine to get a link to the day (suppress a warning with row names):
		suppressWarnings(out <- data.frame(superindDay=superindDay, biomassG=zo$z, LcmMean=LcmMean, condParFactor=condPar[1], condParExponent=condPar[2]))
		out <- cbind(out, temp)
		return(out)
	}
	
	#lonlatmargin <- lonlatmargin * c(diff(range(biomass$X)), diff(range(biomass$Y)))
	
	# The following stratum union matrix is used whcn calculating the mean length to use in the conversion from biomass in grams to NASC:
	#id <- double(length(transects$Input$lonlatSP))
	# Get the union of the SpatialPolygons of the strata:
	#stratumUnion <- maptools::unionSpatialPolygons(transects$Input$lonlatSP, id)
	## Select the fist in case there are holes in the union:
	#stratumUnionMatrix <- Rstox::getMatrixList(stratumUnion)
	#if(is.list(stratumUnionMatrix)){
	#	# This was an attempt to select only the surrounding polygon, but the following matrices might not be holes after all:
	#	stratumUnionMatrix <- stratumUnionMatrix[1]
	#}
	#else{
	#	stratumUnionMatrix <- list(stratumUnionMatrix)
	#}
	stratumMatrix <- transects$Input$lonlat
	
	# Run through the days and interpolate the biomass onto the transect lines, and convert to NASC:
	udays <- unique(logDays)
	#bio <- do.call(rbind, lapply(udays, interpolateOneDayOld, biomass=biomass, superind=superind, transects=transects, logDays=logDays, lonlatmargin=lonlatmargin, stratumUnionMatrix=stratumUnionMatrix, superindFilter=superindFilter, m=m, TS0=TS0))
	bio <- do.call(rbind, lapply(udays, interpolateOneDayOld, biomass=biomass, superind=superind, transects=transects, logDays=logDays, lonlatmargin=lonlatmargin, stratumMatrix=stratumMatrix, superindFilter=superindFilter, m=m, TS0=TS0, union=TRUE, insideRange=insideRange))
	
	# Order by the days:
	bio <- bio[order(bio$superindDay), ]
	# Add the biomass in grams and the NASC, as well as the a and b in the length-weight relationship:
	#transects$Transect <- cbind(transects$Transect, bio[c("biomassG", "NASC", "LcmMean", "condParFactor", "condParExponent")])
	transects$Transect <- cbind(transects$Transect, bio[c("biomassG", "NASC", "LcmMean", "condParFactor", "condParExponent")])
	
	return(transects)

	### # Fill in the biomass values at the indices with days present in the interpolation:
	### out <- double(length(logDays))
	### out[logDays %in% bio[, 1]] <- bio[,2]
	### out
	### #bio[order(bio[,1]), 2]
}


# Small function for adding time which is given as hours from 1950 in the NORWECOM model:
addTimeFromHoursSince1950 <- function(x, timedim="time"){
	# Add time from the dimensions:
	x$data$time <- as.POSIXct("1950-01-01", tz = "UTC") + x$nc$dim[[timedim]]$vals * 60^2
	return(x)
}



#*********************************************
#*********************************************
#' Simulate survey based on NORWECOM data.
#'
#' \code{biomass2tsb} Reads biomass and superind files, draws transects, converts biomass to NASC, draws trawl stations, creates and runs a StoX project, and saves the results. \cr \cr
#' \code{biomass2nasc} Funciton for converting from a biomass density to nautical area scattering coefficient NASC (linear values, not dB) used in \code{biomass2tsb}. \cr \cr
#'
#' @param projectName			The biomass area density in units of g/m^2 (gram per square meter in the horizontal plane). 
#' @param xmlOutputDir			The biomass area density in units of g/m^2 (gram per square meter in the horizontal plane). 
#' @param biomass				The path to a NORWECOM NetCDF4 file with biomass in grams per square meter in a irregular geographical grid. The file should contain the variables "biomass", "longitude" and "latitude".
#' @param superind				The path to a NORWECOM NetCDF4 file with superindividuals. The file should contain the variables "age", "inumb", "length", "weight", "longitude" and "latitude".
#' @param stratum				The path to a file with the stratum polygons given as a two column matrix with stratum name in the first column and a MULTIPOLYGON wkt string in the second conlumn.
#' @param startdate				The start date of the survey, given as "\%d/\%m", e.g., 31 January is "31/1".
#' @param dateshift 			An integer number of days to shift the survey by, negative values shifting to earlier in the year.
#' @param centroid				The centroid of the data, given in longitude, latitude.
#' @param seed					A single integer, or a list of seeds used in the funciton, including the following seeds: ('transect') for drawing transects using surveyPlanner(), ('trawl') for drawing trawl stations along the transects with probability as a funciton of the NASC (see \code{probfact}), and ('bootstrap') for getting the final estimate using runBootstrap() (see \code{nboot}).
#' @param nboot				Number of bootstrap replicates.
#' @param xsd					A named list of xsd versions of the acoustic and biotic file format.
#' @param nTrawl				The number of trawls to draw. Implies a penalty on the total transect time by \code{hoursPerTrawl}.
#' @param hours,type,knots,equalEffort,bearing,distsep,margin	See \code{\link{surveyPlanner}}
#' @param tsn					The tsn code of the species. Default is 161722 (herring). Run taxa <- getNMDinfo("taxa") to get a data frame of species used by the IMR (takes 2 min to download and process).
#' @param m,TS0					The parameters of the target strength-length relationship TS = m * log10(Lcm) + TS0, typically m = 20 and TS0 = -71.9 (herring, from Foote, K. G. 1987. Fish target strengths for use in echo integrator surveys. Journal of the Acoustical Society of America, 82: 981 - 987.)
#' @param condPar				A two element vector representing the condition parameters a and b in the standard allomettric formula W = a * L^b.
#' @param LcmMean				The characteristic fish length representing the mean acoustic backscatter, typically calculated as sqrt(mean(L^2)).
#' @param platform				The platform to use, defaulted to G.O.Sars. Only kept for cosmetic reasons.
#' @param distance,sweepWidth	The trawled distance and the seew width of the simulated trawl.
#' @param probfact				A numeric indicating an addition in the probability of selecting a log distance for trawling. I.e., add probfact / numberOfTransects to each log distance probability NASC / sum(NASC), and normalize to toal probability = 1.
#' @param radius				The radius around the trawl station within which individuals are drawn from the superindividuals. If given as a numeric vector of length 2, the second element is interpreted as the radius inside which superindividuals are selected when converting from biomass in gram/square meter to nautical area scattering coefficient NASC.
#' @param N						The number of individuals to draw for each trawl station.
#' @param cores					The number of cores to use for writing XML files and for the bootstrapping in the StoX project, given either as a named list such as list(biotic=1, acoustic=1, bootstrap=1), og a single numeric setting the cores for all.
#' @param unit					The unit used in the report form the StoX project (see \code{\link{getPlottingUnit}}).
#' @param superindFilter			An optional function acting on the superind data, such as longitude or latitude, e.g., superindFilter <- function(x) x$longitude < 20.
#' @param ...					Arguments passed to \code{\link{surveyPlanner}}.
#' @param Wg					The biomass area density in units of g/m^2 (gram per square meter in the horizontal plane).
#' @param LcmOne				The length in cm representative of the individuals in the survey region. This length should be so that the NASC converted from the biomass in grams per square meters is unbiased. The function \code{\link{getTotalBiomass}} with type="length" can be used to obtain this value.
#' @param a,b					Parameters of the length-weight relationship Wg = a * Lcm^b, where Wg is the weight in grams and Lcm is the length in cm. Typical values are e.g., a = 0.01 and b = 3.
#' 
#' @export
#' @import Rstox
#' @importFrom maptools unionSpatialPolygons
#' @importFrom utils head
#' @importFrom stats median
#' @rdname biomass2tsb
#'
biomass2tsb <- function(
	# For the StoX project and writing XML files:
	projectName, xmlOutputDir, 
	# For the biomass and superindividuals and other global options:
	biomass, superind, stratum, startdate = "31/1", dateshift = 0, centroid = c(0, 68), seed = 0, nboot=10, xsd = list(acoustic="1", biotic="1.4"),
	# For transects and NASC:
	nTrawl = 50, hours = list(240), type = "RectEnclZZ", knots = 10, equalEffort = FALSE, bearing = "along", distsep = 1, margin = 0.1, tsn = 161722, m = 20, TS0 = -71.9, condPar=NULL, LcmMean=NULL, 
	# For trawl:
	platform = 4174, distance = 5, sweepWidth = 25, probfact = 1, radius=10, 
	N=100, cores=list(biotic=1, acoustic=1, bootstrap=1), unit="mt", superindFilter=NULL, stagesToRemove = c(0, 1, 2), 
	...){
		
	# Not used, calculate the pure transect time outside of the function, and do not consider the time used by trawling as a delay of the acoustic logging.
	# hoursPerTrawl = 2, 
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

	### # Read the biomass:
	### ##### NOTE 2: Remove the rename parameter once the files have the requested variable names, also in the readNorwecomSuperind(): #####
	### if(is.character(biomass) && file.exists(biomass)){
	### 	cat("Read biomass file...\n")
	### 	biomass <- readNorwecomBiomass(
	### 		biomass, centroid=centroid, 
	### 		#vars=c("HERbiom", "Long", "Latt"), 
	### 		#rename=list(HERbiom="Biom")
	### 		vars=c("Biom", "Long", "Latt")
	### 	)
	### }
	### 
	### 
	### # Read the superindividuals:
	### ##### NOTE 3: Remove the rename parameter once the files have the requested variable names, also in the readNorwecomSuperind(): #####
	### if(is.character(superind) && file.exists(superind)){
	### 	cat("Read superind file...\n")
	### 	superind <- readNorwecomSuperind(
	### 		superind, centroid=centroid, 
	### 		vars = c("xpos", "ypos", "Long", "Latt", "female", "age", "inumb", "length", "sweight", "pweight"), 
	### 		rename=list(xpos="gridLongInd", ypos="gridLattInd", Long="gridLong", Latt="gridLatt")
	### 	)
	### }
	### 
	### 
	# Read the biomass:
	if(is.character(biomass) && file.exists(biomass)){
		cat("Reading biomass file...\n")
		biomass <- readNorwecomBiomass(biomass, centroid=centroid)
	}
	# Read the superindividuals:
	if(is.character(superind) && file.exists(superind)){
		cat("Reading superind file...\n")
		superind <- readNorwecomSuperind(superind, centroid=centroid, stagesToRemove=stagesToRemove)
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
	print(t0)
	
	transects <- Rstox::surveyPlanner(
		projectName = stratum, # Here 'stratum' is from the getPelfossPaths().
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
	logDays <- findInterval(transects$Transect$time_start, biomass$time)
	daysOfSurvey <- unique(logDays)
	startDayOfSurvey <- min(daysOfSurvey)
	endDayOfSurvey <- max(daysOfSurvey)
	midDayOfSurvey <- round(mean(c(startDayOfSurvey, endDayOfSurvey)))
	startDateOfSurvey <- min(transects$Transect$time_mid)
	endDateOfSurvey <- max(transects$Transect$time_mid)
	
	
	
	
	### # Get condition exponent from the superindividuals:
	### s <- getSuperindOfDay(superind=superind, day=1)
	### Lcm <- s$length
	### Wg <- s$sweight
	### #pW <- s$pweight * 1e-3
	### plot((Lcm^3), Wg)
	### b_true <- 3
	### a_true <- Wg / Lcm^b_true
	### summary(a_true)
	### a_true <- median(a_true)

	### # Get the typical length of the fish, by a weighted average of the superindividuals in the survey region:
	### id <- double(length(transects$Input$lonlatSP))
	### # Get the union of the SpatialPolygons of the strata:
	### stratumUnion <- maptools::unionSpatialPolygons(transects$Input$lonlatSP, id)
	### # Select the fist in case there are holes in the union:
	### stratumUnionMatrix <- Rstox::getMatrixList(stratumUnion)
	### if(is.list(stratumUnionMatrix)){
	### 	stratumUnionMatrix <- stratumUnionMatrix[1]
	### }
	### else{
	### 	stratumUnionMatrix <- list(stratumUnionMatrix)
	### }
	
	## Get the mean length of the fish weighted by the expected backscatter:
	#startDayOfSurvey <- as.numeric(strftime(min(transects$Tra$time_mid), format = "%j"))
	#endDayOfSurvey <- as.numeric(strftime(max(transects$Tra$time_mid), format = "%j"))
	##midDayOfSurvey <- as.numeric(strftime(median(transects$Tra$time_mid), format = "%j"))
	#midDayOfSurvey <- round(mean(c(startDayOfSurvey, endDayOfSurvey)))
	
	
	###	LcmMean <- getTotalBiomass(superind=superind, stratumPolygons=stratumUnionMatrix, day=midDayOfSurvey, type="length")$AllStrata
	# Get the NASC values from the biomass horizontal area density:
	insideRange <- if(length(radius) == 2) radius[2] else NULL
	transects <- interpolateBiomassToTransects(biomass=biomass, superind=superind, transects=transects, centroid=centroid, logDays=logDays, m=m, TS0=TS0, insideRange=insideRange, condPar=condPar, LcmMean=LcmMean)
	# NORWECOM biomass is given in grams
	#Wkg <- Wg * 1e-3
	
	###	NASC <- biomass2nasc(Bg, LcmOne=LcmMean, a=a_true, b=b_true, m=m, TS0=TS0)
	###	if(all(is.na(NASC)) || !isTRUE(any(NASC>0))){
	###		warning("No biomass along the transects")
	###		return(NULL)
	###	}
	
	
	#transects$Transect <- cbind(transects$Transect, NASC=NASC)


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
		# This defines ch as an attribute with value 1:
		sa..ch.1 = round(transects$Transect$NASC, digits=6), 
		# Why did we incldue this?:
		NASC = transects$Transect$NASC, 
		stringsAsFactors=FALSE
	)
	
	############################################################################
	############################################################################


	############################################
	##### 3. Get total theoretical biomass #####
	############################################

	cat("Get total biomass...\n")
	strataInd <- match(transects$Stratum$stratum, names(transects$Input$lonlat))
	totalBiomass <- getTotalBiomass(superind=superind, stratumPolygons=transects$Input$lonlat, strataInd=strataInd, type="biomass", superindFilter=superindFilter)
	
	############################################
	############################################


	################################################################
	##### 4. Simulate trawl stations, return biotic data frame #####
	################################################################

	cat("Draw trawl stations...\n")
	
	# Draw trawl stations (get the indices of the log distances with trawl):
	trawlInd <- drawTrawlStationInd(transects$Transect$NASC, nTrawl, probfact=probfact)
	if((length(trawlInd) == 0)){
		return(NULL)
	}
	
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
	biotic <- simulateTrawl(superind, fishStation, seed=seed$trawl, radius=radius, N=N, superindFilter=superindFilter, LcmMean=LcmMean)
	superindCount <- biotic$superindCount
	biotic <- biotic$biotic
	
	### Add trawl info to the transects$Transect: ###
	# Add superindividual count for the trawl stations:
	transects$Transect$superindCount <- NA
	transects$Transect$superindCount[trawlInd] <- superindCount
	# Add serialno for the trawl stations:
	transects$Transect$serialno <- NA
	transects$Transect$serialno[trawlInd] <- fishStation$serialno
	
	# Discard strata with no biotic data:
	#sumBiotic <- by(transects$Transect$stratum, transects$Transect$superindCount, sum, na.rm=TRUE)
	sumBiotic <- by(transects$Transect$superindCount, transects$Transect$stratum, sum, na.rm=TRUE)
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
	#p1 <- plotPelfoss(
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
	
	
	
	suppressWarnings(dir.create(xmlOutputDir))
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
		#SweptAreaDensity = list(
		#	ProcessData = "ReadProcessData", 
		#	SweptAreaMethod = "LengthDependent", 
		#	BioticData = "ReadBioticXML", 
		#	LengthDist = "TotalLengthDist", 
		#	DistanceMethod = "FullDistance",
		#	sweepWidthMethod = "Constant", 
		#	sweepWidth = sweepWidth
		#), 
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
		#	Superindividuals = "SuperIndAbundance", 
		#	FillVariables = "ImputeByAge", 
		#	Seed = "1", 
		#	FillWeight = "Regression"
		#), 
		#EstimateByPopulationCategory = list(
		#	Superindividuals = "FillMissingData", 
		#	LengthInterval = "2.0", 
		#	Scale = "1000", 
		#	Dim1 = "LenGrp", 
		#	Dim2 = "age", 
		#	Dim3 = "SpecCat"
		#)
		EstimateByPopulationCategory = list(
			Superindividuals = "SuperIndAbundance", 
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
	
	hasFish <- which(transects$Transect$superindCount > 0)
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
	report <- Rstox::getReports(projectName, var="Weight", unit="mt", grp1=NULL, drop.out=FALSE)
	
	# Save project data and close the project:
	Rstox::saveProjectData(projectName)
	Rstox::closeProject(projectName)
	
	#   TSB Ab.Sum.5% Ab.Sum.50% Ab.Sum.95% Ab.Sum.mean Ab.Sum.sd Ab.Sum.cv
	# 1 TSB   4124917    6376444    9903694     6864052   2228430 0.3246522
	#######################################
	#######################################
	

	# Output a list of relevant objects:
	out <- list(report=report, totalBiomass=totalBiomass, superind=superind, transects=transects, biomass=biomass, t0=t0, daysOfSurvey=daysOfSurvey, startDayOfSurvey=startDayOfSurvey, midDayOfSurvey=midDayOfSurvey, endDayOfSurvey=endDayOfSurvey, startDateOfSurvey=startDateOfSurvey, endDateOfSurvey=endDateOfSurvey, projectName=projectName)
	
	# Add total (estimated) stock biomass (TSB) and theoretical stock biomass (ThSB):
	#out$report <- getTSB(out, unit=unit)
	
	# Add function inputs:
	out$input <- biomass2tsbInputs
	
	return(out)
}
#'
#' @export
#' @rdname biomass2tsb
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
	NASC <- 4 * pi * 1852^2 * numFish * sigmabsOne
	
	out <- data.frame(
		NASC = NASC, 
		WgOne = WgOne, 
		numFish = numFish, 
		sigmabsOne = sigmabsOne
		)
	
	return(out)
}


#*********************************************
#*********************************************
#' Plot the output from biomass2tsb() or runPelfoss().
#'
#' @param x				The output from \code{\link{biomass2tsb}} (containing object named report, totalBiomass, superind, transects, biomass).
#' @param ncolour		The number of colours used to plot the biomass field in layers.
#' @param col			A list of hue, saturation and value ranges, between which a sequence of 'ncolour' is used for the biomass layers, or a color vector such as gray(ncolour, 0.2, 0.8).
#' @param logColscale	Logical: If TRUE, use 10 * log10 for the plotted biomass layers.
#' @param firstcol		The color to use for the first biomass level (the default, NA, omits plotting the first layer, which may include 0-values).
#' @param NASCthr		A value below which NASC values are omitted when plotted as black dots with size proportional to NASC^NASCexp.
#' @param NASCmax_size	Maximum size of the NASC dots.
#' @param biomassAlpha	Transparency of the biomass layers.
#' @param trawlSize		The size to use for the red asteriks representing trawl stations.
#' 
#' @export
#' @import ggplot2
#' @importFrom Rstox plotStratum 
#' @importFrom grDevices hsv
#' @rdname plotPelfoss
#'
plotPelfoss <- function(
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
	plot = c("biom", "nasc", "trawl", "super"), 
	daysOfSurvey = NULL, 
	timerange = NULL, 
	#alldays = TRUE, 
	...){
	
	
	# Apply the timerange:
	if(length(timerange) == 2){
		insideTimeRange <- x$transects$Transect$time_mid >= timerange[1] & x$transects$Transect$time_mid <= timerange[2]
		x$transects$Transect <- x$transects$Transect[insideTimeRange, ]
	}
	
	# Run the original surveyPlanner plot:
	if(nrow(x$transects$Transect) == 0){
		p <- ggplot()
		p0 <- p
		return(list(p=p, p0=p0))
	}
	p <- plotStratum(x$transects, ...)
	p0 <- p
	
					### # Get the super individuals:
					### if(any(startsWith(tolower(plot), "biom"))){
					### 	superind <- readNorwecomSuperind(files$superind, centroid=centroid)
					### }
	
	# Get the biomass to plot, either as an average of the days of the survey, or at the mid day of the survey:
	x <- getBiomassOfSurvey(x, daysOfSurvey=daysOfSurvey, fresh=if(length(daysOfSurvey)) TRUE else FALSE)
	
	maxBiomass <- max(x$biomassOfSurvey$biomass, na.rm=TRUE) * biomThr
	if(logColscale){
		biomassSeq <- seq(0, 10*log10(maxBiomass), length.out=ncolour)
		biomassSeq <- 10^(biomassSeq/10)
	}
	else{
		biomassSeq <- seq(0, maxBiomass, length.out=ncolour)
	}
	
	if(is.list(col) && all(c("h", "s", "v") %in% names(col))){
		col <- hsv(h=seq(col$h[1], col$h[2], l=ncolour), s=seq(col$s[1], col$s[2], l=ncolour), v=seq(col$v[1], col$v[2], l=ncolour))
	}
	
	if(length(firstcol)){
		col[1] <- firstcol
	}
	colorInterval <- findInterval(x$biomassOfSurvey$biomass, biomassSeq)
	
	### if("superind" %in%  tolower(plot)){
	### 	x <- getSuperindOfSurvey(x)
	### 	#if(length(x$superindOfSurvey)==0  && length(x$superind)){
	### 	#	x$superindOfSurvey <- as.data.frame(getSuperindOfDay(x$superind, x$midDayOfSurvey)[c("Long", "Latt")])
	### 	#}
	### 	#else if(length(x$superindOfSurvey)==0){
	### 	#	warning("Super individual data missing, and addSuperind cannot be TRUE.")
	### 	#}
	### 	p <- p + geom_point(data=x$superindOfSurvey, aes(x=Long, y=Latt), shape=2, color="black", size=0.2)
	### }
	
	# Add biomass to the transect plot:
	if(any(startsWith(tolower(plot), "biom"))){
		for(i in seq_len(ncolour)){
			p <- p + geom_point(data=x$biomassOfSurvey[colorInterval==i, ], aes_string(x="longitude", y="latitude"), colour=col[i], alpha=biomAlpha, shape=biomShape, size=biomSize)
		}
	}
	
	# Add NASC values larger than NASCthr of the fraction of NASC relative to maxNASC:
	maxNasc <- max(x$transects$Transect$NASC, na.rm=TRUE)
	pointSize <- (x$transects$Transect$NASC/maxNasc)^NASCexp
	
	pointSize[pointSize < NASCthr] <- NA
	x$transects$Transect <- cbind(x$transects$Transect, pointSize=pointSize)
	#if("nasc" %in%  tolower(plot)){
	if(any(startsWith(tolower(plot), "nasc"))){
		p <- p + geom_point(data=x$transects$Transect, aes_string(x="lon_mid", y="lat_mid", size="pointSize"), shape=20)  +  scale_size_area(max_size=NASCmax_size, guide=FALSE)
	}
	
	# Plot the trawl stations
	#if("trawl" %in%  tolower(plot)){
	if(any(startsWith(tolower(plot), "trawl"))){
		p <- p + geom_point(data=x$transects$Transect[ x$transects$Transect$trawl, ], aes_string(x="lon_mid", y="lat_mid"), shape=42, color="red", size=trawlSize)
	}
	
	if(length(stratumcol)){
		p <- p + scale_fill_manual(values = rep(stratumcol, length.out=nrow(x$transects$Stratum)))
	}
	
	suppressWarnings(print(p))
	list(p=p, p0=p0)
}
### #'
### #' @export
### #' @importFrom Rstox getPlottingUnit 
### #' @rdname plotPelfoss
### #'
### plotTSB <- function(x, unit="mt"){
### 	scale <- Rstox::getPlottingUnit(unit, var="weight")$scale
### 	x$totalBiomass$Total <- x$totalBiomass$Total / scale
### 	x$totalBiomass$NonEmptyStrata <- x$totalBiomass$NonEmptyStrata / scale
### 	x$totalBiomass$AllStrata <- x$totalBiomass$AllStrata / scale
### 	TSB <- getTSB(x, unit=unit)$bootstrap$abnd
### 	
### 	plot(x$totalBiomass$Total, ylim=c(0, max(x$totalBiomass$Total)))
### 	points(x$totalBiomass$AllStrata, col=2)
### 	points(x$totalBiomass$NonEmptyStrata, col=3)
### 	abline(h = TSB$Ab.Sum.mean)
### 	abline(h = TSB$'Ab.Sum.50%', lty=2)
### 	abline(v = range(x$daysOfSurvey))
### }
#'
#' @export
#' @keywords internal
#' @importFrom Rstox getPlottingUnit 
#' @rdname plotPelfoss
#'
getTSB <- function(x, unit="mt", digits=4){
	
	# Remove the name column:
	x$report$reportAbundance$bootstrap$abnd <- x$report$reportAbundance$bootstrap$abnd[,-1]
	
	# Scale the Rstox report back to the baseunit, which is grams and the same as used in the input biomass data:
	x$report$reportAbundance$bootstrap$abnd[names(x$report$reportAbundance$bootstrap$abnd) != "Ab.Sum.cv"] <- x$report$reportAbundance$bootstrap$abnd[names(x$report$reportAbundance$bootstrap$abnd) != "Ab.Sum.cv"] * x$report$reportAbundance$bootstrap$scale
	
	# The scale to the requested unit:
	scale <- getPlottingUnit(unit, var="weight")$scale
	x$report$reportAbundance$bootstrap$abnd[names(x$report$reportAbundance$bootstrap$abnd) != "Ab.Sum.cv"] <- x$report$reportAbundance$bootstrap$abnd[names(x$report$reportAbundance$bootstrap$abnd) != "Ab.Sum.cv"] / scale
	#ThSB_inside <- x$totalBiomass$NonEmptyStrata[x$midDayOfSurvey] / scale
	#ThSB_insideAll <- x$totalBiomass$AllStrata[x$midDayOfSurvey] / scale
	#ThSB_total <- x$totalBiomass$Total[x$midDayOfSurvey] / scale
	ThSB_inside <- mean(x$totalBiomass$NonEmptyStrata[x$daysOfSurvey]) / scale
	ThSB_insideAll <- mean(x$totalBiomass$AllStrata[x$daysOfSurvey]) / scale
	ThSB_total <- mean(x$totalBiomass$Total[x$daysOfSurvey]) / scale
	
	x$report$reportAbundance$bootstrap$abnd <- cbind(
		data.frame(
			startDayOfSurvey = x$startDayOfSurvey, 
			midDayOfSurvey = x$midDayOfSurvey, 
			endDayOfSurvey = x$endDayOfSurvey, 
			startDateOfSurvey = x$startDateOfSurvey,
			endDateOfSurvey = x$endDateOfSurvey, stringsAsFactors=FALSE
		),
		x$report$reportAbundance$bootstrap$abnd[, c("Ab.Sum.5%", "Ab.Sum.50%", "Ab.Sum.95%", "Ab.Sum.mean", "Ab.Sum.sd", "Ab.Sum.cv")], 
		data.frame(
			ThSB_inside = ThSB_inside, 
			ThSB_insideAll = ThSB_insideAll, 
			ThSB_total = ThSB_total, 
			# Mean:
			TSB_meanByInside = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.mean"]] / ThSB_inside, 
			TSB_meanByInsideAll = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.mean"]] / ThSB_insideAll, 
			TSB_meanByTotal = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.mean"]] / ThSB_total, 
			# Median:
			TSB_medianByInside = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.50%"]] / ThSB_inside, 
			TSB_medianByInsideAll = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.50%"]] / ThSB_insideAll, 
			TSB_medianByTotal = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.50%"]] / ThSB_total, 
			# 5 percentile:
			TSB_5percentByInside = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.5%"]] / ThSB_inside, 
			TSB_5percentByInsideAll = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.5%"]] / ThSB_insideAll, 
			TSB_5percentByTotal = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.5%"]] / ThSB_total, 
			# 95 percentile:
			TSB_95percentByInside = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.95%"]] / ThSB_inside, 
			TSB_95percentByInsideAll = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.95%"]] / ThSB_insideAll, 
			TSB_95percentByTotal = x$report$reportAbundance$bootstrap$abnd[["Ab.Sum.95%"]] / ThSB_total, 
			stringsAsFactors=FALSE
		)
	)
	
	# Set the significant digits:
	areNumeric <- sapply(x$report$reportAbundance$bootstrap$abnd, is.numeric)
	x$report$reportAbundance$bootstrap$abnd[areNumeric] <- signif(x$report$reportAbundance$bootstrap$abnd[areNumeric], digits)
	x$report$reportAbundance
}
#'
#' @export
#' @importFrom utils head write.table
#' @rdname plotPelfoss
#' 
getAllTSB <- function(dir, run="Run1", outfile=NULL, copy=TRUE){
	# Get PELFOSS files:
	getPelfossFiles <- function(l, lshort, dirs, type="output", copy=TRUE){
		# Define keys to grep by:
		keys <- list(
			output = "PelfossProjectRstox/output/r/report/pelfossOutput.RData", 
			plot = "PelfossProjectRstox/output/r/report/pelfossOutput.png"
		)
		
		# Get the directory in which to put the files:
		allPelfossDir <- file.path(
			dirs$thisrundir, 
			paste("all", type, sep="_"), 
			paste(type, "TSB", sep="_")
		)
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
		direction <- getBetween(x, "_direction_", "_timing_")
		timing <- getBetween(x, "_timing_", "_seed_")
		seed <- getBetween(x, "_seed_", NULL)
		
		data.frame(
			Survey = survey, 
			Resolution = res, 
			Year = year, 
			Direction = direction, 
			Timing = timing, 
			Seed = seed, stringsAsFactors=FALSE
		)
	}
	
	dirs <- getPelfossSkeleton(dir, run=run)
	
	# List all projects:
	projectsDir <- dirs$projectdir
	l <- list.files(projectsDir, recursive=TRUE, full.names=TRUE)
	lshort <- list.files(projectsDir, recursive=TRUE, full.names=FALSE)
	
	# Get and copy files:
	allOutputFiles <- getPelfossFiles(l, lshort, dirs, type="output", copy=copy)
	allPlotFiles <- getPelfossFiles(l, lshort, dirs, type="plot", copy=copy)
	
	# Load all output files into a list:
	data <- lapply(allOutputFiles$new, function(x) mget(load(x)))
	names(data) <- allOutputFiles$filetag
	
	# Get the reports:
	reports <- do.call(rbind, lapply(data, function(x) getTSB(x$pelfossOutput)$bootstrap$abnd))
	#reports <- do.call(rbind, lapply(seq_along(data), function(ind) {print(ind); print(allOutputFiles$filetag[ind]); getTSB(data[[ind]]$pelfossOutput)$bootstrap$abnd}))
	
	reports <- cbind(splitRunName(rownames(reports)), reports, name=rownames(reports), stringsAsFactors=FALSE)
	rownames(reports) <- NULL
	
	reports <- reports[order(reports$Survey, reports$Year),]

	# Write the reports to file:
	if(length(outfile)==0 || !is.character(outfile)){
		outfile <- file.path(dirs$reportsdir, c("allReports.txt", "allReportsSmall.txt"))
	}
	write.table(reports, file=outfile[1], sep="\t", dec=".", row.names=FALSE, fileEncoding="UTF-8")
	if(length(outfile) > 1){
		reports_small <- simplifyAllTSB(reports)
		write.table(reports_small, file=outfile[2], sep="\t", dec=".", row.names=FALSE, fileEncoding="UTF-8")
	}
	
	return(list(reports=reports, data=data, allOutputFiles=allOutputFiles, allPlotFiles=allPlotFiles))
}
simplifyAllTSB <- function(reports, by="InsideAll"){
	columnsToKeep <- c(
		"Survey", 
		"Resolution", 
		"Year", 
		"Direction", 
		"Timing", 
		"Seed", 
		paste0("TSB_meanBy", by), 
		paste0("TSB_5percentBy", by), 
		paste0("TSB_95percentBy", by), 
		"Ab.Sum.cv"
	)
	newNames <- c(
		"Survey", 
		"Resolution", 
		"Year", 
		"Direction", 
		"Timing", 
		"Seed", 
		"TSB_mean", 
		"TSB_5percent", 
		"TSB_95percent", 
		"TSB_CV"
	)
	
	reports_small <- reports[columnsToKeep]
	names(reports_small) <- newNames
	
	return(reports_small)
}
### #'
### #' @export
### #' @importFrom utils head write.table
### #' @rdname plotPelfoss
### #' 
### getPelfossOutput.RData <- function(dir, run="Run1", survey="Herring_IESNS", year=2010, res=4, reversed=FALSE, dateshift=0, seed=0){
### 	# Get the rData file:
### 	files <- getPelfossPaths(dir, run=run, survey=survey, year=year, res=res, reversed=reversed, dateshift=dateshift, seed=seed)
### 	projectName <- file.path(files$projectPath, "PelfossProjectRstox")
### 	projectPaths <- getProjectPaths(projectName)
### 	RData.file <- file.path(projectPaths$RReportDir, "pelfossOutput.RData")
### }

#*********************************************
#*********************************************
#' Simulate trawl stations from superindividuals.
#'
#' \code{simulateTrawl} Simulated trawl data from the superindividuals inside a \code{radius} nautical mile radius around each trawl station,  where the probability of selecting an individual from the superindividual is proportional to its number of (identical) individuals. \cr \cr
#' \code{drawTrawlStationInd} Draws trawl station locations as log-distances (acoustic PSU locations), where the probability of a log-distance to be selected is a mixture of the NASC values and the uniform probability along the track. The mixture is given by the \code{probfact} argument. \cr \cr
#' 
#' @param superind	The superind data read from \code{readNcVarsAndAttrs} (used in \code{readNorwecomBiomass} and \code{readNorwecomSuperind}).
#' @param x	The input NetCDF4 file.
#' @param x	The input NetCDF4 file.
#' @param x	The input NetCDF4 file.
#' @param x	The input NetCDF4 file.
#' @param x	The input NetCDF4 file.
#' @param x	The input NetCDF4 file.
#' 
#' @export
#' @importFrom Rstox getSeedV
#' @importFrom utils head
#' @keywords internal
#' @rdname simulateTrawl
#'
simulateTrawl <- function(superind, fishStation, seed=0, radius=10, nn=NULL, N=100, lengthunit=2, superindFilter=NULL, LcmMean=NULL){
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
	simulateTrawlOne <- function(ind, superind, fishStation, seedV, radius=10, nn=NULL, N=100, validLengthunit=2, superindFilter=NULL, LcmMean=NULL){
		
		# Get the day of the trawl station:
		day <- findInterval(fishStation$starttime[ind], superind$time)
		#thisSuperind <- lapply(superind[c("X", "Y", "inumb", "female", "age", "length", "pweight")], function(x) x[, day])
		thisSuperind <- getSuperindOfDay(superind, superindDay=day, superindFilter=superindFilter)
		
		# First find the superindividuals within the radius around the station:
		d <- edist(
			thisSuperind$x - fishStation$x_mid[ind], 
			thisSuperind$y - fishStation$y_mid[ind]
		)
		
		# Get the superindividuals inside the radius:
		# Also discard invalid superindividuals:
		if(length(nn)){
			inside <- which(thisSuperind$weight < 1e10)
			inside <- order(d[inside])[seq_len(nn)]
		}
		else{
			inside <- which(d <= radius[1] & thisSuperind$weight < 1e10)
		}
		superindCount = length(inside)
		
		# Get the number of fish of the valid superindividuals:
		size <- thisSuperind$inumb[inside]
		
		# We could consider weighting by the distance in the below probability of selecting a fish from the superindividual:
		prob <- size / sum(size)
		set.seed(seedV[ind])
		if(superindCount == 0){
			return(data.frame())
		}
		else if(superindCount == 1){
			s <- rep(inside, N)
		}
		else{
			s <- sample(inside, size=N, prob=prob, replace=TRUE)
		}
		
		# Generate the individual samples:
		individual <- data.frame(
			serialno = fishStation$serialno[ind], 
			superind = s,
			superindCount = superindCount,
			minDist = min(d),
			specimenno = seq_len(N), 
			lengthunit = lengthunit, 
			### OLD: # Use the input LcmMean if present (forcing the lengths to a single value, used for testing):
			### This was removed on 2019-05-22, after discovering that all lengths were equal to LcmMean in the biotic files when this was set in as input to the conversion from biomass to NASC. LcmMean should not be used in the reversed conversion, i.e., from NASC to biomass, in StoX:
			### length = if(length(LcmMean)) LcmMean else thisSuperind$length[s],
			length = thisSuperind$length[s],
			sex = 2 - thisSuperind$female[s], # The definition of biotic xml has male=2 and female=1. In the Norwecom model male=0 and female=1.
			weight.individual = thisSuperind$weight[s],
			age..agedetermination.1 = thisSuperind$age[s]
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
	individual <- lapply(seq_len(nstations), simulateTrawlOne, superind=superind, fishStation=fishStation, seedV=seedV, radius=radius, nn=nn, N=N, validLengthunit=validLengthunit, superindFilter=superindFilter, LcmMean=LcmMean)
	superindCount <- sapply(individual, function(x) if(length(x$superindCount)) head(x$superindCount, 1) else 0)
	minDist <- sapply(individual, function(x) if(length(x$minDist)) head(x$minDist, 1) else NA)
	individual <- do.call("rbind", individual)
	
	# Merge the fish stations and individuals
	out <- list(biotic=merge(fishStation, individual), superindCount=superindCount)
	out
}
#'
#' @export
#' @keywords internal
#' @rdname simulateTrawl
#'
drawTrawlStationInd <- function(NASC, nstations=1, seed=0, probfact=0){
	# Generate the cummulative probability distribution as the cumsum of the NASC, and draw based on this distribution:
	NASC0 <- NASC
	NASC0[is.na(NASC0)] <- 0
	
	if(!any(NASC0 > 0)){
		return(NULL)
	}
	
	prob <- NASC0 / sum(NASC0)
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
#' Extract dayly subsets of biomass and superind.
#'
#' \code{getSuperindOfDay} Gets the superind data of one day. \cr \cr
#' \code{getBiomassOfSurvey} Gets the average biomass data of all days of a survey. \cr \cr
### #' \code{getSuperindOfSurvey} Gets the superind data of the mid day of the survey (only used in \code{plotPelfoss} with "superind" included in the argument \code{plot}). \cr \cr
#' \code{getTotalBiomass} Gets the total biomass/'fish length which would produce an unbiased biomass via acoustic backscattering'. Uses the superind data. \cr \cr
#' 
#' @param biomass,superind	The biomass and superind data read from \code{readNcVarsAndAttrs} (used in \code{readNorwecomBiomass} and \code{readNorwecomSuperind}).
#' @param day				The day on which to extract the superind data.
#' @param NAsByFirst		Logical: If TRUE, identify NAs by the first day of the superind data.
#' @param daysOfSurvey		A vector of the days of the survey over which to average the biomass.
#' @param midDayOfSurvey	The mid day of the survey ar which to extract the superind data.
#' @param stratumPolygons	A named list of stratum polygons as data frames with longitude in the first column and latitude in the second.
#' @param strataInd			A vector of indices of the strata to compute the total biomass over.
#' @param days				A vector of indices of the days to compute the total biomass for.
#' @param type				A string specifying the type of calculation (first element used), one of "biomass" and "length", returning total biomass and the length that would result in unbiased estimated total biomass if bioass was converted to NASC and aggregated throughout the survey region.
#' @param superindFilter		A function such as superindFilter <- function(x){x$Long < 20 & x$Lat < 70}, taking the superind data as input and returning a logical vector discarding certain longitude/latitude combinations.
#' 
#' @export
#' @importFrom sp point.in.polygon
#' @keywords internal
#' @rdname getTotalBiomass
#' 
getSuperindOfDay <- function(superind, superindDay=NULL, superindFilter=NULL, stratumPolygons=NULL, transectsXyOfDay=NULL, insideRange=NULL){
	# Get the number of days:
	numDays <- length(superind$time)
	# Get the indices of variables with days at the second dimension:
	withDays <- which(sapply(superind, function(x) identical(ncol(x), numDays)))
	
	# First subset by superindDay:
	if(length(superindDay)){
		superind[withDays] <- lapply(superind[withDays], "[", , superindDay)
	}
	
	### # Get the indices of superindividuals with valid data of the given day:
	### if(NAsByFirst){
	### 	valid <- which(superind[[withDays[1]]][, 1] < maxValue)
	### }
	### else{
	### 	valid <- which(superind[[withDays[1]]][, day] < maxValue)
	### }
	# Apply also the filter given by the funciton superindFilter:
	if(is.function(superindFilter)){
		insideFilter <- which(superindFilter(superind))
		### valid <- intersect(valid, insideFilter)
		superind[withDays] <- lapply(superind[withDays], "[", insideFilter)
	}
	
	# Get only superindividuals inside the stratumPolygon object, which can be one or more polygons:
	if(length(stratumPolygons)){
		getInside <- function(stratumPolygon, superind){
			inside <- sp::point.in.polygon(
				point.x = superind$longitude, 
				point.y = superind$latitude, 
				pol.x = stratumPolygon[,1], 
				pol.y = stratumPolygon[,2]
			)
		}
		if(is.list(stratumPolygons) && !is.data.frame(stratumPolygons)){
			insidePolygon <- lapply(stratumPolygons, getInside, superind=superind)
			insidePolygon <- do.call(pmax, insidePolygon)
			insidePolygon <- which(insidePolygon > 0)
		}
		else if(length(stratumPolygons)){
			insidePolygon <- getInside(stratumPolygon=stratumPolygons, superind=superind)
			insidePolygon <- which(insidePolygon > 0)
		}
		else{
			insidePolygon <- seq_along(superind$x)
		}
		
		# Subset the superind:
		superind[withDays] <- lapply(superind[withDays], "[", insidePolygon)
	}
	
	
	# Extract inside a radius around all vessel positions:
	if(length(transectsXyOfDay) && length(insideRange)){
		
		# Get the distance between the super individuals and the log distance of the given day
		superindXyOfDay <- cbind(superind$x, superind$y)
		system.time(d <- raster::pointDistance(superindXyOfDay, transectsXyOfDay, longlat=FALSE, allpairs=TRUE))
		
		for(i in seq_len(10)){
			
			# Get superindividuals inside the range:
			thisInsideRange <- insideRange * 2^(i-1)
			insideM <- d <= thisInsideRange
			inside <- rowSums(insideM) > 0
			if(any(inside, na.rm=TRUE)){
				message("'insideRange' = ", thisInsideRange, " for 'superindDay' = ", superindDay)
				break
			}
		}
		
		superind[withDays] <- lapply(superind[withDays], "[", inside)
	}
	
	
	return(superind)
}
#'
#' @export
#' @keywords internal
#' @rdname getTotalBiomass
#' 
getBiomassOfSurvey <- function(x, daysOfSurvey=NULL, fresh=FALSE){
	if(fresh || length(x$biomassOfSurvey)==0){
		if(length(x$daysOfSurvey) && length(daysOfSurvey)==0){
			daysOfSurvey <- x$daysOfSurvey
		}
		
		# Get the average biomass of the days of the survey:
		thisBiomass <- x$biomass$biomass[,,daysOfSurvey]
		dim(thisBiomass) <- c(prod(dim(thisBiomass)[1:2]), length(daysOfSurvey))
		thisBiomass <- rowMeans(thisBiomass)
		x$biomassOfSurvey <- data.frame(longitude = c(x$biomass$longitude), latitude = c(x$biomass$latitude), biomass = thisBiomass)
	}
	return(x)
}
### #'
### #' @export
### #' @keywords internal
### #' @rdname getTotalBiomass
### #' 
### getSuperindOfSurvey <- function(x, midDayOfSurvey=NULL){
### 	if(length(x$superindOfSurvey)==0  && length(x$superind)){
### 		if(length(x$midDayOfSurvey)){
### 			midDayOfSurvey <- x$midDayOfSurvey
### 		}
### 		# Get the supeInd of the mid day of the survey:
### 		x$superindOfSurvey <- as.data.frame(getSuperindOfDay(x$superind, midDayOfSurvey)[c("Long", "Latt")])
### 	}
### 	return(x)
### }
#'
#' @export
#' @importFrom sp point.in.polygon
#' @importFrom stats weighted.mean
#' @keywords internal
#' @rdname getTotalBiomass
#' 
getTotalBiomass <- function(superind, stratumPolygons, strataInd=NULL, days=NULL, type=c("biomass", "length"), superindFilter=NULL, union=FALSE, transectsXyOfDay=NULL, insideRange=NULL){
	
	sumBiomass <- function(d){
		sum(d$weight * d$inumb, na.rm=TRUE)
	}
	
	meanLengthByTS <- function(d){
		sqrt(weighted.mean(d$length^2, d$inumb, na.rm=TRUE))
	}
	
	# Function for getting the biomass of all strata of one day:
	getTotalBiomassStratum <- function(superindDay, superind, stratumPolygons, superindFilter, strataInd, type=c("biomass", "length")){
		# Get the indices in the variable in thisSuperind which are inside the straum:
		getBiomassInside <- function(superindDay, superind, superindFilter, stratumPolygons, strataInd, transectsXyOfDay=NULL, insideRange=NULL, type=c("biomass", "length")){
			
			# Get the superindividualt of the day and polygon:
			thisSuperind <- getSuperindOfDay(superind=superind, superindDay=superindDay, superindFilter=superindFilter, stratumPolygons=stratumPolygons, transectsXyOfDay=transectsXyOfDay, insideRange=insideRange)
			
			# Get the total biomass:
			if(type[1] == "biomass"){
				total <- sumBiomass(thisSuperind)
			}
			else if(type[1] == "length"){
				total <- meanLengthByTS(thisSuperind)
			}
			else{
				stop("Invalid type in getTotalBiomass()")
			}
			#total <- sum(thisSuperind$pweight[inside] * thisSuperind$inumb[inside])
			total
		}
		
		# Get the total biomass of the entire stock:
		total <- getBiomassInside(superindDay=superindDay, superind=superind, superindFilter=superindFilter, stratumPolygons=NULL, transectsXyOfDay=transectsXyOfDay, insideRange=insideRange, type=type)
		
		# Get the biomass of each stratum:
		if(union){
			byStratum <- getBiomassInside(superindDay=superindDay, superind=superind, superindFilter=superindFilter, stratumPolygons=stratumPolygons, transectsXyOfDay=transectsXyOfDay, insideRange=insideRange, type=type)
		}
		else{
			byStratum <- sapply(stratumPolygons, function(thisStratumPolygon) getBiomassInside(superindDay=superindDay, superind=superind, superindFilter=superindFilter, stratumPolygons=thisStratumPolygon, transectsXyOfDay=transectsXyOfDay, insideRange=insideRange, type=type))
		}
		#byStratum <- sapply(stratumPolygons, getBiomassInside, thisSuperind=thisSuperind, type=type)
		
		# Return the biomass in the strata, the sum of all strata, and the total biomass:
		if(length(strataInd)==0){
			strataInd <- seq_along(byStratum)
		}
		c(byStratum, sum(byStratum[strataInd]), sum(byStratum), total)
	}
	
	# Get a sequence of the days:
	if(length(days)==0){
		days <- seq_along(superind$time)
	}
	
	# Get the biomass for all days and combine to a data frame:
	biomass <- lapply(days, getTotalBiomassStratum, superind=superind, stratumPolygons=stratumPolygons, superindFilter=superindFilter, strataInd=strataInd, type=type)
	biomass <- do.call("rbind", biomass)
	biomass <- as.data.frame(biomass)
	
	if(union){
		biomassNames <- c(1, "NonEmptyStrata", "AllStrata", "Total")
	}
	else{
		biomassNames <- c(if(length(names(stratumPolygons))) names(stratumPolygons) else seq_along(stratumPolygons), "NonEmptyStrata", "AllStrata", "Total")
	}
	
	
	names(biomass) <- biomassNames
	# Return:
	biomass
}


#*********************************************
#*********************************************
#' Utility functions of the pelfoss package.
#'
#' \code{pelfossDefaults} Defines default values which can be overridden by the \code{...} option in \code{runPelfoss}. \cr \cr
#' \code{getPelfossPaths} Given a root directory of a pelfoss folder structure (see ?pelfoss) and a survey name, this function returns the paths to biomass and superind NetCDF4-files, stratum polygons, fishery data files (any data format) of the species as interpreted from the survey name. . \cr \cr
#'
#' @param dirs		A list of paths to directories as returned from \code{link{getPelfossSkeleton}}.
#' @param survey	The survey name.
#' @param year		The year.
#' @param res		The resoludion.
#' @param reversed	Logical: If TRUE reverse the direction of the survey as generated by \code{\link{surveyPlanner}}.
#' @param dateshift	The number of days to shift the timing of the survey by.
#' @param seed		The seed used in \code{\link{surveyPlanner}}.
#'
#' @export
#' @keywords internal
#' @rdname pelfossDefaults
#' 
pelfossDefaults <- function(){
	defaults <- list(
		margin = 0.1, # Related to equalEffort
		bearing = "along",
		distsep = 1, # One nmi log distance is standard
		probfact = 1, # This reserves half of the probability of drawing trawls to the NASC and half to chance
		radius = 10, # Increase this to include more superindividuals in the trawl sample
		N = 100, # Draw 100 fish per trawl
		cores = 1, # Use 1 core for both biotic and acoustic xml and bootstrap
		unit = "mt",
		platform = 4174, # G.O.Sars
		distance = 5, # Trawling distance, equal for all stations, thus ineffective in an acousic-trawl survey
		sweepWidth = 25, # Seewp width, see distance above
		condPar=NULL, 
		LcmMean=NULL
	)
	defaults
}
#'
#' @export
#' @keywords internal
#' @rdname pelfossDefaults
#' 
getPelfossSkeleton <- function(dir, run="Run1"){
	# The parameter 'run' can be given as empty, which gives an error (used in getPelfossPaths()):
	if(length(run) == 0 || nchar(run) == 0){
		stop("The argument run must be given as a positive length string")
	}
	
	# The model, stratum and fishery data are located in the 'data' subdir:
	datadir <- file.path(dir, "data")
	modeldir <- file.path(datadir, "model")
	stratumdir <- file.path(datadir, "stratum")
	fisherydir <- file.path(datadir, "fishery")
	
	# Define the direcory of the runs:
	runsdir <- file.path(dir, "runs")
	thisrundir <- file.path(runsdir, run)
	# Define the subdirs of the run:
	reportsdir <- file.path(thisrundir, "reports")
	XMLfilesdir <- file.path(thisrundir, "XMLfiles")
	projectdir <- file.path(thisrundir, "project")
	all_outputdir <- file.path(thisrundir, "all_output")
	all_plotdir <- file.path(thisrundir, "all_plot")
	
	# Return the paths:
	list(
		datadir = datadir, 
		runsdir = runsdir, 
		thisrundir = thisrundir, 
		modeldir = modeldir, 
		stratumdir = stratumdir, 
		fisherydir = fisherydir, 
		reportsdir = reportsdir, 
		XMLfilesdir = XMLfilesdir, 
		projectdir = projectdir, 
		all_outputdir = all_outputdir, 
		all_plotdir = all_plotdir
	)
}
#'
#' @export
#' @keywords internal
#' @importFrom utils head
#' @rdname pelfossDefaults
#' 
getPelfossPaths <- function(dirs, run="Run1", survey="Herring_IESNS", year=2010, res=4, reversed=FALSE, dateshift=0, seed=0){
	# If the dirs is given as a single directory, use the default structure, reuqiring the 'run' parameter to be set:
	if(length(dirs) == 1){
		dirs <- getPelfossSkeleton(dirs, run=run)
	}
	
	# Extract species name from the survey:
	species <- strsplit(survey, "_", fixed=TRUE)[[1]][1]
	# List the files and identify the file type (one of biomass and superind):
	data <- list.files(file.path(dirs$modeldir, species, paste0("res_", res, "km"), paste0("year_", year)), full.names=TRUE)
	type <- sapply(strsplit(basename(data), "_", fixed=TRUE), head, 1)
	
	# Get the biomass and superind files:
	biomass <- data[tolower(type) == "biomass"]
	superind <- data[tolower(type) == "superind"]
	# Get stratum file:
	stratum <- list.files(file.path(dirs$stratumdir, survey), full.names=TRUE)[1]
	
	# Get the directory in which to copy the StoX project containing all files and plots:
	direction <- c("Normal", "Reversed")[reversed + 1]
	projectPath <- file.path(dirs$projectdir, survey, paste0("res_", res, "km"), paste0("year_", year), paste0("direction_", direction, "_timing_", dateshift, "days"), paste0("seed_", seed))
	
	projectPath <- file.path(projectPath, "PelfossProjectRstox")
	projectPaths <- getProjectPaths(projectPath)
	projectRData <- file.path(projectPaths$RReportDir, "pelfossOutput.RData")
	plotSuperIndDir <- file.path(projectPaths$RReportDir, "plotSuperInd")
	
	
	# Get fishery file:
	fishery <- list.files(file.path(dirs$fisherydir, species), full.names=TRUE)[1]
	
	# Return the paths:
	list(biomass=biomass, superind=superind, stratum=stratum, species=species, projectPath=projectPath, projectRData=projectRData, plotSuperIndDir=plotSuperIndDir, survey=survey, fishery=fishery)
}
#'
#' @export
#' @keywords internal
#' @rdname pelfossDefaults
#' 
createPelfossSkeleton <- function(dir, run="Run1"){
	# Define the folders of the pelfoss folder structure:
	dirs <- unlist(getPelfossSkeleton(dir, run=run))
	#dirs <- file.path(dir, c("data", "fishery", "reports", "stratum", "XMLfiles"))
	lapply(dirs, dir.create, recursive=TRUE)
}
#'
#' @export
#' @keywords internal
#' @rdname pelfossDefaults
#' 
fitLengthWeightSuperind <- function(superindDay, superind, count="inumb", superindFilter=NULL, stratumPolygons=NULL, a=NULL, b=NULL, transectsXyOfDay=NULL, insideRange=NULL){
	# Get the superindividual data of the requested day(s):
	superind1 <- getSuperindOfDay(superind=superind, superindDay=superindDay, superindFilter=superindFilter, stratumPolygons=stratumPolygons, transectsXyOfDay=transectsXyOfDay, insideRange=insideRange)
	
	# Define the data frame of weight and length:
	data <- as.data.frame(superind1[c("weight", "length", "inumb")])
	
	fitLengthWeight(data, count=if(length(count)) data[[count]] else NULL, a=a, b=b)
}
#'
#' @export
#' @keywords internal
#' @importFrom stats coef lm
#' @rdname pelfossDefaults
#' 
fitLengthWeight <- function(L, W=NULL, a=NULL, b=NULL, count=NULL, filter=NULL){
	# The standard allomettric formula is 
	#    W = a * L^b, 
	# where W is weight and L is length, and a and b are parameters.
	# Taking the logarithm and rearranging leads to 
	#    log(W) = log(a) + b * log(L)
	# Thus we transform W and L according to this fomula, and transform any estimate of a in the output
	
	# Read the input data, and create a data frame:
	if(length(L)>0 && length(L) == length(W)){
		data <- data.frame(length=L, weight=W)
	}
	else if(is.data.frame(L) && all(c("length", "weight") %in% names(L))){
		data <- L
	}
	else{
		warning("The input must be a data frame with columns 'length' and 'weight'")
		return(NULL)
	}
	# Add the count if present:
	if(length(count) == nrow(data)){
		data$count <- count
	}
	
	# Apply the filter:
	if(is.function(filter)){
		insideFilter <- which(filter(data))
		data <- data[insideFilter, ]
	}
	
	
	# create the formula to use, depending on whether a, b, or both should be estimated:
	useExp <- 0
	if(length(a)==0 && length(b)==1){
		formula <- log(weight) - b * log(length) ~ 1
		fitnames <- "a"
		useExp <- 1
	}
	else if(length(a)==1 && length(b)==0){
		formula <- log(weight) - a ~ log(length) + 0
		fitnames <- "b"
	}
	else if(length(a)==0 && length(b)==0){
		formula <- log(weight) ~ log(length) + 1
		fitnames <- c("a", "b")
		useExp <- 1
	}
	else{
		warning("One of a or b must be empty, otherwise nothing will be estimated")
		return(NULL)
	}
	
	
	fit <- coef(lm(formula, weights=data$count, data=data))
	#fit <- coef(lm(formula, weights=data$count * data$length, data=data))
	
	# Transform a:
	if(useExp==1){
		fit[1] <- exp(fit[1])
	}
	
	# Add the appropriate names:
	names(fit) <- fitnames
	
	
	out <- list(fit=fit, data=data, par=list(a=a, a=b))

	return(out)
}
#'
#' @export
#' @keywords internal
#' @importFrom ggplot2 ggplot aes_string geom_jitter geom_point geom_path
#' @rdname pelfossDefaults
#' 
plotLengthWeight <- function(x, cols=c(1, 2), useJitter=FALSE, ...){
	p <- ggplot(data=x$data, aes_string(x="length", y="weight", if(length(x$data$count)) color=x$data$count))
	
	p <- p + do.call(if(useJitter) geom_jitter else geom_point, list(...))

	rangeL <- range(x$data$length, na.rm=TRUE)
	gridL <- seq(rangeL[1], rangeL[2], l=100)
	p <- p + geom_path(data=data.frame(x=gridL, y=x$fit["a"] * gridL^x$fit["b"]), aes_string(x="x", y="y"), color=cols[2])
	p
}









	
	
### 	
### 	
### 	
### 	# Thus we transform W and L according to this fomula, and transforms any estimate of a in the output:
### 	logW <- log(W)
### 	logL <- log(L)
### 	
### 	useExp <- 0
### 	if(length(a)==0 && length(b)==1){
### 		formula <- logW - b*logL ~ 1
### 		fitnames <- "a"
### 		useExp <- 1
### 	}
### 	else if(length(a)==1 && length(b)==0){
### 		formula <- logW - a ~ logL + 0
### 		fitnames <- "b"
### 	}
### 	else if(length(a)==0 && length(b)==0){
### 		formula <- logW ~ logL + 1
### 		fitnames <- c("a", "b")
### 		useExp <- 1
### 	}
### 	else{
### 		warning("One of a or b must be empty, otherwise nothing will be estimated")
### 		return(NULL)
### 	}
### 	
### 	fit <- coef(lm(formula, weights=weights))
### 	
### 	# Transform a:
### 	if(useExp==1){
### 		fit[1] <- exp(out[1])
### 	}
### 	
### 	names(fit) <- fitnames
### 	
### 	out <- list(fit=fit, W=W, L=L, weights=weights, input_a=a, input_b=b)
### 
### 	return(out)
### }
### 




#*********************************************
#*********************************************
#' Function for getting weights to apply to the strata based on the landings.
#'
#' @param fisheryFile	The path to the csv file containing the catches from the fishery. The file must contain the columns DATE, in the format "yyyy-mm-dd", VEKT, and LON and LAT in decimal degrees.
#' @param stratumFile	A file with strata definitions readable by \code{\link[Rstox]{readStrataPolygons}}.
#' @param strata		a vector of names of the strata to include.
#' @param startdate		A string of the start date of the survey, given as "d/m"
#' @param numdays		An integer giving the number of days to use for the survey.
#'
#' @importFrom sp point.in.polygon
#'
getFisheryWeights <- function(fisheryFile, stratumFile, strata, startdate, numdays){
	# Function for summing the weight inside each stratum:
	insideStratum <- function(stratum, fishery, stratumPolygon, startdate, numdays){
		# Get the rows inside the given date interval:
		dateInt <- paste(fishery$year[1], startdate, sep="/")
		dateInt <- as.Date(dateInt, format="%Y/%d/%m") + c(0, numdays)

		insideDate <- dateInt[1] <= fishery$DATE & fishery$DATE <= dateInt[2]
		allWeight <- sum(fishery$VEKT)
		fishery <- fishery[insideDate, ]
		allWeightDates <- sum(fishery$VEKT)
		
		inside <- sp::point.in.polygon(
			point.x = fishery$LON, 
			point.y = fishery$LAT, 
			pol.x = stratumPolygon$lonlat[[stratum]][,1], 
			pol.y = stratumPolygon$lonlat[[stratum]][,2]
		)
		fishery <- fishery[inside, ]
		
		print(data.frame(allWeightStratum=sum(fishery$VEKT), allWeightDates=allWeightDates, allWeight=allWeight))
	
		sum(fishery$VEKT)
	}
	
	# Read the fishery file:
	fishery <- read.csv(fisheryFile)
	fishery$DATE <- as.Date(fishery$DATE)
	fishery$year <- substr(fishery$DATE, 1,4)

	# Read the strata system of the spawning cruise:
	stratumPolygon <- readStrataPolygons(stratumFile)
	
	
	
	sapply(strata, insideStratum, fishery=fishery, stratumPolygon=stratumPolygon, startdate=startdate, numdays=numdays)
}

#*********************************************
#*********************************************
#' Function for setting the parameters of each survey for paper 1 using the PELFOSS framework.
#'
#' @param run		The name of the simulation run (which is three timings (-30, 0, +30 days), two directions (normal, reversed), and two resolutions (low, high), in total 12 individual runs).
#' @param dir		The directory of the PELFOSS folder structure (see \code{\link{getPelfossSkeleton}}).
#' @param survey	The name of the survey to run, one of "Herring_IESNS", "Herring_NASSHS", "Herring_NASSHS_fishery", "Mackerel_IESSNS", "Mackerel_IESSNS_sept" and "Mackerel_IESSNS_sept_fishery".
#' @param type		The survey design type (see \code{\link[Rstox]{surveyPlanner}}).
#'
#' @export
#'
setSurveyParametersPaper1 <- function(run, dir, survey="Herring_IESNS", type="RectEnclZZ"){
	# Function for getting the effort in each stratum as the total catch by the fishery:
	getFisheryHours <- function(run, dir, survey, strata, startdate, numdays=1){
		files <- getPelfossPaths(run=run, dir=dir, survey=survey)
		fisheryFile <- files$fishery
		stratumFile <- files$stratum
		hours <- getFisheryWeights(fisheryFile=fisheryFile, stratumFile=stratumFile, strata=strata, startdate=startdate, numdays=numdays)
		hours <- hours / sum(hours) * numdays * 24
		hours
	}
	
	if(survey == "Herring_IESNS"){
		nstrata <- 4
		strata <- as.character(c(4, 3, 1, 2))
		rev <- c(FALSE, TRUE, FALSE, FALSE)
		numdays <- 28
		hours <- list(numdays * 24) # Four weeks
		knots <- 15 
		equalEffort <- TRUE 
		centroid <- c(0, 68)
		startdate <- "1/5" # Cruise conducted in May
		trawlDens <- 0.013 # Corresponding to the Test_Rstox project
		tsn <- 161722
		m <- 20
		TS0 <- -71.9
		# For herring, exclude the Barents sea:
		superindFilter <- function(x){
			x$longitude < 20
		}
	}
	else if(survey == "Herring_NASSHS"){
		nstrata <- 13
		strata <- as.character(c(3, 2, 4, 5, 6, 7, 8, 17, 10, 9, 11, 13, 14))
		rev <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
		numdays <- 10
		hours <- list(numdays * 24)
		knots <- 20
		equalEffort <- TRUE 
		centroid <- c(10, 75)
		startdate <- "15/2" # Spring spawning cruise in February-March
		trawlDens <- 0.03 # Corresponding to the Test_Rstox project
		tsn <- 161722
		m <- 20
		TS0 <- -71.9
		# For herring, exclude the Barents sea:
		superindFilter <- function(x){
			x$longitude < 20
		}
	}
	else if(survey == "Herring_NASSHS_fishery"){
		nstrata <- 13
		strata <- as.character(c(3, 2, 4, 5, 6, 7, 8, 17, 10, 9, 11, 13, 14))
		rev <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
		numdays <- 10
		#hours <- list(numdays * 24)
		knots <- 20
		equalEffort <- FALSE 
		centroid <- c(10, 75)
		startdate <- "15/2" # Spring spawning cruise in February-March
		trawlDens <- 0.03 # Corresponding to the Test_Rstox project
		tsn <- 161722
		m <- 20
		TS0 <- -71.9
		
		# Get the distribution of catches in the strata:
		hours <- getFisheryHours(run=run, dir=dir, survey=survey, strata=strata, startdate=startdate, numdays=numdays)
		print(hours)
		
		#files <- getPelfossPaths(survey=survey, dir=dir)
		#fisheryFile <- files$fishery
		#stratumFile <- files$stratum
		#hours <- getFisheryWeights(fisheryFile=fisheryFile, stratumFile=stratumFile, strata=strata, startdate=startdate, numdays=numdays, survey=survey)
		#hours <- hours / sum(hours) * numdays * 24
		
		# For herring, exclude the Barents sea:
		superindFilter <- function(x){
			x$longitude < 20
		}
	}
	else if(survey == "Mackerel_IESSNS"){
		# Old using all strata
		# nstrata <- 10
		# strata <- as.character(c(11, 10, 4, 5, 6, 3, 2, 1, 7, 9))
		# rev <- c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,FALSE, TRUE)
		nstrata <- 5
		strata <- as.character(c(3, 2, 1, 7, 9))
		rev <- c(FALSE, TRUE, FALSE,FALSE, TRUE)
		numdays <- 28
		hours <- list(numdays * 24) # Same time spent as for the Herring_IESNS (approximately 4 weeks)
		knots <- 2 * 15 # Twice the speed of the Herring_IESNS to account for approximately double area.
		equalEffort <- TRUE 
		# Old using all strata:
		#centroid <- c(-10, 65)
		centroid <- c(0, 68)
		startdate <- "1/7" # Conducted in July
		# startdate <- "1/9" # Fictive cruise in September, which is when the fishery is concentrated in 2, 62.5 (longitude, latitude)
		trawlDens <- 0.013 # Corresponding to the Test_Rstox project
		tsn <- 172414
		m <- 20
		TS0 <- -86.5 # Misund and Beltestad 1996
		#TS0 <- -71.9
		superindFilter <- NULL
		# Removed entry 2:5 in multipolygon 1.
		# Removed entry -1.43162167 60 in multipolygon 1.
		# Removed entry 2:30 in multipolygon 4.
		# Removed entry 2:6 in multipolygon 5.
	}
	else if(survey == "Mackerel_IESSNS_sept"){
		nstrata <- 5
		strata <- as.character(c(3, 2, 1, 7, 9))
		rev <- c(FALSE, TRUE, FALSE,FALSE, TRUE)
		numdays <- 28
		hours <- list(numdays * 24) # Same time spent as for the Herring_IESNS (approximately 4 weeks)
		knots <- 2 * 15 # Twice the speed of the Herring_IESNS to account for approximately double area.
		equalEffort <- TRUE 
		centroid <- c(0, 68)
		startdate <- "1/9" # Fictive cruise in September, which is when the fishery is concentrated in 2, 62.5 (longitude, latitude)
		trawlDens <- 0.013 # Corresponding to the Test_Rstox project
		tsn <- 172414
		m <- 20
		TS0 <- -86.5 # Misund and Beltestad 1996
		#TS0 <- -71.9
		superindFilter <- NULL
	}
	else if(survey == "Mackerel_IESSNS_sept_fishery"){
		# Updated these with the fishery stratum system (added one fishery stratum splitting stratum 1 into a southern and a nothern part):
			nstrata <- 7
			strata <- c(3, 2, "1S(11%)", "1F(89%)", "1N", 7, 9)
			rev <- c(FALSE, TRUE, FALSE,FALSE, FALSE,FALSE, TRUE)
		numdays <- 28
		#hours <- numdays * 24 * c(rep(0.02, 3), 0.88, rep(0.02, 3)) # Two percnet in each of the non-fishery strata. Will this work?
		knots <- 2 * 15 # Twice the speed of the Herring_IESNS to account for approximately double area.
		equalEffort <- FALSE 
		centroid <- c(0, 68)
		startdate <- "1/9" # Fictive cruise in September, which is when the fishery is concentrated in 2, 62.5 (longitude, latitude)
		trawlDens <- 0.013 # Corresponding to the Test_Rstox project
		tsn <- 172414
		m <- 20
		TS0 <- -86.5 # Misund and Beltestad 1996
		#TS0 <- -71.9
		superindFilter <- NULL
		
		# Get the distribution of catches in the strata:
		hours <- getFisheryHours(run=run, dir=dir, survey=survey, strata=strata, startdate=startdate, numdays=numdays)
		round(hours / sum(hours) * 100) # 11, 89
	}
	else {
		stop(paste("Survey", survey, "not implemented (use 'Herring_IESNS', 'Herring_NASSHS' or 'Mackerel_IESSNS')."))
	}
	
	out <- list(
		survey = survey, 
		nstrata = nstrata, 
		strata = strata,
		rev = rev,
		hours = hours,
		knots = knots,
		equalEffort = equalEffort, 
		centroid = centroid,
		startdate = startdate, 
		trawlDens = trawlDens, 
		tsn = tsn,
		m = m,
		TS0 = TS0,
		superindFilter = superindFilter, 
		type = type
	)
	
	# Chech that the length of strata matches the nstrata:
	if(out$nstrata != length(out$strata)){
		warning("Mismatch between nstrata and strata.")
	}
	
	# Add number of trawls from the trawling density input:
	sumHours <- if(is.list(out$hours)) unlist(out$hours) else sum(rep(out$hours, length.out=out$nstrata))
	out$nTrawl <- round(out$trawlDens * sumHours * out$knots)
	
	return(out)
}



#*********************************************
#*********************************************
#' Function for running all combinations of parameters of a survey stored in the spec, which is a data frame with one row per run.
#'
#' @param spec		A data frame of the specifications of res and year (these are linked in the current PELFOSS framework), dateshift and reversed.
#' @param run		The name of the simulation run (which is three timings (-30, 0, +30 days), two directions (normal, reversed), and two resolutions (low, high), in total 12 individual runs).
#' @param dir		The directory of the PELFOSS folder structure (see \code{\link{getPelfossSkeleton}}).
#' @param surveyPar	A list of survey parameters as returned from \code{\link{setSurveyParametersPaper1}}.
#' @param seed		A single integer, or a list of seeds used in the funciton, including the following seeds: ('transect') for drawing transects using surveyPlanner(), ('trawl') for drawing trawl stations along the transects with probability as a funciton of the NASC (see \code{probfact}), and ('bootstrap') for getting the final estimate using runBootstrap() (see \code{nboot}).
#' @param nboot		Number of bootstrap replicates.
#' @param radius	A vector of one or two elements giving the radius inside which to sample trawl stations, and the radius inside which to estimate the parameters a, b and L in the length-weight relationship. If only one element is given, the length-weight relationship is estimated from all adult super individuals.
#' @param subset	Either "all" (the default) or an integer vector giving the subset of the simulations to run.
#' @param ...		Additional inputs overriding the defaults returned by pelfossDefaults().
#'
#' @export
#'
runAll <- function(spec, run, dir, surveyPar, seed=1, nboot=100, radius=10, subset="all", ...){
	# Run all the cased in 'spec':
	runs <- nrow(spec)
	times <- double(runs)
	files <- double(runs)
	if(subset=="all"){
		subset <- seq_len(runs)
	}
	else{
		subset <- intersect(subset, seq_len(runs))
	}
	for(i in subset){
		thisSpecList <- spec[i,]
		cat("############################################################\nProject:", surveyPar$survey, "\n")
		print(thisSpecList)
		tempTimes <- system.time(tempFiles <- do.call(runPelfoss, c(list(run=run, dir=dir, surveyPar=surveyPar, seed=seed, nboot=nboot, ...), thisSpecList, list(radius=radius))))
		if(length(tempFiles)){
			files[i] <- tempFiles
		}
		print(tempTimes)
		times[i] <- tempTimes[3]
	}
	# Return the process times:
	print(sum(times))
	cat("############################################################\n\n")
	list(times=times, files=files)
}

#*********************************************
#*********************************************
#' Funciton for running all specified surveys of paper 1 of the PELFOSS project.
#'
#' @param type							The survey design type (see \code{\link[Rstox]{surveyPlanner}}).
#' @param run							The name of the simulation run (which is three timings (-30, 0, +30 days), two directions (normal, reversed), and two resolutions (low, high), in total 12 individual runs).
#' @param survey						The name(s) of the survey(s) to run. Must be one or more of "Herring_IESNS", "Herring_NASSHS", "Herring_NASSHS_fishery", "Mackerel_IESSNS", "Mackerel_IESSNS_sept" and "Mackerel_IESSNS_sept_fishery", which are the surveys used in the Paper 1 of the PELFOSS project.
#' @param dir							The directory of the PELFOSS folder structure (see \code{\link{getPelfossSkeleton}}).
#' @param speedFact						A factor to speed the survey up or down by. Values larger than 1 speed up and values smaller than 1 slow down.
#' @param condPar						A vector of two elements giving the parameters a and b in the length-weight relationship Wg = a * Lcm^b, where Wg is the weight in grams and Lcm is the length in cm. Typical values are e.g., a = 0.01 and b = 3. If not set, this will be estimated from the super individuals, possibly inside a radius given as the second elements of \code{radius}.
#' @param LcmMean						The typical length of the fish, representing the mean acoustic backscatter which refects the square of fish length. This can typically be estimated as sqrt(mean(L^2)). If not set, this will be estimated from the super individuals, possibly inside a radius given as the second elements of \code{radius}.
#' @param year,res,dateshift,reversed	These parameters are linked to the folder structure specifide by \code{\link{getPelfossSkeleton}}, and specifies the year/resolusion, which are linked in the PELFOSS folder structure, the shift in date, which is a vector of integer days to shift the survey by, and a logical stating whether to reverse the direction of the survey. See the default values.
#' @param seed		A single integer, or a list of seeds used in the funciton, including the following seeds: ('transect') for drawing transects using surveyPlanner(), ('trawl') for drawing trawl stations along the transects with probability as a funciton of the NASC (see \code{probfact}), and ('bootstrap') for getting the final estimate using runBootstrap() (see \code{nboot}).
#' @param nboot							Number of bootstrap replicates.
#' @param radius						A vector of one or two elements giving the radius inside which to sample trawl stations, and the radius inside which to estimate the parameters a, b and L in the length-weight relationship. If only one element is given, the length-weight relationship is estimated from all adult super individuals.
#' @param subset						Either "all" (the default) or an integer vector giving the subset of the simulations to run.
#' @param cores					The number of cores to use for writing XML files and for the bootstrapping in the StoX project, given either as a named list such as list(biotic=1, acoustic=1, bootstrap=1), og a single numeric setting the cores for all.
#' @param ...							Additional inputs overriding the defaults returned by pelfossDefaults().
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Script for running experiments for Paper 1 in the PELFOSS project:
#' 
#' # The experiment requires the Rstox and pelfoss package to be installed (see https://github.com/Sea2Data/Rstox and https://github.com/Sea2Data/pelfoss for install instructions):
#' library(pelfoss)
#' library(Rstox)
#' 
#' # Define the directory in which to run the simulation experiments:
#' dir <- "~/Projects/PELFOSS"
#' # Download the simulation experiment folder to the experiment directory:
#' expermientURL <- ""
#' temp <- download.file(expermientURL, dir)
#' unzip(temp)
#' 
#' # Specify to use parallel transect lines:
#' type <- "Parallel"
#' # Specify the number of bootstrap replicates used in the variance estimation of the survey estimates:
#' nboot <- 100
#' # Define the year and resolution, and the shift by date and direction, which will constitute a grid in the function runAllPelfossPaper1():
#' year <- 2012
#' res <- 10
#' dateshift <- c(-30, 0, 30)
#' reversed <- c(FALSE, TRUE)
#' # Define the radius inside which the trawl stations are drawn:
#' radius <- 10
#' # Define the surveys to run:
#' survey = c(
#' 	"Herring_IESNS", 
#' 	"Herring_NASSHS", 
#' 	"Herring_NASSHS_fishery", 
#' 	"Mackerel_IESSNS", 
#' 	"Mackerel_IESSNS_sept", 
#' 	"Mackerel_IESSNS_sept_fishery"
#' )
#' # Define the names of the 
#' runs <- c(
#' 	"Parallel", 
#' 	"Parallel_Fixed_abL", 
#' 	"Parallel_Fixed_abL_SpeedTimes50", 
#' )
#' # Run for seeds 1 and 2:
#' seeds <- c(1, 2)
#' # Write xml files and run boostrapping on multiple cores:
#' cores <- 8
#' 
#' ################
#' ##### Runs #####
#' ################
#' ############################################################
#' ##### 1. Parallel, with with estimation of W0 and L0: ######
#' ############################################################
#' for(seed %in% seeds){
#' 	Parallel_Fixed_abL <- runAllPelfossPaper1(
#' 		type=type, run=runs[1], survey=survey, dir=dir, # Paths and names
#' 	year=year, res=res, dateshift=dateshift, reversed=reversed, # spec
#' 	seed=seed, nboot=nboot, radius=radius, subset="all", cores=cores
#' 	)
#' }
#' ############################################################
#' 
#' ############################################################
#' ##### 2. Parallel, with 'condPar' and 'LcmMean' fixed: #####
#' ############################################################
#' for(seed %in% seeds){
#' 	Parallel_Fixed_abL <- runAllPelfossPaper1(
#' 		type=type, run=runs[2], survey=survey, dir=dir, # Paths and names
#' 		condPar=c(0.3, 2), LcmMean=30, 
#' 		year=year, res=res, dateshift=dateshift, reversed=reversed, # spec
#' 		seed=seed, nboot=nboot, radius=radius, subset="all", cores=cores
#' 	)
#' }
#' ############################################################
#' 
#' ############################################################
#' ##### 3. Parallel, with 'condPar' and 'LcmMean' fixed ######
#' ############# and survey conducted on one day: #############
#' ############################################################
#' for(seed %in% seeds){
#' 	Parallel_Fixed_abL <- runAllPelfossPaper1(
#' 		type=type, run=runs[3], survey=survey, dir=dir, # Paths and names
#' 		speedFact=50, condPar=c(0.3, 2), LcmMean=30, 
#' 		year=year, res=res, dateshift=dateshift, reversed=reversed, # spec
#' 		seed=seed, nboot=nboot, radius=radius, subset="all", cores=cores
#' 	)
#' }
#' ############################################################
#' 	
#' ###################
#' ##### Reports #####
#' ###################
#' # Get all reports (set copy = FALSE to speed up if rerunning this function):
#' reports <- lapply(runs, reportAllPelfoss, dir=dir, copy=TRUE)
#' names(reports) <- runs
#' 
#' # Show the parameters:
#' headPar <- lapply(runs, function(this) head(reports[[this]]$par))
#' names(headPar) <- runs
#' headPar
#' 
#' #################
#' ##### Plots #####
#' #################
#' #  Plots with ylim up ro 2.15:
#' temp <- lapply(runs, plotAllPelfossReports, dir=dir, by="InsideAll", ylim=c(0, 2.15), save="allPlots_max2.15_cv_allSeeds", case_adj=c(0.5, 0.95), add.cv=TRUE)
#' 
#' # Read the reports back in from the experiment referred to in the PELFOSS Paper 1:
#' run <- runs[2]
#' r <- readAllPelfossReports(run, dir=dir)
#' 
#' # Write a table for use in supplementary materials:
#' columnsToKeep <- c(
#' 	"Survey", 
#' 	"Resolution", 
#' 	"Year", 
#' 	"Direction", 
#' 	"Timing", 
#' 	"Seed", 
#' 	"TSB_meanByInsideAll", 
#' 	"TSB_5percentByInsideAll", 
#' 	"TSB_95percentByInsideAll", 
#' 	"Ab.Sum.cv"
#' )
#' newNames <- c(
#' 	"Survey", 
#' 	"Resolution", 
#' 	"Year", 
#' 	"Direction", 
#' 	"Timing", 
#' 	"Seed", 
#' 	"TSB_mean", 
#' 	"TSB_5percent", 
#' 	"TSB_95percent", 
#' 	"TSB_CV"
#' )
#' rsmall <- r[columnsToKeep]
#' names(rsmall) <- newNames
#' head(rsmall)
#' outfile <- file.path(dir, "runs", run, "reports", "allReportsSmall.txt")
#' write.table(rsmall, file=outfile, sep="\t", dec=".", row.names=FALSE, fileEncoding="UTF-8")
#' 
#' #######################
#' ##### Diagnostics #####
#' #######################
#' ##### Absolute relative change from seed 1 to seed 2:
#' summary(abs(dd$TSB_mean[dd$Seed == 2] - dd$TSB_mean[dd$Seed == 1]) /  dd$TSB_mean[dd$Seed == 1])
#' 
#' ##### CVs for the 6 surveys:
#' surveys <- unique(dd$Survey)
#' CV <- sapply(surveys, function(x) summary(dd$TSB_CV[dd$Survey == x]))
#' colnames(CV) <- surveys
#' round(CV, digits=2)
#' #         Herring_IESNS Herring_NASSHS Herring_NASSHS_fishery Mackerel_IESSNS Mackerel_IESSNS_sept Mackerel_IESSNS_sept_fishery
#' # Min.             0.12           0.02                   0.01            0.08                 0.08                         0.05
#' # 1st Qu.          0.13           0.06                   0.02            0.13                 0.10                         0.11
#' # Median           0.18           0.08                   0.03            0.16                 0.12                         0.13
#' # Mean             0.18           0.08                   0.04            0.15                 0.22                         0.18
#' # 3rd Qu.          0.21           0.09                   0.04            0.18                 0.38                         0.29
#' # Max.             0.32           0.14                   0.08            0.19                 0.49                         0.35
#' 
#' ##### Average relative estimates for the Herring_NASSHS shifted by +30 days in the normal and reversed direction:
#' mean(dd$TSB_mean[dd$Survey == "Herring_NASSHS" & dd$Direction == "Normal" & dd$Timing == "30days"])
#' mean(dd$TSB_mean[dd$Survey == "Herring_NASSHS" & dd$Direction == "Reversed" & dd$Timing == "30days"])
#' 
#' ##### The range of the relative estimates for Herring_NASSHS_fishery versus the other herring surveys:
#' # Herring_NASSHS:
#' valid <- dd$Survey == "Herring_NASSHS"
#' range_Herring_NASSHS <- range(dd$TSB_mean[valid])
#' range_Herring_NASSHS
#' diff(range_Herring_NASSHS) / mean(range_Herring_NASSHS)
#' 
#' # Herring_NASSHS_fishery:
#' valid <- dd$Survey == "Herring_NASSHS_fishery" & dd$Timing != "-30days"
#' range_Herring_NASSHS_fishery <- range(dd$TSB_mean[valid])
#' range_Herring_NASSHS_fishery
#' diff(range_Herring_NASSHS_fishery) / mean(range_Herring_NASSHS_fishery)
#' }
runAllPelfossPaper1 <- function(
	type="Parallel", run=type, survey="Herring_IESNS", dir="~/Projects/PELFOSS/delphi", # Paths and names
	speedFact=1, condPar=NULL, LcmMean=NULL, 
	year=c(2010, 2012), res=c(4, 10), dateshift=c(-30, 0, 30), reversed=c(FALSE, TRUE), # spec
	seed=1, nboot=100, radius=10, subset="all", cores=1, # Global settings
	...){
	
	# Create a data frame of the grid of specifications for the runs. Year and resolution are combined, whereas the result is expanded with dateshift and reversed:
	spec <- expand.grid(year=year, dateshift=dateshift, reversed=reversed)
	# Add resolution given year:
	spec$res <- res[match(spec$year, year)]
	spec
	
	# Save all the surveyPar lists:
	out <- list()
	
	# Run through the surveys:
	for(i in seq_along(survey)){
		surveyPar <- setSurveyParametersPaper1(run=run, dir=dir, survey=survey[i], type=type)
		surveyPar <- adjustSpeed(surveyPar, speedFact)
		runtimes <- runAll(spec=spec, run=run, dir=dir, surveyPar=surveyPar, seed=seed, nboot=nboot, radius=radius, subset=subset, condPar=condPar, LcmMean=LcmMean, cores=cores)
		out[[i]] <- list(surveyPar=surveyPar, spec=spec, seed=seed, nboot=nboot, radius=radius, subset=subset, cores=cores, runtimes=runtimes)
	}
	
	return(out)
}

#*********************************************
#*********************************************
#' Funciton for adjusting speed, e.g., to speed up the survey.
#'
#' @param x		The output from \code{\link{biomass2tsb}} or \code{\link{runPelfoss}}.
#' @param fact	A value below which NASC values are omitted when plotted as black dots with size proportional to NASC^NASCexp.
#'
adjustSpeed <- function(x, fact=1){
	x$knots <- x$knots * fact
	if(is.list(x$hours)){
		x$hours <- lapply(x$hours, "/", fact)
	}
	else{
		x$hours <- x$hours / fact
	}
	
	x
}



#*********************************************
#*********************************************
#' Function for reading all reports of one PELFOSS run.
#'
#' @param run							The name of the simulation run (which is three timings (-30, 0, +30 days), two directions (normal, reversed), and two resolutions (low, high), in total 12 individual runs).
#' @param dir							The directory of the PELFOSS folder structure (see \code{\link{getPelfossSkeleton}}).
#' @param seed		A single integer, or a list of seeds used in the funciton, including the following seeds: ('transect') for drawing transects using surveyPlanner(), ('trawl') for drawing trawl stations along the transects with probability as a funciton of the NASC (see \code{probfact}), and ('bootstrap') for getting the final estimate using runBootstrap() (see \code{nboot}).
#'
#' @export
#'
readAllPelfossReports <- function(run, dir, seed=NULL){
	allReports <- file.path(dir, "runs", run, "reports", "allReports.txt")
	out <- read.table(allReports, sep="\t", header=TRUE, stringsAsFactors=FALSE)
	# Select seed:
	if(length(seed)){
		out <- subset(out, Seed==seed)
	}
	out
}

#*********************************************
#*********************************************
#' Function for plotting the PELFOSS report form one PELFOSS run.
#'
#' @param run		The name of the simulation run (which is three timings (-30, 0, +30 days), two directions (normal, reversed), and two resolutions (low, high), in total 12 individual runs).
#' @param dir		The directory of the PELFOSS folder structure (see \code{\link{getPelfossSkeleton}}).
#' @param by		A string specifying which variable to plot. The options are "Inside" (only inside the strata with positive NASC), InsideAll (inside the strata used in the survey) and "Total" (all possible strata in the strata system).
#' @param seed		A single integer, or a list of seeds used in the funciton, including the following seeds: ('transect') for drawing transects using surveyPlanner(), ('trawl') for drawing trawl stations along the transects with probability as a funciton of the NASC (see \code{probfact}), and ('bootstrap') for getting the final estimate using runBootstrap() (see \code{nboot}).
#' @param ylim		The ylim of the plot, useful for setting the same ylim on several related plots.
#' @param save		The folder in which to put the plots in the "reports" folder of \code{dir}. If TRUE this is interpreted as "allPlots".
#' @param case_adj	A two element vector giving the adjustment of the labels of each survey (A, B, ...).
#'
#' @export
#' @import ggplot2
#'
plotAllPelfossReports <- function(run, dir, by="InsideAll", seed=NULL, ylim=NULL, save=TRUE, case_adj=c(0.15, 0.95), add.cv=FALSE, adds=list()){
	
	# Simple funciton for replacing an underscore by line break, used to print the cases in the plot:
	addline_format <- function(x, ...){
	    gsub('_', '\n', x, fixed=TRUE)
	}
	
	# Read the report:
	report <- readAllPelfossReports(run=run, dir=dir, seed=seed)

	# Add case number and convert to factor:
	report$CaseNr <- seq_len(nrow(report))
	report$Direction <- as.factor(report$Direction)
	report$Resolution <- as.factor(report$Resolution)
	report$Timing <- as.factor(report$Timing)
	report$Resolution_Direction <- as.factor(paste0("Res/year: ", report$Resolution, "/", report$Year, " ", "Dir: ", report$Direction))
	
	# Modify the data by the optional 'adds':
	for(i in seq_along(adds)){
		thisname <- names(adds[i])
		reportnames <- names(report)
		if(thisname %in% reportnames){
			if(is.function(adds[[i]])){
				report[[thisname]] <- adds[[i]](report[[thisname]])
			}
			else{
				report[[thisname]] <- adds[[i]]
			}
		}
		
	}

	# Define positions fo labels:
	sep <- which(diff(as.numeric(as.factor(report$Survey))) == 1)
	sep <- c(0, sep, nrow(report))
	sep <- sep + 0.5
	at <- sep[-1] - diff(sep)/2
	surveyNames <- unique(report$Survey)
	
	var5 <- paste0("TSB_5percentBy", by)
	varMean <- paste0("TSB_meanBy", by)
	var95 <- paste0("TSB_95percentBy", by)

	if(length(ylim) == 0){
		ylim <- c(0, max(report[[var95]]))
	}
	
	# With error bars:
	p <- ggplot(data=report, aes_string(x="CaseNr", y=varMean, group="Survey", colour="Timing")) + 
	    geom_errorbar(aes_string(ymin=var5, ymax=var95), width=0.4) +
	    #geom_line() +
	    geom_point(aes(shape=Resolution_Direction), size=4) + 
		geom_abline(intercept=1, slope=0) + 
		geom_vline(xintercept=sep, lty=2) + 
		#geom_text(data=data.frame(x=sep[-length(sep)], y=max(report$TSB_95percentByInsideAll) * 0.9), aes(x=x, y=y, label=paste("(", toupper(letters[seq_len(length(x))]), ")"))) + 
		annotate("text", x=sep[-length(sep)] + mean(diff(sep)) * case_adj[1], y=max(ylim) * case_adj[2], label=paste0("(", toupper(letters[seq_len(length(sep) - 1)]), ")"), size=5) + 
		scale_shape_manual(values=c(17, 19, 2, 1)) + labs(shape='Resolution and direction') + ylab("Estiamted divided by true biomass inside survey region") + xlab("") + 
		scale_x_continuous(breaks=at, labels=addline_format(surveyNames)) + ylim(ylim) + 
		theme(axis.title.y = element_text(size=15), axis.text.x = element_text(size=10, face=c("plain", "plain", "bold", "plain", "plain", "bold")), axis.text.y = element_text(size=15) )
	
	if(add.cv){
		cv_breaks <- unique(
			pretty(ylim), 
			pretty(c(0, max(report[, "Ab.Sum.cv"])))
		)
		p <- p + 
			geom_path(data = report, aes_string(x = "CaseNr", y = "Ab.Sum.cv"), color="black") + 
			geom_point(data = report, aes_string(x = "CaseNr", y = "Ab.Sum.cv", group = "Survey", colour = "Timing")) + 
			scale_y_continuous(sec.axis=sec_axis(~., name="CV", breaks=pretty(c(0, max(report[, "Ab.Sum.cv"])))), limits = ylim)
	}
	#if(length(ylim)){
	#	p <- p + ylim(ylim)
	#}
	#else{
	#	p <- p + ylim(0, NA)
	#}
	
	print(p)
	
	if(isTRUE(save)){
		save <- "allPlots"
	}
	if(length(save) && is.character(save)){
		plotFolder <- file.path(dir, "reports", save)
		pdfFolder <- file.path(plotFolder, "pdf")
		pngFolder <- file.path(plotFolder, "png")
		suppressWarnings(dir.create(pdfFolder, recursive=TRUE))
		suppressWarnings(dir.create(pngFolder, recursive=TRUE))
		ggsave(file.path(pdfFolder, paste0(run, ".pdf")))
		ggsave(file.path(pngFolder, paste0(run, ".png")))
	}
	
	return(p)
}


#*********************************************
#*********************************************
#' Function to generate all reports from a PELFOSS survey.
#'
#' @param run	The name of the simulation run (which is three timings (-30, 0, +30 days), two directions (normal, reversed), and two resolutions (low, high), in total 12 individual runs).
#' @param dir	The directory of the PELFOSS folder structure (see \code{\link{getPelfossSkeleton}}).
#' @param copy	Logical: If TRUE copy the files form the individual runs to the reports directory.
#'
#' @export
#'
reportAllPelfoss <- function(run, dir, copy=TRUE){
	report <- getAllTSB(run=run, dir=dir, copy=copy)
	report$reportsSelection <- report$reports[,c("Survey", "Resolution", "Year", "Direction", "Timing", "midDayOfSurvey", "Seed", "midDayOfSurvey", "Ab.Sum.mean", "Ab.Sum.cv", "ThSB_inside", "ThSB_total", "TSB_meanByInside")]
	report$par <- getRunPar(report)
	report
}
#*********************************************
#*********************************************
#' Add the parameters used in the PELFOSS simulation, stored in the list of reports 'x'.
#'
#' @param x	List of reports as returned from \code{getAllTSB}.
#'
getRunPar <- function(x){
	
	getSummary <- function(x, label=NULL){
		out <- as.data.frame(as.list(c(summary(x))))
		# remove any repots of NA:
		out <- out[!startsWith(names(out), "NA")]
		
		# Add the span:
		out$Span <- out$Max - out$Min
		
		if(length(label)){
			names(out) <- paste(label, names(out), sep="_")
		}
		out
	}
	
	getRunParOneRun <- function(x){
		data.frame(
			getSummary(x$pelfossOutput$transects$Transect$time_mid, label="time_mid"), 
			getSummary(x$pelfossOutput$transects$Transect$LcmMean, label="LcmMean"), 
			getSummary(x$pelfossOutput$transects$Transect$condParFactor, label="condParFactor"), 
			getSummary(x$pelfossOutput$transects$Transect$condParExponent, label="condParExponent")
		)
	}
	
	out <- lapply(x$data, getRunParOneRun)
	out <- do.call(rbind, out)
	out
}


# Function for plotting all days of the modelled fish distribution on top of one of the strata systems:
#'
#' @export
#'
plotAllDays <- function(survey, days=NULL, dir="~/Projects/PELFOSS/delphi", run="Parallel", year=2012, res=10, reversed=FALSE, dateshift=0, seed=1, centroid=c(0, 68), size=1, alpha=0.5, img="png", addTransects=FALSE, ...){

	##### get file paths: #####
	# Get files and dirs:
	files <- getPelfossPaths(dir, run=run, survey=survey, year=year, res=res, reversed=reversed, dateshift=dateshift, seed=seed)
	
	plotDir <- file.path(files$plotSuperIndDir, img)
	suppressWarnings(dir.create(plotDir, recursive=TRUE))
	##########
	
	
	##### Read data: #####
	### # Read the super individuals:	
	### superind <- readNorwecomSuperind(files$superind, centroid=centroid)
	
	# Read the transects:
	if(!isFALSE(addTransects)){
		load(files$projectRData)
	}
	##########
	
	
	##### Get day vector and output file paths: #####
	# Plot each day:
	ndays <- ncol(superind$longitude)
	ind <- sprintf(paste0("%0", nchar(ndays), "d"), seq_len(ndays))
	ncharInd <- nchar(ind[1])
	fileStam <- paste("Superind", run, year, res, sep="_")
	indWildcard <- paste0("%", ncharInd, "d")
	fileNameImages <- paste(paste(fileStam, indWildcard, sep="_"), img, sep=".")
	fileImages <- file.path(plotDir, fileNameImages)
	
	if(length(days) == 0){
		# Get the days correponding to the times of the transects:
		if(!isFALSE(addTransects)){
			day <- findInterval(pelfossOutput$transects$Transect$time_mid, superind$time)
			dayInd <- split(seq_along(day), day)
			days <- seq(min(day), max(day))
		}
		else{
			days <- seq_len(ndays)
		}	
	}
	
	# Get the file to which to save the plot:
	fileName <- paste(paste(fileStam, ind[days], sep="_"), img, sep=".")
	file <- file.path(plotDir, fileName)
	##########
	
	
	##### Plot the data: #####
	### # Run the stratum plot:
	### if(isFALSE(addTransects)){
	### 	p1 <- plotStratum(files$stratum, quiet=TRUE)
	### }
	
	browser()
	message("Here we should implement either biomass field or superindividuals (or both?) in plotPelfoss")
	# Plot each day:
	for(i in seq_along(days)){
		
		if(isTRUE(addTransects)){
			p1 <- Rstox::plotStratum(
				pelfossOutput$transects, 
				quiet = TRUE, 
				transectcol = "black", 
				timerange = c(
					min(pelfossOutput$transects$Transect$time_mid), 
					superind$time[days[i]]
				)
			)
		}
		else if(is.character(addTransects) && identical(tolower(addTransects), "nasc")){
			p1 <- plotPelfoss(
				pelfossOutput, 
				quiet = TRUE, 
				plot = c("nasc", "trawl"), 
				timerange = c(
					min(pelfossOutput$transects$Transect$time_mid), 
					superind$time[days[i]]
				), 
				...
			)$p
		}
			
		# Plot only if there are valid weights:
		valid <- !is.na(superind$weight[, days[i]])
		if(max(superind$weight[valid, days[i]], na.rm=TRUE)>0){
			
			# Get the data frame to plot from, by extracting the current day and the valid superindividuals:
			dat <- data.frame(
				longitude = superind$longitude[valid, days[i]], 
				latitude = superind$latitude[valid, days[i]], 
				age = superind$age[valid, days[i]], 
				weight = superind$weight[valid, days[i]], 
				inumb = superind$inumb[valid, days[i]], 
				mass = superind$inumb[valid, days[i]] * superind$weight[valid, days[i]]
			)
			
			# Plot the super individuals:
			p2 <- p1 + geom_point(data=dat, aes(x=longitude, y=latitude, color=age, fill=mass), alpha=0.5) + 
				#scale_colour_gradientn(colours=rainbow(4, v=c(0.6, 0.7, 0.8, 1))) + 
				scale_fill_gradientn(colours=rainbow(4, v=c(0.6, 0.7, 0.8, 1))) + 
				###scale_size_discrete(range = c(0.1,2)) + 
				#scale_size_continuous(range = c(0.1,2)) + 
				#scale_size(range = c(0.1,2)) + 
				#scale_size_area() + 
				ggtitle(superind$time[days[i]]
				)	
		}
		else{
			p2 <- p1
		}
		
		p2 <- moveGeom(p2, geom="segment", move=100, ind=1)
		suppressWarnings(p2 <- moveGeom(p2, geom="point", move=100, ind=1))
		p2 <- moveGeom(p2, geom="text", move=100, ind=1)
		
		# Save the plot:
		ggsave(file[i], plot=p2, device=img)
	}
	##########

	
	return(file)
}

# Function for generating a movie of all the fish distribution plots:
#'
#' @export
#'
movieAllDays <- function(survey, days=NULL, dir="~/Projects/PELFOSS/delphi", run="Parallel", year=2012, res=10, reversed=FALSE, dateshift=0, seed=1, img="png", mov="mp4", fpsIn=20, fpsOut=fpsIn, bitrate=1e7){

	# Get files and dirs:
	files <- getPelfossPaths(dir, run=run, survey=survey, year=year, res=res, reversed=reversed, dateshift=dateshift, seed=seed)
	
	# Get the plot files:
	plotDir <- file.path(files$plotSuperIndDir, img)
	suppressWarnings(dir.create(plotDir, recursive=TRUE))
	fileImages <- list.files(plotDir, full.names=TRUE)
	if(length(days)){
		fileImages <- fileImages[days]
	}
	
	# Get the path to the movie:
	if(length(mov)){
		movieDir <- file.path(files$plotSuperIndDir, mov)
		suppressWarnings(dir.create(movieDir, recursive=TRUE))
		fileNameVideo <- paste(fileStam, mov, sep=".")
		fileVideo <- file.path(movieDir, fileNameVideo)
	}
	
	
	##### Create movie: #####
	if(length(mov)){
		cmd <- paste(
			"ffmpeg", 
			"-y", 
			"-r", fpsIn, 
			"-start_number", min(days), 
			"-i", fileImages, 
			"-c:v libx264", 
			"-pix_fmt yuv420p", 
			"-r", fpsOut, 
			"-b", bitrate, 
			fileVideo
		)
		
		print(cmd)
		
		system(cmd)
	}
	##########
	
	return(fileVideo)
}

# Function to rearrange layes in a ggplot:
moveGeom <- function(p, geom="segment", move=1, ind=1){
	geoms <- sapply(p$layers, function(x) class(x$geom)[1])
	ord <- seq_along(geoms)
	hit <- grep(geom, geoms, ignore.case=TRUE)
	if(length(hit) > 1){
		warning("There are ", length(hit), " layers matching the geom ", geom, ". Use 'ind' to select the requested layer, and view layers with p$layers")
		hit <- hit[ind]
	}
	ord[hit] <- ord[hit] + move + 0.1
	ord <- order(ord)
	p$layers <- p$layers[ord]
	p
}
