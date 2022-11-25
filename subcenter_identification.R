############################################################################################
#####################################################################################
#####  SUBCENTER IDENTIFICATION ALGORITHM companion to:
####  Ban, Arnott, Macdonald
###   "Identifying Employment Subcenters: The Method of Exponentially
##   Declining Cutoffs"
#
#
# R Code which reads in a shapefile of a metropolitan area with employment data for each
# geography (i.e. Census Tract; TAZ). This data is used in identifying employment subcenters
# based on the accompanying paper’s Density Gradient Cutoff (DGC) methodology. Special cases
# of this function determine subcenters based on the Giuliano--Small (GS) or Exponentially
# Declining Cutoff (EDC) methodologies.
#
# INPUTS:
# [shapefile] -- SpatialPolygonsDataFrame
#  ".shp" file of metropolitan area with employment by zone
# [output]  -- character (default: working directory)
#  Output folder location where to save files
# [employment] -- character
#  Name of the employment variable
# [location] -- character (vector)
#  Vector of the names of ID, County, or Location variables that should be kept with the
#  shapefile (additional variables are then deleted)
# [D]     -- numeric
#  User specified cutoff employment density threshold at metropolitan center.
#  Should be specified in the same units as the algorithm (if units=="metric" then D
#  should be specified as employees/hectare; if units=="imperial" then D should be
#  specified as employees/acre)
# [E]     -- numeric
#  User specified cutoff total employment threshold at metropolitan center.
# [type]   -- ("DGC" (default) or "EDC")
#  If =="EDC" then employs the method of Exponentially Declining Cutoff (EDC)
#  If =="DGC" then employs the method of Density Gradient Cutoff (DGC)
# [alpha]   -- numeric (default: ln2/40)
#  Used only under the "EDC" method. User specified value of the cutoff gradient.
# [theta]   -- numeric in [0,1]
#  Weight used in absolute and relative density and employment cutoffs (theta==0 => all
#  weight on absolute density; equivalent to G.S. methodology)
# [gamma]   -- numeric (optional)
#  Used only under the "DGC" method and can be manually specified.
#  Default calculates employment density gradient estimated by: ln(D_z) = c - gamma∗x_z
#  Where "D" and "x" are zone level employment densities and distance to CBD.
# [units]   -- ("metric" (default) or "imperial")
#  If "metric" then analysis is done with distance in kilometers and area in hectares.
#  If "imperial then analysis done with distance in miles and area in acres.
# [generate] -- TRUE | FALSE (default: TRUE)
#  Set to FALSE if the output files should not be saved - subcenter results will only be
#  available in the R environment.
#
# OUTPUTS:
# subcenters[[1]]
#  ".shp" file with variable "subcenter" which identifies candidate subcenters (= 1)
#  and full subcenters (= 2)
# subcenters[[2]]
#  ".csv" file labeling each subcenter starting from the highest subcenter density
#  with respective employment, area, density and zone information.
# subcenters[[3]]
#  Value of the density cutoff gradient ([alpha] or [gamma]).
 
x <- c("sp", "spdep", "rgdal", "geosphere", "rgeos", "maptools")
uninstalled <- x[!(x %in% installed.packages()[,"Package"])]
if(length(uninstalled)) install.packages(uninstalled)
lapply(x, library, character.only = TRUE)
rm(x, uninstalled)
 
############################################################################################
############################################################################################
subcenters <- function(shapefile, output=getwd(), employment, location, type="DGC", D, E,
 alpha, theta, gamma=NULL, units="metric", generate=TRUE){
 
SESSION <- paste0(unlist(strsplit(as.character(Sys.Date()), "[-]"))[[3]],
 paste0(unlist(strsplit(as.character(Sys.Date()), "[-]"))[[2]],
 substr(unlist(strsplit(as.character(Sys.Date()), "[-]"))[[1]], 3, 4)))
 
# RETREIVE INTERNAL POLYGON ID
shapefile@data$IDpg <- as.factor(sapply(slot(shapefile, "polygons"),
 function(x) slot(x, "ID")))
 
# CHANGE THE NAME OF THE VARIABLE THE USER IDENTIFIES AS EMPLOYMENT
names(shapefile@data)[names(shapefile@data)==employment] <- "employment"
 
shapefile@data$area <- rep(0, length(shapefile@data[,1]))
for(i in 1:length(shapefile@data[,1])){
 shapefile@data$area[i] <- areaPolygon(lapply(slot(shapefile, "polygons"),
  function(x) lapply(slot(x, "Polygons"),
  function(y) slot(y, "coords")))[[i]][[1]])
}
 
if(units=="metric"){
 shapefile@data$area <- shapefile@data$area∗0.0001
} else if(units=="imperial"){
 shapefile@data$area <- shapefile@data$area∗0.000247105
} else {
 shapefile@data$area <- NA
}
 
# CALCULATE EMPLOYMENT DENSITY USED TO DETERMINE THE CBD
shapefile@data$density <- shapefile@data$employment/shapefile@data$area
shapefile@data$density <- ifelse(is.na(shapefile@data$density), 0, shapefile@data$density)
 
# CBD IS DEFINED AS THE CENSUS TRACT WITH THE MAXIMUM EMPLOYMENT DENSITY
# DISTANCE IS CALCULATED BETWEEN CENSUS CENTROIDS AND CALCULATED IN MILES
shapefile@data$distance <- 3963.0∗acos(sin(coordinates(shapefile[shapefile@data$density==
 max(shapefile@data$density),])[2]/57.2958)∗sin(coordinates(shapefile)[,2]/57.2958) +
 cos(coordinates(shapefile[shapefile@data$density==
 max(shapefile@data$density),])[2]/57.2958)∗cos(coordinates(shapefile)[,2]/57.2958)∗
 cos(coordinates(shapefile)[,1]/57.2958 -
 coordinates(shapefile[shapefile@data$density==max(shapefile@data$density),])[1]/57.2958))
 
if(units=="metric"){
 shapefile@data$distance <- shapefile@data$distance∗1.60934
} else if(units=="imperial"){
 shapefile@data$distance <- shapefile@data$distance
} else {
 shapefile@data$distance <- NA
}
 
# GAMMA IS EITHER LOG(2)/40 IN THE NEG. EXPONENTIAL MODEL OR THE EMPLOYMENT GRADIENT IN DGC
# MODELS
 if(type=="EDC") {
  gradient <- alpha
  theta <- 1
 } else if(type=="DGC" & !is.null(gamma)) {
  gradient <- gamma
 } else if(type=="DGC" & is.null(gamma)) {
  gradient <- abs(as.numeric(as.character(lm(log(ifelse(shapefile@data$density > 0,
   shapefile@data$density, 0.5)) ~ shapefile@data$distance)$coefficients[2])))
 }
 
# KEEP ONLY THE VARIABLES NEEDED FOR SC IDENTIFICATION AND IMPORTANT LOCATION AND
# IDENTIFYING VARIABLES
 keeps <- c("IDpg", "employment","area", "distance", "density", location)
 shapefile <- shapefile[ , (colnames(shapefile@data) %in% keeps)]
 rm(keeps)
 
 shapefile@data$longitude <- coordinates(shapefile)[,1]
 shapefile@data$latitude <- coordinates(shapefile)[,2]
 
# FOR EACH INDIVIDUAL CENSUS TRACT DETERMINE WHETHER THE DENSITY MEETS THE EXPONENTIALLY
# DECREASING DENSITY THRESHOLD. IF YES THEN IT CAN BE CONSIDERED A CANDIDATE SUBCENTER (==1)
 shapefile@data$Dcutoff <- D∗exp(-theta∗gradient∗shapefile@data$distance)
 shapefile@data$subcenter <- ifelse(shapefile@data$density > shapefile@data$Dcutoff, 1, 0)
 
# THE TOTAL EMPLOYMENT CUTOFF MUST BE APPLIED TO CONTIGUOUS TRACTS OF CANDIDATE SUBCENTERS
# TO IDENTIFY CONTIGUOUS TRACTS, WE FIRST DEFINE THE ADJACENCY MATRIX WHICH DETERMINES WHAT
# CENSUS TRACTS ARE ADJACENT TO EACH OTHER AND FURTHER DEFINES ALL CONTIGUOUS BLOCKS OF
# CANDIDATE SUBCENTERS.
 candidates <- shapefile[shapefile@data$subcenter==1,]
 adjacency <- gTouches(candidates, returnDense=TRUE, byid=TRUE)
 adjacency[adjacency == FALSE] <- 0
 adjacency[adjacency ==TRUE] <- 1
 
 colnames(adjacency) <- rownames(adjacency)
 n = nrow(adjacency)
 amat <- matrix(0,nrow=n,ncol=n)
 amat[row(amat)==col(amat)] <- 1
 colnames(amat) <- colnames(adjacency)
 rownames(amat) <- rownames(adjacency)
 bmat <- adjacency
 wmat1 <- adjacency
 newnum = sum(bmat)
 cnt = 1
 while (newnum > 0) {
  amat <- amat+bmat
  wmat2 <- wmat1%∗%adjacency
  bmat <- ifelse(wmat2 > 0 & amat==0, 1, 0)
  wmat1 <- wmat2
  newnum = sum(bmat)
  cnt = cnt+1
 }
 
 Ez <- amat%∗%diag(candidates@data$employment)
 Xz <- amat%∗%diag(candidates@data$distance)
 Az <- amat%∗%diag(candidates@data$area)
 EzXz <- Ez∗Xz
 
 IDz <- amat∗as.numeric(colnames(amat))
 IDz[IDz==0] <- 999999999999999999
 IDz <- apply(cbind(apply(IDz, 1, min), apply(IDz, 2, min)), 1, min)
 
# FOR EACH CANDIDATE SUBCENTER WE CALCULATE THE TOTAL EMPLOYMENT, DISTANCE WEIGHTED
# EMPLOYMENT AND AREA FROM THE RESPECTIVE BROADER CONTIGUOUS GROUP OF SUBCENTERS.
 candidates@data$SCemployment <- rowSums(Ez, na.rm = FALSE, dims = 1)
 candidates@data$SCdistance <- rowSums(EzXz, na.rm = FALSE, dims = 1)/candidates@data$SCemployment
 candidates@data$SCarea <- rowSums(Az, na.rm = FALSE, dims = 1)
 candidates@data$SCdensity <- candidates@data$SCemployment/candidates@data$SCarea
 
# EACH CANDIDATE SUBSCENTER IS COMPARED AGAINSED THE TOTAL EMPLOYMENT CUTOFF FOR THE BROADER
# GROUP OF CONTIGUOUS TRACTS.
 candidates@data$Ecutoff <- E∗exp(-theta∗gradient∗candidates@data$SCdistance)
 candidates@data$subcenter <- ifelse(candidates@data$SCemployment > candidates@data$Ecutoff,
  candidates@data$subcenter + 1, candidates@data$subcenter)
 
 candidates@data$SCidz <- IDz
 
 SCs <- candidates@data[,c("subcenter", "SCemployment", "SCdistance", "SCarea",
  "SCdensity", "SCidz")]
 SCs <- merge(SCs, as.data.frame(table(SCs$SCidz)), by.x="SCidz", by.y="Var1", all.x=T)
 colnames(SCs)[colnames(SCs)=="Freq"] <- "Nzones"
 SCs <- SCs[SCs$subcenter==2,]
 SCs <- unique(SCs)
 SCs <- as.data.frame(SCs[order(SCs$SCdensity, decreasing=T),])
 rownames(SCs) <- NULL
 SCs$subcenter <- NULL
 SCs <- as.data.frame(cbind(SCs, tolower(as.roman(1:length(SCs[,1])))))
 colnames(SCs) <- c("SCidz", "Employment", "DistanceCBD", "Area", "Density", "Nzones", "SCID")
 
 candidates@data <- merge(candidates@data, SCs[,c("SCidz", "SCID")], by.x="SCidz",
  by="SCidz", all.x=T, sort=F)
 SCs$SCidz <- NULL
 keeps <- c("IDpg", "subcenter","SCemployment", "SCdistance", "SCarea", "SCdensity",
  "Ecutoff", "SCID")
 candidates <- candidates[ , (colnames(candidates@data) %in% keeps)]
 rm(keeps)
 
 shapefile@data$subcenter <- NULL
 
 shapefile <- sp::merge(shapefile, candidates, by.x="IDpg", by.y="IDpg", all.x=T, sort=F)
 rownames(shapefile) <- rownames(shapefile)
 shapefile@data$subcenter[is.na(shapefile@data$subcenter)] <- 0
 
 SCs <- SCs[,c("SCID", "Density", "Employment", "Area", "DistanceCBD", "Nzones")]
 
# IF GENERATE=TRUE THEN BOTH A SHAPEFILE AND CSV FILE ARE EXPORTED TO THE OUTPUT LOCATION
 if(generate){
  suppressWarnings(writeOGR(shapefile, dsn = output, layer = paste0("subcenter", type,
  ifelse(type=="EDC", strsplit(as.character(as.character(round(alpha, 3))), "[.]")[[1]][2],
  gsub(".", "", as.character(theta), fixed=T)), "_", SESSION), driver="ESRI Shapefile",
  check_exists=TRUE, overwrite_layer=TRUE))
  write.csv(SCs, paste0(output, "/subcenter", type, ifelse(type=="EDC",
  strsplit(as.character(as.character(round(alpha, 3))), "[.]")[[1]][2],
  gsub(".", "", as.character(theta), fixed=T)), "_", SESSION, ".csv"))
 }
 
# THERE ARE THREE OUTPUTS WHICH ARE SAVED IN R: THE SHAPEFILE WITH CANDIDATE SUBCENTERS == 1
# AND FULL SUBCENTERS == 2; A TABLE WITH EACH FULL SUBCENTER IDENTIFIED BY ROMAN NUMERAL WITH
# RESPECTIVE DATA ON TOTAL EMPLOYMENT, DENSITY AND AREA; THE VALUE OF THE GRADIENT ESTIMATED
 
 out <- list(shapefile, SCs, gradient)
 names(out) <- c("SHAPEFILE", "CSV", "Gradient")
 return(out)
}
