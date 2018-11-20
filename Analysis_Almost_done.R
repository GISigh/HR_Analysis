# Functions from wildLifeDI

# ---- roxygen documentation ----
#' @title Volume contour from Raster
#'
#' @description
#'   Compute a percent volume contour polygon from a raster UD.
#' @details
#'   The volras function is a simpler version of the getvolumeUD function from the package \code{adehabitatHR} developed by C. Calenge. It allows the output to be a 'raster looking' polygon (i.e., the cells that are within the UD) or a simplified (smoothed) polygon.
#'   
#' @param x a \code{RasterLayer}
#' @param percent a percent value to get the volume contour, e.g., 95. Note: This is a simple function and only accepts one value at a time.
#' @param simplify (logical; default = TRUE) whether or not to simplify the output home range polygon using \code{gSimplify} from \code{rgeos} with a tolerance value of 1.5 times the spatial resolution of the UD. 
#' 
#' @return
#'   A \code{SpatialPolygonsDataFrame}.
#'
#' @seealso fbtgUD, rspUD, tgkde
#' @examples
#' data(m3)
#' ud <- tgkde(m3,disfun='inv',method='vanderWatt')
#' raster::plot(ud)
#' hr <- volras(ud,95)
#' sp::plot(hr,add=TRUE)
#' 
#' @export
#
# ---- End of roxygen documentation ----
volras <- function(x,percent=95,simplify=TRUE){
  
  x[is.na(x)] <- 0
  pfs <- proj4string(x)
  
  ## standardize it so that the total volume is 1 over the area
  v <- as.vector(values(x))
  index<-1:length(v)
  vord<-v[order(v, decreasing=TRUE)]
  vsu<-cumsum(vord)
  
  cont <- which(vsu > (percent/100)*max(vsu))[1]
  cont.lev <- vord[cont]
  
  #Get all the cells above the cont.lev and make 1
  m <- c(0, cont.lev, 0, cont.lev,cellStats(x,'max'),1)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  x2 <- reclassify(x, rclmat)
  
  #Convert to polygon and simplify if desired
  hr <- rasterToPolygons(x2,fun=function(x){x==1},n=4,na.rm=TRUE,dissolve=TRUE)
  if(simplify){
    hr <- gSimplify(hr,tol=res(x2)[1]*1.5)
  }
  
  #return the volume contour polygon
  return(hr)
}
# ---- roxygen documentation ----
#' @title Dynamic PPA Measure of Animal Space Use
#'
#' @description
#'   The function \code{dynppa} computes the (dynamic) PPA measure of the accessibility space of an animal. 
#'   The PPA method can be thought as an alternative view on the home range; one that explicitly considers the
#'   spatial and temporal constraints on movement given known telemetry fixes, and a (dynamic) measure of maximum
#'   mobility - termed Vmax. The PPA method incorporates dynamic behaviour into the calculation of the vmax parameter 
#'   used to delineate the original version of the PPA method, but the original method is still an option here.
#' @details
#'   The function \code{dyn.ppa} represents an extension to an existing PPA method (Long and Nelson, 2012). 
#'   Dynamic calculation of the PPA method improves upon the original version by flexibly modelling the vmax 
#'   parameter according to wildlife behaviour. See the function \code{dyn.vmax} for more information on how 
#'   to incorporate dynamic behaviour into the vmax parameter estimation. 
#'   
#' @param traj an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{
#'    help(ltraj)}.
#' @param tol parameter used to filter out those segments where the time between fixes is overly 
#'    large (often due to irregular sampling or missing fixes); which leads to an overestimation of the 
#'    activity space via the PPA method. Default is the maximum sampling interval from \code{traj}.
#' @param dissolve (logical) whether or not to dissolve output elliplse polygons to create a single
#'    output polygon, or keep the individual segment PPA ellipses. Default = TRUE.
#' @param proj4string a string object containing the projection information to be passed included in the output 
#'    \code{SpatialPolygonsDataFrame} object. For more information see the \code{CRS-class} in the packages
#'    \code{sp} and \code{rgdal}. Default is \code{NA}.
#' @param ePoints number of vertices used to construct each PPA ellipse. More points will necessarily provide
#'    a more detailed ellipse shape, but will slow computation; default is 360.
#' @param ... additional parameters to be passed to the function \code{dynvmax}. For example, should include
#'    \code{method} and/or \code{dynamic} parameters, see the documentation for \code{dynvmax} for more detailed
#'    information on what to include here.
#'    
#' @return
#'   This function returns a \code{SpatialPolygonsDataFrame} representing the dynamic PPA measure of the accessibility
#'   space of an individual animal.
#'
#' @references
#'   Long, JA, Nelson, TA. (2012) Time geography and wildlife home range delineation. \emph{Journal of Wildlife
#'   Management}, 76(2):407-413.\cr \cr
#'   Long, JA, Nelson, TA. (2014) Home range and habitat analysis using dynamic time geography. \emph{Journal of 
#'   Wildlife Management}. Accepted: 2014-12-03.\cr
#'   
# @keywords 
#' @seealso dynvmax
#' @examples 
#' data(m3)
#' ppa1 <- dynppa(m3,method='vanderWatt')
#' ppa2 <- dynppa(m3,method='vanderWatt',dynamic='focal')
#' sp::plot(ppa1)
#' sp::plot(ppa2,add=TRUE,border='red')
#' 
#' @export
#
# ---- End of roxygen documentation ----

dynppa <- function(traj, tol=max(ld(traj)$dt,na.rm=TRUE),dissolve=TRUE,proj4string=CRS(as.character(NA)),ePoints=360, ...){
  #check to see if ltraj is not a type II
  if (attributes(traj)$typeII == FALSE) {stop("The trajectory object is not a TypeII ltraj object")}
  
  
  #Append the local Vmax values to the trajectory using the LocalVmax function
  traj <- dynvmax(traj, ...)
  
  #convert to dataframe
  trDF <- ld(traj)
  n <- dim(trDF)[1]
  
  #loop through trajectory dataset and calculate PPA
  polyList <- vector('list',length=(n-1))
  for (i in 1:(n-1)){ #replace with function and use lapply
    #check to see if the local Vmax value exists
    if (!is.na(trDF$dynVmax[i])){      
      #check to see if the time interval is below the tolerance level
      if (trDF$dt[i] <= tol){
        cpX <- (trDF$x[i] + trDF$x[i+1])/2
        cpY <- (trDF$y[i] + trDF$y[i+1])/2
        #check to see if the angle is NA, if it is make it zero (no movement)
        if (is.na(trDF$abs.angle[i]) == TRUE) {
          thetaRot <- 0
        } else {thetaRot <- trDF$abs.angle[i]}
        
        c <- trDF$dist[i]/2
        a <- (trDF$dt[i]*trDF$dynVmax[i])/2
        b <- sqrt((a^2)-(c^2))
        
        #Compute the ellipse
        polyList[[i]] <- ppaEllipse(cpX,cpY,a,b,thetaRot,ePoints)
      }
    }
  }
  
  #---------------------------------------------
  # Create spatial polygons from the output
  #---------------------------------------------
  ind.poly <- which(lapply(polyList, is.null) == FALSE)
  polyList <- polyList[ind.poly]
  polyList <- lapply(polyList,list)
  tempPoly <- mapply(Polygons, polyList, ind.poly)
  data <- data.frame(ind=ind.poly, date=trDF$date[ind.poly])
  spPoly <- SpatialPolygonsDataFrame(SpatialPolygons(tempPoly,proj4string=proj4string),data,match.ID=FALSE)
  
  #"union" multiple spatial polygons stored in single sp object
  #(dissolve in GIS) 
  if (dissolve == TRUE){
    spUnion <- gUnaryUnion(spPoly)
    spUnion <- SpatialPolygonsDataFrame(spUnion,data=data.frame(id=1:length(spUnion)))
    return(spUnion)
  } else {
    return(spPoly)
  }
}
#End of Function
#===============================================================================

# ---- roxygen documentation ----
#' @title Dynamic Calculation of the Vmax Parameter
#'
#' @description
#'   The function \code{dynvmax} computes a dynamic version of the Vmax parameter for the PPA method. It can be used to incorporate changes in animal movement behaviour into the PPA method caluculation to better model that area accessible to an individual animal given the set of known telemetry locations in space and time.
#'   
#' @details
#'   The function \code{dynvmax} represents an intermediary function used to extend and improve upon an existing PPA home range method (Long and Nelson, 2012) as described in the paper (Long and Nelson, 2014). Four options are available for computing the vmax parameter dynamically and are passed into the \code{dynvmax} function using \code{dynamic} option.\cr 
#'  \cr 1) \code{NA} -- if \code{dynamic} = \code{'NA'} (the default) the function estimates the original, 
#'          non-dynamic estimate of Vmax which is a global estimate, as per Long & Nelson (2012). 
#'  \cr 2) \code{focal} -- a moving window approach whereby a window of size \code{w} is moved along the 
#'        trajectory and vmax computed dynamically within each window and assigned to the central segment.
#'  \cr 3) \code{cumulative} -- A moving window of size \code{w} is again used, only in this case the value 
#'        is assigned to the end segment. This represents the vmax calculation of the previous \code{w} segments.
#'  \cr 4) \code{class} -- A priori analysis (e.g., obtained via state-space models, or from expert knowledge) 
#'        is used to identify discrete behavioural states in the telemetry data and these stored in a column 
#'        which is then passed into the function.\cr\cr 
#'The \code{class} method is the preferred choice, as it allows the use of more sophisticated models for identifying behavioural shifts in telemetry data where we would expect to see clear differences in the Vmax parameter based on changing movement behaviour.\cr\cr 
#'The use of the \code{'focal'} or \code{'cumulative'} dynamic methods uses a moving window approach, which is sensitive to edge effects at the initial and ending times of the trajectory. Thus, the dynamic Vmax parameter is only computed for those segments that have a valid window and the dataset is shrunk by \code{w-1} segments.
#'   
#' @param traj an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{
#'    help(ltraj)}.
#' @param dynamic one of \code{'NA'}, \code{'focal'}, \code{'cumulative'}, or \code{'class'}; which signifies 
#'    whether or how to dynamically compute the Vmax parameter. See \bold{Details} for more information on 
#'    each of the choices.
#' @param method method for computing the Vmax parameter dynamically; can be one of several options:
#'    -- \code{"Robson"} for the Robson & Whitlock (1964) method,\cr
#'    -- \code{"RobsonLL"} for the R & W (1964) lower \eqn{(1-\alpha)*100\%} C.I. limit, \cr
#'    -- \code{"RobsonUL"} for the R & W (1964) upper \eqn{(1-\alpha)*100\%} C.I. limit, \cr
#'    -- \code{"vanderWatt"} for the van der Watt (1980) method, \cr
#'    -- \code{"vanderWattLL"} for the van der Watt (1980) lower \eqn{(1-\alpha)*100\%} C.I. limit, \cr
#'    -- \code{"vanderWattUL"} for the van der Watt (1980) upper \eqn{(1-\alpha)*100\%} C.I. limit. \cr
#' @param w (optional) window size (only used with \code{dynamic = 'focal'} or \code{'cumulative'}). 
#' @param class.col (optional) character indicating the name of the column in the \code{infolocs} dataframe
#'    of \code{traj} containing the categorized behavioural states of the animal (which can be stored as
#'    a character or numeric column).
#' @param k (optional) value for the \emph{k} parameter in the van der watt (1980) method; default is 5.
#' @param alpha (optional) value for the \eqn{\alpha} parameter if using upper or lower C.I. methods; default is 0.05.
#' @param manualVmax (optional) Character name of column in \code{traj} storing user input column of vmax values (typically call the column dynVmax).
#' @param vmaxtrunc (optional) due to irregular sampling intervals, or errors in GPS location, or other
#'    effects, the calculation of the vmax parameter through the statistical methods outlined above can be
#'    heavily influenced by high outliers. Thus, it may be useful to exclude those segments from calculation 
#'    of the dynamic Vmax parameter. Default is \code{NA}.
#' 
#' @return
#'   This function returns the original \code{traj} object with a new column -- \code{dynVmax} in the \code{infolocs} dataframe
#'   containing the dynamic vmax parameter for each trajectory segment.
#'
#' @references
#'   Long, JA, Nelson, TA. (2012) Time geography and wildlife home range delineation. \emph{Journal 
#'   of Wildlife Management}. 76(2):407-413.\cr\cr
#'   Long, JA, Nelson, TA. (2015) Home range and habitat analysis using dynamic time geography. \emph{Journal of -
#'   Wildlife Management}. 79(3):481-490.\cr\cr
#'   Robson, DS, Whitlock, JH. (1964) Estimation of a truncation point. \emph{Biometrika}
#'   51:33-39.\cr\cr
#'   van der Watt, P. (1980) A note on estimation bounds of random variables. \emph{Biometrika}
#'   67(3):712-714.\cr
# @keywords 
#' @seealso dynppa
#' @examples
#' data(m3)
#' m3R <- dynvmax(m3,dynamic='focal',method='Robson')
#' m3V <- dynvmax(m3,dynamic='focal',method='vanderWatt')
#' m3c <- dynvmax(m3,dynamic='cumulative')
#' 
#' @export
#
# ---- End of roxygen documentation ----

dynvmax <- function(traj,dynamic="NA",w=9,class.col="dt",method="Robson",k=5,alpha=0.05,manualVmax=NA,vmaxtrunc=NA){
  #check to see if ltraj is not a type II
  if (attributes(traj)$typeII == FALSE) {stop("The trajectory object is not a TypeII ltraj object")}
  
  #check to see if k < w iff the vanderWatt methods are used
  #if (method == "vanderWatt" || method == "vanderWattLL" || method == "vanderWattUL"){
  #  if(k > w) {stop("The k parameter is greater than the w parameter.")}}
  
  #store trajectory data as a dataframe for indexing
  trDF <- ld(traj)
  n <- dim(trDF)[1]
  
  #compute segment velocities
  trDF$vi <- c(trDF$dist[1:(n-1)] / trDF$dt[1:(n-1)],NA)
  
  if (is.na(vmaxtrunc)){
    vmaxtrunc <- max(trDF$vi,na.rm=TRUE)
  }  
  
  trDF$dynVmax <- NA
  
  #original vmax calculation
  if (dynamic == "NA"){
    v.temp <- trDF$vi
    v.temp <- v.temp[which(v.temp <= vmaxtrunc)]
    v <- sort(v.temp,decreasing=T)
    trDF$dynVmax <- switch(method,
                           #Robson & Whitlock (1964) method --- main estimate and UL come
                           # directly from the paper, however the LL was inferred
                           Robson = v[1] + (v[1]-v[2]),
                           RobsonLL = v[1] + ((1-(1-alpha))/(1-alpha))*(v[1]-v[2]),
                           RobsonUL = v[1] + ((1-alpha)/alpha)*(v[1]-v[2]),
                           #van der Watt (1980) method --- all come directly from the paper.
                           vanderWatt = (((k+2)/(k+1))*v[1]) - ((1/(k+1))*v[k]),
                           vanderWattLL = v[1] + ((1-(1-alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
                           vanderWattUL = v[1] + ((1-(alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
                           stop(paste("The vMax method is spelt incorrectly: ",method)))
  }
  
  #Focal local Vmax
  if (dynamic == "focal"){
    w1 <- trunc(w/2)
    for (i in (w1+1):(n-1-w1)){
      v.temp <- trDF$vi[(i-w1):(i+w1)]
      v.temp <- v.temp[which(v.temp <= vmaxtrunc)]
      v <- sort(v.temp,decreasing=T)
      
      trDF$dynVmax[i] = switch(method,
                               #Robson & Whitlock (1964) method --- main estimate and UL come
                               # directly from the paper, however the LL was inferred
                               Robson = v[1] + (v[1]-v[2]),
                               RobsonLL = v[1] + ((1-(1-alpha))/(1-alpha))*(v[1]-v[2]),
                               RobsonUL = v[1] + ((1-alpha)/alpha)*(v[1]-v[2]),
                               #van der Watt (1980) method --- all come directly from the paper.
                               vanderWatt = (((k+2)/(k+1))*v[1]) - ((1/(k+1))*v[k]),
                               vanderWattLL = v[1] + ((1-(1-alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
                               vanderWattUL = v[1] + ((1-(alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
                               stop(paste("The vMax method is spelt incorrectly: ",method)))
    }
  }
  
  #cumulative local Vmax
  if (dynamic == "cumulative"){
    for (i in w:n){
      v.temp <- trDF$vi[(i-w+1):i]
      v.temp <- v.temp[which(v.temp <= vmaxtrunc)]
      v <- sort(v.temp,decreasing=T)
      
      trDF$dynVmax[i] = switch(method,
                               #Robson & Whitlock (1964) method --- main estimate and UL come
                               # directly from the paper, however the LL was inferred
                               Robson = v[1] + (v[1]-v[2]),
                               RobsonLL = v[1] + ((1-(1-alpha))/(1-alpha))*(v[1]-v[2]),
                               RobsonUL = v[1] + ((1-alpha)/alpha)*(v[1]-v[2]),
                               #van der Watt (1980) method --- all come directly from the paper.
                               vanderWatt = (((k+2)/(k+1))*v[1]) - ((1/(k+1))*v[k]),
                               vanderWattLL = v[1] + ((1-(1-alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
                               vanderWattUL = v[1] + ((1-(alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
                               stop(paste("The vMax method is spelt incorrectly: ",method)))
    }
  }
  
  #based on behaviour classes
  if (dynamic == "class"){
    #Check to see the input field matches one from the column
    colID <- which(names(trDF) == class.col)[1]
    if (is.na(colID)==TRUE) {stop(paste("The column name used for the class.col argument does not exist."))}
    #get the unique values of from class column
    classes <- unique(trDF[,colID])
    for (class in classes){
      ind <- which(trDF[,colID] == class)
      v.temp <- trDF$vi[ind]
      v.temp <- v.temp[which(v.temp <= vmaxtrunc)]
      v <- sort(v.temp,decreasing=T)
      
      trDF$dynVmax[ind] = switch(method,
                                 #Robson & Whitlock (1964) method --- main estimate and UL come
                                 # directly from the paper, however the LL was inferred
                                 Robson = v[1] + (v[1]-v[2]),
                                 RobsonLL = v[1] + ((1-(1-alpha))/(1-alpha))*(v[1]-v[2]),
                                 RobsonUL = v[1] + ((1-alpha)/alpha)*(v[1]-v[2]),
                                 #van der Watt (1980) method --- all come directly from the paper.
                                 vanderWatt = (((k+2)/(k+1))*v[1]) - ((1/(k+1))*v[k]),
                                 vanderWattLL = v[1] + ((1-(1-alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
                                 vanderWattUL = v[1] + ((1-(alpha)^(1/k))^(-1)-1)^(-1)*(v[1]-v[k]),
                                 stop(paste("The vMax method is spelt incorrectly: ",method)))
    }
  }
  #User-defined vector of Vmax values as column in ltraj object.
  if (!is.na(manualVmax)){
    col <- which(names(trDF)==manualVmax)
    trDF$dynVmax <- trDF[,col]
  }
  
  #-----------------------------------------------------
  ## Process vmaxtrunc
  vmax.inflate <- 1.05   #could be an optional variable to pass in
  ind <- which(trDF$vi > vmaxtrunc)
  trDF$dynVmax[ind] <- trDF$vi[ind]*vmax.inflate
  
  #-----------------------------------------------------
  traj <- dl(trDF)
  return(traj)
}

#============= End =============================================================


# ---- roxygen documentation ----
#' @title Time geographic kernel density estimation
#'
#' @description
#'   Compute the time geographic kernel density estimate of an animals utilization distribution following the methods described in Downs et al. (2011). 
#' @details
#'   The function \code{tgkde} can be used to delineate an animals home range using the time geographic kernel density estimation procedure described by Downs et al. (2011). Specifically, it modifies the shape of the traditional kernel to consider the time geographic limits of movement opportunity - termed the geoellipse by Downs, which is analagous to the potential path area concept described by Long & Nelson (2012,2015). Several basic functions -- including inverse distance, inverse distance squared, exponential, and normal -- can be used to quantify movement probabilities within the geoellipse. The output is then the utilization distribution of the animal, confined to the accessibility space defined by the potential path area home range as in Long & Nelson (2012,2015).
#'   The function \code{volras} can be used to extract volume contours (e.g., 95% volume contour) which are commonly used to delineate home ranges based on the utilization distribution.
#'
#' @param traj an object of the class \code{ltraj} which contains the time-stamped movement fixes of the first object. Note this object must be a \code{type II ltraj} object. For more information on objects of this type see \code{ help(ltraj)}.
#' @param disfun one of { \code{'inv', 'inv2', 'exp', 'norm'} } representing the shape of the distance decay function used to model movement probabilities within the geoellipse (see Downs et al. (2011) for more information). 
#' @param c2 parameter used to control shape of different distance decay functions, default = 1.
#' @param grid spatial resolution (pixel size) of output utilization raster in appropriate units. Default is chosen based on the x and y range of the input telemetry data. Alternatively, a \code{RasterLayer} can be passed in upon which the UD is computed.
#' @param ... additional parameters to be passed to the function \code{dynvmax}. For example, should include options for \code{dynamic} and \code{method}; see the documentation for \code{dynvmax} for more detailed information on what to include here.
#'    
#' @return
#'   This function returns a \code{RasterLayer} representing the utilization distribution of the animal
#'
#' @references
#' Downs, J.A., Horner, M.W., Tucker, A.D. (2011) Time-geographic density estimation for home range analysis. Annals of GIS. 17(3): 163-171.
# @keywords 
#' @seealso dynvmax, dynppa, volras
#' @examples
#' data(m3)
#' ud <- tgkde(m3,disfun='inv',method='vanderWatt')
#' raster::plot(ud)
#' 
#' @export
#
# ---- End of roxygen documentation ----

tgkde <- function(traj,disfun='inv',c2=1,grid=NA,...){
  #check to see if ltraj is not a type II
  if (attributes(traj)$typeII == FALSE) {stop("The trajectory object is not a TypeII ltraj object")}
  
  #Append the local Vmax values to the trajectory using the dynvmax function
  traj <- dynvmax(traj, ...)
  
  #convert to dataframe
  trDF <- ld(traj)
  n <- dim(trDF)[1]
  
  #If grid is rasterLayer
  if (class(grid) == 'RasterLayer'){
    xy <- data.frame(coordinates(grid))
  } else {
    #-------------------------------------------
    # Create the output grid
    #-------------------------------------------
    x.range <- max(trDF$x) - min(trDF$x)
    y.range <- max(trDF$y) - min(trDF$y)
    
    if (is.na(grid)){
      grid <- 0.01*min(c(x.range,y.range)) 
    }
    xx <- seq(min(trDF$x)-0.1*x.range,max(trDF$x)+0.1*x.range,by=grid)
    yy <- seq(min(trDF$y)-0.1*y.range,max(trDF$y)+0.1*y.range,by=grid)
    xy <- expand.grid(x=xx,y=yy)
    #-------------------------------------------
  }
  xy$z <- 0
  
  #loop through Traj object and compute the tgkde estimate
  x <- xy$x
  y <- xy$y
  for (i in 1:(n-1)){
    sx <- trDF$x[i]
    sy <- trDF$y[i]
    ex <- trDF$x[i+1]
    ey <- trDF$y[i+1]
    dt <- trDF$dt[i]
    vmax <- trDF$dynVmax[i]
    if (is.na(vmax)){next}
    dd <- sqrt((sx-x)^2 + (sy-y)^2) + sqrt((x-ex)^2 + (y-ey)^2)
    #dp <- sqrt((sx-ex)^2 + (sy-ey)^2)
    dmax <- vmax*dt
    ind <- which(dd > dmax)
    #insert function here perhaps using switch?
    d <- dd/dmax
    g <- switch(disfun,
                inv = 1/(c2*d),              #inverse distance
                inv2 = 1/(c2*(d^2)),         #inverse distance ^2        
                exp = exp(-c2*d),        #exponential function 
                norm = exp(-c2*(d^2)),   #normal function
                stop(paste('The distance decay function',disfun,'does not exist.'))
    )
    
    g[ind] <- 0
    #In Theory, need to normalize here, because each ellipse should sum to dt
    #g <- dt*g/sum(g)
    
    xy$z <- xy$z + g
  }
  
  #---- This all seems rather arbitrary? -----
  ##Normalize following the eqn. in Downs et al. (2011)
  ##get the 'average' vmax
  vmax. <- mean(trDF$dynVmax,na.rm=TRUE)
  ##get the 'overall' time difference
  DT <- as.numeric(difftime(trDF$date[n], trDF$date[1], units='secs'))
  ##Normalize the values
  xy$z <- xy$z * (1 / ((n-1)*(DT*vmax.)^2))
  
  # Downs et al. 2011 method does not sum to 1... 
  # based on normalization of each ellipse, just dividing by n-1 should work
  #xy$z <- xy$z / (n-1)
  
  #----------------------------
  #  Format output to Raster
  #----------------------------
  coordinates(xy) = ~x+y
  gridded(xy) <- TRUE
  ras <- raster(xy)
  #----------------------------
  return(ras)
}
# ---- roxygen documentation ----
#' @title PPA Ellipse
#'
#' @description
#'   Internal ellipse calculation function.
#' @details
#'   Internal function for calculating ellipses in time geographic analysis.
#' @param x first coordinate
#' @param y second coordinate
#' @param a semi-major axis
#' @param b semi-minor axis
#' @param theta rotation angle of the ellipse (in radians)
#' @param steps number of segments, from ePoints parameter in \code{dyn.ppa.hr}

#' @return
#'   This function returns a polygon ellipse.
#'
# @references
#' @keywords internal 
#' @seealso dynppa
# @examples
# @export
#
# ---- End of roxygen documentation ----
#===============================================================================
# function: ppaEllipse
# purpose: Calculate Time geography ellipses
#---------------------------------------------------------------
#x,y are center points
#a,b are semi-major and semi-minor axis
#theta is rotation angle (in radians)
#steps is number of segments to draw (default: 360 from ePoints parameter)
# This formulation uses the general parametric form of an ellipse
ppaEllipse <- function(x,y,a,b,theta,steps){
  X=rep(0,steps);Y=rep(0,steps)
  sinTheta <- sin(theta)
  cosTheta <- cos(theta)
  
  for (i in 1:steps){
    alpha <- i*(360/steps)*(pi/180)
    sinAlpha <- sin(alpha)
    cosAlpha <- cos(alpha)
    X[i] = x + (a*cosAlpha*cosTheta - b*sinAlpha*sinTheta)
    Y[i] = y + (a*cosAlpha*sinTheta + b*sinAlpha*cosTheta)
  }
  #save X,Y as Points Data frame
  ptDF <- data.frame(x=X,y=Y)
  ptDF[steps+1,] <- ptDF[1,]     #add first point to end to close ellipse
  #turn ellipses from points to polygons
  ellipsePoly <- Polygon(ptDF,hole=F)
  return(ellipsePoly)
}
#End of Function  
#===============================================================================
pacman::p_load(spatstat)
pacman::p_load(SDMTools)
pacman::p_load(adehabitatHR)
pacman::p_load(sp)
pacman::p_load(rgeos)
pacman::p_load(sf)
pacman::p_load(maptools)
pacman::p_load(purrr)
pacman::p_load(igraph)
pacman::p_load(tlocoh)
pacman::p_load(caTools)
pacman::p_load(ggplot2)
pacman::p_load(geoR)
pacman::p_load(plyr)
pacman::p_load(MASS)
pacman::p_load(SDraw)

#T-locoh
# install.packages("tlocoh", dependencies=TRUE, repos=c("http://R-Forge.R-project.org", "http://cran.cnr.berkeley.edu"))
require(tlocoh)
pacman::p_load(pbapply)
pacman::p_load(FNN)
pacman::p_load(rgdal)
pacman::p_load(raster)
pacman::p_load(png)
pacman::p_load(move)
pacman::p_load(XML)
pacman::p_load(transport)
pacman::p_load(xlsx)

################################################################################
###Perforated Home Ranges

  ##################setting up the entire point pattern ##############################
#Cluster Function
nclust <- function(x0, y0, radius, n) {
    return(runifdisc(n, radius, centre = c(x0, y0)))
}

#Create Dataframe to save area calculations for MCP
mcp.df <- data.frame(matrix(ncol = 4, nrow = 100))
mcp.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("MCP_Area50%", "MCP_Area 25%", "MCP_Area 10%", "MCP_Area 05%")
colnames(mcp.df) <- x
x <- c("Obs", "MCP_AUC 100%","MCP_AUC 50%", "MCP_AUC 25%", "MCP_AUC 10%", "MCP_AUC 05%")
colnames(mcp.df.auc) <- x
mcp.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("MCP_EMD 100%","MCP_EMD 50%", "MCP_EMD 25%", "MCP_EMD 10%", "MCP_EMD 05%")
colnames(mcp.df.emd) <- x
mcp.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("MCP_edge 100%", "100% Patches","MCP_edge 50%", "50% Patches","MCP_edge 25%", "25% Patches", "MCP_edge 10%", "10% Patches", "MCP_edge 05%", "05% Patches")
colnames(mcp.df.edge) <- x

#Create Dataframe to save area calculations for KDE
KDE.df <- data.frame(matrix(ncol = 4, nrow = 100))
KDE.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("KDE_Area50%", "KDE_Area 25%", "KDE_Area 10%", "KDE_Area 05%")
colnames(KDE.df) <- x
x <- c("Obs", "KDE_AUC 100%","KDE_AUC 50%", "KDE_AUC 25%", "KDE_AUC 10%", "KDE_AUC 05%")
colnames(KDE.df.auc) <- x
KDE.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("KDE_EMD 100%","KDE_EMD 50%", "KDE_EMD 25%", "KDE_EMD 10%", "KDE_EMD 05%")
colnames(KDE.df.emd) <- x
KDE.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("KDE_edge 100%", "100% Patches","KDE_edge 50%", "50% Patches","KDE_edge 25%", "25% Patches", "KDE_edge 10%", "10% Patches", "KDE_edge 05%", "05% Patches")
colnames(KDE.df.edge) <- x

#Create Dataframe to save area calculations for PPA
PPA.df <- data.frame(matrix(ncol = 4, nrow = 100))
PPA.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("PPA_Area50%", "PPA_Area 25%", "PPA_Area 10%", "PPA_Area 05%")
colnames(PPA.df) <- x
x <- c("Obs", "PPA_AUC 100%","PPA_AUC 50%", "PPA_AUC 25%", "PPA_AUC 10%", "PPA_AUC 05%")
colnames(PPA.df.auc) <- x
PPA.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("PPA_EMD 100%","PPA_EMD 50%", "PPA_EMD 25%", "PPA_EMD 10%", "PPA_EMD 05%")
colnames(PPA.df.emd) <- x
PPA.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("PPA_edge 100%", "100% Patches","PPA_edge 50%", "50% Patches","PPA_edge 25%", "25% Patches", "PPA_edge 10%", "10% Patches", "PPA_edge 05%", "05% Patches")
colnames(PPA.df.edge) <- x

#Create Dataframe to save area calculations for TDGE
tdge.df <- data.frame(matrix(ncol = 4, nrow = 100))
tdge.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("TDGE_Area50%", "TDGE_Area 25%", "TDGE_Area 10%", "TDGE_Area 05%")
colnames(tdge.df) <- x
x <- c("Obs", "TDGE_AUC 100%","TDGE_AUC 50%", "TDGE_AUC 25%", "TDGE_AUC 10%", "TDGE_AUC 05%")
colnames(tdge.df.auc) <- x
tdge.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("tdge_EMD 100%","tdge_EMD 50%", "tdge_EMD 25%", "tdge_EMD 10%", "tdge_EMD 05%")
colnames(tdge.df.emd) <- x
tdge.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("tdge_edge 100%", "100% Patches","tdge_edge 50%", "50% Patches","tdge_edge 25%", "25% Patches", "tdge_edge 10%", "10% Patches", "tdge_edge 05%", "05% Patches")
colnames(tdge.df.edge) <- x

#Create Dataframe to save area calculations for BRB
BRB.df <- data.frame(matrix(ncol = 4, nrow = 100))
BRB.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("BRB_Area50%", "BRB_Area 25%", "BRB_Area 10%", "BRB_Area 05%")
colnames(BRB.df) <- x
x <- c("Obs", "BRB_AUC 100%","BRB_AUC 50%", "BRB_AUC 25%", "BRB_AUC 10%", "BRB_AUC 05%")
colnames(BRB.df.auc) <- x
BRB.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("BRB_EMD 100%","BRB_EMD 50%", "BRB_EMD 25%", "BRB_EMD 10%", "BRB_EMD 05%")
colnames(BRB.df.emd) <- x
BRB.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("BRB_edge 100%", "100% Patches","BRB_edge 50%", "50% Patches","BRB_edge 25%", "25% Patches", "BRB_edge 10%", "10% Patches", "BRB_edge 05%", "05% Patches")
colnames(BRB.df.edge) <- x

#Create Dataframe to save area calculations for T-LoCoHO
TLoCoHo.df <- data.frame(matrix(ncol = 4, nrow = 100))
TLoCoHo.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("T-LoCoHO_Area50%", "T-LoCoHO_Area 25%", "T-LoCoHO_Area 10%", "T-LoCoHO_Area 05%")
colnames(TLoCoHo.df) <- x
x <- c("Obs", "T-LoCoH_AUC 100%","T-LoCoH_AUC 50%", "T-LoCoH_AUC 25%", "T-LoCoH_AUC 10%", "T-LoCoH_AUC 05%")
colnames(TLoCoHo.df.auc) <- x
TLoCoHo.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("TLoCoHo_EMD 100%","TLoCoHo_EMD 50%", "TLoCoHo_EMD 25%", "TLoCoHo_EMD 10%", "TLoCoHo_EMD 05%")
colnames(TLoCoHo.df.emd) <- x
TLoCoHo.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("TLoCoHo_edge 100%", "100% Patches","TLoCoHo_edge 50%", "50% Patches","TLoCoHo_edge 25%", "25% Patches", "TLoCoHo_edge 10%", "10% Patches", "TLoCoHo_edge 05%", "05% Patches")
colnames(TLoCoHo.df.edge) <- x

i = 1
###############Run loop for sampling###############################################
for (i in 1:4) {
  pc = rPoissonCluster(100, 0.05, nclust, radius = 0.25, n = 40)
  pc.df = as.data.frame(pc)
  # dfOrder(pc.df, c(1,2))
  sp = SpatialPoints(pc.df)
  #WGS 84
  projection(sp) = CRS("+init=epsg:4326")
  #To UTM 31 N
  sp<- spTransform(sp, CRS("+init=epsg:32631"))
  #Create sample for MCP
  s <- sample(sp, 10)
  #MCP
  cp <- mcp(s, percent = 95)
  #Clip points by MCP
  allPts<- gDifference(sp, cp)
  ################################################################
  #Set Up df for Traj
  ptNumber = length(allPts)
  time = seq(as.Date("2000/1/1"), by = "day", length.out = ptNumber)
  allPts$Date = time
  allPts$id = "Brock"
  allPts.df = as.data.frame(allPts)
  #Set up data for AUC
  #Create Hex grid for present/absence counts
  ptsBuffer = gBuffer(sp, width = 4000)
  hexPts <-spsample(ptsBuffer,type="hexagonal",cellsize=5000)
  hexPols <- HexPoints2SpatialPolygons(hexPts)
  #Define true for AUC
  pts2 <- allPts.df[, -c(1:2)] # delete columns 5 through 7
  pts2 <- SpatialPoints(pts2)
  projection(pts2) = CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pts2, hexPols)))
  #MCP
  pred.df.MCP <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
  names(pred.df.MCP ) <- "GridNumber"
  pred.df.MCP$GridNumber <- 1:length(hexPols)
  pred.df.MCP$Obs <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
  pred.df.MCP$Obs[is.na(pred.df.MCP$Obs)] <- 0
  pred.df.MCP$Presence <- pred.df.MCP$Obs > 0
  #PPA
  pred.df.PPA <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
  names(pred.df.PPA ) <- "GridNumber"
  pred.df.PPA$GridNumber <- 1:length(hexPols)
  pred.df.PPA$Obs <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
  pred.df.PPA$Obs[is.na(pred.df.PPA$Obs)] <- 0
  pred.df.PPA$Presence <- pred.df.PPA$Obs > 0
  #tdge
  pred.df.tdge <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
  names(pred.df.tdge ) <- "GridNumber"
  pred.df.tdge$GridNumber <- 1:length(hexPols)
  pred.df.tdge$Obs <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
  pred.df.tdge$Obs[is.na(pred.df.tdge$Obs)] <- 0
  pred.df.tdge$Presence <- pred.df.tdge$Obs > 0
  #BRB
  pred.df.BRB <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
  names(pred.df.BRB ) <- "GridNumber"
  pred.df.BRB$GridNumber <- 1:length(hexPols)
  pred.df.BRB$Obs <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
  pred.df.BRB$Obs[is.na(pred.df.BRB$Obs)] <- 0
  pred.df.BRB$Presence <- pred.df.BRB$Obs > 0
  #T-Locoh
  pred.df.Tlocoh <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
  names(pred.df.Tlocoh ) <- "GridNumber"
  pred.df.Tlocoh$GridNumber <- 1:length(hexPols)
  pred.df.Tlocoh$Obs <- as.numeric(res$Freq[match(pred.df.Tlocoh$GridNumber, res$Var1)])
  pred.df.Tlocoh$Obs[is.na(pred.df.Tlocoh$Obs)] <- 0
  pred.df.Tlocoh$Presence <- pred.df.Tlocoh$Obs > 0
  ###################################################################
 
  #########################End Point Pattern Set UP ##############################
  # Sample the Points at differ levels
  #Sample 50 %
  sampleSize = length(allPts) * .50
  pts.50 = spsample(allPts,n = sampleSize, "random")
  #Sample 25 %
  sampleSize = length(allPts) * .25
  pts.25 = spsample(allPts,n = sampleSize, "random")
  #Sample 10 %
  sampleSize = length(allPts) * .10
  pts.10 = spsample(allPts,n = sampleSize, "random")
  #Sample 5 %
  sampleSize = length(allPts) * .05
  pts.05 = spsample(allPts,n = sampleSize, "random")
  #################################Run Home Ranges###############################
  #MCP Home Ranges
  mcp.100 = mcp(allPts, percent = 95)
  #MCP Core Area
  # mcp.100.perforated.50 = mcp(allPts, percent = 50)
  #MCP Home Range 50%
  mcp.50 = mcp(pts.50, percent = 95)
  #MCP Core Area 50%
  # mcp.50.perforated.50 = mcp(pts.50, percent = 50)
  #MCP Home Range 50%
  mcp.25= mcp(pts.25, percent = 95)
  #MCP Core Area 50%
  # mcp.25.perforated.50 = mcp(pts.25, percent = 50)
  #MCP Home Range 50%
  mcp.10 = mcp(pts.10, percent = 95)
  #MCP Core Area 50%
  # mcp.10.perforated.50 = mcp(pts.10, percent = 50)
  #MCP Home Range 50%
  mcp.05 = mcp(pts.05, percent = 95)
  #MCP Core Area 50%
  # mcp.05.perforated.50 = mcp(pts.05, percent = 50)
  #Perferate the MCP to get correct area
  area.100 <- gArea(gDifference(mcp.100,cp))
  area.50 <- gArea(mcp.50)
  area.25 <- gArea(mcp.25)
  area.10 <- gArea(mcp.10)
  area.05 <- gArea(mcp.05)
  mcp.df[i,1] <- area.50/area.100
  mcp.df[i,2] <- area.25/area.100
  mcp.df[i,3] <- area.10/area.100
  mcp.df[i,4] <- area.05/area.100
  mcp.df.edge[i,1] <- (lineLength(mcp.100, byid = FALSE)/area.100) * 1000
  mcp.df.edge[i,3] <- (lineLength(mcp.50, byid = FALSE)/area.50) * 1000
  mcp.df.edge[i,5] <- (lineLength(mcp.25, byid = FALSE)/area.25) * 1000
  mcp.df.edge[i,7] <- (lineLength(mcp.10, byid = FALSE)/area.10) * 1000
  mcp.df.edge[i,9] <- (lineLength(mcp.05, byid = FALSE)/area.05) * 1000
  mcp.df.edge[i,2] <- length(mcp.100)
  mcp.df.edge[i,4] <- length(mcp.50)
  mcp.df.edge[i,6] <- length(mcp.25)
  mcp.df.edge[i,8] <- length(mcp.10)
  mcp.df.edge[i,10] <- length(mcp.05)
  ###################################################################
  #AUC MCP
  pred.100.mcp <- gIntersection(mcp.100, sp)
  projection(pred.100.mcp) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.100.mcp, hexPols)))
  pred.df.MCP$pred.100 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
  pred.df.MCP$pred.100[is.na(pred.df.MCP$pred.100)] <- 0
  #reorder dataframe
  pred.df.MCP <- pred.df.MCP[c(1,3,2,4)]
  
  pred.50.mcp <- gIntersection(mcp.50, sp)
  projection(pred.50.mcp) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.50.mcp, hexPols)))
  pred.df.MCP$pred.50 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
  pred.df.MCP$pred.50[is.na(pred.df.MCP$pred.50)] <- 0
  
  pred.25.mcp <- gIntersection(mcp.25, sp)
  projection(pred.25.mcp) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.25.mcp, hexPols)))
  pred.df.MCP$pred.25 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
  pred.df.MCP$pred.25[is.na(pred.df.MCP$pred.25)] <- 0
  
  pred.10.mcp <- gIntersection(mcp.10, sp)
  projection(pred.10.mcp) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.10.mcp, hexPols)))
  pred.df.MCP$pred.10 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
  pred.df.MCP$pred.10[is.na(pred.df.MCP$pred.10)] <- 0
  
  pred.05.mcp <- gIntersection(mcp.05, sp)
  projection(pred.05.mcp) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.05.mcp, hexPols)))
  pred.df.MCP$pred.05 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
  pred.df.MCP$pred.05[is.na(pred.df.MCP$pred.05)] <- 0
  
  tAUC <- colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
  colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
  # tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
  mcp.df.auc[i,1] <- tAUC[1,1]
  mcp.df.auc[i,2] <- tAUC[1,2]
  mcp.df.auc[i,3] <- tAUC[1,3]
  mcp.df.auc[i,4] <- tAUC[1,4]
  mcp.df.auc[i,5] <- tAUC[1,5]
  mcp.df.auc[i,6] <- tAUC[1,6]
  
  ##################################################################
  #KDE Homerange
  KDE.100 <- kernelUD(allPts, h="href")
  KDE.100 <- getverticeshr(KDE.100,percent = 95)
  KDE.50 <- kernelUD(pts.50, h="href")
  KDE.50 <- getverticeshr(KDE.50,percent = 95)
  KDE.25 <- kernelUD(pts.25, h="href")
  KDE.25 <- getverticeshr(KDE.25,percent = 95)
  KDE.10 <- kernelUD(pts.10, h="href")
  KDE.10 <- getverticeshr(KDE.10,percent = 95)
  KDE.05 <- kernelUD(pts.05, h="href")
  KDE.05 <- getverticeshr(KDE.05,percent = 95)
  #KDE Area
  area.100 <- KDE.100$area
  area.50 <- KDE.50$area
  area.25 <- KDE.25$area
  area.10 <- KDE.10$area
  area.05 <- KDE.05$area
  #Assign Area percents to datafrmae
  KDE.df[i,1] <- area.50/area.100
  KDE.df[i,2] <- area.25/area.100
  KDE.df[i,3] <- area.10/area.100
  KDE.df[i,4] <- area.05/area.100
  KDE.df.edge[i,1] <- (lineLength(KDE.100, byid = FALSE)/area.100) * 1000
  KDE.df.edge[i,3] <- (lineLength(KDE.50, byid = FALSE)/area.50) * 1000
  KDE.df.edge[i,5] <- (lineLength(KDE.25, byid = FALSE)/area.25) * 1000
  KDE.df.edge[i,7] <- (lineLength(KDE.10, byid = FALSE)/area.10) * 1000
  KDE.df.edge[i,9] <- (lineLength(KDE.05, byid = FALSE)/area.05) * 1000
  KDE.df.edge[i,2] <- length(KDE.100)
  KDE.df.edge[i,4] <- length(KDE.50)
  KDE.df.edge[i,6] <- length(KDE.25)
  KDE.df.edge[i,8] <- length(KDE.10)
  KDE.df.edge[i,10] <- length(KDE.05)
  ###################################################################
  #AUC KDE
  #Define true AUC
  pts2 <- allPts.df[, -c(1:2)] # delete columns 5 through 7
  pts2 <- SpatialPoints(pts2)
  projection(pts2) = CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pts2, hexPols)))
  #KDE
  pred.df.KDE <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
  names(pred.df.KDE) <- "GridNumber"
  pred.df.KDE$GridNumber <- 1:length(hexPols)
  pred.df.KDE$Obs <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
  pred.df.KDE$Obs[is.na(pred.df.KDE$Obs)] <- 0
  pred.df.KDE$Presence <- pred.df.KDE$Obs > 0
  pred.df.KDE <- pred.df.KDE[c(1,3,2)]
  
  #100 %
  pred.100.KDE <- gIntersection(KDE.100, sp)
  projection(pred.100.KDE) <- CRS("+init=epsg:32631")
  plot(pred.100.KDE, add = TRUE)
  res <- as.data.frame(table(over(pred.100.KDE, hexPols)))
  pred.df.KDE$pred.100 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
  pred.df.KDE$pred.100[is.na(pred.df.KDE$pred.100)] <- 0
  
  
  #50%
  pred.50.KDE <- gIntersection(KDE.50, sp)
  projection(pred.50.KDE) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.50.KDE, hexPols)))
  pred.df.KDE$pred.50 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
  pred.df.KDE$pred.50[is.na(pred.df.KDE$pred.50)] <- 0
  
  
  #25%
  pred.25.KDE <- gIntersection(KDE.25, sp)
  projection(pred.25.KDE) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.25.KDE, hexPols)))
  pred.df.KDE$pred.25 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
  pred.df.KDE$pred.25[is.na(pred.df.KDE$pred.25)] <- 0
  
  
  #10%
  pred.10.KDE <- gIntersection(KDE.10, sp)
  projection(pred.10.KDE) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.10.KDE, hexPols)))
  pred.df.KDE$pred.10 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
  pred.df.KDE$pred.10[is.na(pred.df.KDE$pred.10)] <- 0
  
  #05%
  pred.05.KDE <- gIntersection(KDE.05, sp)
  projection(pred.05.KDE) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.05.KDE, hexPols)))
  pred.df.KDE$pred.05 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
  pred.df.KDE$pred.05[is.na(pred.df.KDE$pred.05)] <- 0
  
  #auc
  tAUC <- colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
  colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
  # tAUC <- colAUC(sp,KDE.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
  KDE.df.auc[i,1] <- tAUC[1,1]
  KDE.df.auc[i,2] <- tAUC[1,2]
  KDE.df.auc[i,3] <- tAUC[1,3]
  KDE.df.auc[i,4] <- tAUC[1,4]
  KDE.df.auc[i,5] <- tAUC[1,5]
  KDE.df.auc[i,6] <- tAUC[1,6]
  
  ######################################################################################################
  #PPA, TDGE, BRB, AND t-LOCOCH  Homerange
  #Converting to ltraj for time dynamics
  PPA.100 <- as.ltraj(xy = allPts.df[,c("x","y")], date = as.POSIXct(strptime(allPts.df$Date,"%Y-%m-%d"))
                    , id = allPts.df$id, proj4string = CRS("+init=epsg:32631"))
  tdge.100 <- tgkde(PPA.100)
  tdge.100 <- volras(tdge.100,95)
  projection(tdge.100) = CRS("+init=epsg:32631")
  vv <- BRB.D(PPA.100, Tmax=180*60, Lmin=60)
  #Hmin is super important here. 
  BRB.100 <- BRB(PPA.100, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
  BRB.100 <- getverticeshr(BRB.100,percent = 95)
  projection(BRB.100 ) = CRS("+init=epsg:32631")
  PPA.100 <- dynppa(PPA.100)
  projection(PPA.100) = CRS("+init=epsg:32631")
  #T-locoh
  T.100 <- xyt.lxy(xy = allPts.df[,c("x","y")], dt= as.POSIXct(strptime(allPts.df$Date,"%Y-%m-%d"))
                   , proj4string= CRS("+init=epsg:32631"), id=allPts$id)
  toni.lxy <- lxy.ptsh.add(T.100)
  # S value of 60 percent of the data "half way between 40-80%"
  s.value = toni.lxy$ptsh
  s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
  T.100 <- lxy.nn.add(T.100, s=s.value, k=25)
  T.100 <- lxy.lhs(T.100, k=3*2:8, s=s.value)
  T.100 <- lhs.iso.add(T.100)
  T.100 <- isopleths(T.100)
  T.100.shp <- T.100[[1]][T.100[[1]]$iso.level == 0.95,]
  #####################################################################
  #AUC
  #T-locoh
  pred.100 <- gIntersection(T.100.shp, sp)
  projection(pred.100) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.100, hexPols)))
  pred.df.Tlocoh$pred.100 <- as.numeric(res$Freq[match(pred.df.Tlocoh$GridNumber, res$Var1)])
  pred.df.Tlocoh$pred.100[is.na(pred.df.Tlocoh$pred.100)] <- 0
  #Tdge
  pred.100 <- gIntersection(tdge.100, sp)
  projection(pred.100) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.100, hexPols)))
  pred.df.tdge$pred.100 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
  pred.df.tdge$pred.100[is.na(pred.df.tdge$pred.100)] <- 0
  #BRB
  pred.100 <- gIntersection(BRB.100, sp)
  projection(pred.100) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.100, hexPols)))
  pred.df.BRB$pred.100 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
  pred.df.BRB$pred.100[is.na(pred.df.BRB$pred.100)] <- 0
  #PPA
  pred.100 <- gIntersection(PPA.100, sp)
  projection(pred.100) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.100, hexPols)))
  pred.df.PPA$pred.100 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
  pred.df.PPA$pred.100[is.na(pred.df.PPA$pred.100)] <- 0
  
  
  #####################################################################
  #Even Sample
  PPA.50 <- allPts.df[seq(1, NROW(allPts.df), by = 2),]
  #Need tp run xyt before converting to PPA
  T.50 <- xyt.lxy(xy = PPA.50[,c("x","y")], dt= as.POSIXct(strptime(PPA.50$Date,"%Y-%m-%d"))
                  , proj4string= CRS("+init=epsg:32631"), id=PPA.50$id)
  PPA.50 <- as.ltraj(xy = PPA.50[,c("x","y")], date = as.POSIXct(strptime(PPA.50$Date,"%Y-%m-%d"))
                      , id = PPA.50$id, proj4string = CRS("+init=epsg:32631"))
  tdge.50 <- tgkde(PPA.50)
  tdge.50 <- volras(tdge.50,95)
  projection(tdge.50) = CRS("+init=epsg:32631")
  vv <- BRB.D(PPA.50, Tmax=180*60, Lmin=60)
  BRB.50 <- BRB(PPA.50, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
  BRB.50 <- getverticeshr(BRB.50,percent = 95)
  projection(BRB.50) = CRS("+init=epsg:32631")
  PPA.50 <- dynppa(PPA.50)
  projection(PPA.50 ) = CRS("+init=epsg:32631")
  #T-locoh
  toni.lxy <- lxy.ptsh.add(T.50)
  # S value of 60 percent of the data "half way between 40-80%"
  s.value = toni.lxy$ptsh
  s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
  T.50 <- lxy.nn.add(T.50, s=s.value, k=25)
  T.50 <- lxy.lhs(T.50, k=3*2:8, s=s.value)
  T.50 <- lhs.iso.add(T.50)
  T.50 <- isopleths(T.50)
  #Get isopleth
  T.50.shp <- T.50[[1]][T.50[[1]]$iso.level == 0.95,]
  #######################################################################
  #AUC
  #T-locoh
  pred.50 <- gIntersection(T.50.shp, sp)
  projection(pred.50) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.50, hexPols)))
  pred.df.Tlocoh$pred.50 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
  pred.df.Tlocoh$pred.50[is.na(pred.df.Tlocoh$pred.50)] <- 0
  #Tdge
  pred.50 <- gIntersection(tdge.50, sp)
  projection(pred.50) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.50, hexPols)))
  pred.df.tdge$pred.50 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
  pred.df.tdge$pred.50[is.na(pred.df.tdge$pred.50)] <- 0
  #BRB
  pred.50 <- gIntersection(BRB.50, sp)
  projection(pred.50) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.50, hexPols)))
  pred.df.BRB$pred.50 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
  pred.df.BRB$pred.50[is.na(pred.df.BRB$pred.50)] <- 0
  #PPA
  pred.50 <- gIntersection(PPA.50, sp)
  projection(pred.50) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.50, hexPols)))
  pred.df.PPA$pred.50 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
  pred.df.PPA$pred.50[is.na(pred.df.PPA$pred.50)] <- 0
  ########################################################################
  PPA.25 <- allPts.df[seq(1, NROW(allPts.df), by = 4),]
  T.25 <- xyt.lxy(xy = PPA.25[,c("x","y")], dt= as.POSIXct(strptime(PPA.25$Date,"%Y-%m-%d"))
                  , proj4string= CRS("+init=epsg:32631"), id=PPA.25$id)
  PPA.25 <- as.ltraj(xy = PPA.25[,c("x","y")], date = as.POSIXct(strptime(PPA.25$Date,"%Y-%m-%d"))
                     , id = PPA.25$id, proj4string = CRS("+init=epsg:32631"))
  tdge.25 <- tgkde(PPA.25)
  tdge.25 <- volras(tdge.25,95)
  projection(tdge.25) = CRS("+init=epsg:32631")
  vv <- BRB.D(PPA.25, Tmax=180*60, Lmin=60)
  BRB.25 <- BRB(PPA.25, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
  BRB.25 <- getverticeshr(BRB.25,percent = 95)
  projection(BRB.25) = CRS("+init=epsg:32631")
  PPA.25 <- dynppa(PPA.25)
  projection(PPA.25) = CRS("+init=epsg:32631")
  #T-locoh
  toni.lxy <- lxy.ptsh.add(T.25)
  # S value of 60 percent of the data "half way between 40-80%"
  s.value = toni.lxy$ptsh
  s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
  T.25 <- lxy.nn.add(T.25, s=s.value, k=25)
  T.25 <- lxy.lhs(T.25, k=3*2:8, s=s.value)
  T.25 <- lhs.iso.add(T.25)
  T.25 <- isopleths(T.25)
  T.25.shp <- T.25[[1]][T.25[[1]]$iso.level == 0.95,]
  #######################################################################
  pred.25 <- gIntersection(T.25.shp, sp)
  projection(pred.25) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.25, hexPols)))
  pred.df.Tlocoh$pred.25 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
  pred.df.Tlocoh$pred.25[is.na(pred.df.Tlocoh$pred.25)] <- 0
  #Tdge
  pred.25 <- gIntersection(tdge.25, sp)
  projection(pred.25) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.25, hexPols)))
  pred.df.tdge$pred.25 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
  pred.df.tdge$pred.25[is.na(pred.df.tdge$pred.25)] <- 0
  #BRB
  pred.25 <- gIntersection(BRB.25, sp)
  projection(pred.25) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.25, hexPols)))
  pred.df.BRB$pred.25 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
  pred.df.BRB$pred.25[is.na(pred.df.BRB$pred.25)] <- 0
  #PPA
  pred.25 <- gIntersection(PPA.25, sp)
  projection(pred.25) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.25, hexPols)))
  pred.df.PPA$pred.25 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
  pred.df.PPA$pred.25[is.na(pred.df.PPA$pred.25)] <- 0
  ########################################################################
  PPA.10 <- allPts.df[seq(1, NROW(allPts.df), by = 10),]
  T.10 <- xyt.lxy(xy = PPA.10[,c("x","y")], dt= as.POSIXct(strptime(PPA.10$Date,"%Y-%m-%d"))
                  , proj4string= CRS("+init=epsg:32631"), id=PPA.10$id)
  PPA.10<- as.ltraj(xy = PPA.10[,c("x","y")], date = as.POSIXct(strptime(PPA.10$Date,"%Y-%m-%d"))
                     , id = PPA.10$id, proj4string = CRS("+init=epsg:32631"))
  tdge.10 <- tgkde(PPA.10)
  tdge.10 <- volras(tdge.10,95)
  projection(tdge.10 ) = CRS("+init=epsg:32631")
  vv <- BRB.D(PPA.10, Tmax=180*60, Lmin=60)
  BRB.10 <- BRB(PPA.10, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
  BRB.10 <- getverticeshr(BRB.10,percent = 95)
  projection(BRB.10) = CRS("+init=epsg:32631")
  PPA.10 <- dynppa(PPA.10)
  projection(PPA.10 ) = CRS("+init=epsg:32631")
  #T-locoh
  toni.lxy <- lxy.ptsh.add(T.10)
  # S value of 60 percent of the data "half way between 40-80%"
  s.value = toni.lxy$ptsh
  s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
  T.10 <- lxy.nn.add(T.10, s=s.value, k=25)
  T.10 <- lxy.lhs(T.10, k=3*2:8, s=s.value)
  T.10 <- lhs.iso.add(T.10)
  T.10 <- isopleths(T.10)
  T.10.shp <- T.10[[1]][T.10[[1]]$iso.level == 0.95,]
  #######################################################################
  pred.10 <- gIntersection(T.10.shp, sp)
  projection(pred.10) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.10, hexPols)))
  pred.df.Tlocoh$pred.10 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
  pred.df.Tlocoh$pred.10[is.na(pred.df.Tlocoh$pred.10)] <- 0
  #Tdge
  pred.10 <- gIntersection(tdge.10, sp)
  projection(pred.10) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.10, hexPols)))
  pred.df.tdge$pred.10 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
  pred.df.tdge$pred.10[is.na(pred.df.tdge$pred.10)] <- 0
  #BRB
  pred.10 <- gIntersection(BRB.10, sp)
  projection(pred.10) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.10, hexPols)))
  pred.df.BRB$pred.10 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
  pred.df.BRB$pred.10[is.na(pred.df.BRB$pred.10)] <- 0
  #PPA
  pred.10 <- gIntersection(PPA.10, sp)
  projection(pred.10) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.10, hexPols)))
  pred.df.PPA$pred.10 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
  pred.df.PPA$pred.10[is.na(pred.df.PPA$pred.10)] <- 0
  ########################################################################
  PPA.05 <- allPts.df[seq(1, NROW(allPts.df), by = 4),]
  T.05 <- xyt.lxy(xy = PPA.05[,c("x","y")], dt= as.POSIXct(strptime(PPA.05$Date,"%Y-%m-%d"))
                  , proj4string= CRS("+init=epsg:32631"), id=PPA.05$id)
  PPA.05 <- as.ltraj(xy = PPA.05[,c("x","y")], date = as.POSIXct(strptime(PPA.05$Date,"%Y-%m-%d"))
                     , id = PPA.05$id, proj4string = CRS("+init=epsg:32631"))
  tdge.05 <- tgkde(PPA.05)
  tdge.05 <- volras(tdge.05,95)
  projection(tdge.05) = CRS("+init=epsg:32631")
  vv <- BRB.D(PPA.05, Tmax=180*60, Lmin=60)
  BRB.05 <- BRB(PPA.05, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
  BRB.05 <- getverticeshr(BRB.05,percent = 95)
  projection(BRB.05) = CRS("+init=epsg:32631")
  PPA.05 <- dynppa(PPA.05)
  projection(PPA.05) = CRS("+init=epsg:32631")
  #T-locoh
  toni.lxy <- lxy.ptsh.add(T.05)
  # S value of 60 percent of the data "half way between 40-80%"
  s.value = toni.lxy$ptsh
  s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
  T.05 <- lxy.nn.add(T.05, s=s.value, k=25)
  T.05 <- lxy.lhs(T.05, k=3*2:8, s=s.value)
  T.05 <- lhs.iso.add(T.05)
  T.05 <- isopleths(T.05)
  T.05.shp <- T.05[[1]][T.05[[1]]$iso.level == 0.95,]
  #######################################################################
  pred.05 <- gIntersection(T.05.shp, sp)
  projection(pred.05) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.05, hexPols)))
  pred.df.Tlocoh$pred.05 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
  pred.df.Tlocoh$pred.05[is.na(pred.df.Tlocoh$pred.05)] <- 0
  #Tdge
  pred.05 <- gIntersection(tdge.05, sp)
  projection(pred.05) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.05, hexPols)))
  pred.df.tdge$pred.05 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
  pred.df.tdge$pred.05[is.na(pred.df.tdge$pred.05)] <- 0
  #BRB
  pred.05 <- gIntersection(BRB.05, sp)
  projection(pred.05) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.05, hexPols)))
  pred.df.BRB$pred.05 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
  pred.df.BRB$pred.05[is.na(pred.df.BRB$pred.05)] <- 0
  #PPA
  pred.05 <- gIntersection(PPA.100, sp)
  projection(pred.05) <- CRS("+init=epsg:32631")
  res <- as.data.frame(table(over(pred.05, hexPols)))
  pred.df.PPA$pred.05 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
  pred.df.PPA$pred.05[is.na(pred.df.PPA$pred.05)] <- 0
  ########################################################################
  #PPA Area
  area.100 <- gArea(gDifference(PPA.100,cp))
  area.50 <- gArea(PPA.50)
  area.25 <- gArea(PPA.25)
  area.10 <- gArea(PPA.10)
  area.05 <- gArea(PPA.05)
  #Edge Complexity 
  PPA.df.edge[i,1] <- (lineLength(PPA.100, byid = FALSE)/area.100) * 1000
  PPA.df.edge[i,3] <- (lineLength(PPA.50, byid = FALSE)/area.50) * 1000
  PPA.df.edge[i,5] <- (lineLength(PPA.25, byid = FALSE)/area.25) * 1000
  PPA.df.edge[i,7] <- (lineLength(PPA.10, byid = FALSE)/area.10) * 1000
  PPA.df.edge[i,9] <- (lineLength(PPA.05, byid = FALSE)/area.05) * 1000
  PPA.df.edge[i,2] <- length(PPA.100)
  PPA.df.edge[i,4] <- length(PPA.50)
  PPA.df.edge[i,6] <- length(PPA.25)
  PPA.df.edge[i,8] <- length(PPA.10)
  PPA.df.edge[i,10] <- length(PPA.05)
  #Assign Area percents to datafrmae
  PPA.df[i,1] <- area.50/area.100
  PPA.df[i,2] <- area.25/area.100
  PPA.df[i,3] <- area.10/area.100
  PPA.df[i,4] <- area.05/area.100
  #TDGE Area
  area.100 <- gArea(gDifference(tdge.100,cp))
  area.50 <- gArea(tdge.50)
  area.25 <- gArea(tdge.25)
  area.10 <- gArea(tdge.10)
  area.05 <- gArea(tdge.05)
  #Assign Area percents to datafrmae
  tdge.df[i,1] <- area.50/area.100
  tdge.df[i,2] <- area.25/area.100
  tdge.df[i,3] <- area.10/area.100
  tdge.df[i,4] <- area.05/area.100
  #Edge Complexity 
  tdge.df.edge[i,1] <- (lineLength(tdge.100, byid = FALSE)/area.100) * 1000
  tdge.df.edge[i,3] <- (lineLength(tdge.50, byid = FALSE)/area.50) * 1000
  tdge.df.edge[i,5] <- (lineLength(tdge.25, byid = FALSE)/area.25) * 1000
  tdge.df.edge[i,7] <- (lineLength(tdge.10, byid = FALSE)/area.10) * 1000
  tdge.df.edge[i,9] <- (lineLength(tdge.05, byid = FALSE)/area.05) * 1000
  tdge.df.edge[i,2] <- length(tdge.100)
  tdge.df.edge[i,4] <- length(tdge.50)
  tdge.df.edge[i,6] <- length(tdge.25)
  tdge.df.edge[i,8] <- length(tdge.10)
  tdge.df.edge[i,10] <- length(tdge.05)
  #BRB Area
  area.100 <- gArea(BRB.100)
  area.50 <- gArea(BRB.50)
  area.25 <- gArea(BRB.25)
  area.10 <- gArea(BRB.10)
  area.05 <- gArea(BRB.05)
  #Assign Area percents to datafrmae
  BRB.df[i,1] <- area.50/area.100
  BRB.df[i,2] <- area.25/area.100
  BRB.df[i,3] <- area.10/area.100
  BRB.df[i,4] <- area.05/area.100
  #Edge Complexity 
  BRB.df.edge[i,1] <- (lineLength(BRB.100, byid = FALSE)/area.100) * 1000
  BRB.df.edge[i,3] <- (lineLength(BRB.50, byid = FALSE)/area.50) * 1000
  BRB.df.edge[i,5] <- (lineLength(BRB.25, byid = FALSE)/area.25) * 1000
  BRB.df.edge[i,7] <- (lineLength(BRB.10, byid = FALSE)/area.10) * 1000
  BRB.df.edge[i,9] <- (lineLength(BRB.05, byid = FALSE)/area.05) * 1000
  BRB.df.edge[i,2] <- length(BRB.100)
  BRB.df.edge[i,4] <- length(BRB.50)
  BRB.df.edge[i,6] <- length(BRB.25)
  BRB.df.edge[i,8] <- length(BRB.10)
  BRB.df.edge[i,10] <- length(BRB.05)
  #T-locoh area
  area.100 <- T.100[[1]]@data[5,2]
  area.50 <- T.50[[1]]@data[5,2]
  area.25 <- T.25[[1]]@data[5,2]
  area.10 <- T.10[[1]]@data[5,2]
  area.05 <- T.05[[1]]@data[5,2]
  #Assign Area percents to datafrmae
  TLoCoHo.df[i,1] <- area.50/area.100
  TLoCoHo.df[i,2] <- area.25/area.100
  TLoCoHo.df[i,3] <- area.10/area.100
  TLoCoHo.df[i,4] <- area.05/area.100
  #Edge Complexity 
  TLoCoHo.df.edge[i,1] <- (lineLength(T.100.shp, byid = FALSE)/area.100) * 1000
  TLoCoHo.df.edge[i,3] <- (lineLength(T.50.shp, byid = FALSE)/area.50) * 1000
  TLoCoHo.df.edge[i,5] <- (lineLength(T.25.shp, byid = FALSE)/area.25) * 1000
  TLoCoHo.df.edge[i,7] <- (lineLength(T.10.shp, byid = FALSE)/area.10) * 1000
  TLoCoHo.df.edge[i,9] <- (lineLength(T.05.shp, byid = FALSE)/area.05) * 1000
  TLoCoHo.df.edge[i,2] <- length(T.100.shp)
  TLoCoHo.df.edge[i,4] <- length(T.50.shp)
  TLoCoHo.df.edge[i,6] <- length(T.25.shp)
  TLoCoHo.df.edge[i,8] <- length(T.10.shp)
  TLoCoHo.df.edge[i,10] <- length(T.05.shp)
  #####AUC#############################################
  pred.df.Tlocoh <-pred.df.Tlocoh[c(1,3,2,4,5,6,7,8)]
  pred.df.tdge <-pred.df.tdge[c(1,3,2,4,5,6,7,8)]
  pred.df.BRB <-pred.df.BRB[c(1,3,2,4,5,6,7,8)]
  pred.df.PPA <-pred.df.PPA[c(1,3,2,4,5,6,7,8)]
  tAUC <- colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
  colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
  # tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
  KDE.df.auc[i,1] <- tAUC[1,1]
  KDE.df.auc[i,2] <- tAUC[1,2]
  KDE.df.auc[i,3] <- tAUC[1,3]
  KDE.df.auc[i,4] <- tAUC[1,4]
  KDE.df.auc[i,5] <- tAUC[1,5]
  KDE.df.auc[i,6] <- tAUC[1,6]
  tAUC <- colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
  colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
  # tAUC <- colAUC(allPts,MCP.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
  mcp.df.auc[i,1] <- tAUC[1,1]
  mcp.df.auc[i,2] <- tAUC[1,2]
  mcp.df.auc[i,3] <- tAUC[1,3]
  mcp.df.auc[i,4] <- tAUC[1,4]
  mcp.df.auc[i,5] <- tAUC[1,5]
  mcp.df.auc[i,6] <- tAUC[1,6]
  tAUC <- colAUC(pred.df.Tlocoh[,3:8], pred.df.Tlocoh[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
  colAUC(pred.df.Tlocoh[,3:8], pred.df.Tlocoh[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
  # tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
  TLoCoHo.df.auc[i,1] <- tAUC[1,1]
  TLoCoHo.df.auc[i,2] <- tAUC[1,2]
  TLoCoHo.df.auc[i,3] <- tAUC[1,3]
  TLoCoHo.df.auc[i,4] <- tAUC[1,4]
  TLoCoHo.df.auc[i,5] <- tAUC[1,5]
  TLoCoHo.df.auc[i,6] <- tAUC[1,6]
  tAUC <- colAUC(pred.df.tdge[,3:8], pred.df.tdge[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
  colAUC(pred.df.tdge[,3:8], pred.df.tdge[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
  # tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
  tdge.df.auc[i,1] <- tAUC[1,1]
  tdge.df.auc[i,2] <- tAUC[1,2]
  tdge.df.auc[i,3] <- tAUC[1,3]
  tdge.df.auc[i,4] <- tAUC[1,4]
  tdge.df.auc[i,5] <- tAUC[1,5]
  tdge.df.auc[i,6] <- tAUC[1,6]
  tAUC <- colAUC(pred.df.BRB[,3:8], pred.df.BRB[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
  colAUC(pred.df.BRB[,3:8], pred.df.BRB[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
  # tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
  BRB.df.auc[i,1] <- tAUC[1,1]
  BRB.df.auc[i,2] <- tAUC[1,2]
  BRB.df.auc[i,3] <- tAUC[1,3]
  BRB.df.auc[i,4] <- tAUC[1,4]
  BRB.df.auc[i,5] <- tAUC[1,5]
  BRB.df.auc[i,6] <- tAUC[1,6]
  tAUC <- colAUC(pred.df.PPA[,3:8], pred.df.PPA[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
  colAUC(pred.df.PPA[,3:8], pred.df.PPA[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
  # tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
  PPA.df.auc[i,1] <- tAUC[1,1]
  PPA.df.auc[i,2] <- tAUC[1,2]
  PPA.df.auc[i,3] <- tAUC[1,3]
  PPA.df.auc[i,4] <- tAUC[1,4]
  PPA.df.auc[i,5] <- tAUC[1,5]
  PPA.df.auc[i,6] <- tAUC[1,6]
  
  ##################################################################
  #EMD
  r.real <- raster(ncol=200, nrow=200, res=10)
  extent(r.real) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.real[] <- 0
  tab <- table(cellFromXY(r.real, pts2))
  r.real[as.numeric(names(tab))] <- tab
  r.real <- raster::as.matrix(r.real)
  r.real <- pp(r.real)
  
  #T-locoh
  r.100 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(T.100.shp, sp)
  extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.100[] <- 0
  tab <- table(cellFromXY(r.100, e))
  r.100[as.numeric(names(tab))] <- tab
  r.100 <- raster::as.matrix(r.100)
  r.100 <- pp(r.100)
  
  r.50 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(T.50.shp, sp)
  extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.50[] <- 0
  tab <- table(cellFromXY(r.50, e))
  r.50[as.numeric(names(tab))] <- tab
  r.50 <- raster::as.matrix(r.50)
  r.50 <- pp(r.50)
  
  r.25 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(T.25.shp, sp)
  extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.25[] <- 0
  tab <- table(cellFromXY(r.25, e))
  r.25[as.numeric(names(tab))] <- tab
  r.25 <- raster::as.matrix(r.25)
  r.25 <- pp(r.25)
  
  r.10 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(T.10.shp, sp)
  extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.10[] <- 0
  tab <- table(cellFromXY(r.10, e))
  r.10[as.numeric(names(tab))] <- tab
  r.10 <- raster::as.matrix(r.10)
  r.10 <- pp(r.10)
  
  r.05 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(T.05.shp, sp)
  extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.05[] <- 0
  tab <- table(cellFromXY(r.05, e))
  r.05[as.numeric(names(tab))] <- tab
  r.05 <- raster::as.matrix(r.05)
  r.05 <- pp(r.05)

  TLoCoHo.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
  TLoCoHo.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
  TLoCoHo.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
  TLoCoHo.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
  TLoCoHo.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)
  
  #tdge
  r.100 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(tdge.100, sp)
  extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.100[] <- 0
  tab <- table(cellFromXY(r.100, e))
  r.100[as.numeric(names(tab))] <- tab
  r.100 <- raster::as.matrix(r.100)
  r.100 <- pp(r.100)
  
  r.50 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(tdge.50, sp)
  extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.50[] <- 0
  tab <- table(cellFromXY(r.50, e))
  r.50[as.numeric(names(tab))] <- tab
  r.50 <- raster::as.matrix(r.50)
  r.50 <- pp(r.50)
  
  r.25 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(tdge.25, sp)
  extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.25[] <- 0
  tab <- table(cellFromXY(r.25, e))
  r.25[as.numeric(names(tab))] <- tab
  r.25 <- raster::as.matrix(r.25)
  r.25 <- pp(r.25)
  
  r.10 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(tdge.10, sp)
  extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.10[] <- 0
  tab <- table(cellFromXY(r.10, e))
  r.10[as.numeric(names(tab))] <- tab
  r.10 <- raster::as.matrix(r.10)
  r.10 <- pp(r.10)
  
  r.05 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(tdge.05, sp)
  extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.05[] <- 0
  tab <- table(cellFromXY(r.05, e))
  r.05[as.numeric(names(tab))] <- tab
  r.05 <- raster::as.matrix(r.05)
  r.05 <- pp(r.05)
  
  tdge.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
  tdge.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
  tdge.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
  tdge.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
  tdge.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)
  
  #BRB
  r.100 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(BRB.100, sp)
  extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.100[] <- 0
  tab <- table(cellFromXY(r.100, e))
  r.100[as.numeric(names(tab))] <- tab
  r.100 <- raster::as.matrix(r.100)
  r.100 <- pp(r.100)
  
  r.50 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(BRB.50, sp)
  extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.50[] <- 0
  tab <- table(cellFromXY(r.50, e))
  r.50[as.numeric(names(tab))] <- tab
  r.50 <- raster::as.matrix(r.50)
  r.50 <- pp(r.50)
  
  r.25 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(BRB.25, sp)
  extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.25[] <- 0
  tab <- table(cellFromXY(r.25, e))
  r.25[as.numeric(names(tab))] <- tab
  r.25 <- raster::as.matrix(r.25)
  r.25 <- pp(r.25)
  
  r.10 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(BRB.10, sp)
  extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.10[] <- 0
  tab <- table(cellFromXY(r.10, e))
  r.10[as.numeric(names(tab))] <- tab
  r.10 <- raster::as.matrix(r.10)
  r.10 <- pp(r.10)
  
  r.05 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(BRB.05, sp)
  extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.05[] <- 0
  tab <- table(cellFromXY(r.05, e))
  r.05[as.numeric(names(tab))] <- tab
  r.05 <- raster::as.matrix(r.05)
  r.05 <- pp(r.05)
  
  BRB.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
  BRB.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
  BRB.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
  BRB.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
  BRB.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)
  
  #PPA
  r.100 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(PPA.100, sp)
  extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.100[] <- 0
  tab <- table(cellFromXY(r.100, e))
  r.100[as.numeric(names(tab))] <- tab
  r.100 <- raster::as.matrix(r.100)
  r.100 <- pp(r.100)
  
  r.50 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(PPA.50, sp)
  extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.50[] <- 0
  tab <- table(cellFromXY(r.50, e))
  r.50[as.numeric(names(tab))] <- tab
  r.50 <- raster::as.matrix(r.50)
  r.50 <- pp(r.50)
  
  r.25 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(PPA.25, sp)
  extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.25[] <- 0
  tab <- table(cellFromXY(r.25, e))
  r.25[as.numeric(names(tab))] <- tab
  r.25 <- raster::as.matrix(r.25)
  r.25 <- pp(r.25)
  
  r.10 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(PPA.10, sp)
  extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.10[] <- 0
  tab <- table(cellFromXY(r.10, e))
  r.10[as.numeric(names(tab))] <- tab
  r.10 <- raster::as.matrix(r.10)
  r.10 <- pp(r.10)
  
  r.05 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(PPA.05, sp)
  extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.05[] <- 0
  tab <- table(cellFromXY(r.05, e))
  r.05[as.numeric(names(tab))] <- tab
  r.05 <- raster::as.matrix(r.05)
  r.05 <- pp(r.05)
  
  PPA.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
  PPA.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
  PPA.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
  PPA.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
  PPA.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)
  
  #MCP
  r.100 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(mcp.100, sp)
  extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.100[] <- 0
  tab <- table(cellFromXY(r.100, e))
  r.100[as.numeric(names(tab))] <- tab
  r.100 <- raster::as.matrix(r.100)
  r.100 <- pp(r.100)
  
  r.50 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(mcp.50, sp)
  extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.50[] <- 0
  tab <- table(cellFromXY(r.50, e))
  r.50[as.numeric(names(tab))] <- tab
  r.50 <- raster::as.matrix(r.50)
  r.50 <- pp(r.50)
  
  r.25 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(mcp.25, sp)
  extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.25[] <- 0
  tab <- table(cellFromXY(r.25, e))
  r.25[as.numeric(names(tab))] <- tab
  r.25 <- raster::as.matrix(r.25)
  r.25 <- pp(r.25)
  
  r.10 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(mcp.10, sp)
  extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.10[] <- 0
  tab <- table(cellFromXY(r.10, e))
  r.10[as.numeric(names(tab))] <- tab
  r.10 <- raster::as.matrix(r.10)
  r.10 <- pp(r.10)
  
  r.05 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(mcp.05, sp)
  extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.05[] <- 0
  tab <- table(cellFromXY(r.05, e))
  r.05[as.numeric(names(tab))] <- tab
  r.05 <- raster::as.matrix(r.05)
  r.05 <- pp(r.05)
  
  mcp.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
  mcp.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
  mcp.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
  mcp.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
  mcp.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)
  
  #KDE
  r.100 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(KDE.100, sp)
  extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.100[] <- 0
  tab <- table(cellFromXY(r.100, e))
  r.100[as.numeric(names(tab))] <- tab
  r.100 <- raster::as.matrix(r.100)
  r.100 <- pp(r.100)
  
  r.50 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(KDE.50, sp)
  extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.50[] <- 0
  tab <- table(cellFromXY(r.50, e))
  r.50[as.numeric(names(tab))] <- tab
  r.50 <- raster::as.matrix(r.50)
  r.50 <- pp(r.50)
  
  r.25 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(KDE.25, sp)
  extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.25[] <- 0
  tab <- table(cellFromXY(r.25, e))
  r.25[as.numeric(names(tab))] <- tab
  r.25 <- raster::as.matrix(r.25)
  r.25 <- pp(r.25)
  
  r.10 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(KDE.10, sp)
  extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.10[] <- 0
  tab <- table(cellFromXY(r.10, e))
  r.10[as.numeric(names(tab))] <- tab
  r.10 <- raster::as.matrix(r.10)
  r.10 <- pp(r.10)
  
  r.05 <- raster(ncol=200, nrow=200, res=10)
  e <- gIntersection(KDE.05, sp)
  extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
  r.05[] <- 0
  tab <- table(cellFromXY(r.05, e))
  r.05[as.numeric(names(tab))] <- tab
  r.05 <- raster::as.matrix(r.05)
  r.05 <- pp(r.05)
  
  KDE.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
  KDE.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
  KDE.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
  KDE.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
  KDE.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)
  
  
}
write.xlsx(mcp.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\mcp_area.xlsx') 
write.xlsx(mcp.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\mcp_edge.xlsx') 
write.xlsx(mcp.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\mcp_auc.xlsx') 
write.xlsx(mcp.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\mcp_emd.xlsx') 
write.xlsx(KDE.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\KDE_area.xlsx') 
write.xlsx(KDE.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\KDE_edge.xlsx') 
write.xlsx(KDE.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\KDE_auc.xlsx') 
write.xlsx(KDE.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\KDE_emd.xlsx') 
write.xlsx(PPA.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\PPA_area.xlsx') 
write.xlsx(PPA.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\PPA_edge.xlsx') 
write.xlsx(PPA.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\PPA_auc.xlsx') 
write.xlsx(PPA.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\PPA_emd.xlsx') 
write.xlsx(tdge.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\tdge_area.xlsx') 
write.xlsx(tdge.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\tdge_edge.xlsx') 
write.xlsx(tdge.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\tdge_auc.xlsx') 
write.xlsx(tdge.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\tdge_emd.xlsx')
write.xlsx(TLoCoHo.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\TLoCoHo_area.xlsx') 
write.xlsx(TLoCoHo.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\TLoCoHo_edge.xlsx') 
write.xlsx(TLoCoHo.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\TLoCoHo_auc.xlsx') 
write.xlsx(TLoCoHo.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\TLoCoHo_emd.xlsx') 
write.xlsx(BRB.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\BRB_area.xlsx') 
write.xlsx(BRB.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\BRB_edge.xlsx') 
write.xlsx(BRB.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\BRB_auc.xlsx') 
write.xlsx(BRB.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\perforated\\BRB_emd.xlsx') 
################################################################################
###Concave HomeRanges
nclust <- function(x0, y0, radius, n) {
  return(runifdisc(n, radius, centre=c(x0, y0)))
}

#Create Dataframe to save area calculations for MCP
mcp.df <- data.frame(matrix(ncol = 4, nrow = 100))
mcp.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("MCP_Area50%", "MCP_Area 25%", "MCP_Area 10%", "MCP_Area 05%")
colnames(mcp.df) <- x
x <- c("Obs", "MCP_AUC 100%","MCP_AUC 50%", "MCP_AUC 25%", "MCP_AUC 10%", "MCP_AUC 05%")
colnames(mcp.df.auc) <- x
mcp.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("MCP_EMD 100%","MCP_EMD 50%", "MCP_EMD 25%", "MCP_EMD 10%", "MCP_EMD 05%")
colnames(mcp.df.emd) <- x
mcp.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("MCP_edge 100%", "100% Patches","MCP_edge 50%", "50% Patches","MCP_edge 25%", "25% Patches", "MCP_edge 10%", "10% Patches", "MCP_edge 05%", "05% Patches")
colnames(mcp.df.edge) <- x

#Create Dataframe to save area calculations for KDE
KDE.df <- data.frame(matrix(ncol = 4, nrow = 100))
KDE.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("KDE_Area50%", "KDE_Area 25%", "KDE_Area 10%", "KDE_Area 05%")
colnames(KDE.df) <- x
x <- c("Obs", "KDE_AUC 100%","KDE_AUC 50%", "KDE_AUC 25%", "KDE_AUC 10%", "KDE_AUC 05%")
colnames(KDE.df.auc) <- x
KDE.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("KDE_EMD 100%","KDE_EMD 50%", "KDE_EMD 25%", "KDE_EMD 10%", "KDE_EMD 05%")
colnames(KDE.df.emd) <- x
KDE.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("KDE_edge 100%", "100% Patches","KDE_edge 50%", "50% Patches","KDE_edge 25%", "25% Patches", "KDE_edge 10%", "10% Patches", "KDE_edge 05%", "05% Patches")
colnames(KDE.df.edge) <- x

#Create Dataframe to save area calculations for PPA
PPA.df <- data.frame(matrix(ncol = 4, nrow = 100))
PPA.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("PPA_Area50%", "PPA_Area 25%", "PPA_Area 10%", "PPA_Area 05%")
colnames(PPA.df) <- x
x <- c("Obs", "PPA_AUC 100%","PPA_AUC 50%", "PPA_AUC 25%", "PPA_AUC 10%", "PPA_AUC 05%")
colnames(PPA.df.auc) <- x
PPA.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("PPA_EMD 100%","PPA_EMD 50%", "PPA_EMD 25%", "PPA_EMD 10%", "PPA_EMD 05%")
colnames(PPA.df.emd) <- x
PPA.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("PPA_edge 100%", "100% Patches","PPA_edge 50%", "50% Patches","PPA_edge 25%", "25% Patches", "PPA_edge 10%", "10% Patches", "PPA_edge 05%", "05% Patches")
colnames(PPA.df.edge) <- x

#Create Dataframe to save area calculations for TDGE
tdge.df <- data.frame(matrix(ncol = 4, nrow = 100))
tdge.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("TDGE_Area50%", "TDGE_Area 25%", "TDGE_Area 10%", "TDGE_Area 05%")
colnames(tdge.df) <- x
x <- c("Obs", "TDGE_AUC 100%","TDGE_AUC 50%", "TDGE_AUC 25%", "TDGE_AUC 10%", "TDGE_AUC 05%")
colnames(tdge.df.auc) <- x
tdge.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("tdge_EMD 100%","tdge_EMD 50%", "tdge_EMD 25%", "tdge_EMD 10%", "tdge_EMD 05%")
colnames(tdge.df.emd) <- x
tdge.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("tdge_edge 100%", "100% Patches","tdge_edge 50%", "50% Patches","tdge_edge 25%", "25% Patches", "tdge_edge 10%", "10% Patches", "tdge_edge 05%", "05% Patches")
colnames(tdge.df.edge) <- x

#Create Dataframe to save area calculations for BRB
BRB.df <- data.frame(matrix(ncol = 4, nrow = 100))
BRB.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("BRB_Area50%", "BRB_Area 25%", "BRB_Area 10%", "BRB_Area 05%")
colnames(BRB.df) <- x
x <- c("Obs", "BRB_AUC 100%","BRB_AUC 50%", "BRB_AUC 25%", "BRB_AUC 10%", "BRB_AUC 05%")
colnames(BRB.df.auc) <- x
BRB.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("BRB_EMD 100%","BRB_EMD 50%", "BRB_EMD 25%", "BRB_EMD 10%", "BRB_EMD 05%")
colnames(BRB.df.emd) <- x
BRB.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("BRB_edge 100%", "100% Patches","BRB_edge 50%", "50% Patches","BRB_edge 25%", "25% Patches", "BRB_edge 10%", "10% Patches", "BRB_edge 05%", "05% Patches")
colnames(BRB.df.edge) <- x

#Create Dataframe to save area calculations for T-LoCoHO
TLoCoHo.df <- data.frame(matrix(ncol = 4, nrow = 100))
TLoCoHo.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("T-LoCoHO_Area50%", "T-LoCoHO_Area 25%", "T-LoCoHO_Area 10%", "T-LoCoHO_Area 05%")
colnames(TLoCoHo.df) <- x
x <- c("Obs", "T-LoCoH_AUC 100%","T-LoCoH_AUC 50%", "T-LoCoH_AUC 25%", "T-LoCoH_AUC 10%", "T-LoCoH_AUC 05%")
colnames(TLoCoHo.df.auc) <- x
TLoCoHo.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("TLoCoHo_EMD 100%","TLoCoHo_EMD 50%", "TLoCoHo_EMD 25%", "TLoCoHo_EMD 10%", "TLoCoHo_EMD 05%")
colnames(TLoCoHo.df.emd) <- x
TLoCoHo.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("TLoCoHo_edge 100%", "100% Patches","TLoCoHo_edge 50%", "50% Patches","TLoCoHo_edge 25%", "25% Patches", "TLoCoHo_edge 10%", "10% Patches", "TLoCoHo_edge 05%", "05% Patches")
colnames(TLoCoHo.df.edge) <- x

i = 1
for (i in 1:4) {

#pc = point cluster
pc = rPoissonCluster(10, 0.05, nclust, radius=0.2, n=300)
# number of points
ptNumber = pc$n
#Convert to dataframe
pc.df = as.data.frame(pc)
#Convert to spatial Points
sp = SpatialPoints(pc.df )
#WGS 84
projection(sp) = CRS("+init=epsg:4326")
#Plots points
sp<- spTransform(sp, CRS("+init=epsg:32631"))
#Create sample for MCP
s <- sample(sp, 10)
#MCP
cp <- mcp(s, percent = 95)
#Clip points by MCP
allPts<- gDifference(sp, cp)

################################################################
#Set Up df for Traj
ptNumber = length(allPts)
time = seq(as.Date("2000/1/1"), by = "day", length.out = ptNumber)
allPts$Date = time
allPts$id = "Brock"
allPts.df = as.data.frame(allPts)
#Set up data for AUC
#Create Hex grid for present/absence counts
ptsBuffer = gBuffer(sp, width = 4000)
hexPts <-spsample(ptsBuffer,type="hexagonal",cellsize=5000)
hexPols <- HexPoints2SpatialPolygons(hexPts)
#Define true for AUC
pts2 <- allPts.df[, -c(1:2)] # delete columns 5 through 7
pts2 <- SpatialPoints(pts2)
projection(pts2) = CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pts2, hexPols)))
#MCP
pred.df.MCP <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.MCP ) <- "GridNumber"
pred.df.MCP$GridNumber <- 1:length(hexPols)
pred.df.MCP$Obs <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$Obs[is.na(pred.df.MCP$Obs)] <- 0
pred.df.MCP$Presence <- pred.df.MCP$Obs > 0
#PPA
pred.df.PPA <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.PPA ) <- "GridNumber"
pred.df.PPA$GridNumber <- 1:length(hexPols)
pred.df.PPA$Obs <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$Obs[is.na(pred.df.PPA$Obs)] <- 0
pred.df.PPA$Presence <- pred.df.PPA$Obs > 0
#tdge
pred.df.tdge <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.tdge ) <- "GridNumber"
pred.df.tdge$GridNumber <- 1:length(hexPols)
pred.df.tdge$Obs <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$Obs[is.na(pred.df.tdge$Obs)] <- 0
pred.df.tdge$Presence <- pred.df.tdge$Obs > 0
#BRB
pred.df.BRB <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.BRB ) <- "GridNumber"
pred.df.BRB$GridNumber <- 1:length(hexPols)
pred.df.BRB$Obs <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$Obs[is.na(pred.df.BRB$Obs)] <- 0
pred.df.BRB$Presence <- pred.df.BRB$Obs > 0
#T-Locoh
pred.df.Tlocoh <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.Tlocoh ) <- "GridNumber"
pred.df.Tlocoh$GridNumber <- 1:length(hexPols)
pred.df.Tlocoh$Obs <- as.numeric(res$Freq[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$Obs[is.na(pred.df.Tlocoh$Obs)] <- 0
pred.df.Tlocoh$Presence <- pred.df.Tlocoh$Obs > 0
###################################################################

#########################End Point Pattern Set UP ##############################
# Sample the Points at differ levels
#Sample 50 %
sampleSize = length(allPts) * .50
pts.50 = spsample(allPts,n = sampleSize, "random")
#Sample 25 %
sampleSize = length(allPts) * .25
pts.25 = spsample(allPts,n = sampleSize, "random")
#Sample 10 %
sampleSize = length(allPts) * .10
pts.10 = spsample(allPts,n = sampleSize, "random")
#Sample 5 %
sampleSize = length(allPts) * .05
pts.05 = spsample(allPts,n = sampleSize, "random")
#################################Run Home Ranges###############################
#MCP Home Ranges
mcp.100 = mcp(allPts, percent = 95)
#MCP Core Area
# mcp.100.perforated.50 = mcp(allPts, percent = 50)
#MCP Home Range 50%
mcp.50 = mcp(pts.50, percent = 95)
#MCP Core Area 50%
# mcp.50.perforated.50 = mcp(pts.50, percent = 50)
#MCP Home Range 50%
mcp.25= mcp(pts.25, percent = 95)
#MCP Core Area 50%
# mcp.25.perforated.50 = mcp(pts.25, percent = 50)
#MCP Home Range 50%
mcp.10 = mcp(pts.10, percent = 95)
#MCP Core Area 50%
# mcp.10.perforated.50 = mcp(pts.10, percent = 50)
#MCP Home Range 50%
mcp.05 = mcp(pts.05, percent = 95)
#MCP Core Area 50%
# mcp.05.perforated.50 = mcp(pts.05, percent = 50)
#Perferate the MCP to get correct area
area.100 <- gArea(gDifference(mcp.100,cp))
area.50 <- gArea(mcp.50)
area.25 <- gArea(mcp.25)
area.10 <- gArea(mcp.10)
area.05 <- gArea(mcp.05)
mcp.df[i,1] <- area.50/area.100
mcp.df[i,2] <- area.25/area.100
mcp.df[i,3] <- area.10/area.100
mcp.df[i,4] <- area.05/area.100
mcp.df.edge[i,1] <- (lineLength(mcp.100, byid = FALSE)/area.100) * 1000
mcp.df.edge[i,3] <- (lineLength(mcp.50, byid = FALSE)/area.50) * 1000
mcp.df.edge[i,5] <- (lineLength(mcp.25, byid = FALSE)/area.25) * 1000
mcp.df.edge[i,7] <- (lineLength(mcp.10, byid = FALSE)/area.10) * 1000
mcp.df.edge[i,9] <- (lineLength(mcp.05, byid = FALSE)/area.05) * 1000
mcp.df.edge[i,2] <- length(mcp.100)
mcp.df.edge[i,4] <- length(mcp.50)
mcp.df.edge[i,6] <- length(mcp.25)
mcp.df.edge[i,8] <- length(mcp.10)
mcp.df.edge[i,10] <- length(mcp.05)
###################################################################
#AUC MCP
pred.100.mcp <- gIntersection(mcp.100, sp)
projection(pred.100.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100.mcp, hexPols)))
pred.df.MCP$pred.100 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.100[is.na(pred.df.MCP$pred.100)] <- 0
#reorder dataframe
pred.df.MCP <- pred.df.MCP[c(1,3,2,4)]

pred.50.mcp <- gIntersection(mcp.50, sp)
projection(pred.50.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50.mcp, hexPols)))
pred.df.MCP$pred.50 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.50[is.na(pred.df.MCP$pred.50)] <- 0

pred.25.mcp <- gIntersection(mcp.25, sp)
projection(pred.25.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25.mcp, hexPols)))
pred.df.MCP$pred.25 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.25[is.na(pred.df.MCP$pred.25)] <- 0

pred.10.mcp <- gIntersection(mcp.10, sp)
projection(pred.10.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10.mcp, hexPols)))
pred.df.MCP$pred.10 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.10[is.na(pred.df.MCP$pred.10)] <- 0

pred.05.mcp <- gIntersection(mcp.05, sp)
projection(pred.05.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05.mcp, hexPols)))
pred.df.MCP$pred.05 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.05[is.na(pred.df.MCP$pred.05)] <- 0

tAUC <- colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
mcp.df.auc[i,1] <- tAUC[1,1]
mcp.df.auc[i,2] <- tAUC[1,2]
mcp.df.auc[i,3] <- tAUC[1,3]
mcp.df.auc[i,4] <- tAUC[1,4]
mcp.df.auc[i,5] <- tAUC[1,5]
mcp.df.auc[i,6] <- tAUC[1,6]

##################################################################
#KDE Homerange
KDE.100 <- kernelUD(allPts, h="href")
KDE.100 <- getverticeshr(KDE.100,percent = 95)
KDE.50 <- kernelUD(pts.50, h="href")
KDE.50 <- getverticeshr(KDE.50,percent = 95)
KDE.25 <- kernelUD(pts.25, h="href")
KDE.25 <- getverticeshr(KDE.25,percent = 95)
KDE.10 <- kernelUD(pts.10, h="href")
KDE.10 <- getverticeshr(KDE.10,percent = 95)
KDE.05 <- kernelUD(pts.05, h="href")
KDE.05 <- getverticeshr(KDE.05,percent = 95)
#KDE Area
area.100 <- KDE.100$area
area.50 <- KDE.50$area
area.25 <- KDE.25$area
area.10 <- KDE.10$area
area.05 <- KDE.05$area
#Assign Area percents to datafrmae
KDE.df[i,1] <- area.50/area.100
KDE.df[i,2] <- area.25/area.100
KDE.df[i,3] <- area.10/area.100
KDE.df[i,4] <- area.05/area.100
KDE.df.edge[i,1] <- (lineLength(KDE.100, byid = FALSE)/area.100) * 1000
KDE.df.edge[i,3] <- (lineLength(KDE.50, byid = FALSE)/area.50) * 1000
KDE.df.edge[i,5] <- (lineLength(KDE.25, byid = FALSE)/area.25) * 1000
KDE.df.edge[i,7] <- (lineLength(KDE.10, byid = FALSE)/area.10) * 1000
KDE.df.edge[i,9] <- (lineLength(KDE.05, byid = FALSE)/area.05) * 1000
KDE.df.edge[i,2] <- length(KDE.100)
KDE.df.edge[i,4] <- length(KDE.50)
KDE.df.edge[i,6] <- length(KDE.25)
KDE.df.edge[i,8] <- length(KDE.10)
KDE.df.edge[i,10] <- length(KDE.05)
###################################################################
#AUC KDE
#Define true AUC
pts2 <- allPts.df[, -c(1:2)] # delete columns 5 through 7
pts2 <- SpatialPoints(pts2)
projection(pts2) = CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pts2, hexPols)))
#KDE
pred.df.KDE <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.KDE) <- "GridNumber"
pred.df.KDE$GridNumber <- 1:length(hexPols)
pred.df.KDE$Obs <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$Obs[is.na(pred.df.KDE$Obs)] <- 0
pred.df.KDE$Presence <- pred.df.KDE$Obs > 0
pred.df.KDE <- pred.df.KDE[c(1,3,2)]

#100 %
pred.100.KDE <- gIntersection(KDE.100, sp)
projection(pred.100.KDE) <- CRS("+init=epsg:32631")
plot(pred.100.KDE, add = TRUE)
res <- as.data.frame(table(over(pred.100.KDE, hexPols)))
pred.df.KDE$pred.100 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.100[is.na(pred.df.KDE$pred.100)] <- 0


#50%
pred.50.KDE <- gIntersection(KDE.50, sp)
projection(pred.50.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50.KDE, hexPols)))
pred.df.KDE$pred.50 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.50[is.na(pred.df.KDE$pred.50)] <- 0


#25%
pred.25.KDE <- gIntersection(KDE.25, sp)
projection(pred.25.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25.KDE, hexPols)))
pred.df.KDE$pred.25 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.25[is.na(pred.df.KDE$pred.25)] <- 0


#10%
pred.10.KDE <- gIntersection(KDE.10, sp)
projection(pred.10.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10.KDE, hexPols)))
pred.df.KDE$pred.10 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.10[is.na(pred.df.KDE$pred.10)] <- 0

#05%
pred.05.KDE <- gIntersection(KDE.05, sp)
projection(pred.05.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05.KDE, hexPols)))
pred.df.KDE$pred.05 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.05[is.na(pred.df.KDE$pred.05)] <- 0

#auc
tAUC <- colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(sp,KDE.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
KDE.df.auc[i,1] <- tAUC[1,1]
KDE.df.auc[i,2] <- tAUC[1,2]
KDE.df.auc[i,3] <- tAUC[1,3]
KDE.df.auc[i,4] <- tAUC[1,4]
KDE.df.auc[i,5] <- tAUC[1,5]
KDE.df.auc[i,6] <- tAUC[1,6]

######################################################################################################
#PPA, TDGE, BRB, AND t-LOCOCH  Homerange
#Converting to ltraj for time dynamics
PPA.100 <- as.ltraj(xy = allPts.df[,c("x","y")], date = as.POSIXct(strptime(allPts.df$Date,"%Y-%m-%d"))
                    , id = allPts.df$id, proj4string = CRS("+init=epsg:32631"))
tdge.100 <- tgkde(PPA.100)
tdge.100 <- volras(tdge.100,95)
projection(tdge.100) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.100, Tmax=180*60, Lmin=60)
#Hmin is super important here. 
BRB.100 <- BRB(PPA.100, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.100 <- getverticeshr(BRB.100,percent = 95)
projection(BRB.100 ) = CRS("+init=epsg:32631")
PPA.100 <- dynppa(PPA.100)
projection(PPA.100) = CRS("+init=epsg:32631")
#T-locoh
T.100 <- xyt.lxy(xy = allPts.df[,c("x","y")], dt= as.POSIXct(strptime(allPts.df$Date,"%Y-%m-%d"))
                 , proj4string= CRS("+init=epsg:32631"), id=allPts$id)
toni.lxy <- lxy.ptsh.add(T.100)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.100 <- lxy.nn.add(T.100, s=s.value, k=25)
T.100 <- lxy.lhs(T.100, k=3*2:8, s=s.value)
T.100 <- lhs.iso.add(T.100)
T.100 <- isopleths(T.100)
T.100.shp <- T.100[[1]][T.100[[1]]$iso.level == 0.95,]
#####################################################################
#AUC
#T-locoh
pred.100 <- gIntersection(T.100.shp, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.Tlocoh$pred.100 <- as.numeric(res$Freq[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.100[is.na(pred.df.Tlocoh$pred.100)] <- 0
#Tdge
pred.100 <- gIntersection(tdge.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.tdge$pred.100 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.100[is.na(pred.df.tdge$pred.100)] <- 0
#BRB
pred.100 <- gIntersection(BRB.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.BRB$pred.100 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.100[is.na(pred.df.BRB$pred.100)] <- 0
#PPA
pred.100 <- gIntersection(PPA.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.PPA$pred.100 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.100[is.na(pred.df.PPA$pred.100)] <- 0


#####################################################################
#Even Sample
PPA.50 <- allPts.df[seq(1, NROW(allPts.df), by = 2),]
#Need tp run xyt before converting to PPA
T.50 <- xyt.lxy(xy = PPA.50[,c("x","y")], dt= as.POSIXct(strptime(PPA.50$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.50$id)
PPA.50 <- as.ltraj(xy = PPA.50[,c("x","y")], date = as.POSIXct(strptime(PPA.50$Date,"%Y-%m-%d"))
                   , id = PPA.50$id, proj4string = CRS("+init=epsg:32631"))
tdge.50 <- tgkde(PPA.50)
tdge.50 <- volras(tdge.50,95)
projection(tdge.50) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.50, Tmax=180*60, Lmin=60)
BRB.50 <- BRB(PPA.50, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.50 <- getverticeshr(BRB.50,percent = 95)
projection(BRB.50) = CRS("+init=epsg:32631")
PPA.50 <- dynppa(PPA.50)
projection(PPA.50 ) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.50)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.50 <- lxy.nn.add(T.50, s=s.value, k=25)
T.50 <- lxy.lhs(T.50, k=3*2:8, s=s.value)
T.50 <- lhs.iso.add(T.50)
T.50 <- isopleths(T.50)
#Get isopleth
T.50.shp <- T.50[[1]][T.50[[1]]$iso.level == 0.95,]
#######################################################################
#AUC
#T-locoh
pred.50 <- gIntersection(T.50.shp, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.Tlocoh$pred.50 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.50[is.na(pred.df.Tlocoh$pred.50)] <- 0
#Tdge
pred.50 <- gIntersection(tdge.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.tdge$pred.50 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.50[is.na(pred.df.tdge$pred.50)] <- 0
#BRB
pred.50 <- gIntersection(BRB.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.BRB$pred.50 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.50[is.na(pred.df.BRB$pred.50)] <- 0
#PPA
pred.50 <- gIntersection(PPA.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.PPA$pred.50 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.50[is.na(pred.df.PPA$pred.50)] <- 0
########################################################################
PPA.25 <- allPts.df[seq(1, NROW(allPts.df), by = 4),]
T.25 <- xyt.lxy(xy = PPA.25[,c("x","y")], dt= as.POSIXct(strptime(PPA.25$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.25$id)
PPA.25 <- as.ltraj(xy = PPA.25[,c("x","y")], date = as.POSIXct(strptime(PPA.25$Date,"%Y-%m-%d"))
                   , id = PPA.25$id, proj4string = CRS("+init=epsg:32631"))
tdge.25 <- tgkde(PPA.25)
tdge.25 <- volras(tdge.25,95)
projection(tdge.25) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.25, Tmax=180*60, Lmin=60)
BRB.25 <- BRB(PPA.25, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.25 <- getverticeshr(BRB.25,percent = 95)
projection(BRB.25) = CRS("+init=epsg:32631")
PPA.25 <- dynppa(PPA.25)
projection(PPA.25) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.25)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.25 <- lxy.nn.add(T.25, s=s.value, k=25)
T.25 <- lxy.lhs(T.25, k=3*2:8, s=s.value)
T.25 <- lhs.iso.add(T.25)
T.25 <- isopleths(T.25)
T.25.shp <- T.25[[1]][T.25[[1]]$iso.level == 0.95,]
#######################################################################
pred.25 <- gIntersection(T.25.shp, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.Tlocoh$pred.25 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.25[is.na(pred.df.Tlocoh$pred.25)] <- 0
#Tdge
pred.25 <- gIntersection(tdge.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.tdge$pred.25 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.25[is.na(pred.df.tdge$pred.25)] <- 0
#BRB
pred.25 <- gIntersection(BRB.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.BRB$pred.25 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.25[is.na(pred.df.BRB$pred.25)] <- 0
#PPA
pred.25 <- gIntersection(PPA.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.PPA$pred.25 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.25[is.na(pred.df.PPA$pred.25)] <- 0
########################################################################
PPA.10 <- allPts.df[seq(1, NROW(allPts.df), by = 10),]
T.10 <- xyt.lxy(xy = PPA.10[,c("x","y")], dt= as.POSIXct(strptime(PPA.10$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.10$id)
PPA.10<- as.ltraj(xy = PPA.10[,c("x","y")], date = as.POSIXct(strptime(PPA.10$Date,"%Y-%m-%d"))
                  , id = PPA.10$id, proj4string = CRS("+init=epsg:32631"))
tdge.10 <- tgkde(PPA.10)
tdge.10 <- volras(tdge.10,95)
projection(tdge.10 ) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.10, Tmax=180*60, Lmin=60)
BRB.10 <- BRB(PPA.10, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.10 <- getverticeshr(BRB.10,percent = 95)
projection(BRB.10) = CRS("+init=epsg:32631")
PPA.10 <- dynppa(PPA.10)
projection(PPA.10 ) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.10)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.10 <- lxy.nn.add(T.10, s=s.value, k=25)
T.10 <- lxy.lhs(T.10, k=3*2:8, s=s.value)
T.10 <- lhs.iso.add(T.10)
T.10 <- isopleths(T.10)
T.10.shp <- T.10[[1]][T.10[[1]]$iso.level == 0.95,]
#######################################################################
pred.10 <- gIntersection(T.10.shp, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.Tlocoh$pred.10 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.10[is.na(pred.df.Tlocoh$pred.10)] <- 0
#Tdge
pred.10 <- gIntersection(tdge.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.tdge$pred.10 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.10[is.na(pred.df.tdge$pred.10)] <- 0
#BRB
pred.10 <- gIntersection(BRB.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.BRB$pred.10 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.10[is.na(pred.df.BRB$pred.10)] <- 0
#PPA
pred.10 <- gIntersection(PPA.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.PPA$pred.10 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.10[is.na(pred.df.PPA$pred.10)] <- 0
########################################################################
PPA.05 <- allPts.df[seq(1, NROW(allPts.df), by = 4),]
T.05 <- xyt.lxy(xy = PPA.05[,c("x","y")], dt= as.POSIXct(strptime(PPA.05$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.05$id)
PPA.05 <- as.ltraj(xy = PPA.05[,c("x","y")], date = as.POSIXct(strptime(PPA.05$Date,"%Y-%m-%d"))
                   , id = PPA.05$id, proj4string = CRS("+init=epsg:32631"))
tdge.05 <- tgkde(PPA.05)
tdge.05 <- volras(tdge.05,95)
projection(tdge.05) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.05, Tmax=180*60, Lmin=60)
BRB.05 <- BRB(PPA.05, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.05 <- getverticeshr(BRB.05,percent = 95)
projection(BRB.05) = CRS("+init=epsg:32631")
PPA.05 <- dynppa(PPA.05)
projection(PPA.05) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.05)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.05 <- lxy.nn.add(T.05, s=s.value, k=25)
T.05 <- lxy.lhs(T.05, k=3*2:8, s=s.value)
T.05 <- lhs.iso.add(T.05)
T.05 <- isopleths(T.05)
T.05.shp <- T.05[[1]][T.05[[1]]$iso.level == 0.95,]
#######################################################################
pred.05 <- gIntersection(T.05.shp, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.Tlocoh$pred.05 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.05[is.na(pred.df.Tlocoh$pred.05)] <- 0
#Tdge
pred.05 <- gIntersection(tdge.05, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.tdge$pred.05 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.05[is.na(pred.df.tdge$pred.05)] <- 0
#BRB
pred.05 <- gIntersection(BRB.05, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.BRB$pred.05 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.05[is.na(pred.df.BRB$pred.05)] <- 0
#PPA
pred.05 <- gIntersection(PPA.100, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.PPA$pred.05 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.05[is.na(pred.df.PPA$pred.05)] <- 0
########################################################################
#PPA Area
area.100 <- gArea(gDifference(PPA.100,cp))
area.50 <- gArea(PPA.50)
area.25 <- gArea(PPA.25)
area.10 <- gArea(PPA.10)
area.05 <- gArea(PPA.05)
#Edge Complexity 
PPA.df.edge[i,1] <- (lineLength(PPA.100, byid = FALSE)/area.100) * 1000
PPA.df.edge[i,3] <- (lineLength(PPA.50, byid = FALSE)/area.50) * 1000
PPA.df.edge[i,5] <- (lineLength(PPA.25, byid = FALSE)/area.25) * 1000
PPA.df.edge[i,7] <- (lineLength(PPA.10, byid = FALSE)/area.10) * 1000
PPA.df.edge[i,9] <- (lineLength(PPA.05, byid = FALSE)/area.05) * 1000
PPA.df.edge[i,2] <- length(PPA.100)
PPA.df.edge[i,4] <- length(PPA.50)
PPA.df.edge[i,6] <- length(PPA.25)
PPA.df.edge[i,8] <- length(PPA.10)
PPA.df.edge[i,10] <- length(PPA.05)
#Assign Area percents to datafrmae
PPA.df[i,1] <- area.50/area.100
PPA.df[i,2] <- area.25/area.100
PPA.df[i,3] <- area.10/area.100
PPA.df[i,4] <- area.05/area.100
#TDGE Area
area.100 <- gArea(gDifference(tdge.100,cp))
area.50 <- gArea(tdge.50)
area.25 <- gArea(tdge.25)
area.10 <- gArea(tdge.10)
area.05 <- gArea(tdge.05)
#Assign Area percents to datafrmae
tdge.df[i,1] <- area.50/area.100
tdge.df[i,2] <- area.25/area.100
tdge.df[i,3] <- area.10/area.100
tdge.df[i,4] <- area.05/area.100
#Edge Complexity 
tdge.df.edge[i,1] <- (lineLength(tdge.100, byid = FALSE)/area.100) * 1000
tdge.df.edge[i,3] <- (lineLength(tdge.50, byid = FALSE)/area.50) * 1000
tdge.df.edge[i,5] <- (lineLength(tdge.25, byid = FALSE)/area.25) * 1000
tdge.df.edge[i,7] <- (lineLength(tdge.10, byid = FALSE)/area.10) * 1000
tdge.df.edge[i,9] <- (lineLength(tdge.05, byid = FALSE)/area.05) * 1000
tdge.df.edge[i,2] <- length(tdge.100)
tdge.df.edge[i,4] <- length(tdge.50)
tdge.df.edge[i,6] <- length(tdge.25)
tdge.df.edge[i,8] <- length(tdge.10)
tdge.df.edge[i,10] <- length(tdge.05)
#BRB Area
area.100 <- gArea(BRB.100)
area.50 <- gArea(BRB.50)
area.25 <- gArea(BRB.25)
area.10 <- gArea(BRB.10)
area.05 <- gArea(BRB.05)
#Assign Area percents to datafrmae
BRB.df[i,1] <- area.50/area.100
BRB.df[i,2] <- area.25/area.100
BRB.df[i,3] <- area.10/area.100
BRB.df[i,4] <- area.05/area.100
#Edge Complexity 
BRB.df.edge[i,1] <- (lineLength(BRB.100, byid = FALSE)/area.100) * 1000
BRB.df.edge[i,3] <- (lineLength(BRB.50, byid = FALSE)/area.50) * 1000
BRB.df.edge[i,5] <- (lineLength(BRB.25, byid = FALSE)/area.25) * 1000
BRB.df.edge[i,7] <- (lineLength(BRB.10, byid = FALSE)/area.10) * 1000
BRB.df.edge[i,9] <- (lineLength(BRB.05, byid = FALSE)/area.05) * 1000
BRB.df.edge[i,2] <- length(BRB.100)
BRB.df.edge[i,4] <- length(BRB.50)
BRB.df.edge[i,6] <- length(BRB.25)
BRB.df.edge[i,8] <- length(BRB.10)
BRB.df.edge[i,10] <- length(BRB.05)
#T-locoh area
area.100 <- T.100[[1]]@data[5,2]
area.50 <- T.50[[1]]@data[5,2]
area.25 <- T.25[[1]]@data[5,2]
area.10 <- T.10[[1]]@data[5,2]
area.05 <- T.05[[1]]@data[5,2]
#Assign Area percents to datafrmae
TLoCoHo.df[i,1] <- area.50/area.100
TLoCoHo.df[i,2] <- area.25/area.100
TLoCoHo.df[i,3] <- area.10/area.100
TLoCoHo.df[i,4] <- area.05/area.100
#Edge Complexity 
TLoCoHo.df.edge[i,1] <- (lineLength(T.100.shp, byid = FALSE)/area.100) * 1000
TLoCoHo.df.edge[i,3] <- (lineLength(T.50.shp, byid = FALSE)/area.50) * 1000
TLoCoHo.df.edge[i,5] <- (lineLength(T.25.shp, byid = FALSE)/area.25) * 1000
TLoCoHo.df.edge[i,7] <- (lineLength(T.10.shp, byid = FALSE)/area.10) * 1000
TLoCoHo.df.edge[i,9] <- (lineLength(T.05.shp, byid = FALSE)/area.05) * 1000
TLoCoHo.df.edge[i,2] <- length(T.100.shp)
TLoCoHo.df.edge[i,4] <- length(T.50.shp)
TLoCoHo.df.edge[i,6] <- length(T.25.shp)
TLoCoHo.df.edge[i,8] <- length(T.10.shp)
TLoCoHo.df.edge[i,10] <- length(T.05.shp)
#####AUC#############################################
pred.df.Tlocoh <-pred.df.Tlocoh[c(1,3,2,4,5,6,7,8)]
pred.df.tdge <-pred.df.tdge[c(1,3,2,4,5,6,7,8)]
pred.df.BRB <-pred.df.BRB[c(1,3,2,4,5,6,7,8)]
pred.df.PPA <-pred.df.PPA[c(1,3,2,4,5,6,7,8)]
tAUC <- colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
KDE.df.auc[i,1] <- tAUC[1,1]
KDE.df.auc[i,2] <- tAUC[1,2]
KDE.df.auc[i,3] <- tAUC[1,3]
KDE.df.auc[i,4] <- tAUC[1,4]
KDE.df.auc[i,5] <- tAUC[1,5]
KDE.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,MCP.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
mcp.df.auc[i,1] <- tAUC[1,1]
mcp.df.auc[i,2] <- tAUC[1,2]
mcp.df.auc[i,3] <- tAUC[1,3]
mcp.df.auc[i,4] <- tAUC[1,4]
mcp.df.auc[i,5] <- tAUC[1,5]
mcp.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.Tlocoh[,3:8], pred.df.Tlocoh[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.Tlocoh[,3:8], pred.df.Tlocoh[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
TLoCoHo.df.auc[i,1] <- tAUC[1,1]
TLoCoHo.df.auc[i,2] <- tAUC[1,2]
TLoCoHo.df.auc[i,3] <- tAUC[1,3]
TLoCoHo.df.auc[i,4] <- tAUC[1,4]
TLoCoHo.df.auc[i,5] <- tAUC[1,5]
TLoCoHo.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.tdge[,3:8], pred.df.tdge[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.tdge[,3:8], pred.df.tdge[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
tdge.df.auc[i,1] <- tAUC[1,1]
tdge.df.auc[i,2] <- tAUC[1,2]
tdge.df.auc[i,3] <- tAUC[1,3]
tdge.df.auc[i,4] <- tAUC[1,4]
tdge.df.auc[i,5] <- tAUC[1,5]
tdge.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.BRB[,3:8], pred.df.BRB[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.BRB[,3:8], pred.df.BRB[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
BRB.df.auc[i,1] <- tAUC[1,1]
BRB.df.auc[i,2] <- tAUC[1,2]
BRB.df.auc[i,3] <- tAUC[1,3]
BRB.df.auc[i,4] <- tAUC[1,4]
BRB.df.auc[i,5] <- tAUC[1,5]
BRB.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.PPA[,3:8], pred.df.PPA[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.PPA[,3:8], pred.df.PPA[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
PPA.df.auc[i,1] <- tAUC[1,1]
PPA.df.auc[i,2] <- tAUC[1,2]
PPA.df.auc[i,3] <- tAUC[1,3]
PPA.df.auc[i,4] <- tAUC[1,4]
PPA.df.auc[i,5] <- tAUC[1,5]
PPA.df.auc[i,6] <- tAUC[1,6]

##################################################################
#EMD
r.real <- raster(ncol=200, nrow=200, res=10)
extent(r.real) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.real[] <- 0
tab <- table(cellFromXY(r.real, pts2))
r.real[as.numeric(names(tab))] <- tab
r.real <- raster::as.matrix(r.real)
r.real <- pp(r.real)

#T-locoh
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.100.shp, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.50.shp, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.25.shp, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.10.shp, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.05.shp, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

TLoCoHo.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
TLoCoHo.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
TLoCoHo.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
TLoCoHo.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
TLoCoHo.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#tdge
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

tdge.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
tdge.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
tdge.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
tdge.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
tdge.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#BRB
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

BRB.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
BRB.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
BRB.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
BRB.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
BRB.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#PPA
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

PPA.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
PPA.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
PPA.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
PPA.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
PPA.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#MCP
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

mcp.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
mcp.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
mcp.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
mcp.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
mcp.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#KDE
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

KDE.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
KDE.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
KDE.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
KDE.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
KDE.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)


}
write.xlsx(mcp.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\mcp_area.xlsx') 
write.xlsx(mcp.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\mcp_edge.xlsx') 
write.xlsx(mcp.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\mcp_auc.xlsx') 
write.xlsx(mcp.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\mcp_emd.xlsx') 
write.xlsx(KDE.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\KDE_area.xlsx') 
write.xlsx(KDE.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\KDE_edge.xlsx') 
write.xlsx(KDE.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\KDE_auc.xlsx') 
write.xlsx(KDE.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\KDE_emd.xlsx') 
write.xlsx(PPA.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\PPA_area.xlsx') 
write.xlsx(PPA.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\PPA_edge.xlsx') 
write.xlsx(PPA.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\PPA_auc.xlsx') 
write.xlsx(PPA.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\PPA_emd.xlsx') 
write.xlsx(tdge.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\tdge_area.xlsx') 
write.xlsx(tdge.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\tdge_edge.xlsx') 
write.xlsx(tdge.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\tdge_auc.xlsx') 
write.xlsx(tdge.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\tdge_emd.xlsx')
write.xlsx(TLoCoHo.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\TLoCoHo_area.xlsx') 
write.xlsx(TLoCoHo.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\TLoCoHo_edge.xlsx') 
write.xlsx(TLoCoHo.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\TLoCoHo_auc.xlsx') 
write.xlsx(TLoCoHo.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\TLoCoHo_emd.xlsx') 
write.xlsx(BRB.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\BRB_area.xlsx') 
write.xlsx(BRB.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\BRB_edge.xlsx') 
write.xlsx(BRB.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\BRB_auc.xlsx') 
write.xlsx(BRB.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\concave\\BRB_emd.xlsx') 
################################################################################
###Disjoint HomeRanges
nclust <- function(x0, y0, radius, n) {
  return(runifdisc(n, radius, centre=c(x0, y0)))
}

#Create Dataframe to save area calculations for MCP
mcp.df <- data.frame(matrix(ncol = 4, nrow = 100))
mcp.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("MCP_Area50%", "MCP_Area 25%", "MCP_Area 10%", "MCP_Area 05%")
colnames(mcp.df) <- x
x <- c("Obs", "MCP_AUC 100%","MCP_AUC 50%", "MCP_AUC 25%", "MCP_AUC 10%", "MCP_AUC 05%")
colnames(mcp.df.auc) <- x
mcp.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("MCP_EMD 100%","MCP_EMD 50%", "MCP_EMD 25%", "MCP_EMD 10%", "MCP_EMD 05%")
colnames(mcp.df.emd) <- x
mcp.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("MCP_edge 100%", "100% Patches","MCP_edge 50%", "50% Patches","MCP_edge 25%", "25% Patches", "MCP_edge 10%", "10% Patches", "MCP_edge 05%", "05% Patches")
colnames(mcp.df.edge) <- x

#Create Dataframe to save area calculations for KDE
KDE.df <- data.frame(matrix(ncol = 4, nrow = 100))
KDE.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("KDE_Area50%", "KDE_Area 25%", "KDE_Area 10%", "KDE_Area 05%")
colnames(KDE.df) <- x
x <- c("Obs", "KDE_AUC 100%","KDE_AUC 50%", "KDE_AUC 25%", "KDE_AUC 10%", "KDE_AUC 05%")
colnames(KDE.df.auc) <- x
KDE.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("KDE_EMD 100%","KDE_EMD 50%", "KDE_EMD 25%", "KDE_EMD 10%", "KDE_EMD 05%")
colnames(KDE.df.emd) <- x
KDE.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("KDE_edge 100%", "100% Patches","KDE_edge 50%", "50% Patches","KDE_edge 25%", "25% Patches", "KDE_edge 10%", "10% Patches", "KDE_edge 05%", "05% Patches")
colnames(KDE.df.edge) <- x

#Create Dataframe to save area calculations for PPA
PPA.df <- data.frame(matrix(ncol = 4, nrow = 100))
PPA.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("PPA_Area50%", "PPA_Area 25%", "PPA_Area 10%", "PPA_Area 05%")
colnames(PPA.df) <- x
x <- c("Obs", "PPA_AUC 100%","PPA_AUC 50%", "PPA_AUC 25%", "PPA_AUC 10%", "PPA_AUC 05%")
colnames(PPA.df.auc) <- x
PPA.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("PPA_EMD 100%","PPA_EMD 50%", "PPA_EMD 25%", "PPA_EMD 10%", "PPA_EMD 05%")
colnames(PPA.df.emd) <- x
PPA.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("PPA_edge 100%", "100% Patches","PPA_edge 50%", "50% Patches","PPA_edge 25%", "25% Patches", "PPA_edge 10%", "10% Patches", "PPA_edge 05%", "05% Patches")
colnames(PPA.df.edge) <- x

#Create Dataframe to save area calculations for TDGE
tdge.df <- data.frame(matrix(ncol = 4, nrow = 100))
tdge.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("TDGE_Area50%", "TDGE_Area 25%", "TDGE_Area 10%", "TDGE_Area 05%")
colnames(tdge.df) <- x
x <- c("Obs", "TDGE_AUC 100%","TDGE_AUC 50%", "TDGE_AUC 25%", "TDGE_AUC 10%", "TDGE_AUC 05%")
colnames(tdge.df.auc) <- x
tdge.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("tdge_EMD 100%","tdge_EMD 50%", "tdge_EMD 25%", "tdge_EMD 10%", "tdge_EMD 05%")
colnames(tdge.df.emd) <- x
tdge.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("tdge_edge 100%", "100% Patches","tdge_edge 50%", "50% Patches","tdge_edge 25%", "25% Patches", "tdge_edge 10%", "10% Patches", "tdge_edge 05%", "05% Patches")
colnames(tdge.df.edge) <- x

#Create Dataframe to save area calculations for BRB
BRB.df <- data.frame(matrix(ncol = 4, nrow = 100))
BRB.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("BRB_Area50%", "BRB_Area 25%", "BRB_Area 10%", "BRB_Area 05%")
colnames(BRB.df) <- x
x <- c("Obs", "BRB_AUC 100%","BRB_AUC 50%", "BRB_AUC 25%", "BRB_AUC 10%", "BRB_AUC 05%")
colnames(BRB.df.auc) <- x
BRB.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("BRB_EMD 100%","BRB_EMD 50%", "BRB_EMD 25%", "BRB_EMD 10%", "BRB_EMD 05%")
colnames(BRB.df.emd) <- x
BRB.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("BRB_edge 100%", "100% Patches","BRB_edge 50%", "50% Patches","BRB_edge 25%", "25% Patches", "BRB_edge 10%", "10% Patches", "BRB_edge 05%", "05% Patches")
colnames(BRB.df.edge) <- x

#Create Dataframe to save area calculations for T-LoCoHO
TLoCoHo.df <- data.frame(matrix(ncol = 4, nrow = 100))
TLoCoHo.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("T-LoCoHO_Area50%", "T-LoCoHO_Area 25%", "T-LoCoHO_Area 10%", "T-LoCoHO_Area 05%")
colnames(TLoCoHo.df) <- x
x <- c("Obs", "T-LoCoH_AUC 100%","T-LoCoH_AUC 50%", "T-LoCoH_AUC 25%", "T-LoCoH_AUC 10%", "T-LoCoH_AUC 05%")
colnames(TLoCoHo.df.auc) <- x
TLoCoHo.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("TLoCoHo_EMD 100%","TLoCoHo_EMD 50%", "TLoCoHo_EMD 25%", "TLoCoHo_EMD 10%", "TLoCoHo_EMD 05%")
colnames(TLoCoHo.df.emd) <- x
TLoCoHo.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("TLoCoHo_edge 100%", "100% Patches","TLoCoHo_edge 50%", "50% Patches","TLoCoHo_edge 25%", "25% Patches", "TLoCoHo_edge 10%", "10% Patches", "TLoCoHo_edge 05%", "05% Patches")
colnames(TLoCoHo.df.edge) <- x

i = 1
for(i in 1:4){
#WGS 84
pc = rPoissonCluster(5, 0.3, nclust, radius=0.25, n=500)
pc.df = as.data.frame(pc)
#Convert to spatial Points
sp = SpatialPoints(pc.df )
projection(sp) = CRS("+init=epsg:4326")
#Plots points
sp<- spTransform(sp, CRS("+init=epsg:32631"))
#Create sample for MCP
s <- sample(sp, 10)
#MCP
cp <- mcp(s, percent = 95)
#Clip points by MCP
allPts<- gDifference(sp, cp)
################################################################
#Set Up df for Traj
ptNumber = length(allPts)
time = seq(as.Date("2000/1/1"), by = "day", length.out = ptNumber)
allPts$Date = time
allPts$id = "Brock"
allPts.df = as.data.frame(allPts)
#Set up data for AUC
#Create Hex grid for present/absence counts
ptsBuffer = gBuffer(sp, width = 4000)
hexPts <-spsample(ptsBuffer,type="hexagonal",cellsize=5000)
hexPols <- HexPoints2SpatialPolygons(hexPts)
#Define true for AUC
pts2 <- allPts.df[, -c(1:2)] # delete columns 5 through 7
pts2 <- SpatialPoints(pts2)
projection(pts2) = CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pts2, hexPols)))
#MCP
pred.df.MCP <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.MCP ) <- "GridNumber"
pred.df.MCP$GridNumber <- 1:length(hexPols)
pred.df.MCP$Obs <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$Obs[is.na(pred.df.MCP$Obs)] <- 0
pred.df.MCP$Presence <- pred.df.MCP$Obs > 0
#PPA
pred.df.PPA <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.PPA ) <- "GridNumber"
pred.df.PPA$GridNumber <- 1:length(hexPols)
pred.df.PPA$Obs <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$Obs[is.na(pred.df.PPA$Obs)] <- 0
pred.df.PPA$Presence <- pred.df.PPA$Obs > 0
#tdge
pred.df.tdge <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.tdge ) <- "GridNumber"
pred.df.tdge$GridNumber <- 1:length(hexPols)
pred.df.tdge$Obs <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$Obs[is.na(pred.df.tdge$Obs)] <- 0
pred.df.tdge$Presence <- pred.df.tdge$Obs > 0
#BRB
pred.df.BRB <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.BRB ) <- "GridNumber"
pred.df.BRB$GridNumber <- 1:length(hexPols)
pred.df.BRB$Obs <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$Obs[is.na(pred.df.BRB$Obs)] <- 0
pred.df.BRB$Presence <- pred.df.BRB$Obs > 0
#T-Locoh
pred.df.Tlocoh <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.Tlocoh ) <- "GridNumber"
pred.df.Tlocoh$GridNumber <- 1:length(hexPols)
pred.df.Tlocoh$Obs <- as.numeric(res$Freq[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$Obs[is.na(pred.df.Tlocoh$Obs)] <- 0
pred.df.Tlocoh$Presence <- pred.df.Tlocoh$Obs > 0
###################################################################

#########################End Point Pattern Set UP ##############################
# Sample the Points at differ levels
#Sample 50 %
sampleSize = length(allPts) * .50
pts.50 = spsample(allPts,n = sampleSize, "random")
#Sample 25 %
sampleSize = length(allPts) * .25
pts.25 = spsample(allPts,n = sampleSize, "random")
#Sample 10 %
sampleSize = length(allPts) * .10
pts.10 = spsample(allPts,n = sampleSize, "random")
#Sample 5 %
sampleSize = length(allPts) * .05
pts.05 = spsample(allPts,n = sampleSize, "random")
#################################Run Home Ranges###############################
#MCP Home Ranges
mcp.100 = mcp(allPts, percent = 95)
#MCP Core Area
# mcp.100.perforated.50 = mcp(allPts, percent = 50)
#MCP Home Range 50%
mcp.50 = mcp(pts.50, percent = 95)
#MCP Core Area 50%
# mcp.50.perforated.50 = mcp(pts.50, percent = 50)
#MCP Home Range 50%
mcp.25= mcp(pts.25, percent = 95)
#MCP Core Area 50%
# mcp.25.perforated.50 = mcp(pts.25, percent = 50)
#MCP Home Range 50%
mcp.10 = mcp(pts.10, percent = 95)
#MCP Core Area 50%
# mcp.10.perforated.50 = mcp(pts.10, percent = 50)
#MCP Home Range 50%
mcp.05 = mcp(pts.05, percent = 95)
#MCP Core Area 50%
# mcp.05.perforated.50 = mcp(pts.05, percent = 50)
#Perferate the MCP to get correct area
area.100 <- gArea(gDifference(mcp.100,cp))
area.50 <- gArea(mcp.50)
area.25 <- gArea(mcp.25)
area.10 <- gArea(mcp.10)
area.05 <- gArea(mcp.05)
mcp.df[i,1] <- area.50/area.100
mcp.df[i,2] <- area.25/area.100
mcp.df[i,3] <- area.10/area.100
mcp.df[i,4] <- area.05/area.100
mcp.df.edge[i,1] <- (lineLength(mcp.100, byid = FALSE)/area.100) * 1000
mcp.df.edge[i,3] <- (lineLength(mcp.50, byid = FALSE)/area.50) * 1000
mcp.df.edge[i,5] <- (lineLength(mcp.25, byid = FALSE)/area.25) * 1000
mcp.df.edge[i,7] <- (lineLength(mcp.10, byid = FALSE)/area.10) * 1000
mcp.df.edge[i,9] <- (lineLength(mcp.05, byid = FALSE)/area.05) * 1000
mcp.df.edge[i,2] <- length(mcp.100)
mcp.df.edge[i,4] <- length(mcp.50)
mcp.df.edge[i,6] <- length(mcp.25)
mcp.df.edge[i,8] <- length(mcp.10)
mcp.df.edge[i,10] <- length(mcp.05)
###################################################################
#AUC MCP
pred.100.mcp <- gIntersection(mcp.100, sp)
projection(pred.100.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100.mcp, hexPols)))
pred.df.MCP$pred.100 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.100[is.na(pred.df.MCP$pred.100)] <- 0
#reorder dataframe
pred.df.MCP <- pred.df.MCP[c(1,3,2,4)]

pred.50.mcp <- gIntersection(mcp.50, sp)
projection(pred.50.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50.mcp, hexPols)))
pred.df.MCP$pred.50 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.50[is.na(pred.df.MCP$pred.50)] <- 0

pred.25.mcp <- gIntersection(mcp.25, sp)
projection(pred.25.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25.mcp, hexPols)))
pred.df.MCP$pred.25 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.25[is.na(pred.df.MCP$pred.25)] <- 0

pred.10.mcp <- gIntersection(mcp.10, sp)
projection(pred.10.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10.mcp, hexPols)))
pred.df.MCP$pred.10 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.10[is.na(pred.df.MCP$pred.10)] <- 0

pred.05.mcp <- gIntersection(mcp.05, sp)
projection(pred.05.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05.mcp, hexPols)))
pred.df.MCP$pred.05 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.05[is.na(pred.df.MCP$pred.05)] <- 0

tAUC <- colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
mcp.df.auc[i,1] <- tAUC[1,1]
mcp.df.auc[i,2] <- tAUC[1,2]
mcp.df.auc[i,3] <- tAUC[1,3]
mcp.df.auc[i,4] <- tAUC[1,4]
mcp.df.auc[i,5] <- tAUC[1,5]
mcp.df.auc[i,6] <- tAUC[1,6]

##################################################################
#KDE Homerange
KDE.100 <- kernelUD(allPts, h="href")
KDE.100 <- getverticeshr(KDE.100,percent = 95)
KDE.50 <- kernelUD(pts.50, h="href")
KDE.50 <- getverticeshr(KDE.50,percent = 95)
KDE.25 <- kernelUD(pts.25, h="href")
KDE.25 <- getverticeshr(KDE.25,percent = 95)
KDE.10 <- kernelUD(pts.10, h="href")
KDE.10 <- getverticeshr(KDE.10,percent = 95)
KDE.05 <- kernelUD(pts.05, h="href")
KDE.05 <- getverticeshr(KDE.05,percent = 95)
#KDE Area
area.100 <- KDE.100$area
area.50 <- KDE.50$area
area.25 <- KDE.25$area
area.10 <- KDE.10$area
area.05 <- KDE.05$area
#Assign Area percents to datafrmae
KDE.df[i,1] <- area.50/area.100
KDE.df[i,2] <- area.25/area.100
KDE.df[i,3] <- area.10/area.100
KDE.df[i,4] <- area.05/area.100
KDE.df.edge[i,1] <- (lineLength(KDE.100, byid = FALSE)/area.100) * 1000
KDE.df.edge[i,3] <- (lineLength(KDE.50, byid = FALSE)/area.50) * 1000
KDE.df.edge[i,5] <- (lineLength(KDE.25, byid = FALSE)/area.25) * 1000
KDE.df.edge[i,7] <- (lineLength(KDE.10, byid = FALSE)/area.10) * 1000
KDE.df.edge[i,9] <- (lineLength(KDE.05, byid = FALSE)/area.05) * 1000
KDE.df.edge[i,2] <- length(KDE.100)
KDE.df.edge[i,4] <- length(KDE.50)
KDE.df.edge[i,6] <- length(KDE.25)
KDE.df.edge[i,8] <- length(KDE.10)
KDE.df.edge[i,10] <- length(KDE.05)
###################################################################
#AUC KDE
#Define true AUC
pts2 <- allPts.df[, -c(1:2)] # delete columns 5 through 7
pts2 <- SpatialPoints(pts2)
projection(pts2) = CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pts2, hexPols)))
#KDE
pred.df.KDE <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.KDE) <- "GridNumber"
pred.df.KDE$GridNumber <- 1:length(hexPols)
pred.df.KDE$Obs <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$Obs[is.na(pred.df.KDE$Obs)] <- 0
pred.df.KDE$Presence <- pred.df.KDE$Obs > 0
pred.df.KDE <- pred.df.KDE[c(1,3,2)]

#100 %
pred.100.KDE <- gIntersection(KDE.100, sp)
projection(pred.100.KDE) <- CRS("+init=epsg:32631")
plot(pred.100.KDE, add = TRUE)
res <- as.data.frame(table(over(pred.100.KDE, hexPols)))
pred.df.KDE$pred.100 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.100[is.na(pred.df.KDE$pred.100)] <- 0


#50%
pred.50.KDE <- gIntersection(KDE.50, sp)
projection(pred.50.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50.KDE, hexPols)))
pred.df.KDE$pred.50 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.50[is.na(pred.df.KDE$pred.50)] <- 0


#25%
pred.25.KDE <- gIntersection(KDE.25, sp)
projection(pred.25.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25.KDE, hexPols)))
pred.df.KDE$pred.25 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.25[is.na(pred.df.KDE$pred.25)] <- 0


#10%
pred.10.KDE <- gIntersection(KDE.10, sp)
projection(pred.10.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10.KDE, hexPols)))
pred.df.KDE$pred.10 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.10[is.na(pred.df.KDE$pred.10)] <- 0

#05%
pred.05.KDE <- gIntersection(KDE.05, sp)
projection(pred.05.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05.KDE, hexPols)))
pred.df.KDE$pred.05 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.05[is.na(pred.df.KDE$pred.05)] <- 0

#auc
tAUC <- colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(sp,KDE.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
KDE.df.auc[i,1] <- tAUC[1,1]
KDE.df.auc[i,2] <- tAUC[1,2]
KDE.df.auc[i,3] <- tAUC[1,3]
KDE.df.auc[i,4] <- tAUC[1,4]
KDE.df.auc[i,5] <- tAUC[1,5]
KDE.df.auc[i,6] <- tAUC[1,6]

######################################################################################################
#PPA, TDGE, BRB, AND t-LOCOCH  Homerange
#Converting to ltraj for time dynamics
PPA.100 <- as.ltraj(xy = allPts.df[,c("x","y")], date = as.POSIXct(strptime(allPts.df$Date,"%Y-%m-%d"))
                    , id = allPts.df$id, proj4string = CRS("+init=epsg:32631"))
tdge.100 <- tgkde(PPA.100)
tdge.100 <- volras(tdge.100,95)
projection(tdge.100) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.100, Tmax=180*60, Lmin=60)
#Hmin is super important here. 
BRB.100 <- BRB(PPA.100, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.100 <- getverticeshr(BRB.100,percent = 95)
projection(BRB.100 ) = CRS("+init=epsg:32631")
PPA.100 <- dynppa(PPA.100)
projection(PPA.100) = CRS("+init=epsg:32631")
#T-locoh
T.100 <- xyt.lxy(xy = allPts.df[,c("x","y")], dt= as.POSIXct(strptime(allPts.df$Date,"%Y-%m-%d"))
                 , proj4string= CRS("+init=epsg:32631"), id=allPts$id)
toni.lxy <- lxy.ptsh.add(T.100)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.100 <- lxy.nn.add(T.100, s=s.value, k=25)
T.100 <- lxy.lhs(T.100, k=3*2:8, s=s.value)
T.100 <- lhs.iso.add(T.100)
T.100 <- isopleths(T.100)
T.100.shp <- T.100[[1]][T.100[[1]]$iso.level == 0.95,]
#####################################################################
#AUC
#T-locoh
pred.100 <- gIntersection(T.100.shp, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.Tlocoh$pred.100 <- as.numeric(res$Freq[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.100[is.na(pred.df.Tlocoh$pred.100)] <- 0
#Tdge
pred.100 <- gIntersection(tdge.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.tdge$pred.100 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.100[is.na(pred.df.tdge$pred.100)] <- 0
#BRB
pred.100 <- gIntersection(BRB.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.BRB$pred.100 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.100[is.na(pred.df.BRB$pred.100)] <- 0
#PPA
pred.100 <- gIntersection(PPA.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.PPA$pred.100 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.100[is.na(pred.df.PPA$pred.100)] <- 0


#####################################################################
#Even Sample
PPA.50 <- allPts.df[seq(1, NROW(allPts.df), by = 2),]
#Need tp run xyt before converting to PPA
T.50 <- xyt.lxy(xy = PPA.50[,c("x","y")], dt= as.POSIXct(strptime(PPA.50$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.50$id)
PPA.50 <- as.ltraj(xy = PPA.50[,c("x","y")], date = as.POSIXct(strptime(PPA.50$Date,"%Y-%m-%d"))
                   , id = PPA.50$id, proj4string = CRS("+init=epsg:32631"))
tdge.50 <- tgkde(PPA.50)
tdge.50 <- volras(tdge.50,95)
projection(tdge.50) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.50, Tmax=180*60, Lmin=60)
BRB.50 <- BRB(PPA.50, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.50 <- getverticeshr(BRB.50,percent = 95)
projection(BRB.50) = CRS("+init=epsg:32631")
PPA.50 <- dynppa(PPA.50)
projection(PPA.50 ) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.50)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.50 <- lxy.nn.add(T.50, s=s.value, k=25)
T.50 <- lxy.lhs(T.50, k=3*2:8, s=s.value)
T.50 <- lhs.iso.add(T.50)
T.50 <- isopleths(T.50)
#Get isopleth
T.50.shp <- T.50[[1]][T.50[[1]]$iso.level == 0.95,]
#######################################################################
#AUC
#T-locoh
pred.50 <- gIntersection(T.50.shp, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.Tlocoh$pred.50 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.50[is.na(pred.df.Tlocoh$pred.50)] <- 0
#Tdge
pred.50 <- gIntersection(tdge.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.tdge$pred.50 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.50[is.na(pred.df.tdge$pred.50)] <- 0
#BRB
pred.50 <- gIntersection(BRB.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.BRB$pred.50 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.50[is.na(pred.df.BRB$pred.50)] <- 0
#PPA
pred.50 <- gIntersection(PPA.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.PPA$pred.50 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.50[is.na(pred.df.PPA$pred.50)] <- 0
########################################################################
PPA.25 <- allPts.df[seq(1, NROW(allPts.df), by = 4),]
T.25 <- xyt.lxy(xy = PPA.25[,c("x","y")], dt= as.POSIXct(strptime(PPA.25$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.25$id)
PPA.25 <- as.ltraj(xy = PPA.25[,c("x","y")], date = as.POSIXct(strptime(PPA.25$Date,"%Y-%m-%d"))
                   , id = PPA.25$id, proj4string = CRS("+init=epsg:32631"))
tdge.25 <- tgkde(PPA.25)
tdge.25 <- volras(tdge.25,95)
projection(tdge.25) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.25, Tmax=180*60, Lmin=60)
BRB.25 <- BRB(PPA.25, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.25 <- getverticeshr(BRB.25,percent = 95)
projection(BRB.25) = CRS("+init=epsg:32631")
PPA.25 <- dynppa(PPA.25)
projection(PPA.25) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.25)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.25 <- lxy.nn.add(T.25, s=s.value, k=25)
T.25 <- lxy.lhs(T.25, k=3*2:8, s=s.value)
T.25 <- lhs.iso.add(T.25)
T.25 <- isopleths(T.25)
T.25.shp <- T.25[[1]][T.25[[1]]$iso.level == 0.95,]
#######################################################################
pred.25 <- gIntersection(T.25.shp, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.Tlocoh$pred.25 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.25[is.na(pred.df.Tlocoh$pred.25)] <- 0
#Tdge
pred.25 <- gIntersection(tdge.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.tdge$pred.25 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.25[is.na(pred.df.tdge$pred.25)] <- 0
#BRB
pred.25 <- gIntersection(BRB.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.BRB$pred.25 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.25[is.na(pred.df.BRB$pred.25)] <- 0
#PPA
pred.25 <- gIntersection(PPA.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.PPA$pred.25 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.25[is.na(pred.df.PPA$pred.25)] <- 0
########################################################################
PPA.10 <- allPts.df[seq(1, NROW(allPts.df), by = 10),]
T.10 <- xyt.lxy(xy = PPA.10[,c("x","y")], dt= as.POSIXct(strptime(PPA.10$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.10$id)
PPA.10<- as.ltraj(xy = PPA.10[,c("x","y")], date = as.POSIXct(strptime(PPA.10$Date,"%Y-%m-%d"))
                  , id = PPA.10$id, proj4string = CRS("+init=epsg:32631"))
tdge.10 <- tgkde(PPA.10)
tdge.10 <- volras(tdge.10,95)
projection(tdge.10 ) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.10, Tmax=180*60, Lmin=60)
BRB.10 <- BRB(PPA.10, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.10 <- getverticeshr(BRB.10,percent = 95)
projection(BRB.10) = CRS("+init=epsg:32631")
PPA.10 <- dynppa(PPA.10)
projection(PPA.10 ) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.10)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.10 <- lxy.nn.add(T.10, s=s.value, k=25)
T.10 <- lxy.lhs(T.10, k=3*2:8, s=s.value)
T.10 <- lhs.iso.add(T.10)
T.10 <- isopleths(T.10)
T.10.shp <- T.10[[1]][T.10[[1]]$iso.level == 0.95,]
#######################################################################
pred.10 <- gIntersection(T.10.shp, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.Tlocoh$pred.10 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.10[is.na(pred.df.Tlocoh$pred.10)] <- 0
#Tdge
pred.10 <- gIntersection(tdge.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.tdge$pred.10 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.10[is.na(pred.df.tdge$pred.10)] <- 0
#BRB
pred.10 <- gIntersection(BRB.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.BRB$pred.10 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.10[is.na(pred.df.BRB$pred.10)] <- 0
#PPA
pred.10 <- gIntersection(PPA.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.PPA$pred.10 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.10[is.na(pred.df.PPA$pred.10)] <- 0
########################################################################
PPA.05 <- allPts.df[seq(1, NROW(allPts.df), by = 4),]
T.05 <- xyt.lxy(xy = PPA.05[,c("x","y")], dt= as.POSIXct(strptime(PPA.05$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.05$id)
PPA.05 <- as.ltraj(xy = PPA.05[,c("x","y")], date = as.POSIXct(strptime(PPA.05$Date,"%Y-%m-%d"))
                   , id = PPA.05$id, proj4string = CRS("+init=epsg:32631"))
tdge.05 <- tgkde(PPA.05)
tdge.05 <- volras(tdge.05,95)
projection(tdge.05) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.05, Tmax=180*60, Lmin=60)
BRB.05 <- BRB(PPA.05, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.05 <- getverticeshr(BRB.05,percent = 95)
projection(BRB.05) = CRS("+init=epsg:32631")
PPA.05 <- dynppa(PPA.05)
projection(PPA.05) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.05)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.05 <- lxy.nn.add(T.05, s=s.value, k=25)
T.05 <- lxy.lhs(T.05, k=3*2:8, s=s.value)
T.05 <- lhs.iso.add(T.05)
T.05 <- isopleths(T.05)
T.05.shp <- T.05[[1]][T.05[[1]]$iso.level == 0.95,]
#######################################################################
pred.05 <- gIntersection(T.05.shp, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.Tlocoh$pred.05 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.05[is.na(pred.df.Tlocoh$pred.05)] <- 0
#Tdge
pred.05 <- gIntersection(tdge.05, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.tdge$pred.05 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.05[is.na(pred.df.tdge$pred.05)] <- 0
#BRB
pred.05 <- gIntersection(BRB.05, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.BRB$pred.05 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.05[is.na(pred.df.BRB$pred.05)] <- 0
#PPA
pred.05 <- gIntersection(PPA.100, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.PPA$pred.05 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.05[is.na(pred.df.PPA$pred.05)] <- 0
########################################################################
#PPA Area
area.100 <- gArea(gDifference(PPA.100,cp))
area.50 <- gArea(PPA.50)
area.25 <- gArea(PPA.25)
area.10 <- gArea(PPA.10)
area.05 <- gArea(PPA.05)
#Edge Complexity 
PPA.df.edge[i,1] <- (lineLength(PPA.100, byid = FALSE)/area.100) * 1000
PPA.df.edge[i,3] <- (lineLength(PPA.50, byid = FALSE)/area.50) * 1000
PPA.df.edge[i,5] <- (lineLength(PPA.25, byid = FALSE)/area.25) * 1000
PPA.df.edge[i,7] <- (lineLength(PPA.10, byid = FALSE)/area.10) * 1000
PPA.df.edge[i,9] <- (lineLength(PPA.05, byid = FALSE)/area.05) * 1000
PPA.df.edge[i,2] <- length(PPA.100)
PPA.df.edge[i,4] <- length(PPA.50)
PPA.df.edge[i,6] <- length(PPA.25)
PPA.df.edge[i,8] <- length(PPA.10)
PPA.df.edge[i,10] <- length(PPA.05)
#Assign Area percents to datafrmae
PPA.df[i,1] <- area.50/area.100
PPA.df[i,2] <- area.25/area.100
PPA.df[i,3] <- area.10/area.100
PPA.df[i,4] <- area.05/area.100
#TDGE Area
area.100 <- gArea(gDifference(tdge.100,cp))
area.50 <- gArea(tdge.50)
area.25 <- gArea(tdge.25)
area.10 <- gArea(tdge.10)
area.05 <- gArea(tdge.05)
#Assign Area percents to datafrmae
tdge.df[i,1] <- area.50/area.100
tdge.df[i,2] <- area.25/area.100
tdge.df[i,3] <- area.10/area.100
tdge.df[i,4] <- area.05/area.100
#Edge Complexity 
tdge.df.edge[i,1] <- (lineLength(tdge.100, byid = FALSE)/area.100) * 1000
tdge.df.edge[i,3] <- (lineLength(tdge.50, byid = FALSE)/area.50) * 1000
tdge.df.edge[i,5] <- (lineLength(tdge.25, byid = FALSE)/area.25) * 1000
tdge.df.edge[i,7] <- (lineLength(tdge.10, byid = FALSE)/area.10) * 1000
tdge.df.edge[i,9] <- (lineLength(tdge.05, byid = FALSE)/area.05) * 1000
tdge.df.edge[i,2] <- length(tdge.100)
tdge.df.edge[i,4] <- length(tdge.50)
tdge.df.edge[i,6] <- length(tdge.25)
tdge.df.edge[i,8] <- length(tdge.10)
tdge.df.edge[i,10] <- length(tdge.05)
#BRB Area
area.100 <- gArea(BRB.100)
area.50 <- gArea(BRB.50)
area.25 <- gArea(BRB.25)
area.10 <- gArea(BRB.10)
area.05 <- gArea(BRB.05)
#Assign Area percents to datafrmae
BRB.df[i,1] <- area.50/area.100
BRB.df[i,2] <- area.25/area.100
BRB.df[i,3] <- area.10/area.100
BRB.df[i,4] <- area.05/area.100
#Edge Complexity 
BRB.df.edge[i,1] <- (lineLength(BRB.100, byid = FALSE)/area.100) * 1000
BRB.df.edge[i,3] <- (lineLength(BRB.50, byid = FALSE)/area.50) * 1000
BRB.df.edge[i,5] <- (lineLength(BRB.25, byid = FALSE)/area.25) * 1000
BRB.df.edge[i,7] <- (lineLength(BRB.10, byid = FALSE)/area.10) * 1000
BRB.df.edge[i,9] <- (lineLength(BRB.05, byid = FALSE)/area.05) * 1000
BRB.df.edge[i,2] <- length(BRB.100)
BRB.df.edge[i,4] <- length(BRB.50)
BRB.df.edge[i,6] <- length(BRB.25)
BRB.df.edge[i,8] <- length(BRB.10)
BRB.df.edge[i,10] <- length(BRB.05)
#T-locoh area
area.100 <- T.100[[1]]@data[5,2]
area.50 <- T.50[[1]]@data[5,2]
area.25 <- T.25[[1]]@data[5,2]
area.10 <- T.10[[1]]@data[5,2]
area.05 <- T.05[[1]]@data[5,2]
#Assign Area percents to datafrmae
TLoCoHo.df[i,1] <- area.50/area.100
TLoCoHo.df[i,2] <- area.25/area.100
TLoCoHo.df[i,3] <- area.10/area.100
TLoCoHo.df[i,4] <- area.05/area.100
#Edge Complexity 
TLoCoHo.df.edge[i,1] <- (lineLength(T.100.shp, byid = FALSE)/area.100) * 1000
TLoCoHo.df.edge[i,3] <- (lineLength(T.50.shp, byid = FALSE)/area.50) * 1000
TLoCoHo.df.edge[i,5] <- (lineLength(T.25.shp, byid = FALSE)/area.25) * 1000
TLoCoHo.df.edge[i,7] <- (lineLength(T.10.shp, byid = FALSE)/area.10) * 1000
TLoCoHo.df.edge[i,9] <- (lineLength(T.05.shp, byid = FALSE)/area.05) * 1000
TLoCoHo.df.edge[i,2] <- length(T.100.shp)
TLoCoHo.df.edge[i,4] <- length(T.50.shp)
TLoCoHo.df.edge[i,6] <- length(T.25.shp)
TLoCoHo.df.edge[i,8] <- length(T.10.shp)
TLoCoHo.df.edge[i,10] <- length(T.05.shp)
#####AUC#############################################
pred.df.Tlocoh <-pred.df.Tlocoh[c(1,3,2,4,5,6,7,8)]
pred.df.tdge <-pred.df.tdge[c(1,3,2,4,5,6,7,8)]
pred.df.BRB <-pred.df.BRB[c(1,3,2,4,5,6,7,8)]
pred.df.PPA <-pred.df.PPA[c(1,3,2,4,5,6,7,8)]
tAUC <- colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
KDE.df.auc[i,1] <- tAUC[1,1]
KDE.df.auc[i,2] <- tAUC[1,2]
KDE.df.auc[i,3] <- tAUC[1,3]
KDE.df.auc[i,4] <- tAUC[1,4]
KDE.df.auc[i,5] <- tAUC[1,5]
KDE.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,MCP.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
mcp.df.auc[i,1] <- tAUC[1,1]
mcp.df.auc[i,2] <- tAUC[1,2]
mcp.df.auc[i,3] <- tAUC[1,3]
mcp.df.auc[i,4] <- tAUC[1,4]
mcp.df.auc[i,5] <- tAUC[1,5]
mcp.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.Tlocoh[,3:8], pred.df.Tlocoh[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.Tlocoh[,3:8], pred.df.Tlocoh[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
TLoCoHo.df.auc[i,1] <- tAUC[1,1]
TLoCoHo.df.auc[i,2] <- tAUC[1,2]
TLoCoHo.df.auc[i,3] <- tAUC[1,3]
TLoCoHo.df.auc[i,4] <- tAUC[1,4]
TLoCoHo.df.auc[i,5] <- tAUC[1,5]
TLoCoHo.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.tdge[,3:8], pred.df.tdge[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.tdge[,3:8], pred.df.tdge[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
tdge.df.auc[i,1] <- tAUC[1,1]
tdge.df.auc[i,2] <- tAUC[1,2]
tdge.df.auc[i,3] <- tAUC[1,3]
tdge.df.auc[i,4] <- tAUC[1,4]
tdge.df.auc[i,5] <- tAUC[1,5]
tdge.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.BRB[,3:8], pred.df.BRB[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.BRB[,3:8], pred.df.BRB[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
BRB.df.auc[i,1] <- tAUC[1,1]
BRB.df.auc[i,2] <- tAUC[1,2]
BRB.df.auc[i,3] <- tAUC[1,3]
BRB.df.auc[i,4] <- tAUC[1,4]
BRB.df.auc[i,5] <- tAUC[1,5]
BRB.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.PPA[,3:8], pred.df.PPA[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.PPA[,3:8], pred.df.PPA[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
PPA.df.auc[i,1] <- tAUC[1,1]
PPA.df.auc[i,2] <- tAUC[1,2]
PPA.df.auc[i,3] <- tAUC[1,3]
PPA.df.auc[i,4] <- tAUC[1,4]
PPA.df.auc[i,5] <- tAUC[1,5]
PPA.df.auc[i,6] <- tAUC[1,6]

##################################################################
#EMD
r.real <- raster(ncol=200, nrow=200, res=10)
extent(r.real) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.real[] <- 0
tab <- table(cellFromXY(r.real, pts2))
r.real[as.numeric(names(tab))] <- tab
r.real <- raster::as.matrix(r.real)
r.real <- pp(r.real)

#T-locoh
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.100.shp, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.50.shp, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.25.shp, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.10.shp, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.05.shp, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

TLoCoHo.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
TLoCoHo.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
TLoCoHo.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
TLoCoHo.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
TLoCoHo.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#tdge
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

tdge.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
tdge.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
tdge.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
tdge.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
tdge.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#BRB
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

BRB.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
BRB.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
BRB.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
BRB.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
BRB.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#PPA
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

PPA.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
PPA.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
PPA.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
PPA.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
PPA.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#MCP
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

mcp.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
mcp.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
mcp.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
mcp.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
mcp.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#KDE
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

KDE.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
KDE.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
KDE.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
KDE.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
KDE.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

}
write.xlsx(mcp.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\mcp_area.xlsx') 
write.xlsx(mcp.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\mcp_edge.xlsx') 
write.xlsx(mcp.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\mcp_auc.xlsx') 
write.xlsx(mcp.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\mcp_emd.xlsx') 
write.xlsx(KDE.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\KDE_area.xlsx') 
write.xlsx(KDE.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\KDE_edge.xlsx') 
write.xlsx(KDE.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\KDE_auc.xlsx') 
write.xlsx(KDE.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\KDE_emd.xlsx') 
write.xlsx(PPA.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\PPA_area.xlsx') 
write.xlsx(PPA.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\PPA_edge.xlsx') 
write.xlsx(PPA.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\PPA_auc.xlsx') 
write.xlsx(PPA.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\PPA_emd.xlsx') 
write.xlsx(tdge.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\tdge_area.xlsx') 
write.xlsx(tdge.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\tdge_edge.xlsx') 
write.xlsx(tdge.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\tdge_auc.xlsx') 
write.xlsx(tdge.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\tdge_emd.xlsx')
write.xlsx(TLoCoHo.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\TLoCoHo_area.xlsx') 
write.xlsx(TLoCoHo.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\TLoCoHo_edge.xlsx') 
write.xlsx(TLoCoHo.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\TLoCoHo_auc.xlsx') 
write.xlsx(TLoCoHo.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\TLoCoHo_emd.xlsx') 
write.xlsx(BRB.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\BRB_area.xlsx') 
write.xlsx(BRB.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\BRB_edge.xlsx') 
write.xlsx(BRB.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\BRB_auc.xlsx') 
write.xlsx(BRB.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\disjoint\\BRB_emd.xlsx') 


################################################################################
###Convex HomeRanges
nclust <- function(x0, y0, radius, n) {
  return(runifdisc(n, radius, centre=c(x0, y0)))
}

#Create Dataframe to save area calculations for MCP
mcp.df <- data.frame(matrix(ncol = 4, nrow = 100))
mcp.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("MCP_Area50%", "MCP_Area 25%", "MCP_Area 10%", "MCP_Area 05%")
colnames(mcp.df) <- x
x <- c("Obs", "MCP_AUC 100%","MCP_AUC 50%", "MCP_AUC 25%", "MCP_AUC 10%", "MCP_AUC 05%")
colnames(mcp.df.auc) <- x
mcp.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("MCP_EMD 100%","MCP_EMD 50%", "MCP_EMD 25%", "MCP_EMD 10%", "MCP_EMD 05%")
colnames(mcp.df.emd) <- x
mcp.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("MCP_edge 100%", "100% Patches","MCP_edge 50%", "50% Patches","MCP_edge 25%", "25% Patches", "MCP_edge 10%", "10% Patches", "MCP_edge 05%", "05% Patches")
colnames(mcp.df.edge) <- x

#Create Dataframe to save area calculations for KDE
KDE.df <- data.frame(matrix(ncol = 4, nrow = 100))
KDE.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("KDE_Area50%", "KDE_Area 25%", "KDE_Area 10%", "KDE_Area 05%")
colnames(KDE.df) <- x
x <- c("Obs", "KDE_AUC 100%","KDE_AUC 50%", "KDE_AUC 25%", "KDE_AUC 10%", "KDE_AUC 05%")
colnames(KDE.df.auc) <- x
KDE.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("KDE_EMD 100%","KDE_EMD 50%", "KDE_EMD 25%", "KDE_EMD 10%", "KDE_EMD 05%")
colnames(KDE.df.emd) <- x
KDE.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("KDE_edge 100%", "100% Patches","KDE_edge 50%", "50% Patches","KDE_edge 25%", "25% Patches", "KDE_edge 10%", "10% Patches", "KDE_edge 05%", "05% Patches")
colnames(KDE.df.edge) <- x

#Create Dataframe to save area calculations for PPA
PPA.df <- data.frame(matrix(ncol = 4, nrow = 100))
PPA.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("PPA_Area50%", "PPA_Area 25%", "PPA_Area 10%", "PPA_Area 05%")
colnames(PPA.df) <- x
x <- c("Obs", "PPA_AUC 100%","PPA_AUC 50%", "PPA_AUC 25%", "PPA_AUC 10%", "PPA_AUC 05%")
colnames(PPA.df.auc) <- x
PPA.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("PPA_EMD 100%","PPA_EMD 50%", "PPA_EMD 25%", "PPA_EMD 10%", "PPA_EMD 05%")
colnames(PPA.df.emd) <- x
PPA.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("PPA_edge 100%", "100% Patches","PPA_edge 50%", "50% Patches","PPA_edge 25%", "25% Patches", "PPA_edge 10%", "10% Patches", "PPA_edge 05%", "05% Patches")
colnames(PPA.df.edge) <- x

#Create Dataframe to save area calculations for TDGE
tdge.df <- data.frame(matrix(ncol = 4, nrow = 100))
tdge.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("TDGE_Area50%", "TDGE_Area 25%", "TDGE_Area 10%", "TDGE_Area 05%")
colnames(tdge.df) <- x
x <- c("Obs", "TDGE_AUC 100%","TDGE_AUC 50%", "TDGE_AUC 25%", "TDGE_AUC 10%", "TDGE_AUC 05%")
colnames(tdge.df.auc) <- x
tdge.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("tdge_EMD 100%","tdge_EMD 50%", "tdge_EMD 25%", "tdge_EMD 10%", "tdge_EMD 05%")
colnames(tdge.df.emd) <- x
tdge.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("tdge_edge 100%", "100% Patches","tdge_edge 50%", "50% Patches","tdge_edge 25%", "25% Patches", "tdge_edge 10%", "10% Patches", "tdge_edge 05%", "05% Patches")
colnames(tdge.df.edge) <- x

#Create Dataframe to save area calculations for BRB
BRB.df <- data.frame(matrix(ncol = 4, nrow = 100))
BRB.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("BRB_Area50%", "BRB_Area 25%", "BRB_Area 10%", "BRB_Area 05%")
colnames(BRB.df) <- x
x <- c("Obs", "BRB_AUC 100%","BRB_AUC 50%", "BRB_AUC 25%", "BRB_AUC 10%", "BRB_AUC 05%")
colnames(BRB.df.auc) <- x
BRB.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("BRB_EMD 100%","BRB_EMD 50%", "BRB_EMD 25%", "BRB_EMD 10%", "BRB_EMD 05%")
colnames(BRB.df.emd) <- x
BRB.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("BRB_edge 100%", "100% Patches","BRB_edge 50%", "50% Patches","BRB_edge 25%", "25% Patches", "BRB_edge 10%", "10% Patches", "BRB_edge 05%", "05% Patches")
colnames(BRB.df.edge) <- x

#Create Dataframe to save area calculations for T-LoCoHO
TLoCoHo.df <- data.frame(matrix(ncol = 4, nrow = 100))
TLoCoHo.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("T-LoCoHO_Area50%", "T-LoCoHO_Area 25%", "T-LoCoHO_Area 10%", "T-LoCoHO_Area 05%")
colnames(TLoCoHo.df) <- x
x <- c("Obs", "T-LoCoH_AUC 100%","T-LoCoH_AUC 50%", "T-LoCoH_AUC 25%", "T-LoCoH_AUC 10%", "T-LoCoH_AUC 05%")
colnames(TLoCoHo.df.auc) <- x
TLoCoHo.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("TLoCoHo_EMD 100%","TLoCoHo_EMD 50%", "TLoCoHo_EMD 25%", "TLoCoHo_EMD 10%", "TLoCoHo_EMD 05%")
colnames(TLoCoHo.df.emd) <- x
TLoCoHo.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("TLoCoHo_edge 100%", "100% Patches","TLoCoHo_edge 50%", "50% Patches","TLoCoHo_edge 25%", "25% Patches", "TLoCoHo_edge 10%", "10% Patches", "TLoCoHo_edge 05%", "05% Patches")
colnames(TLoCoHo.df.edge) <- x

i = 1
for(i in 1:4){
pc = rPoissonCluster(100, 0.05, nclust, radius=0.25, n=80)
pc.df = as.data.frame(pc)
sp = SpatialPoints(pc.df )

#Create sample for MCP
s <- sample(sp, 12)
#MCP
cp <- mcp(s, percent=95)
#Clip points by MCP
sp <- sp[cp, ]

#WGS 84
projection(sp) = CRS("+init=epsg:4326")
#Plots points
sp<- spTransform(sp, CRS("+init=epsg:32631"))
#Create sample for MCP
s <- sample(sp, 10)
#MCP
cp <- mcp(s, percent = 95)
#Clip points by MCP
allPts<- gDifference(sp, cp)

################################################################
#Set Up df for Traj
ptNumber = length(allPts)
time = seq(as.Date("2000/1/1"), by = "day", length.out = ptNumber)
allPts$Date = time
allPts$id = "Brock"
allPts.df = as.data.frame(allPts)
#Set up data for AUC
#Create Hex grid for present/absence counts
ptsBuffer = gBuffer(sp, width = 4000)
hexPts <-spsample(ptsBuffer,type="hexagonal",cellsize=5000)
hexPols <- HexPoints2SpatialPolygons(hexPts)
#Define true for AUC
pts2 <- allPts.df[, -c(1:2)] # delete columns 5 through 7
pts2 <- SpatialPoints(pts2)
projection(pts2) = CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pts2, hexPols)))
#MCP
pred.df.MCP <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.MCP ) <- "GridNumber"
pred.df.MCP$GridNumber <- 1:length(hexPols)
pred.df.MCP$Obs <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$Obs[is.na(pred.df.MCP$Obs)] <- 0
pred.df.MCP$Presence <- pred.df.MCP$Obs > 0
#PPA
pred.df.PPA <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.PPA ) <- "GridNumber"
pred.df.PPA$GridNumber <- 1:length(hexPols)
pred.df.PPA$Obs <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$Obs[is.na(pred.df.PPA$Obs)] <- 0
pred.df.PPA$Presence <- pred.df.PPA$Obs > 0
#tdge
pred.df.tdge <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.tdge ) <- "GridNumber"
pred.df.tdge$GridNumber <- 1:length(hexPols)
pred.df.tdge$Obs <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$Obs[is.na(pred.df.tdge$Obs)] <- 0
pred.df.tdge$Presence <- pred.df.tdge$Obs > 0
#BRB
pred.df.BRB <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.BRB ) <- "GridNumber"
pred.df.BRB$GridNumber <- 1:length(hexPols)
pred.df.BRB$Obs <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$Obs[is.na(pred.df.BRB$Obs)] <- 0
pred.df.BRB$Presence <- pred.df.BRB$Obs > 0
#T-Locoh
pred.df.Tlocoh <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.Tlocoh ) <- "GridNumber"
pred.df.Tlocoh$GridNumber <- 1:length(hexPols)
pred.df.Tlocoh$Obs <- as.numeric(res$Freq[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$Obs[is.na(pred.df.Tlocoh$Obs)] <- 0
pred.df.Tlocoh$Presence <- pred.df.Tlocoh$Obs > 0
###################################################################

#########################End Point Pattern Set UP ##############################
# Sample the Points at differ levels
#Sample 50 %
sampleSize = length(allPts) * .50
pts.50 = spsample(allPts,n = sampleSize, "random")
#Sample 25 %
sampleSize = length(allPts) * .25
pts.25 = spsample(allPts,n = sampleSize, "random")
#Sample 10 %
sampleSize = length(allPts) * .10
pts.10 = spsample(allPts,n = sampleSize, "random")
#Sample 5 %
sampleSize = length(allPts) * .05
pts.05 = spsample(allPts,n = sampleSize, "random")
#################################Run Home Ranges###############################
#MCP Home Ranges
mcp.100 = mcp(allPts, percent = 95)
#MCP Core Area
# mcp.100.perforated.50 = mcp(allPts, percent = 50)
#MCP Home Range 50%
mcp.50 = mcp(pts.50, percent = 95)
#MCP Core Area 50%
# mcp.50.perforated.50 = mcp(pts.50, percent = 50)
#MCP Home Range 50%
mcp.25= mcp(pts.25, percent = 95)
#MCP Core Area 50%
# mcp.25.perforated.50 = mcp(pts.25, percent = 50)
#MCP Home Range 50%
mcp.10 = mcp(pts.10, percent = 95)
#MCP Core Area 50%
# mcp.10.perforated.50 = mcp(pts.10, percent = 50)
#MCP Home Range 50%
mcp.05 = mcp(pts.05, percent = 95)
#MCP Core Area 50%
# mcp.05.perforated.50 = mcp(pts.05, percent = 50)
#Perferate the MCP to get correct area
area.100 <- gArea(gDifference(mcp.100,cp))
area.50 <- gArea(mcp.50)
area.25 <- gArea(mcp.25)
area.10 <- gArea(mcp.10)
area.05 <- gArea(mcp.05)
mcp.df[i,1] <- area.50/area.100
mcp.df[i,2] <- area.25/area.100
mcp.df[i,3] <- area.10/area.100
mcp.df[i,4] <- area.05/area.100
mcp.df.edge[i,1] <- (lineLength(mcp.100, byid = FALSE)/area.100) * 1000
mcp.df.edge[i,3] <- (lineLength(mcp.50, byid = FALSE)/area.50) * 1000
mcp.df.edge[i,5] <- (lineLength(mcp.25, byid = FALSE)/area.25) * 1000
mcp.df.edge[i,7] <- (lineLength(mcp.10, byid = FALSE)/area.10) * 1000
mcp.df.edge[i,9] <- (lineLength(mcp.05, byid = FALSE)/area.05) * 1000
mcp.df.edge[i,2] <- length(mcp.100)
mcp.df.edge[i,4] <- length(mcp.50)
mcp.df.edge[i,6] <- length(mcp.25)
mcp.df.edge[i,8] <- length(mcp.10)
mcp.df.edge[i,10] <- length(mcp.05)
###################################################################
#AUC MCP
pred.100.mcp <- gIntersection(mcp.100, sp)
projection(pred.100.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100.mcp, hexPols)))
pred.df.MCP$pred.100 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.100[is.na(pred.df.MCP$pred.100)] <- 0
#reorder dataframe
pred.df.MCP <- pred.df.MCP[c(1,3,2,4)]

pred.50.mcp <- gIntersection(mcp.50, sp)
projection(pred.50.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50.mcp, hexPols)))
pred.df.MCP$pred.50 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.50[is.na(pred.df.MCP$pred.50)] <- 0

pred.25.mcp <- gIntersection(mcp.25, sp)
projection(pred.25.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25.mcp, hexPols)))
pred.df.MCP$pred.25 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.25[is.na(pred.df.MCP$pred.25)] <- 0

pred.10.mcp <- gIntersection(mcp.10, sp)
projection(pred.10.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10.mcp, hexPols)))
pred.df.MCP$pred.10 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.10[is.na(pred.df.MCP$pred.10)] <- 0

pred.05.mcp <- gIntersection(mcp.05, sp)
projection(pred.05.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05.mcp, hexPols)))
pred.df.MCP$pred.05 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.05[is.na(pred.df.MCP$pred.05)] <- 0

tAUC <- colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
mcp.df.auc[i,1] <- tAUC[1,1]
mcp.df.auc[i,2] <- tAUC[1,2]
mcp.df.auc[i,3] <- tAUC[1,3]
mcp.df.auc[i,4] <- tAUC[1,4]
mcp.df.auc[i,5] <- tAUC[1,5]
mcp.df.auc[i,6] <- tAUC[1,6]

##################################################################
#KDE Homerange
KDE.100 <- kernelUD(allPts, h="href")
KDE.100 <- getverticeshr(KDE.100,percent = 95)
KDE.50 <- kernelUD(pts.50, h="href")
KDE.50 <- getverticeshr(KDE.50,percent = 95)
KDE.25 <- kernelUD(pts.25, h="href")
KDE.25 <- getverticeshr(KDE.25,percent = 95)
KDE.10 <- kernelUD(pts.10, h="href")
KDE.10 <- getverticeshr(KDE.10,percent = 95)
KDE.05 <- kernelUD(pts.05, h="href")
KDE.05 <- getverticeshr(KDE.05,percent = 95)
#KDE Area
area.100 <- KDE.100$area
area.50 <- KDE.50$area
area.25 <- KDE.25$area
area.10 <- KDE.10$area
area.05 <- KDE.05$area
#Assign Area percents to datafrmae
KDE.df[i,1] <- area.50/area.100
KDE.df[i,2] <- area.25/area.100
KDE.df[i,3] <- area.10/area.100
KDE.df[i,4] <- area.05/area.100
KDE.df.edge[i,1] <- (lineLength(KDE.100, byid = FALSE)/area.100) * 1000
KDE.df.edge[i,3] <- (lineLength(KDE.50, byid = FALSE)/area.50) * 1000
KDE.df.edge[i,5] <- (lineLength(KDE.25, byid = FALSE)/area.25) * 1000
KDE.df.edge[i,7] <- (lineLength(KDE.10, byid = FALSE)/area.10) * 1000
KDE.df.edge[i,9] <- (lineLength(KDE.05, byid = FALSE)/area.05) * 1000
KDE.df.edge[i,2] <- length(KDE.100)
KDE.df.edge[i,4] <- length(KDE.50)
KDE.df.edge[i,6] <- length(KDE.25)
KDE.df.edge[i,8] <- length(KDE.10)
KDE.df.edge[i,10] <- length(KDE.05)
###################################################################
#AUC KDE
#Define true AUC
pts2 <- allPts.df[, -c(1:2)] # delete columns 5 through 7
pts2 <- SpatialPoints(pts2)
projection(pts2) = CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pts2, hexPols)))
#KDE
pred.df.KDE <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.KDE) <- "GridNumber"
pred.df.KDE$GridNumber <- 1:length(hexPols)
pred.df.KDE$Obs <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$Obs[is.na(pred.df.KDE$Obs)] <- 0
pred.df.KDE$Presence <- pred.df.KDE$Obs > 0
pred.df.KDE <- pred.df.KDE[c(1,3,2)]

#100 %
pred.100.KDE <- gIntersection(KDE.100, sp)
projection(pred.100.KDE) <- CRS("+init=epsg:32631")
plot(pred.100.KDE, add = TRUE)
res <- as.data.frame(table(over(pred.100.KDE, hexPols)))
pred.df.KDE$pred.100 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.100[is.na(pred.df.KDE$pred.100)] <- 0


#50%
pred.50.KDE <- gIntersection(KDE.50, sp)
projection(pred.50.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50.KDE, hexPols)))
pred.df.KDE$pred.50 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.50[is.na(pred.df.KDE$pred.50)] <- 0


#25%
pred.25.KDE <- gIntersection(KDE.25, sp)
projection(pred.25.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25.KDE, hexPols)))
pred.df.KDE$pred.25 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.25[is.na(pred.df.KDE$pred.25)] <- 0


#10%
pred.10.KDE <- gIntersection(KDE.10, sp)
projection(pred.10.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10.KDE, hexPols)))
pred.df.KDE$pred.10 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.10[is.na(pred.df.KDE$pred.10)] <- 0

#05%
pred.05.KDE <- gIntersection(KDE.05, sp)
projection(pred.05.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05.KDE, hexPols)))
pred.df.KDE$pred.05 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.05[is.na(pred.df.KDE$pred.05)] <- 0

#auc
tAUC <- colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(sp,KDE.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
KDE.df.auc[i,1] <- tAUC[1,1]
KDE.df.auc[i,2] <- tAUC[1,2]
KDE.df.auc[i,3] <- tAUC[1,3]
KDE.df.auc[i,4] <- tAUC[1,4]
KDE.df.auc[i,5] <- tAUC[1,5]
KDE.df.auc[i,6] <- tAUC[1,6]

######################################################################################################
#PPA, TDGE, BRB, AND t-LOCOCH  Homerange
#Converting to ltraj for time dynamics
PPA.100 <- as.ltraj(xy = allPts.df[,c("x","y")], date = as.POSIXct(strptime(allPts.df$Date,"%Y-%m-%d"))
                    , id = allPts.df$id, proj4string = CRS("+init=epsg:32631"))
tdge.100 <- tgkde(PPA.100)
tdge.100 <- volras(tdge.100,95)
projection(tdge.100) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.100, Tmax=180*60, Lmin=60)
#Hmin is super important here. 
BRB.100 <- BRB(PPA.100, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.100 <- getverticeshr(BRB.100,percent = 95)
projection(BRB.100 ) = CRS("+init=epsg:32631")
PPA.100 <- dynppa(PPA.100)
projection(PPA.100) = CRS("+init=epsg:32631")
#T-locoh
T.100 <- xyt.lxy(xy = allPts.df[,c("x","y")], dt= as.POSIXct(strptime(allPts.df$Date,"%Y-%m-%d"))
                 , proj4string= CRS("+init=epsg:32631"), id=allPts$id)
toni.lxy <- lxy.ptsh.add(T.100)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.100 <- lxy.nn.add(T.100, s=s.value, k=25)
T.100 <- lxy.lhs(T.100, k=3*2:8, s=s.value)
T.100 <- lhs.iso.add(T.100)
T.100 <- isopleths(T.100)
T.100.shp <- T.100[[1]][T.100[[1]]$iso.level == 0.95,]
#####################################################################
#AUC
#T-locoh
pred.100 <- gIntersection(T.100.shp, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.Tlocoh$pred.100 <- as.numeric(res$Freq[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.100[is.na(pred.df.Tlocoh$pred.100)] <- 0
#Tdge
pred.100 <- gIntersection(tdge.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.tdge$pred.100 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.100[is.na(pred.df.tdge$pred.100)] <- 0
#BRB
pred.100 <- gIntersection(BRB.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.BRB$pred.100 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.100[is.na(pred.df.BRB$pred.100)] <- 0
#PPA
pred.100 <- gIntersection(PPA.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.PPA$pred.100 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.100[is.na(pred.df.PPA$pred.100)] <- 0


#####################################################################
#Even Sample
PPA.50 <- allPts.df[seq(1, NROW(allPts.df), by = 2),]
#Need tp run xyt before converting to PPA
T.50 <- xyt.lxy(xy = PPA.50[,c("x","y")], dt= as.POSIXct(strptime(PPA.50$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.50$id)
PPA.50 <- as.ltraj(xy = PPA.50[,c("x","y")], date = as.POSIXct(strptime(PPA.50$Date,"%Y-%m-%d"))
                   , id = PPA.50$id, proj4string = CRS("+init=epsg:32631"))
tdge.50 <- tgkde(PPA.50)
tdge.50 <- volras(tdge.50,95)
projection(tdge.50) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.50, Tmax=180*60, Lmin=60)
BRB.50 <- BRB(PPA.50, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.50 <- getverticeshr(BRB.50,percent = 95)
projection(BRB.50) = CRS("+init=epsg:32631")
PPA.50 <- dynppa(PPA.50)
projection(PPA.50 ) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.50)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.50 <- lxy.nn.add(T.50, s=s.value, k=25)
T.50 <- lxy.lhs(T.50, k=3*2:8, s=s.value)
T.50 <- lhs.iso.add(T.50)
T.50 <- isopleths(T.50)
#Get isopleth
T.50.shp <- T.50[[1]][T.50[[1]]$iso.level == 0.95,]
#######################################################################
#AUC
#T-locoh
pred.50 <- gIntersection(T.50.shp, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.Tlocoh$pred.50 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.50[is.na(pred.df.Tlocoh$pred.50)] <- 0
#Tdge
pred.50 <- gIntersection(tdge.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.tdge$pred.50 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.50[is.na(pred.df.tdge$pred.50)] <- 0
#BRB
pred.50 <- gIntersection(BRB.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.BRB$pred.50 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.50[is.na(pred.df.BRB$pred.50)] <- 0
#PPA
pred.50 <- gIntersection(PPA.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.PPA$pred.50 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.50[is.na(pred.df.PPA$pred.50)] <- 0
########################################################################
PPA.25 <- allPts.df[seq(1, NROW(allPts.df), by = 4),]
T.25 <- xyt.lxy(xy = PPA.25[,c("x","y")], dt= as.POSIXct(strptime(PPA.25$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.25$id)
PPA.25 <- as.ltraj(xy = PPA.25[,c("x","y")], date = as.POSIXct(strptime(PPA.25$Date,"%Y-%m-%d"))
                   , id = PPA.25$id, proj4string = CRS("+init=epsg:32631"))
tdge.25 <- tgkde(PPA.25)
tdge.25 <- volras(tdge.25,95)
projection(tdge.25) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.25, Tmax=180*60, Lmin=60)
BRB.25 <- BRB(PPA.25, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.25 <- getverticeshr(BRB.25,percent = 95)
projection(BRB.25) = CRS("+init=epsg:32631")
PPA.25 <- dynppa(PPA.25)
projection(PPA.25) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.25)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.25 <- lxy.nn.add(T.25, s=s.value, k=25)
T.25 <- lxy.lhs(T.25, k=3*2:8, s=s.value)
T.25 <- lhs.iso.add(T.25)
T.25 <- isopleths(T.25)
T.25.shp <- T.25[[1]][T.25[[1]]$iso.level == 0.95,]
#######################################################################
pred.25 <- gIntersection(T.25.shp, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.Tlocoh$pred.25 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.25[is.na(pred.df.Tlocoh$pred.25)] <- 0
#Tdge
pred.25 <- gIntersection(tdge.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.tdge$pred.25 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.25[is.na(pred.df.tdge$pred.25)] <- 0
#BRB
pred.25 <- gIntersection(BRB.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.BRB$pred.25 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.25[is.na(pred.df.BRB$pred.25)] <- 0
#PPA
pred.25 <- gIntersection(PPA.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.PPA$pred.25 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.25[is.na(pred.df.PPA$pred.25)] <- 0
########################################################################
PPA.10 <- allPts.df[seq(1, NROW(allPts.df), by = 10),]
T.10 <- xyt.lxy(xy = PPA.10[,c("x","y")], dt= as.POSIXct(strptime(PPA.10$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.10$id)
PPA.10<- as.ltraj(xy = PPA.10[,c("x","y")], date = as.POSIXct(strptime(PPA.10$Date,"%Y-%m-%d"))
                  , id = PPA.10$id, proj4string = CRS("+init=epsg:32631"))
tdge.10 <- tgkde(PPA.10)
tdge.10 <- volras(tdge.10,95)
projection(tdge.10 ) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.10, Tmax=180*60, Lmin=60)
BRB.10 <- BRB(PPA.10, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.10 <- getverticeshr(BRB.10,percent = 95)
projection(BRB.10) = CRS("+init=epsg:32631")
PPA.10 <- dynppa(PPA.10)
projection(PPA.10 ) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.10)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.10 <- lxy.nn.add(T.10, s=s.value, k=25)
T.10 <- lxy.lhs(T.10, k=3*2:8, s=s.value)
T.10 <- lhs.iso.add(T.10)
T.10 <- isopleths(T.10)
T.10.shp <- T.10[[1]][T.10[[1]]$iso.level == 0.95,]
#######################################################################
pred.10 <- gIntersection(T.10.shp, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.Tlocoh$pred.10 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.10[is.na(pred.df.Tlocoh$pred.10)] <- 0
#Tdge
pred.10 <- gIntersection(tdge.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.tdge$pred.10 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.10[is.na(pred.df.tdge$pred.10)] <- 0
#BRB
pred.10 <- gIntersection(BRB.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.BRB$pred.10 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.10[is.na(pred.df.BRB$pred.10)] <- 0
#PPA
pred.10 <- gIntersection(PPA.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.PPA$pred.10 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.10[is.na(pred.df.PPA$pred.10)] <- 0
########################################################################
PPA.05 <- allPts.df[seq(1, NROW(allPts.df), by = 4),]
T.05 <- xyt.lxy(xy = PPA.05[,c("x","y")], dt= as.POSIXct(strptime(PPA.05$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.05$id)
PPA.05 <- as.ltraj(xy = PPA.05[,c("x","y")], date = as.POSIXct(strptime(PPA.05$Date,"%Y-%m-%d"))
                   , id = PPA.05$id, proj4string = CRS("+init=epsg:32631"))
tdge.05 <- tgkde(PPA.05)
tdge.05 <- volras(tdge.05,95)
projection(tdge.05) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.05, Tmax=180*60, Lmin=60)
BRB.05 <- BRB(PPA.05, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.05 <- getverticeshr(BRB.05,percent = 95)
projection(BRB.05) = CRS("+init=epsg:32631")
PPA.05 <- dynppa(PPA.05)
projection(PPA.05) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.05)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.05 <- lxy.nn.add(T.05, s=s.value, k=25)
T.05 <- lxy.lhs(T.05, k=3*2:8, s=s.value)
T.05 <- lhs.iso.add(T.05)
T.05 <- isopleths(T.05)
T.05.shp <- T.05[[1]][T.05[[1]]$iso.level == 0.95,]
#######################################################################
pred.05 <- gIntersection(T.05.shp, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.Tlocoh$pred.05 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.05[is.na(pred.df.Tlocoh$pred.05)] <- 0
#Tdge
pred.05 <- gIntersection(tdge.05, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.tdge$pred.05 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.05[is.na(pred.df.tdge$pred.05)] <- 0
#BRB
pred.05 <- gIntersection(BRB.05, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.BRB$pred.05 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.05[is.na(pred.df.BRB$pred.05)] <- 0
#PPA
pred.05 <- gIntersection(PPA.100, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.PPA$pred.05 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.05[is.na(pred.df.PPA$pred.05)] <- 0
########################################################################
#PPA Area
area.100 <- gArea(gDifference(PPA.100,cp))
area.50 <- gArea(PPA.50)
area.25 <- gArea(PPA.25)
area.10 <- gArea(PPA.10)
area.05 <- gArea(PPA.05)
#Edge Complexity 
PPA.df.edge[i,1] <- (lineLength(PPA.100, byid = FALSE)/area.100) * 1000
PPA.df.edge[i,3] <- (lineLength(PPA.50, byid = FALSE)/area.50) * 1000
PPA.df.edge[i,5] <- (lineLength(PPA.25, byid = FALSE)/area.25) * 1000
PPA.df.edge[i,7] <- (lineLength(PPA.10, byid = FALSE)/area.10) * 1000
PPA.df.edge[i,9] <- (lineLength(PPA.05, byid = FALSE)/area.05) * 1000
PPA.df.edge[i,2] <- length(PPA.100)
PPA.df.edge[i,4] <- length(PPA.50)
PPA.df.edge[i,6] <- length(PPA.25)
PPA.df.edge[i,8] <- length(PPA.10)
PPA.df.edge[i,10] <- length(PPA.05)
#Assign Area percents to datafrmae
PPA.df[i,1] <- area.50/area.100
PPA.df[i,2] <- area.25/area.100
PPA.df[i,3] <- area.10/area.100
PPA.df[i,4] <- area.05/area.100
#TDGE Area
area.100 <- gArea(gDifference(tdge.100,cp))
area.50 <- gArea(tdge.50)
area.25 <- gArea(tdge.25)
area.10 <- gArea(tdge.10)
area.05 <- gArea(tdge.05)
#Assign Area percents to datafrmae
tdge.df[i,1] <- area.50/area.100
tdge.df[i,2] <- area.25/area.100
tdge.df[i,3] <- area.10/area.100
tdge.df[i,4] <- area.05/area.100
#Edge Complexity 
tdge.df.edge[i,1] <- (lineLength(tdge.100, byid = FALSE)/area.100) * 1000
tdge.df.edge[i,3] <- (lineLength(tdge.50, byid = FALSE)/area.50) * 1000
tdge.df.edge[i,5] <- (lineLength(tdge.25, byid = FALSE)/area.25) * 1000
tdge.df.edge[i,7] <- (lineLength(tdge.10, byid = FALSE)/area.10) * 1000
tdge.df.edge[i,9] <- (lineLength(tdge.05, byid = FALSE)/area.05) * 1000
tdge.df.edge[i,2] <- length(tdge.100)
tdge.df.edge[i,4] <- length(tdge.50)
tdge.df.edge[i,6] <- length(tdge.25)
tdge.df.edge[i,8] <- length(tdge.10)
tdge.df.edge[i,10] <- length(tdge.05)
#BRB Area
area.100 <- gArea(BRB.100)
area.50 <- gArea(BRB.50)
area.25 <- gArea(BRB.25)
area.10 <- gArea(BRB.10)
area.05 <- gArea(BRB.05)
#Assign Area percents to datafrmae
BRB.df[i,1] <- area.50/area.100
BRB.df[i,2] <- area.25/area.100
BRB.df[i,3] <- area.10/area.100
BRB.df[i,4] <- area.05/area.100
#Edge Complexity 
BRB.df.edge[i,1] <- (lineLength(BRB.100, byid = FALSE)/area.100) * 1000
BRB.df.edge[i,3] <- (lineLength(BRB.50, byid = FALSE)/area.50) * 1000
BRB.df.edge[i,5] <- (lineLength(BRB.25, byid = FALSE)/area.25) * 1000
BRB.df.edge[i,7] <- (lineLength(BRB.10, byid = FALSE)/area.10) * 1000
BRB.df.edge[i,9] <- (lineLength(BRB.05, byid = FALSE)/area.05) * 1000
BRB.df.edge[i,2] <- length(BRB.100)
BRB.df.edge[i,4] <- length(BRB.50)
BRB.df.edge[i,6] <- length(BRB.25)
BRB.df.edge[i,8] <- length(BRB.10)
BRB.df.edge[i,10] <- length(BRB.05)
#T-locoh area
area.100 <- T.100[[1]]@data[5,2]
area.50 <- T.50[[1]]@data[5,2]
area.25 <- T.25[[1]]@data[5,2]
area.10 <- T.10[[1]]@data[5,2]
area.05 <- T.05[[1]]@data[5,2]
#Assign Area percents to datafrmae
TLoCoHo.df[i,1] <- area.50/area.100
TLoCoHo.df[i,2] <- area.25/area.100
TLoCoHo.df[i,3] <- area.10/area.100
TLoCoHo.df[i,4] <- area.05/area.100
#Edge Complexity 
TLoCoHo.df.edge[i,1] <- (lineLength(T.100.shp, byid = FALSE)/area.100) * 1000
TLoCoHo.df.edge[i,3] <- (lineLength(T.50.shp, byid = FALSE)/area.50) * 1000
TLoCoHo.df.edge[i,5] <- (lineLength(T.25.shp, byid = FALSE)/area.25) * 1000
TLoCoHo.df.edge[i,7] <- (lineLength(T.10.shp, byid = FALSE)/area.10) * 1000
TLoCoHo.df.edge[i,9] <- (lineLength(T.05.shp, byid = FALSE)/area.05) * 1000
TLoCoHo.df.edge[i,2] <- length(T.100.shp)
TLoCoHo.df.edge[i,4] <- length(T.50.shp)
TLoCoHo.df.edge[i,6] <- length(T.25.shp)
TLoCoHo.df.edge[i,8] <- length(T.10.shp)
TLoCoHo.df.edge[i,10] <- length(T.05.shp)
#####AUC#############################################
pred.df.Tlocoh <-pred.df.Tlocoh[c(1,3,2,4,5,6,7,8)]
pred.df.tdge <-pred.df.tdge[c(1,3,2,4,5,6,7,8)]
pred.df.BRB <-pred.df.BRB[c(1,3,2,4,5,6,7,8)]
pred.df.PPA <-pred.df.PPA[c(1,3,2,4,5,6,7,8)]
tAUC <- colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
KDE.df.auc[i,1] <- tAUC[1,1]
KDE.df.auc[i,2] <- tAUC[1,2]
KDE.df.auc[i,3] <- tAUC[1,3]
KDE.df.auc[i,4] <- tAUC[1,4]
KDE.df.auc[i,5] <- tAUC[1,5]
KDE.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,MCP.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
mcp.df.auc[i,1] <- tAUC[1,1]
mcp.df.auc[i,2] <- tAUC[1,2]
mcp.df.auc[i,3] <- tAUC[1,3]
mcp.df.auc[i,4] <- tAUC[1,4]
mcp.df.auc[i,5] <- tAUC[1,5]
mcp.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.Tlocoh[,3:8], pred.df.Tlocoh[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.Tlocoh[,3:8], pred.df.Tlocoh[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
TLoCoHo.df.auc[i,1] <- tAUC[1,1]
TLoCoHo.df.auc[i,2] <- tAUC[1,2]
TLoCoHo.df.auc[i,3] <- tAUC[1,3]
TLoCoHo.df.auc[i,4] <- tAUC[1,4]
TLoCoHo.df.auc[i,5] <- tAUC[1,5]
TLoCoHo.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.tdge[,3:8], pred.df.tdge[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.tdge[,3:8], pred.df.tdge[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
tdge.df.auc[i,1] <- tAUC[1,1]
tdge.df.auc[i,2] <- tAUC[1,2]
tdge.df.auc[i,3] <- tAUC[1,3]
tdge.df.auc[i,4] <- tAUC[1,4]
tdge.df.auc[i,5] <- tAUC[1,5]
tdge.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.BRB[,3:8], pred.df.BRB[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.BRB[,3:8], pred.df.BRB[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
BRB.df.auc[i,1] <- tAUC[1,1]
BRB.df.auc[i,2] <- tAUC[1,2]
BRB.df.auc[i,3] <- tAUC[1,3]
BRB.df.auc[i,4] <- tAUC[1,4]
BRB.df.auc[i,5] <- tAUC[1,5]
BRB.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.PPA[,3:8], pred.df.PPA[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.PPA[,3:8], pred.df.PPA[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
PPA.df.auc[i,1] <- tAUC[1,1]
PPA.df.auc[i,2] <- tAUC[1,2]
PPA.df.auc[i,3] <- tAUC[1,3]
PPA.df.auc[i,4] <- tAUC[1,4]
PPA.df.auc[i,5] <- tAUC[1,5]
PPA.df.auc[i,6] <- tAUC[1,6]

##################################################################
#EMD
r.real <- raster(ncol=200, nrow=200, res=10)
extent(r.real) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.real[] <- 0
tab <- table(cellFromXY(r.real, pts2))
r.real[as.numeric(names(tab))] <- tab
r.real <- raster::as.matrix(r.real)
r.real <- pp(r.real)

#T-locoh
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.100.shp, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.50.shp, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.25.shp, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.10.shp, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.05.shp, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

TLoCoHo.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
TLoCoHo.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
TLoCoHo.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
TLoCoHo.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
TLoCoHo.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#tdge
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

tdge.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
tdge.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
tdge.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
tdge.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
tdge.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#BRB
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

BRB.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
BRB.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
BRB.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
BRB.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
BRB.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#PPA
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

PPA.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
PPA.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
PPA.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
PPA.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
PPA.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#MCP
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

mcp.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
mcp.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
mcp.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
mcp.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
mcp.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#KDE
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

KDE.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
KDE.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
KDE.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
KDE.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
KDE.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

}
write.xlsx(mcp.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\mcp_area.xlsx') 
write.xlsx(mcp.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\mcp_edge.xlsx') 
write.xlsx(mcp.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\mcp_auc.xlsx') 
write.xlsx(mcp.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\mcp_emd.xlsx') 
write.xlsx(KDE.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\KDE_area.xlsx') 
write.xlsx(KDE.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\KDE_edge.xlsx') 
write.xlsx(KDE.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\KDE_auc.xlsx') 
write.xlsx(KDE.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\KDE_emd.xlsx') 
write.xlsx(PPA.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\PPA_area.xlsx') 
write.xlsx(PPA.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\PPA_edge.xlsx') 
write.xlsx(PPA.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\PPA_auc.xlsx') 
write.xlsx(PPA.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\PPA_emd.xlsx') 
write.xlsx(tdge.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\tdge_area.xlsx') 
write.xlsx(tdge.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\tdge_edge.xlsx') 
write.xlsx(tdge.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\tdge_auc.xlsx') 
write.xlsx(tdge.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\tdge_emd.xlsx')
write.xlsx(TLoCoHo.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\TLoCoHo_area.xlsx') 
write.xlsx(TLoCoHo.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\TLoCoHo_edge.xlsx') 
write.xlsx(TLoCoHo.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\TLoCoHo_auc.xlsx') 
write.xlsx(TLoCoHo.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\TLoCoHo_emd.xlsx') 
write.xlsx(BRB.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\BRB_area.xlsx') 
write.xlsx(BRB.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\BRB_edge.xlsx') 
write.xlsx(BRB.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\BRB_auc.xlsx') 
write.xlsx(BRB.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\convex\\BRB_emd.xlsx') 

################################################################################
###Linear HomeRanges
nclust <- function(x0, y0, radius, n) {
  return(runifdisc(n, radius, centre=c(x0, y0)))
}

#Create Dataframe to save area calculations for MCP
mcp.df <- data.frame(matrix(ncol = 4, nrow = 100))
mcp.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("MCP_Area50%", "MCP_Area 25%", "MCP_Area 10%", "MCP_Area 05%")
colnames(mcp.df) <- x
x <- c("Obs", "MCP_AUC 100%","MCP_AUC 50%", "MCP_AUC 25%", "MCP_AUC 10%", "MCP_AUC 05%")
colnames(mcp.df.auc) <- x
mcp.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("MCP_EMD 100%","MCP_EMD 50%", "MCP_EMD 25%", "MCP_EMD 10%", "MCP_EMD 05%")
colnames(mcp.df.emd) <- x
mcp.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("MCP_edge 100%", "100% Patches","MCP_edge 50%", "50% Patches","MCP_edge 25%", "25% Patches", "MCP_edge 10%", "10% Patches", "MCP_edge 05%", "05% Patches")
colnames(mcp.df.edge) <- x

#Create Dataframe to save area calculations for KDE
KDE.df <- data.frame(matrix(ncol = 4, nrow = 100))
KDE.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("KDE_Area50%", "KDE_Area 25%", "KDE_Area 10%", "KDE_Area 05%")
colnames(KDE.df) <- x
x <- c("Obs", "KDE_AUC 100%","KDE_AUC 50%", "KDE_AUC 25%", "KDE_AUC 10%", "KDE_AUC 05%")
colnames(KDE.df.auc) <- x
KDE.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("KDE_EMD 100%","KDE_EMD 50%", "KDE_EMD 25%", "KDE_EMD 10%", "KDE_EMD 05%")
colnames(KDE.df.emd) <- x
KDE.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("KDE_edge 100%", "100% Patches","KDE_edge 50%", "50% Patches","KDE_edge 25%", "25% Patches", "KDE_edge 10%", "10% Patches", "KDE_edge 05%", "05% Patches")
colnames(KDE.df.edge) <- x

#Create Dataframe to save area calculations for PPA
PPA.df <- data.frame(matrix(ncol = 4, nrow = 100))
PPA.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("PPA_Area50%", "PPA_Area 25%", "PPA_Area 10%", "PPA_Area 05%")
colnames(PPA.df) <- x
x <- c("Obs", "PPA_AUC 100%","PPA_AUC 50%", "PPA_AUC 25%", "PPA_AUC 10%", "PPA_AUC 05%")
colnames(PPA.df.auc) <- x
PPA.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("PPA_EMD 100%","PPA_EMD 50%", "PPA_EMD 25%", "PPA_EMD 10%", "PPA_EMD 05%")
colnames(PPA.df.emd) <- x
PPA.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("PPA_edge 100%", "100% Patches","PPA_edge 50%", "50% Patches","PPA_edge 25%", "25% Patches", "PPA_edge 10%", "10% Patches", "PPA_edge 05%", "05% Patches")
colnames(PPA.df.edge) <- x

#Create Dataframe to save area calculations for TDGE
tdge.df <- data.frame(matrix(ncol = 4, nrow = 100))
tdge.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("TDGE_Area50%", "TDGE_Area 25%", "TDGE_Area 10%", "TDGE_Area 05%")
colnames(tdge.df) <- x
x <- c("Obs", "TDGE_AUC 100%","TDGE_AUC 50%", "TDGE_AUC 25%", "TDGE_AUC 10%", "TDGE_AUC 05%")
colnames(tdge.df.auc) <- x
tdge.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("tdge_EMD 100%","tdge_EMD 50%", "tdge_EMD 25%", "tdge_EMD 10%", "tdge_EMD 05%")
colnames(tdge.df.emd) <- x
tdge.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("tdge_edge 100%", "100% Patches","tdge_edge 50%", "50% Patches","tdge_edge 25%", "25% Patches", "tdge_edge 10%", "10% Patches", "tdge_edge 05%", "05% Patches")
colnames(tdge.df.edge) <- x

#Create Dataframe to save area calculations for BRB
BRB.df <- data.frame(matrix(ncol = 4, nrow = 100))
BRB.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("BRB_Area50%", "BRB_Area 25%", "BRB_Area 10%", "BRB_Area 05%")
colnames(BRB.df) <- x
x <- c("Obs", "BRB_AUC 100%","BRB_AUC 50%", "BRB_AUC 25%", "BRB_AUC 10%", "BRB_AUC 05%")
colnames(BRB.df.auc) <- x
BRB.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("BRB_EMD 100%","BRB_EMD 50%", "BRB_EMD 25%", "BRB_EMD 10%", "BRB_EMD 05%")
colnames(BRB.df.emd) <- x
BRB.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("BRB_edge 100%", "100% Patches","BRB_edge 50%", "50% Patches","BRB_edge 25%", "25% Patches", "BRB_edge 10%", "10% Patches", "BRB_edge 05%", "05% Patches")
colnames(BRB.df.edge) <- x

#Create Dataframe to save area calculations for T-LoCoHO
TLoCoHo.df <- data.frame(matrix(ncol = 4, nrow = 100))
TLoCoHo.df.auc <- data.frame(matrix(ncol = 6, nrow = 100))
x <- c("T-LoCoHO_Area50%", "T-LoCoHO_Area 25%", "T-LoCoHO_Area 10%", "T-LoCoHO_Area 05%")
colnames(TLoCoHo.df) <- x
x <- c("Obs", "T-LoCoH_AUC 100%","T-LoCoH_AUC 50%", "T-LoCoH_AUC 25%", "T-LoCoH_AUC 10%", "T-LoCoH_AUC 05%")
colnames(TLoCoHo.df.auc) <- x
TLoCoHo.df.emd <- data.frame(matrix(ncol = 5, nrow = 100))
x <- c("TLoCoHo_EMD 100%","TLoCoHo_EMD 50%", "TLoCoHo_EMD 25%", "TLoCoHo_EMD 10%", "TLoCoHo_EMD 05%")
colnames(TLoCoHo.df.emd) <- x
TLoCoHo.df.edge <- data.frame(matrix(ncol = 10, nrow = 100))
x <- c("TLoCoHo_edge 100%", "100% Patches","TLoCoHo_edge 50%", "50% Patches","TLoCoHo_edge 25%", "25% Patches", "TLoCoHo_edge 10%", "10% Patches", "TLoCoHo_edge 05%", "05% Patches")
colnames(TLoCoHo.df.edge) <- x

i = 1

for(i in 1:4){
pc = rPoissonCluster(100, 0.05, nclust, radius=0.25, n=80)
pc.df = as.data.frame(pc)
sp = SpatialPoints(pc.df)
# Sample lines
s <- sample(sp, 5)
L = SpatialLines(list(Lines(list(Line(coordinates(s))),"X"))) 
buffer = (gBuffer(L,width=0.1,  joinStyle="MITRE", capStyle="SQUARE"))
sp <- sp[buffer, ]
#WGS 84
projection(sp) = CRS("+init=epsg:4326")
#Plots points
sp<- spTransform(sp, CRS("+init=epsg:32631"))
#Create sample for MCP
s <- sample(sp, 10)
#MCP
cp <- mcp(s, percent = 95)
#Clip points by MCP
allPts<- gDifference(sp, cp)

################################################################
#Set Up df for Traj
ptNumber = length(allPts)
time = seq(as.Date("2000/1/1"), by = "day", length.out = ptNumber)
allPts$Date = time
allPts$id = "Brock"
allPts.df = as.data.frame(allPts)
#Set up data for AUC
#Create Hex grid for present/absence counts
ptsBuffer = gBuffer(sp, width = 4000)
hexPts <-spsample(ptsBuffer,type="hexagonal",cellsize=5000)
hexPols <- HexPoints2SpatialPolygons(hexPts)
#Define true for AUC
pts2 <- allPts.df[, -c(1:2)] # delete columns 5 through 7
pts2 <- SpatialPoints(pts2)
projection(pts2) = CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pts2, hexPols)))
#MCP
pred.df.MCP <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.MCP ) <- "GridNumber"
pred.df.MCP$GridNumber <- 1:length(hexPols)
pred.df.MCP$Obs <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$Obs[is.na(pred.df.MCP$Obs)] <- 0
pred.df.MCP$Presence <- pred.df.MCP$Obs > 0
#PPA
pred.df.PPA <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.PPA ) <- "GridNumber"
pred.df.PPA$GridNumber <- 1:length(hexPols)
pred.df.PPA$Obs <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$Obs[is.na(pred.df.PPA$Obs)] <- 0
pred.df.PPA$Presence <- pred.df.PPA$Obs > 0
#tdge
pred.df.tdge <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.tdge ) <- "GridNumber"
pred.df.tdge$GridNumber <- 1:length(hexPols)
pred.df.tdge$Obs <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$Obs[is.na(pred.df.tdge$Obs)] <- 0
pred.df.tdge$Presence <- pred.df.tdge$Obs > 0
#BRB
pred.df.BRB <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.BRB ) <- "GridNumber"
pred.df.BRB$GridNumber <- 1:length(hexPols)
pred.df.BRB$Obs <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$Obs[is.na(pred.df.BRB$Obs)] <- 0
pred.df.BRB$Presence <- pred.df.BRB$Obs > 0
#T-Locoh
pred.df.Tlocoh <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.Tlocoh ) <- "GridNumber"
pred.df.Tlocoh$GridNumber <- 1:length(hexPols)
pred.df.Tlocoh$Obs <- as.numeric(res$Freq[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$Obs[is.na(pred.df.Tlocoh$Obs)] <- 0
pred.df.Tlocoh$Presence <- pred.df.Tlocoh$Obs > 0
###################################################################

#########################End Point Pattern Set UP ##############################
# Sample the Points at differ levels
#Sample 50 %
sampleSize = length(allPts) * .50
pts.50 = spsample(allPts,n = sampleSize, "random")
#Sample 25 %
sampleSize = length(allPts) * .25
pts.25 = spsample(allPts,n = sampleSize, "random")
#Sample 10 %
sampleSize = length(allPts) * .10
pts.10 = spsample(allPts,n = sampleSize, "random")
#Sample 5 %
sampleSize = length(allPts) * .05
pts.05 = spsample(allPts,n = sampleSize, "random")
#################################Run Home Ranges###############################
#MCP Home Ranges
mcp.100 = mcp(allPts, percent = 95)
#MCP Core Area
# mcp.100.perforated.50 = mcp(allPts, percent = 50)
#MCP Home Range 50%
mcp.50 = mcp(pts.50, percent = 95)
#MCP Core Area 50%
# mcp.50.perforated.50 = mcp(pts.50, percent = 50)
#MCP Home Range 50%
mcp.25= mcp(pts.25, percent = 95)
#MCP Core Area 50%
# mcp.25.perforated.50 = mcp(pts.25, percent = 50)
#MCP Home Range 50%
mcp.10 = mcp(pts.10, percent = 95)
#MCP Core Area 50%
# mcp.10.perforated.50 = mcp(pts.10, percent = 50)
#MCP Home Range 50%
mcp.05 = mcp(pts.05, percent = 95)
#MCP Core Area 50%
# mcp.05.perforated.50 = mcp(pts.05, percent = 50)
#Perferate the MCP to get correct area
area.100 <- gArea(gDifference(mcp.100,cp))
area.50 <- gArea(mcp.50)
area.25 <- gArea(mcp.25)
area.10 <- gArea(mcp.10)
area.05 <- gArea(mcp.05)
mcp.df[i,1] <- area.50/area.100
mcp.df[i,2] <- area.25/area.100
mcp.df[i,3] <- area.10/area.100
mcp.df[i,4] <- area.05/area.100
mcp.df.edge[i,1] <- (lineLength(mcp.100, byid = FALSE)/area.100) * 1000
mcp.df.edge[i,3] <- (lineLength(mcp.50, byid = FALSE)/area.50) * 1000
mcp.df.edge[i,5] <- (lineLength(mcp.25, byid = FALSE)/area.25) * 1000
mcp.df.edge[i,7] <- (lineLength(mcp.10, byid = FALSE)/area.10) * 1000
mcp.df.edge[i,9] <- (lineLength(mcp.05, byid = FALSE)/area.05) * 1000
mcp.df.edge[i,2] <- length(mcp.100)
mcp.df.edge[i,4] <- length(mcp.50)
mcp.df.edge[i,6] <- length(mcp.25)
mcp.df.edge[i,8] <- length(mcp.10)
mcp.df.edge[i,10] <- length(mcp.05)
###################################################################
#AUC MCP
pred.100.mcp <- gIntersection(mcp.100, sp)
projection(pred.100.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100.mcp, hexPols)))
pred.df.MCP$pred.100 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.100[is.na(pred.df.MCP$pred.100)] <- 0
#reorder dataframe
pred.df.MCP <- pred.df.MCP[c(1,3,2,4)]

pred.50.mcp <- gIntersection(mcp.50, sp)
projection(pred.50.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50.mcp, hexPols)))
pred.df.MCP$pred.50 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.50[is.na(pred.df.MCP$pred.50)] <- 0

pred.25.mcp <- gIntersection(mcp.25, sp)
projection(pred.25.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25.mcp, hexPols)))
pred.df.MCP$pred.25 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.25[is.na(pred.df.MCP$pred.25)] <- 0

pred.10.mcp <- gIntersection(mcp.10, sp)
projection(pred.10.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10.mcp, hexPols)))
pred.df.MCP$pred.10 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.10[is.na(pred.df.MCP$pred.10)] <- 0

pred.05.mcp <- gIntersection(mcp.05, sp)
projection(pred.05.mcp) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05.mcp, hexPols)))
pred.df.MCP$pred.05 <- as.numeric(res$Freq[match(pred.df.MCP$GridNumber, res$Var1)])
pred.df.MCP$pred.05[is.na(pred.df.MCP$pred.05)] <- 0

tAUC <- colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
mcp.df.auc[i,1] <- tAUC[1,1]
mcp.df.auc[i,2] <- tAUC[1,2]
mcp.df.auc[i,3] <- tAUC[1,3]
mcp.df.auc[i,4] <- tAUC[1,4]
mcp.df.auc[i,5] <- tAUC[1,5]
mcp.df.auc[i,6] <- tAUC[1,6]

##################################################################
#KDE Homerange
KDE.100 <- kernelUD(allPts, h="href")
KDE.100 <- getverticeshr(KDE.100,percent = 95)
KDE.50 <- kernelUD(pts.50, h="href")
KDE.50 <- getverticeshr(KDE.50,percent = 95)
KDE.25 <- kernelUD(pts.25, h="href")
KDE.25 <- getverticeshr(KDE.25,percent = 95)
KDE.10 <- kernelUD(pts.10, h="href")
KDE.10 <- getverticeshr(KDE.10,percent = 95)
KDE.05 <- kernelUD(pts.05, h="href")
KDE.05 <- getverticeshr(KDE.05,percent = 95)
#KDE Area
area.100 <- KDE.100$area
area.50 <- KDE.50$area
area.25 <- KDE.25$area
area.10 <- KDE.10$area
area.05 <- KDE.05$area
#Assign Area percents to datafrmae
KDE.df[i,1] <- area.50/area.100
KDE.df[i,2] <- area.25/area.100
KDE.df[i,3] <- area.10/area.100
KDE.df[i,4] <- area.05/area.100
KDE.df.edge[i,1] <- (lineLength(KDE.100, byid = FALSE)/area.100) * 1000
KDE.df.edge[i,3] <- (lineLength(KDE.50, byid = FALSE)/area.50) * 1000
KDE.df.edge[i,5] <- (lineLength(KDE.25, byid = FALSE)/area.25) * 1000
KDE.df.edge[i,7] <- (lineLength(KDE.10, byid = FALSE)/area.10) * 1000
KDE.df.edge[i,9] <- (lineLength(KDE.05, byid = FALSE)/area.05) * 1000
KDE.df.edge[i,2] <- length(KDE.100)
KDE.df.edge[i,4] <- length(KDE.50)
KDE.df.edge[i,6] <- length(KDE.25)
KDE.df.edge[i,8] <- length(KDE.10)
KDE.df.edge[i,10] <- length(KDE.05)
###################################################################
#AUC KDE
#Define true AUC
pts2 <- allPts.df[, -c(1:2)] # delete columns 5 through 7
pts2 <- SpatialPoints(pts2)
projection(pts2) = CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pts2, hexPols)))
#KDE
pred.df.KDE <- data.frame(matrix(ncol = 1, nrow = length(hexPols)))
names(pred.df.KDE) <- "GridNumber"
pred.df.KDE$GridNumber <- 1:length(hexPols)
pred.df.KDE$Obs <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$Obs[is.na(pred.df.KDE$Obs)] <- 0
pred.df.KDE$Presence <- pred.df.KDE$Obs > 0
pred.df.KDE <- pred.df.KDE[c(1,3,2)]

#100 %
pred.100.KDE <- gIntersection(KDE.100, sp)
projection(pred.100.KDE) <- CRS("+init=epsg:32631")
plot(pred.100.KDE, add = TRUE)
res <- as.data.frame(table(over(pred.100.KDE, hexPols)))
pred.df.KDE$pred.100 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.100[is.na(pred.df.KDE$pred.100)] <- 0


#50%
pred.50.KDE <- gIntersection(KDE.50, sp)
projection(pred.50.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50.KDE, hexPols)))
pred.df.KDE$pred.50 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.50[is.na(pred.df.KDE$pred.50)] <- 0


#25%
pred.25.KDE <- gIntersection(KDE.25, sp)
projection(pred.25.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25.KDE, hexPols)))
pred.df.KDE$pred.25 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.25[is.na(pred.df.KDE$pred.25)] <- 0


#10%
pred.10.KDE <- gIntersection(KDE.10, sp)
projection(pred.10.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10.KDE, hexPols)))
pred.df.KDE$pred.10 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.10[is.na(pred.df.KDE$pred.10)] <- 0

#05%
pred.05.KDE <- gIntersection(KDE.05, sp)
projection(pred.05.KDE) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05.KDE, hexPols)))
pred.df.KDE$pred.05 <- as.numeric(res$Freq[match(pred.df.KDE$GridNumber, res$Var1)])
pred.df.KDE$pred.05[is.na(pred.df.KDE$pred.05)] <- 0

#auc
tAUC <- colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(sp,KDE.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
KDE.df.auc[i,1] <- tAUC[1,1]
KDE.df.auc[i,2] <- tAUC[1,2]
KDE.df.auc[i,3] <- tAUC[1,3]
KDE.df.auc[i,4] <- tAUC[1,4]
KDE.df.auc[i,5] <- tAUC[1,5]
KDE.df.auc[i,6] <- tAUC[1,6]

######################################################################################################
#PPA, TDGE, BRB, AND t-LOCOCH  Homerange
#Converting to ltraj for time dynamics
PPA.100 <- as.ltraj(xy = allPts.df[,c("x","y")], date = as.POSIXct(strptime(allPts.df$Date,"%Y-%m-%d"))
                    , id = allPts.df$id, proj4string = CRS("+init=epsg:32631"))
tdge.100 <- tgkde(PPA.100)
tdge.100 <- volras(tdge.100,95)
projection(tdge.100) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.100, Tmax=180*60, Lmin=60)
#Hmin is super important here. 
BRB.100 <- BRB(PPA.100, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.100 <- getverticeshr(BRB.100,percent = 95)
projection(BRB.100 ) = CRS("+init=epsg:32631")
PPA.100 <- dynppa(PPA.100)
projection(PPA.100) = CRS("+init=epsg:32631")
#T-locoh
T.100 <- xyt.lxy(xy = allPts.df[,c("x","y")], dt= as.POSIXct(strptime(allPts.df$Date,"%Y-%m-%d"))
                 , proj4string= CRS("+init=epsg:32631"), id=allPts$id)
toni.lxy <- lxy.ptsh.add(T.100)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.100 <- lxy.nn.add(T.100, s=s.value, k=25)
T.100 <- lxy.lhs(T.100, k=3*2:8, s=s.value)
T.100 <- lhs.iso.add(T.100)
T.100 <- isopleths(T.100)
T.100.shp <- T.100[[1]][T.100[[1]]$iso.level == 0.95,]
#####################################################################
#AUC
#T-locoh
pred.100 <- gIntersection(T.100.shp, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.Tlocoh$pred.100 <- as.numeric(res$Freq[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.100[is.na(pred.df.Tlocoh$pred.100)] <- 0
#Tdge
pred.100 <- gIntersection(tdge.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.tdge$pred.100 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.100[is.na(pred.df.tdge$pred.100)] <- 0
#BRB
pred.100 <- gIntersection(BRB.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.BRB$pred.100 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.100[is.na(pred.df.BRB$pred.100)] <- 0
#PPA
pred.100 <- gIntersection(PPA.100, sp)
projection(pred.100) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.100, hexPols)))
pred.df.PPA$pred.100 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.100[is.na(pred.df.PPA$pred.100)] <- 0


#####################################################################
#Even Sample
PPA.50 <- allPts.df[seq(1, NROW(allPts.df), by = 2),]
#Need tp run xyt before converting to PPA
T.50 <- xyt.lxy(xy = PPA.50[,c("x","y")], dt= as.POSIXct(strptime(PPA.50$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.50$id)
PPA.50 <- as.ltraj(xy = PPA.50[,c("x","y")], date = as.POSIXct(strptime(PPA.50$Date,"%Y-%m-%d"))
                   , id = PPA.50$id, proj4string = CRS("+init=epsg:32631"))
tdge.50 <- tgkde(PPA.50)
tdge.50 <- volras(tdge.50,95)
projection(tdge.50) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.50, Tmax=180*60, Lmin=60)
BRB.50 <- BRB(PPA.50, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.50 <- getverticeshr(BRB.50,percent = 95)
projection(BRB.50) = CRS("+init=epsg:32631")
PPA.50 <- dynppa(PPA.50)
projection(PPA.50 ) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.50)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.50 <- lxy.nn.add(T.50, s=s.value, k=25)
T.50 <- lxy.lhs(T.50, k=3*2:8, s=s.value)
T.50 <- lhs.iso.add(T.50)
T.50 <- isopleths(T.50)
#Get isopleth
T.50.shp <- T.50[[1]][T.50[[1]]$iso.level == 0.95,]
#######################################################################
#AUC
#T-locoh
pred.50 <- gIntersection(T.50.shp, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.Tlocoh$pred.50 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.50[is.na(pred.df.Tlocoh$pred.50)] <- 0
#Tdge
pred.50 <- gIntersection(tdge.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.tdge$pred.50 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.50[is.na(pred.df.tdge$pred.50)] <- 0
#BRB
pred.50 <- gIntersection(BRB.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.BRB$pred.50 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.50[is.na(pred.df.BRB$pred.50)] <- 0
#PPA
pred.50 <- gIntersection(PPA.50, sp)
projection(pred.50) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.50, hexPols)))
pred.df.PPA$pred.50 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.50[is.na(pred.df.PPA$pred.50)] <- 0
########################################################################
PPA.25 <- allPts.df[seq(1, NROW(allPts.df), by = 4),]
T.25 <- xyt.lxy(xy = PPA.25[,c("x","y")], dt= as.POSIXct(strptime(PPA.25$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.25$id)
PPA.25 <- as.ltraj(xy = PPA.25[,c("x","y")], date = as.POSIXct(strptime(PPA.25$Date,"%Y-%m-%d"))
                   , id = PPA.25$id, proj4string = CRS("+init=epsg:32631"))
tdge.25 <- tgkde(PPA.25)
tdge.25 <- volras(tdge.25,95)
projection(tdge.25) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.25, Tmax=180*60, Lmin=60)
BRB.25 <- BRB(PPA.25, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.25 <- getverticeshr(BRB.25,percent = 95)
projection(BRB.25) = CRS("+init=epsg:32631")
PPA.25 <- dynppa(PPA.25)
projection(PPA.25) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.25)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.25 <- lxy.nn.add(T.25, s=s.value, k=25)
T.25 <- lxy.lhs(T.25, k=3*2:8, s=s.value)
T.25 <- lhs.iso.add(T.25)
T.25 <- isopleths(T.25)
T.25.shp <- T.25[[1]][T.25[[1]]$iso.level == 0.95,]
#######################################################################
pred.25 <- gIntersection(T.25.shp, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.Tlocoh$pred.25 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.25[is.na(pred.df.Tlocoh$pred.25)] <- 0
#Tdge
pred.25 <- gIntersection(tdge.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.tdge$pred.25 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.25[is.na(pred.df.tdge$pred.25)] <- 0
#BRB
pred.25 <- gIntersection(BRB.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.BRB$pred.25 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.25[is.na(pred.df.BRB$pred.25)] <- 0
#PPA
pred.25 <- gIntersection(PPA.25, sp)
projection(pred.25) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.25, hexPols)))
pred.df.PPA$pred.25 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.25[is.na(pred.df.PPA$pred.25)] <- 0
########################################################################
PPA.10 <- allPts.df[seq(1, NROW(allPts.df), by = 10),]
T.10 <- xyt.lxy(xy = PPA.10[,c("x","y")], dt= as.POSIXct(strptime(PPA.10$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.10$id)
PPA.10<- as.ltraj(xy = PPA.10[,c("x","y")], date = as.POSIXct(strptime(PPA.10$Date,"%Y-%m-%d"))
                  , id = PPA.10$id, proj4string = CRS("+init=epsg:32631"))
tdge.10 <- tgkde(PPA.10)
tdge.10 <- volras(tdge.10,95)
projection(tdge.10 ) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.10, Tmax=180*60, Lmin=60)
BRB.10 <- BRB(PPA.10, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.10 <- getverticeshr(BRB.10,percent = 95)
projection(BRB.10) = CRS("+init=epsg:32631")
PPA.10 <- dynppa(PPA.10)
projection(PPA.10 ) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.10)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.10 <- lxy.nn.add(T.10, s=s.value, k=25)
T.10 <- lxy.lhs(T.10, k=3*2:8, s=s.value)
T.10 <- lhs.iso.add(T.10)
T.10 <- isopleths(T.10)
T.10.shp <- T.10[[1]][T.10[[1]]$iso.level == 0.95,]
#######################################################################
pred.10 <- gIntersection(T.10.shp, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.Tlocoh$pred.10 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.10[is.na(pred.df.Tlocoh$pred.10)] <- 0
#Tdge
pred.10 <- gIntersection(tdge.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.tdge$pred.10 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.10[is.na(pred.df.tdge$pred.10)] <- 0
#BRB
pred.10 <- gIntersection(BRB.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.BRB$pred.10 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.10[is.na(pred.df.BRB$pred.10)] <- 0
#PPA
pred.10 <- gIntersection(PPA.10, sp)
projection(pred.10) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.10, hexPols)))
pred.df.PPA$pred.10 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.10[is.na(pred.df.PPA$pred.10)] <- 0
########################################################################
PPA.05 <- allPts.df[seq(1, NROW(allPts.df), by = 4),]
T.05 <- xyt.lxy(xy = PPA.05[,c("x","y")], dt= as.POSIXct(strptime(PPA.05$Date,"%Y-%m-%d"))
                , proj4string= CRS("+init=epsg:32631"), id=PPA.05$id)
PPA.05 <- as.ltraj(xy = PPA.05[,c("x","y")], date = as.POSIXct(strptime(PPA.05$Date,"%Y-%m-%d"))
                   , id = PPA.05$id, proj4string = CRS("+init=epsg:32631"))
tdge.05 <- tgkde(PPA.05)
tdge.05 <- volras(tdge.05,95)
projection(tdge.05) = CRS("+init=epsg:32631")
vv <- BRB.D(PPA.05, Tmax=180*60, Lmin=60)
BRB.05 <- BRB(PPA.05, D = vv, Tmax = 180*60, Lmin = 60*48, hmin=5000)
BRB.05 <- getverticeshr(BRB.05,percent = 95)
projection(BRB.05) = CRS("+init=epsg:32631")
PPA.05 <- dynppa(PPA.05)
projection(PPA.05) = CRS("+init=epsg:32631")
#T-locoh
toni.lxy <- lxy.ptsh.add(T.05)
# S value of 60 percent of the data "half way between 40-80%"
s.value = toni.lxy$ptsh
s.value = s.value$`Brock`[[1]]$s.ptsh[20,1]
T.05 <- lxy.nn.add(T.05, s=s.value, k=25)
T.05 <- lxy.lhs(T.05, k=3*2:8, s=s.value)
T.05 <- lhs.iso.add(T.05)
T.05 <- isopleths(T.05)
T.05.shp <- T.05[[1]][T.05[[1]]$iso.level == 0.95,]
#######################################################################
pred.05 <- gIntersection(T.05.shp, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.Tlocoh$pred.05 <- as.numeric(res$Var1[match(pred.df.Tlocoh$GridNumber, res$Var1)])
pred.df.Tlocoh$pred.05[is.na(pred.df.Tlocoh$pred.05)] <- 0
#Tdge
pred.05 <- gIntersection(tdge.05, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.tdge$pred.05 <- as.numeric(res$Freq[match(pred.df.tdge$GridNumber, res$Var1)])
pred.df.tdge$pred.05[is.na(pred.df.tdge$pred.05)] <- 0
#BRB
pred.05 <- gIntersection(BRB.05, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.BRB$pred.05 <- as.numeric(res$Freq[match(pred.df.BRB$GridNumber, res$Var1)])
pred.df.BRB$pred.05[is.na(pred.df.BRB$pred.05)] <- 0
#PPA
pred.05 <- gIntersection(PPA.100, sp)
projection(pred.05) <- CRS("+init=epsg:32631")
res <- as.data.frame(table(over(pred.05, hexPols)))
pred.df.PPA$pred.05 <- as.numeric(res$Freq[match(pred.df.PPA$GridNumber, res$Var1)])
pred.df.PPA$pred.05[is.na(pred.df.PPA$pred.05)] <- 0
########################################################################
#PPA Area
area.100 <- gArea(gDifference(PPA.100,cp))
area.50 <- gArea(PPA.50)
area.25 <- gArea(PPA.25)
area.10 <- gArea(PPA.10)
area.05 <- gArea(PPA.05)
#Edge Complexity 
PPA.df.edge[i,1] <- (lineLength(PPA.100, byid = FALSE)/area.100) * 1000
PPA.df.edge[i,3] <- (lineLength(PPA.50, byid = FALSE)/area.50) * 1000
PPA.df.edge[i,5] <- (lineLength(PPA.25, byid = FALSE)/area.25) * 1000
PPA.df.edge[i,7] <- (lineLength(PPA.10, byid = FALSE)/area.10) * 1000
PPA.df.edge[i,9] <- (lineLength(PPA.05, byid = FALSE)/area.05) * 1000
PPA.df.edge[i,2] <- length(PPA.100)
PPA.df.edge[i,4] <- length(PPA.50)
PPA.df.edge[i,6] <- length(PPA.25)
PPA.df.edge[i,8] <- length(PPA.10)
PPA.df.edge[i,10] <- length(PPA.05)
#Assign Area percents to datafrmae
PPA.df[i,1] <- area.50/area.100
PPA.df[i,2] <- area.25/area.100
PPA.df[i,3] <- area.10/area.100
PPA.df[i,4] <- area.05/area.100
#TDGE Area
area.100 <- gArea(gDifference(tdge.100,cp))
area.50 <- gArea(tdge.50)
area.25 <- gArea(tdge.25)
area.10 <- gArea(tdge.10)
area.05 <- gArea(tdge.05)
#Assign Area percents to datafrmae
tdge.df[i,1] <- area.50/area.100
tdge.df[i,2] <- area.25/area.100
tdge.df[i,3] <- area.10/area.100
tdge.df[i,4] <- area.05/area.100
#Edge Complexity 
tdge.df.edge[i,1] <- (lineLength(tdge.100, byid = FALSE)/area.100) * 1000
tdge.df.edge[i,3] <- (lineLength(tdge.50, byid = FALSE)/area.50) * 1000
tdge.df.edge[i,5] <- (lineLength(tdge.25, byid = FALSE)/area.25) * 1000
tdge.df.edge[i,7] <- (lineLength(tdge.10, byid = FALSE)/area.10) * 1000
tdge.df.edge[i,9] <- (lineLength(tdge.05, byid = FALSE)/area.05) * 1000
tdge.df.edge[i,2] <- length(tdge.100)
tdge.df.edge[i,4] <- length(tdge.50)
tdge.df.edge[i,6] <- length(tdge.25)
tdge.df.edge[i,8] <- length(tdge.10)
tdge.df.edge[i,10] <- length(tdge.05)
#BRB Area
area.100 <- gArea(BRB.100)
area.50 <- gArea(BRB.50)
area.25 <- gArea(BRB.25)
area.10 <- gArea(BRB.10)
area.05 <- gArea(BRB.05)
#Assign Area percents to datafrmae
BRB.df[i,1] <- area.50/area.100
BRB.df[i,2] <- area.25/area.100
BRB.df[i,3] <- area.10/area.100
BRB.df[i,4] <- area.05/area.100
#Edge Complexity 
BRB.df.edge[i,1] <- (lineLength(BRB.100, byid = FALSE)/area.100) * 1000
BRB.df.edge[i,3] <- (lineLength(BRB.50, byid = FALSE)/area.50) * 1000
BRB.df.edge[i,5] <- (lineLength(BRB.25, byid = FALSE)/area.25) * 1000
BRB.df.edge[i,7] <- (lineLength(BRB.10, byid = FALSE)/area.10) * 1000
BRB.df.edge[i,9] <- (lineLength(BRB.05, byid = FALSE)/area.05) * 1000
BRB.df.edge[i,2] <- length(BRB.100)
BRB.df.edge[i,4] <- length(BRB.50)
BRB.df.edge[i,6] <- length(BRB.25)
BRB.df.edge[i,8] <- length(BRB.10)
BRB.df.edge[i,10] <- length(BRB.05)
#T-locoh area
area.100 <- T.100[[1]]@data[5,2]
area.50 <- T.50[[1]]@data[5,2]
area.25 <- T.25[[1]]@data[5,2]
area.10 <- T.10[[1]]@data[5,2]
area.05 <- T.05[[1]]@data[5,2]
#Assign Area percents to datafrmae
TLoCoHo.df[i,1] <- area.50/area.100
TLoCoHo.df[i,2] <- area.25/area.100
TLoCoHo.df[i,3] <- area.10/area.100
TLoCoHo.df[i,4] <- area.05/area.100
#Edge Complexity 
TLoCoHo.df.edge[i,1] <- (lineLength(T.100.shp, byid = FALSE)/area.100) * 1000
TLoCoHo.df.edge[i,3] <- (lineLength(T.50.shp, byid = FALSE)/area.50) * 1000
TLoCoHo.df.edge[i,5] <- (lineLength(T.25.shp, byid = FALSE)/area.25) * 1000
TLoCoHo.df.edge[i,7] <- (lineLength(T.10.shp, byid = FALSE)/area.10) * 1000
TLoCoHo.df.edge[i,9] <- (lineLength(T.05.shp, byid = FALSE)/area.05) * 1000
TLoCoHo.df.edge[i,2] <- length(T.100.shp)
TLoCoHo.df.edge[i,4] <- length(T.50.shp)
TLoCoHo.df.edge[i,6] <- length(T.25.shp)
TLoCoHo.df.edge[i,8] <- length(T.10.shp)
TLoCoHo.df.edge[i,10] <- length(T.05.shp)
#####AUC#############################################
pred.df.Tlocoh <-pred.df.Tlocoh[c(1,3,2,4,5,6,7,8)]
pred.df.tdge <-pred.df.tdge[c(1,3,2,4,5,6,7,8)]
pred.df.BRB <-pred.df.BRB[c(1,3,2,4,5,6,7,8)]
pred.df.PPA <-pred.df.PPA[c(1,3,2,4,5,6,7,8)]
tAUC <- colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.KDE[,3:8], pred.df.KDE[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
KDE.df.auc[i,1] <- tAUC[1,1]
KDE.df.auc[i,2] <- tAUC[1,2]
KDE.df.auc[i,3] <- tAUC[1,3]
KDE.df.auc[i,4] <- tAUC[1,4]
KDE.df.auc[i,5] <- tAUC[1,5]
KDE.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.MCP[,3:8], pred.df.MCP[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,MCP.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
mcp.df.auc[i,1] <- tAUC[1,1]
mcp.df.auc[i,2] <- tAUC[1,2]
mcp.df.auc[i,3] <- tAUC[1,3]
mcp.df.auc[i,4] <- tAUC[1,4]
mcp.df.auc[i,5] <- tAUC[1,5]
mcp.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.Tlocoh[,3:8], pred.df.Tlocoh[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.Tlocoh[,3:8], pred.df.Tlocoh[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
TLoCoHo.df.auc[i,1] <- tAUC[1,1]
TLoCoHo.df.auc[i,2] <- tAUC[1,2]
TLoCoHo.df.auc[i,3] <- tAUC[1,3]
TLoCoHo.df.auc[i,4] <- tAUC[1,4]
TLoCoHo.df.auc[i,5] <- tAUC[1,5]
TLoCoHo.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.tdge[,3:8], pred.df.tdge[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.tdge[,3:8], pred.df.tdge[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
tdge.df.auc[i,1] <- tAUC[1,1]
tdge.df.auc[i,2] <- tAUC[1,2]
tdge.df.auc[i,3] <- tAUC[1,3]
tdge.df.auc[i,4] <- tAUC[1,4]
tdge.df.auc[i,5] <- tAUC[1,5]
tdge.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.BRB[,3:8], pred.df.BRB[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.BRB[,3:8], pred.df.BRB[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
BRB.df.auc[i,1] <- tAUC[1,1]
BRB.df.auc[i,2] <- tAUC[1,2]
BRB.df.auc[i,3] <- tAUC[1,3]
BRB.df.auc[i,4] <- tAUC[1,4]
BRB.df.auc[i,5] <- tAUC[1,5]
BRB.df.auc[i,6] <- tAUC[1,6]
tAUC <- colAUC(pred.df.PPA[,3:8], pred.df.PPA[,2], plotROC=FALSE, alg=c("Wilcoxon","ROC")) 
colAUC(pred.df.PPA[,3:8], pred.df.PPA[,2], plotROC=TRUE, alg=c("Wilcoxon","ROC"))
# tAUC <- colAUC(allPts,mcp.05, plotROC=FALSE, alg=c("Wilcoxon","ROC"))
PPA.df.auc[i,1] <- tAUC[1,1]
PPA.df.auc[i,2] <- tAUC[1,2]
PPA.df.auc[i,3] <- tAUC[1,3]
PPA.df.auc[i,4] <- tAUC[1,4]
PPA.df.auc[i,5] <- tAUC[1,5]
PPA.df.auc[i,6] <- tAUC[1,6]

##################################################################
#EMD
r.real <- raster(ncol=200, nrow=200, res=10)
extent(r.real) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.real[] <- 0
tab <- table(cellFromXY(r.real, pts2))
r.real[as.numeric(names(tab))] <- tab
r.real <- raster::as.matrix(r.real)
r.real <- pp(r.real)

#T-locoh
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.100.shp, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.50.shp, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.25.shp, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.10.shp, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(T.05.shp, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

TLoCoHo.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
TLoCoHo.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
TLoCoHo.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
TLoCoHo.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
TLoCoHo.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#tdge
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(tdge.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

tdge.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
tdge.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
tdge.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
tdge.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
tdge.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#BRB
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(BRB.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

BRB.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
BRB.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
BRB.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
BRB.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
BRB.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#PPA
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(PPA.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

PPA.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
PPA.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
PPA.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
PPA.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
PPA.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#MCP
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(mcp.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

mcp.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
mcp.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
mcp.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
mcp.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
mcp.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)

#KDE
r.100 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.100, sp)
extent(r.100) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.100[] <- 0
tab <- table(cellFromXY(r.100, e))
r.100[as.numeric(names(tab))] <- tab
r.100 <- raster::as.matrix(r.100)
r.100 <- pp(r.100)

r.50 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.50, sp)
extent(r.50) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.50[] <- 0
tab <- table(cellFromXY(r.50, e))
r.50[as.numeric(names(tab))] <- tab
r.50 <- raster::as.matrix(r.50)
r.50 <- pp(r.50)

r.25 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.25, sp)
extent(r.25) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.25[] <- 0
tab <- table(cellFromXY(r.25, e))
r.25[as.numeric(names(tab))] <- tab
r.25 <- raster::as.matrix(r.25)
r.25 <- pp(r.25)

r.10 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.10, sp)
extent(r.10) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.10[] <- 0
tab <- table(cellFromXY(r.10, e))
r.10[as.numeric(names(tab))] <- tab
r.10 <- raster::as.matrix(r.10)
r.10 <- pp(r.10)

r.05 <- raster(ncol=200, nrow=200, res=10)
e <- gIntersection(KDE.05, sp)
extent(r.05) <- extent(bbox(pts2)[1,1], bbox(pts2)[1,2],bbox(pts2)[2,1],bbox(pts2)[2,2])
r.05[] <- 0
tab <- table(cellFromXY(r.05, e))
r.05[as.numeric(names(tab))] <- tab
r.05 <- raster::as.matrix(r.05)
r.05 <- pp(r.05)

KDE.df.emd[i,1] <- wasserstein(r.real,r.100,p=1)
KDE.df.emd[i,2] <- wasserstein(r.real,r.50,p=1)
KDE.df.emd[i,3] <- wasserstein(r.real,r.25,p=1)
KDE.df.emd[i,4] <- wasserstein(r.real,r.10,p=1)
KDE.df.emd[i,5] <- wasserstein(r.real,r.05,p=1)
}
write.xlsx(mcp.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\mcp_area.xlsx') 
write.xlsx(mcp.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\mcp_edge.xlsx') 
write.xlsx(mcp.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\mcp_auc.xlsx') 
write.xlsx(mcp.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\mcp_emd.xlsx') 
write.xlsx(KDE.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\KDE_area.xlsx') 
write.xlsx(KDE.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\KDE_edge.xlsx') 
write.xlsx(KDE.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\KDE_auc.xlsx') 
write.xlsx(KDE.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\KDE_emd.xlsx') 
write.xlsx(PPA.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\PPA_area.xlsx') 
write.xlsx(PPA.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\PPA_edge.xlsx') 
write.xlsx(PPA.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\PPA_auc.xlsx') 
write.xlsx(PPA.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\PPA_emd.xlsx') 
write.xlsx(tdge.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\tdge_area.xlsx') 
write.xlsx(tdge.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\tdge_edge.xlsx') 
write.xlsx(tdge.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\tdge_auc.xlsx') 
write.xlsx(tdge.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\tdge_emd.xlsx')
write.xlsx(TLoCoHo.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\TLoCoHo_area.xlsx') 
write.xlsx(TLoCoHo.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\TLoCoHo_edge.xlsx') 
write.xlsx(TLoCoHo.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\TLoCoHo_auc.xlsx') 
write.xlsx(TLoCoHo.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\TLoCoHo_emd.xlsx') 
write.xlsx(BRB.df, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\BRB_area.xlsx') 
write.xlsx(BRB.df.edge, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\BRB_edge.xlsx') 
write.xlsx(BRB.df.auc, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\BRB_auc.xlsx') 
write.xlsx(BRB.df.emd, 'C:\\Users\\311808\\Documents\\HR_Analysis\\outputs\\linear\\BRB_emd.xlsx') 
################################################################################
