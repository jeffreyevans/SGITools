#' @title Sage grouse fence collision model
#' @description Calculates sage grouse fence collision risk following methods presented in Stevens et al., (2013) 
#'
#' @param x            SpatialPointsDataFrame sp class object of lek locations
#' @param y            RasterLayer raster class object of elevation (DEM)
#' @param d            Buffer distance (default to 3000m as defined in model) 
#' @param pad          pad d by raster resoultion ^2 to account for edge effect (TRUE/FALSE)
#' @param file.name    Output file name, if specified results will be mosaiced into single raster 
#' @param cleanup      Delete lek-level rasters (TRUE/FALSE) 
#'
#' @return raster files in tiff format. If "file.name" is specified individual lek rasters are mosaiced into a single raster.
#'
#' @details
#' Uses point locations of sage-grouse lek and a 30m Digital Elevation Model with same extent and projection as lek data to build model. 
#' Expected estimate range 0-3 where; low risk 0-0.49, moderate risk 0.50-0.99, high risk >= 1
#'
#' \itemize{ 
#' \item  Create a point layer representing sage-grouse lek locations.
#' \item  Create dissolved lek buffers 
#' \item  Subset elevation raster to each buffer
#' \item  Calculate terrain ruggedness index raster (Riley et al. 1999) 
#' \item  Calculate Euclidean distance to nearest lek(s)
#' \item  Use TRI and distance to lek and the following algebraic equation to predict collision risk (See Stevens et al., 2013 for details) 
#' \item  regression equation: y = 78 * exp(beta0 + beta1 * TRI + beta2 * distance), 
#' \item    where beta0 = -3.3254377759, beta1 = -0.2504710567, beta2 = -0.0006119843.
#' }
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#' 
#' @references Riley, S.J., S.D. DeGloria and R. Elliot (1999) A terrain ruggedness index that quantifies topographic heterogeneity, Intermountain Journal of Sciences 5(1-4):23-27.
#' @references Stevens, B. S., D.E. Naugle, B. Dennis, J.W. Connelly, T. Griffiths, and K. P. Reese (2013) Mapping sage-grouse fence-collision risk: Spatially explicit models for targeting conservation implementation. Wildlife Society Bulletin, 37(2):409-415
#'
#' @examples 
#' \dontrun{
#'  library(raster)
#'  library(rgeos)
#'  library(rgdal)
#'  
#'  data(leks)  
#'  data(elev)
#'
#'  ext <- as(extent(elev), "SpatialPolygons")
#'  ext <- gBuffer(ext, width = -3000)
#'    proj4string(ext) <- proj4string(leks)
#'  leks <- intersect(leks, ext)
#'  
#'  fence.collision(leks, elev, file.name = "fence_collision.tif", cleanup = TRUE)  
#'    fc <- raster("fence_collision.tif")
#'    plot(fc)
#'   
#'  # Reclassify to risk ranges and plot 
#'  rc <- function(x) {
#'    ifelse(x <= 0.49, 1, 
#'      ifelse(x > 0.40 & x <= 0.99, 2, 
#'  	  ifelse(x > 1, 3, NA)))
#'  }
#'  fc.class <- calc(fc, fun = rc, filename = "fence_collision_class.tif", overwrite = TRUE) 
#'    plot(fc.class, col=c("green","yellow","red"), legend = FALSE)
#'      points(leks, pch = 20, cex = 0.75)
#'      legend("topright", legend=c("Low risk","Moderate risk","High risk"),
#'  	       fill=c("green","yellow","red"))
#'  	title("Sage grouse fence collision risk - Pinedale anticline, WY")
#' }	
#' @export
fence.collision <- function(x, y, d = 3000, pad = TRUE, file.name = NULL, 
                            cleanup = FALSE) { 
  if (!class(x) == "SpatialPointsDataFrame" & !class(x) == "SpatialPoints") 
      stop(deparse(substitute(x)), " must be a sp points object")
  if (!class(y) == "RasterLayer") 
      stop(deparse(substitute(y)), " must be a RasterLayer class object")							
  if(pad) d = d + raster::res(y)[1] * 2
    b <- rgeos::gBuffer(x, byid = FALSE, id = NULL, width = d, quadsegs = 10)
      b <- sp::disaggregate(b)
  for(i in 1:length(b)) { 
	psub <- b[i, ]
	  sp::proj4string(psub) <- sp::proj4string(x)
	    x.sub <- x[which(rgeos::gContains(psub, x, byid = TRUE) == TRUE),] 
		  cr <- raster::crop(y, raster::extent(psub), snap = "out")
            y.sub <- raster::mask(x = cr, mask = psub)   
          rough <- tri(y.sub)
	    edist <- raster::distanceFromPoints(y.sub, sp::SpatialPoints(x.sub))
	  edist <- raster::mask(x = edist, mask = psub)
 	y.hat = (78 * exp(-3.3254377759 + -0.2504710567 * rough + -0.0006119843 * edist))
	  sp::proj4string(y.hat) <- sp::proj4string(y) 	
	  if(pad) {
	    d <- (d - raster::res(y)[1] * 2)
        b.sub <- rgeos::gBuffer(psub, byid = FALSE, id = NULL,  
		                        width = -raster::res(y)[1], quadsegs = 10)
        # b.sub <- sp::disaggregate(b.sub)
        y.hat <- raster::mask(x = y.hat, mask = b.sub)
          sp::proj4string(y.hat) <- sp::proj4string(y) 			
 	  }
	raster::writeRaster(y.hat, filename = paste0("fc", i, ".tif"), overwrite = TRUE)
	  cat("lek model =", i, "of", length(b), "\n")
    }
  if(!is.null(file.name)) {
    cat("Mosaicing results to: ", file.name, "\n")
	rlist <- list()
	  for(i in 1:length(b)) { rlist[[i]] <- raster::raster(paste0("fc", i, ".tif")) }
	    rlist$fun <- "mean"
	    rlist$overwrite = TRUE
	    rlist$filename = file.name
    	  do.call(raster::mosaic, rlist)
        cat("fence collision model written to:", file.name, "\n")	
	  if(cleanup) unlink(paste0("fc", 1:length(b), ".tif"))						  
    }
}  
