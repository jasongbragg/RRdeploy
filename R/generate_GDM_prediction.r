#' Generates a prediction for a genetically local area
#' for a selected species and location
#' 
#' @param  longitude         -- selected longitude
#' @param  latitude          -- selected latitude
#' @param  gdm_model         -- an R list bundling everything required to predict under a GDM 
#' @param  gdm_model_file    -- an R object file containing a model list (called "md")
#' @param  species_name      -- name of species (as it appears in a gdm model file, which ends "_genetic_model.Rd"
#' @param  model_dir         -- directory in which gdm model can be found 
#' @param  out_dir           -- location to place raster of results (if NULL, return raster object)
#' @return a raster containing the prediction of genetically local area
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' generate_GDM_prediction <- function( longitude=150.45, latitude=-33.66, gdm_model_file="CallSerr_genetic_model.Rd", out_dir )
#' }

generate_GDM_prediction <- function(longitude, latitude, gdm_model_file=NULL, species_name=NULL, model_dir=NULL, gdm_model=NULL, out_dir=NULL) {

   s1_point <- c(longitude, latitude) 
   
   ### load required packages
   require(gdm); require(raster); require(fields); require(sp);

   md <- NULL

   ### if model is supplied, use it
   if ( !is.null(gdm_model) ) {
      md <- gdm_model
   } else {

      # if model is not supplied as list
      # first see if file is supplied
      if (!is.null(gdm_model_file)) {

         gmf <- gdm_model_file

      } else {

         # if file name not supplied
         # attempt to deduce file name with species name and directory
         # if still not found, error
         if ( is.null(species_name) | is.null(model_dir) ) {
            cat("  No model provided. Must specify either:                       \n")
            cat("                             model list object                  \n")
            cat("                             a file containing object           \n")
            cat("                             or taxon name plus model directory \n")
            stop()
         }
         gfm <- paste(model_dir, species_name, "_genetic_model.Rd", sep="")

      }
      load(gmf)
   }

   if (is.null(md)) {
      cat( " Something went wrong. GDM model object not found by generate_GDM_prediction() function \n" )
   }

   sdata <- md$sdata
   r.pts <- rasterToPoints(sdata, spatial=TRUE)

   r.pts@data <- data.frame(long=coordinates(r.pts)[,1],
                      lat=coordinates(r.pts)[,2], r.pts@data)                         

   s2_ll      <- r.pts@data[,1:2]
   null_dist  <- rep(1, nrow(s2_ll))
   null_wght  <- rep(1, nrow(s2_ll))

   s1_ll     <- cbind(rep(s1_point[1], nrow(s2_ll)), rep(s1_point[2], nrow(s2_ll)))

   gdm_prd   <- cbind(null_dist, null_wght, s1_ll, s2_ll)
   colnames(gdm_prd)  <- c("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord")
   gdm_prd   <- as.data.frame(gdm_prd)

   if(md$Q) {
      qdata <- md$qdata
      s1_Q  <- extract(qdata, s1_ll)
      s2_Q  <- extract(qdata, s2_ll)

      colnames(s1_Q) <- paste("s1.Q", 1:ncol(s1_Q), sep="")
      colnames(s2_Q) <- paste("s2.Q", 1:ncol(s2_Q), sep="")

      gdm_prd   <- cbind(gdm_prd, s1_Q, s2_Q)

   }

   if(md$E) {
      edata  <- md$edata
      enames <- names(edata)
      s1_E  <- extract(edata, s1_ll)
      s2_E  <- extract(edata, s2_ll)

      colnames(s1_E) <- paste("s1.", enames, 1:ncol(s1_E), sep="")
      colnames(s2_E) <- paste("s2.", enames, 1:ncol(s2_E), sep="")

      gdm_prd   <- cbind(gdm_prd, s1_E, s2_E)
 
      # in case environmental variable undefined in parts of srast
      gdm_prd <- gdm_prd[ rowSums(is.na(gdm_prd)) <= 0, ]
      s2_ll   <- gdm_prd[,5:6]
   }

   gdm_model <- md$model
   gdm_prediction  <- predict(gdm_model, gdm_prd)
   point2all_rast  <- raster(sdata)
   point2all_rast  <- rasterize(s2_ll,
                               point2all_rast,
                               field=gdm_prediction)

   clip_prediction <- mask(point2all_rast, md$confidence_polygon) 

   thresh <- md$threshold
   genetically_local_raster <-  clip_prediction
   values(genetically_local_raster)[values(genetically_local_raster) > thresh ] <- NA
   values(genetically_local_raster)[ !is.na(values(genetically_local_raster)) ] <- 1

   if ( is.null(out_dir) ) {
      return(genetically_local_raster)
   } else {
      # get taxon name from the model file name
      tmpGenLocalRasterFile <- paste(sessionDir,"/geneticLocal.tif",sep="")
      writeRaster(genetically_local_raster, tmpGenLocalRasterFile, format="GTiff", overwrite = TRUE)
      return(tmpGenLocalRasterFile)
   }
   
}

