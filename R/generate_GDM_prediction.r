#' Generates a prediction for a genetically local area
#' for a selected species and location
#' 
#' @param  selected_location -- a selected point [Required]
#' @param  gdm_model_file    -- an R object bundling everything needed for model [Required]

#' @return a raster containing the prediction of genetically local area
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export
#' @examples
#' \dontrun{
#' generate_GDM_prediction <- function( c(150.45, -33.66), CallSerr_genetic_model.Rd )
#' }

generate_GDM_prediction <- function(selected_location, gdm_model_file) {

   s1_point <- selected_location
   
   ### load required packages
   require(gdm); require(raster); require(fields); require(sp);

   ### load object containing model
   ### and all data required to run it
   load(gdm_model_file)

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

   return(genetically_local_raster)
}

