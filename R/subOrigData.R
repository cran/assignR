subOrigData = function(marker = "d2H", taxon = NULL, group = NULL, dataset = NULL,
                       age_code = NULL, mask = NULL, ref_scale = "VSMOW_H",
                       niter = 5000, genplot = TRUE) {
  
  #load data in function environment
  knownOrig = knownOrig
  knownOrig_samples = knownOrig$samples
  knownOrig_sites = knownOrig$sites
  knownOrig_sources = knownOrig$sources
  
  result = knownOrig_samples

  if(length(marker) > 1){
    stop("only one marker currently allowed")
  }
  if(!marker %in% colnames(knownOrig_samples)){
    stop("marker must be column name for isotope data field")
  }
  
  if(!is.null(taxon)){
    if(!all(taxon %in% unique(knownOrig_samples$Taxon))){
      warning("One or more taxa not present in database")
    }
    result = result[result$Taxon %in% taxon,]
  }
  
  if(!is.null(group)){
    if(!all(group %in% unique(knownOrig_samples$Group))){
      warning("One or more groups not present in database")
    }
    result = result[result$Group %in% group,]
  }
  
  if(!is.null(dataset)){
    if(!is.numeric(dataset)){
      warning("dataset format should now be a numeric dataset ID, see knownOrig_sources.rda")
    } else if(!all(dataset %in% unique(knownOrig_sources$Dataset_ID))){
      warning("One or more datasets not present in database")
    }
    result = result[result$Dataset_ID %in% dataset,] 
  }
  
  if(!is.null(age_code)){
    if(!all(age_code %in% unique(knownOrig_samples$Age_class))){
      warning("One or more age codes not present in database")
    }
    result = result[result$Age_class %in% age_code,]
  }
  
  result = result[!is.na(result[,marker]),]
  if(nrow(result) == 0){
    stop("No samples match query")
  }

  mask = check_mask(mask, knownOrig_sites) 
  if(!is.null(mask)){
    result_sites = knownOrig_sites[mask,]

    if(length(result_sites) > 0){
      result = result[result$Site_ID %in% result_sites$Site_ID,]
      if(nrow(result) > 0) {
        result_sites = result_sites[result_sites$Site_ID %in% 
                                      result$Site_ID,]
      } else{
        stop("No samples found in mask\n")
      }
    } else{
      stop("No sites found in mask\n")
    }
    
  } else{
    result_sites = knownOrig_sites[knownOrig_sites$Site_ID %in%
                                     result$Site_ID,]
  }
  
  message(paste(length(result[,1]),"samples are found from", 
          length(result_sites), "sites"))

  result_sources = knownOrig_sources[knownOrig_sources$Dataset_ID %in%
                                       result$Dataset_ID,]
  
  if(!is.null(ref_scale)){
    if(marker == "d2H"){
      result = merge(result, result_sources[,c("Dataset_ID", "H_cal")], 
            by = "Dataset_ID", all.x = TRUE)
    } else{
      result = merge(result, result_sources[,c("Dataset_ID", "O_cal")], 
            by = "Dataset_ID", all.x = TRUE)
    }
    class(result) = "SOD"
    trans_out = refTrans(result, marker, ref_scale, niter)
    result_data = merge(result_sites, trans_out$data, by = "Site_ID", 
                        all.x = FALSE, duplicateGeoms = TRUE)

    return_obj = list("data" = result_data, "sources" =
                        result_sources, "chains" = trans_out$chains,
                      "marker" = marker)
    class(return_obj) = "subOrigData"
    message(paste(length(result_data$Sample_ID), "samples from", 
                  length(unique(result_data$Site_ID)), 
                  "sites in the transformed dataset"))
  } else{
    result_data = merge(result_sites, result, by = "Site_ID", 
                        all.x = FALSE, duplicateGeoms = TRUE)
    
    return_obj = list("data" = result_data, "sources" = result_sources,
                      "chains" = NULL, "marker" = marker)
    class(return_obj) = "subOrigData"
  }
  
  if(!inherits(genplot, "logical")) {
    message("genplot should be logical (T or F), using default = T")
    genplot = TRUE
  }
  
  if(genplot){
    if(is.null(mask)){
      wrld_simpl = wrld_simpl
      plot(wrld_simpl, axes = TRUE)
      points(result_data, col = "red", cex = 0.5)
    } else{
      plot(mask, axes = TRUE)
      points(result_data, col = "red")
    }
  }
  
  return(return_obj)
}
