library(PharmacoGx)
library(openxlsx)

options(stringsAsFactors=FALSE)

##create sensitivity slices
    
    getGRAYrawData <-
      function(result.type=c("array", "list")){
        
        gray.raw.drug.sensitivity <- read.csv(file="/pfs/downloadGRAY2017SensRaw/Gray_data_raw_dose_response.csv", header=T, row.names=1, stringsAsFactors=FALSE)
        gray.raw.drug.sensitivity.list <- do.call(c, apply(gray.raw.drug.sensitivity, 1, list))
        
        gray.conc <- read.csv(file="/pfs/downloadGRAY2017SensRaw/Gray_drug_conc.csv", header=T, row.names=1, stringsAsFactors=FALSE)
        
        concentrations.no <- 9
        
        if(result.type == "array"){
          ## create the gray.drug.response object including information viablilities and concentrations for each cell/drug pair
          obj <- array(NA, dim=c(length(unique(gray.raw.drug.sensitivity[ , "cellline"])), length(unique(gray.raw.drug.sensitivity[ , "drug"])), 2, concentrations.no), dimnames=list(unique(gray.raw.drug.sensitivity[ , "cellline"]), unique(gray.raw.drug.sensitivity[ , "drug"]), c("concentration", "viability"), 1:concentrations.no))
        }
        fnexperiment <-
          function(values){
            cellline <- values["cellline"]
            drug <- values["drug"]
            values["drug_group_id"] <- gsub(" ", "", values["drug_group_id"]) #removes space from drug-id values -i.e. (" 50") - which was causing some values AUC/IC50 values to not be matched
            doses <- as.numeric(gray.conc[which(gray.conc[,"drug_group_id"] == values["drug_group_id"] & gray.conc[,"drug"] == values["drug"]), grep("^c", colnames(gray.conc))]) * 10 ^ 6 # micro molar
            
            if(concentrations.no > length(doses)) {doses <- c(doses, rep(NA, concentrations.no - length(doses)))}
            
            #responses <- as.numeric(unlist(strsplit(input.matrix["Activity Data\n(raw median data)"], split=",")))  #nature paper raw data
            
            responses <- NULL
            background <- median(as.numeric(values[grep("^background_od", names(values))]))
            ctrl <- median(as.numeric(values[sprintf("od%s.%s",0,1:3)])) - background
            ctrl <- ifelse(ctrl < 0, 1, ctrl)
            for( i in 1:concentrations.no)
            {
              res <- median(as.numeric(values[sprintf("od%s.%s",i, 1:3)])) - background
              res <- ifelse(res < 0, 0, res)
              responses <- c(responses, res/ctrl)
            }
            responses <- responses * 100
            if(result.type == "array"){
              obj[cellline,drug, "concentration", 1:length(doses)] <<- doses
              obj[cellline,drug, "viability", 1:length(responses)] <<- responses
            }else{
              return(list(cell=cellline, drug=drug, doses=doses, responses=responses))#paste(doses, collapse = ","), responses=paste(responses, collapse = ",")))
            }
          }
        
        gray.raw.drug.sensitivity.res <- mapply(fnexperiment, values=gray.raw.drug.sensitivity.list)
        if(result.type == "array"){
          return(list("data"=obj, "concentrations.no"=concentrations.no))
        }else{
          return(list("data"=gray.raw.drug.sensitivity.res, "concentrations.no"=concentrations.no))
        }
      }
    raw.sensitivity <- getGRAYrawData(result.type="list")
    
    
    con_tested <- raw.sensitivity$concentrations.no
    raw.sensitivity <- t(raw.sensitivity$data)
    raw.sensitivity <- t(apply(raw.sensitivity,1, function(x){unlist(x)}))
    
    rownames(raw.sensitivity)  <- sprintf("drugid_%s_%s",as.character(raw.sensitivity[ ,2]),as.character(raw.sensitivity[ ,1]))
    ## handle replicates
    tt <- rownames(raw.sensitivity)
    for(i in 1:length(tt)) {
      xx <- which(tt == tt[i])
      if(length(xx) > 1) {
        for(j in 1:length(xx)) {
          tt[xx[j]] <- paste(tt[xx[j]], j, sep="_")
        }
      }
    }
    rownames(raw.sensitivity) <- tt
    
    tt <- paste0("doses", con_tested)
    sensitivity.info <- raw.sensitivity[ , c(1, 2, 3, grep(tt, colnames(raw.sensitivity)))]
    colnames(sensitivity.info) <- c("cellid", "drugid", "min.Dose.uM", "max.Dose.uM")
    sensitivity.info <- cbind(sensitivity.info, "nbr.conc.tested"=con_tested)
    raw.sensitivity <- raw.sensitivity[ ,-c(1,2)]
    raw.sensitivity <- array(c(as.matrix(raw.sensitivity[ ,1:con_tested]), as.matrix(raw.sensitivity[ ,(con_tested+1):(2*con_tested)])), c(nrow(raw.sensitivity), con_tested, 2),
                             dimnames=list(rownames(raw.sensitivity), colnames(raw.sensitivity[ ,1:con_tested]), c("Dose", "Viability")))
    
    

    save(raw.sensitivity, sensitivity.info, tt, con_tested, file="/pfs/out/drug_norm_post.RData")


    raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity)/1000))

    dir.create("/pfs/out/slices/")

    for(i in seq_along(raw.sensitivity.x)){

      slce <- raw.sensitivity[raw.sensitivity.x[[i]],,]
      saveRDS(slce, file=paste0("/pfs/out/slices/gray2017_raw_sens_", i, ".rds"))

    }
