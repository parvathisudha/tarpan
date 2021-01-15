# author = Jonathan Huang
library(data.table)
library(dplyr)
library(stringr)
## notes: counts are in x,y, where x is tier 1 and y is tier 2. Tier 2 is more permissive. Going with tier 1 cause DP is in tier 1

####
## Helper Functions

get_base_count <- function(info, format, base){
  # info is normal/tumor. Format is the format column
  target_index <- match(paste0(base, "U"), unlist(str_split(format, ":")))
  count <- unlist(str_split(info, ":"))[target_index] %>% str_split(",") %>% unlist
  return(as.numeric(count[1]))
}

get_snv_vaf <- function(info, format, base, ref){
  # info is normal/tumor. Format is the format column
  # is_alt is for if you're not looking at the vaf for the alt for some reason....
  format <- str_split(format, ":") %>% unlist
  target_index <- match(paste0(base, "U"), format)
  ref_index <- match(paste0(ref, "U"), format)
  info <- unlist(str_split(info, ":"))
  ref <- info[ref_index] %>% str_split(",") %>% unlist
  ref <- as.numeric(ref[1])
  count <- info[target_index] %>% str_split(",") %>% unlist
  count <- as.numeric(count[1])
  ## for tier 1 use index of 1, for 2 use index of 2
  return( signif(count/(ref + count), digits = 3))
}

## need to edit
get_indel_vaf <- function(info, format){
  # info is normal/tumor. Format is the format column
  format <- str_split(format, ":") %>% unlist
  alt_index <- match("TIR", format)
  DP_index <- match("TAR", format)
  info <- unlist(str_split(info, ":"))
  DP <- info[DP_index] %>% str_split(",") %>% unlist  
  DP <- DP[1] %>% as.numeric
  count <- info[alt_index] %>% str_split(",") %>% unlist
  count <- as.numeric(count[1])
  ## for tier 1 use index of 1, for 2 use index of 2
  return(signif(count/(DP + count), digits=3))
}


####
## Main function which runs the helper functions

get_all_mutation_info <- function(dbhandle, sample) {
  query <- "SELECT * from GENOM_MUTATIONS where sample_id = '"
  query <- paste0(query, sample,"'")
  genome_muts <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)
  genome_muts$mut_id <- NULL
  genome_muts$sample_id <- NULL
  genome_muts <- genome_muts %>% as.data.table %>% get_vaf
  genome_muts
}

## it seems he doesn't keep mut_tools in the table so I made edits to pretty much grab it from the format column
get_vaf <- function(dt, tool=FALSE){
  ## dt is the data.table, tool is the tool i.e strelka_snvs/indel etc
  s_snvs <- function(dt){
    # dt[grepl(":AU:", format), c("normal_vaf", "tumor_vaf") := .(1, 0)]
    dt[grepl(":AU:", format), c("normal_vaf", "tumor_vaf") :=
         .(mapply(get_snv_vaf, normal, format, alt, ref),
           mapply(get_snv_vaf, tumor, format, alt, ref))]
    return(dt)
  }
  s_indels <- function(dt){
    dt[grepl(":TIR:", format), c("normal_vaf", "tumor_vaf") :=
        .(mapply(get_indel_vaf, normal, format),
          mapply(get_indel_vaf, tumor, format))]
    return(dt)
  }
  
  if (tool == FALSE) {
    ## if false just run both i guess
    dt <- dt %>% s_snvs %>% s_indels
  } else if (tool == "strelka_snvs") { dt <- dt %>% s_snvs}
  else if (tool == "strelka_indels") { dt <- dt %>% s_indels}
  
  ### change all NAs to char "NA"
  # dt[is.na(normal_vaf), normal_vaf := "N/A"]
  # dt[is.na(tumor_vaf), normal_vaf := "N/A"]
  
  
  return(dt)
}

