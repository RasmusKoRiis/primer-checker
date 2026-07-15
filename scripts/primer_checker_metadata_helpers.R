primer_checker_metadata_columns <- c(
  "SampleID",
  "Sample_Date",
  "CT_H1",
  "CT_H3",
  "CT_INFA",
  "CT_INFB",
  "Triplex-InfA_CT",
  "Triplex-InfB_CT",
  "CT_RSVA",
  "CT_RSVB",
  "Triplex-SC2_CT"
)

first_existing_column <- function(data, candidates) {
  matches <- candidates[candidates %in% names(data)]
  if (length(matches) == 0) {
    return(NULL)
  }
  matches[[1]]
}

metadata_value <- function(data, candidates, default = "") {
  column <- first_existing_column(data, candidates)
  if (is.null(column)) {
    return(rep(default, nrow(data)))
  }
  value <- data[[column]]
  value <- as.character(value)
  value[is.na(value)] <- default
  value
}

format_metadata_date <- function(value) {
  if (inherits(value, "Date")) {
    return(format(value, "%Y-%m-%d"))
  }
  parsed <- suppressWarnings(as.Date(value))
  out <- ifelse(is.na(parsed), as.character(value), format(parsed, "%Y-%m-%d"))
  out[is.na(out)] <- ""
  out
}

empty_primer_checker_metadata <- function(n) {
  data.frame(
    SampleID = rep("", n),
    Sample_Date = rep("", n),
    CT_H1 = rep("", n),
    CT_H3 = rep("", n),
    CT_INFA = rep("", n),
    CT_INFB = rep("", n),
    `Triplex-InfA_CT` = rep("", n),
    `Triplex-InfB_CT` = rep("", n),
    CT_RSVA = rep("", n),
    CT_RSVB = rep("", n),
    `Triplex-SC2_CT` = rep("", n),
    check.names = FALSE
  )
}

normalise_primer_checker_metadata <- function(metadata) {
  for (column in primer_checker_metadata_columns) {
    if (!column %in% names(metadata)) {
      metadata[[column]] <- ""
    }
  }
  metadata <- metadata[, primer_checker_metadata_columns, drop = FALSE]
  metadata[] <- lapply(metadata, function(value) {
    value <- as.character(value)
    value[is.na(value)] <- ""
    value
  })
  metadata
}

write_primer_checker_metadata <- function(metadata, output_dir, filename_prefix) {
  metadata <- normalise_primer_checker_metadata(metadata)
  output_path <- file.path(
    output_dir,
    paste0(filename_prefix, " PRIMER CHECKER METADATA - ", format(Sys.Date(), "%U-%Y"), ".csv")
  )
  write.csv(metadata, output_path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  message("Primer checker metadata written to: ", output_path)
  invisible(output_path)
}

build_influenza_primer_checker_metadata <- function(data) {
  metadata <- empty_primer_checker_metadata(nrow(data))
  metadata$SampleID <- metadata_value(data, c("key", "Submitting_Sample_Id", "sample_id", "SampleID"))
  metadata$Sample_Date <- format_metadata_date(metadata_value(data, c("prove_tatt", "Collection_Date", "Sample_Date")))

  subtype <- metadata_value(data, c("Subtype", "ngs_sekvens_resultat"))
  inf_type <- metadata_value(data, c("INFType"))

  ct_h1 <- metadata_value(data, c("CT_H1", "ct_h1", "Ct_H1", "ct_H1", "h1_ct", "H1_CT", "pcr_h1_ct", "PCR_H1_CT"))
  ct_h3 <- metadata_value(data, c("CT_H3", "ct_h3", "Ct_H3", "ct_H3", "h3_ct", "H3_CT", "pcr_h3_ct", "PCR_H3_CT"))
  ct_infa <- metadata_value(data, c("CT_INFA", "ct_infa", "Ct_INFA", "ct_inf_a", "infa_ct", "InfA_CT"))
  ct_bvic <- metadata_value(data, c("pcr_bvic_ct", "PCR_BVIC_CT"))
  ct_byam <- metadata_value(data, c("pcr_byam_ct", "PCR_BYAM_CT"))
  ct_infb <- metadata_value(data, c("CT_INFB", "ct_infb", "Ct_INFB", "ct_inf_b", "infb_ct", "InfB_CT"))
  ct_infb <- ifelse(ct_infb != "", ct_infb, ifelse(ct_bvic != "", ct_bvic, ct_byam))
  triplex_infa <- metadata_value(data, c("Triplex-InfA_CT", "Triplex_InfA_CT", "triplex_infa_ct", "triplex_InfA_ct", "InfA_Triplex_CT"))
  triplex_infb <- metadata_value(data, c("Triplex-InfB_CT", "Triplex_InfB_CT", "triplex_infb_ct", "triplex_InfB_ct", "InfB_Triplex_CT"))

  metadata$CT_H1 <- ifelse(grepl("H1", subtype, ignore.case = TRUE), ct_h1, "")
  metadata$CT_H3 <- ifelse(grepl("H3", subtype, ignore.case = TRUE), ct_h3, "")
  metadata$CT_INFA <- ifelse(inf_type == "A" | grepl("^A/", subtype), ct_infa, "")
  metadata$CT_INFB <- ifelse(inf_type == "B" | grepl("^B|VICTORIA|YAMAGATA", subtype, ignore.case = TRUE), ct_infb, "")
  metadata$`Triplex-InfA_CT` <- ifelse(inf_type == "A" | grepl("^A/", subtype), triplex_infa, "")
  metadata$`Triplex-InfB_CT` <- ifelse(inf_type == "B" | grepl("^B|VICTORIA|YAMAGATA", subtype, ignore.case = TRUE), triplex_infb, "")
  normalise_primer_checker_metadata(metadata)
}

build_rsv_primer_checker_metadata <- function(data) {
  metadata <- empty_primer_checker_metadata(nrow(data))
  metadata$SampleID <- metadata_value(data, c("key", "Submitting_Sample_Id", "sample_id", "SampleID"))
  metadata$Sample_Date <- format_metadata_date(metadata_value(data, c("prove_tatt", "collection_date", "Sample_Date")))

  subtype <- metadata_value(data, c("Subtype", "ngs_sekvens_resultat"))
  ct_rsva <- metadata_value(data, c("CT_RSVA", "ct_rsva", "Ct_RSVA", "rsva_ct", "RSVA_CT"))
  ct_rsvb <- metadata_value(data, c("CT_RSVB", "ct_rsvb", "Ct_RSVB", "rsvb_ct", "RSVB_CT"))

  metadata$CT_RSVA <- ifelse(grepl("RSVA|^A$", subtype, ignore.case = TRUE), ct_rsva, "")
  metadata$CT_RSVB <- ifelse(grepl("RSVB|^B$", subtype, ignore.case = TRUE), ct_rsvb, "")
  normalise_primer_checker_metadata(metadata)
}

build_sc2_primer_checker_metadata <- function(data) {
  metadata <- empty_primer_checker_metadata(nrow(data))
  metadata$SampleID <- metadata_value(data, c("key", "Submitting_Sample_Id", "sample_id", "SampleID"))
  metadata$Sample_Date <- format_metadata_date(metadata_value(data, c("prove_tatt", "covv_collection_date", "Sample_Date")))
  metadata$`Triplex-SC2_CT` <- metadata_value(
    data,
    c("Triplex-SC2_CT", "Triplex_SC2_CT", "triplex_sc2_ct", "SC2_CT", "ct_sc2", "Ct_SC2", "sars_cov_2_ct", "SARS_CoV_2_CT")
  )
  normalise_primer_checker_metadata(metadata)
}


