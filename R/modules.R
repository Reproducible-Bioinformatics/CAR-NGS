#' Modules
#'
#' The `modules` function runs an RNA-seq data processing pipeline inside a
#' Docker container.
#' It analyzes a gene expression count matrix and associated metadata to
#' identify gene modules.
#' The pipeline generates heatmaps, gene ontology (GO) analysis for biological
#' processes, and outputs TPM-normalized counts for each gene module.
#' The results are saved in the specified output directory.
#'
#' @param input_dir_path A character string specifying the path to the
#' directory containing input files (count matrix and metadata).
#' @param organism A character string specifying the organism. Supported values
#' are "Homo sapiens", "Mus musculus", or "Drosophila melanogaster".
#' @param count_matrix A character string specifying the name of the count
#' matrix file.
#' @param metadata_file A character string specifying the name of the metadata
#' file.
#' @param count_sep (Optional) A character string indicating the separator for
#' the count matrix file (e.g., "tab", ",").
#' @param meta_sep (Optional) A character string indicating the separator for
#' the metadata file (e.g., "tab", ",").
#'
#' @return The function does not return any values directly. It saves results
#' (heatmaps, GO analysis, and normalized counts) as files in the output
#' directory.
#' @examples
#' \dontrun{
#' modules(
#'   "/path/to/input", "Homo sapiens", "counts.csv", "metadata.csv",
#'   count_sep = ",", meta_sep = ","
#' )
#' }
#'
#' @export
modules <- function(input_dir_path,
                    organism,
                    count_matrix,
                    metadata_file,
                    count_sep = NULL,
                    meta_sep = NULL) {
  # Type checking.
  if (typeof(input_dir_path) != "character") {
    stop(
      paste(
        "input_dir_path type is",
        paste0(typeof(input_dir_path), "."),
        "It should be \"character\""
      )
    )
  }
  if (typeof(organism) != "character") {
    stop(
      paste(
        "organism type is",
        paste0(typeof(organism), "."),
        "It should be \"character\""
      )
    )
  }
  if (typeof(count_matrix) != "character") {
    stop(paste(
      "count_matrix type is",
      paste0(typeof(count_matrix), "."),
      "It should be \"character\""
    ))
  }
  if (typeof(metadata_file) != "character") {
    stop(paste(
      "metadata_file type is",
      paste0(typeof(metadata_file), "."),
      "It should be \"character\""
    ))
  }
  if (!is.null(count_sep) && typeof(count_sep) != "character") {
    stop(paste(
      "count_sep type is",
      paste0(typeof(count_sep), "."),
      "It should be \"character\""
    ))
  }
  if (!is.null(meta_sep) && typeof(meta_sep) != "character") {
    stop(paste(
      "meta_sep type is",
      paste0(typeof(meta_sep), "."),
      "It should be \"character\""
    ))
  }

  # Check if organisms is in supported organisms.
  supported_organisms <- list(
    "Homo sapiens" = "Hs",
    "Mus musculus" = "Mm",
    "Drosophila melanogaster" = "Dm"
  )
  if (!(organism %in% names(supported_organisms))) {
    stop(
      paste(
        "Invalid organism:",
        organism,
        "Supported values are:",
        paste(names(supported_organisms), collapse = ", ")
      )
    )
  }

  # Check if input_dir_path exists
  if (!rrundocker::is_running_in_docker()) {
    if (!dir.exists(input_dir_path)) {
      stop(paste("input_dir_path:", input_dir_path, "does not exist."))
    }
  }

  # Setting the separators as NULL
  if (is.null(count_sep)) {
    count_sep <- "NULL"
  }
  if (is.null(meta_sep)) {
    meta_sep <- "NULL"
  }

  # Executing the docker job
  rrundocker::run_in_docker(
    image_name = paste0("hedgelab/rnaseq_modules:image5"),
    volumes = list(
      c(input_dir_path, "/scratch")
    ),
    additional_arguments = c(
      "Rscript /home/modules.R",
      organism,
      count_matrix,
      metadata_file,
      count_sep,
      meta_sep
    )
  )
}
