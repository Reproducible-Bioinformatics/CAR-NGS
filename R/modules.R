
modules <- function(input_dir_path, organism, count_matrix, metadata_file, count_sep=NULL, meta_sep=NULL){

 # Type checking.
  if (typeof(input_dir_path) != "character") {
    stop(paste("input_dir_path type is", paste0(typeof(input_dir_path), "."), "It should be \"character\""))
  }
  if (typeof(organism) != "character") {
    stop(paste("organism type is", paste0(typeof(organism), "."), "It should be \"character\""))
  }
  if (typeof(count_matrix) != "character") {
    stop(paste("count_matrix type is", paste0(typeof(count_matrix), "."), "It should be \"character\""))
  }
  if (typeof(metadata_file) != "character") {
    stop(paste("metadata_file type is", paste0(typeof(metadata_file), "."), "It should be \"character\""))
  }
  if (!is.null(count_sep) && typeof(count_sep) != "character") {
    stop(paste("count_sep type is", paste0(typeof(count_sep), "."),"It should be \"character\""))
  }
  if (!is.null(meta_sep) && typeof(meta_sep) != "character") {
    stop(paste("meta_sep type is", paste0(typeof(meta_sep), "."),"It should be \"character\""))
  }

  # Check if input_dir_path exists
  if (!is_running_in_docker()) {
    if (!dir.exists(input_dir_path)) {
      stop(paste("input_dir_path:", input_dir_path, "does not exist."))
    }
  }

  # Setting the separators as NULL
  if (is.null(count_sep)) { count_sep = "NULL" }
  if (is.null(meta_sep)) { meta_sep = "NULL" }

  # Executing the docker job
  run_in_docker(
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
