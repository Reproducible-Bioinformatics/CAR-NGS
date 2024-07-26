dimensions <- function(result_dir_path, input_file_path, separator) {
  # Type checking.
  if (typeof(result_dir_path) != "character") {
    stop(
      paste(
        "result_dir_path type is",
        paste0(typeof(result_dir_path), "."),
        "It should be \"character\""
      )
    )
  }
  if (typeof(input_file_path) != "character") {
    stop(
      paste(
        "input_file_path type is",
        paste0(typeof(input_file_path), "."),
        "It should be \"character\""
      )
    )
  }
  if (typeof(separator) != "character") {
    stop(
      paste(
        "separator type is",
        paste0(typeof(separator), "."),
        "It should be \"character\""
      )
    )
  }

  if (!is_running_in_docker()) {
    if (!file.exists(input_file_path)) {
      stop(paste("input_file_path:", input_file_path, "does not exist."))
    }
    if (!dir.exists(result_dir_path)) {
      stop(paste("result_dir_path:", result_dir_path, "does not exist."))
    }
  }

  # Obtain parent folder from input_file_path.
  parent_folder <- dirname(input_file_path)

  # Obtain matrix name and file format from input_file_path.
  input_file_path_parts <- strsplit(basename(input_file_path), "\\.")[[1]]
  matrix_name <- input_file_path_parts[1]
  format <- input_file_path_parts[2] # Uses extension as an heuristic.

  run_in_docker(
    image_name = "docker.io/repbioinfo/r332.2017.01:latest",
    volumes = list(
      c(result_dir_path, "/data"),
      c(parent_folder, "/scratch")
    ),
    additional_arguments = c(
      "Rscript /home/main.R",
      matrix_name,
      format,
      paste0('"', separator, '"')
    )
  )
}

#' Maps a path to host volumes.
#'
#' @param path a normalized, absolute path.
#' @return The absolute host path, if it exists.
docker_mount_mapper <- function(path) {
  if (!is_running_in_docker()) {
    return(path)
  }

  # Since we are running in docker, we can use our hostname as an heuristic
  # to obtain our id.
  # TODO: Check if this assumption is right for other container engines and
  #       generalize this implementation.
  hostname <- Sys.info()["nodename"]
  output <- system2("docker",
    args = paste("inspect -f '{{ json .Mounts }}'", hostname),
    stdout = TRUE,
  )
  parsed_output <- jsonlite::fromJSON(output)

  # Iterate over mounts, return first match.
  for (i in seq_len(nrow(parsed_output))) {
    destination <- parsed_output[i, ]$Destination
    if (startsWith(path, destination)) {
      source <- parsed_output[i, ]$Source
      path <- sub(destination, source, path)
      return(path)
    }
  }
  return(path)
}

#' Gets the absolute path of a file.
#'
#' @param path a normalized, absolute path.
#' @return The absolute host path, if it exists.
absolute_path_mapper <- function(path) {
  return(normalizePath(path, mustWork = FALSE))
}

normalize_path <- function(path, path_mappers = c()) {
  path_mappers <- c(
    absolute_path_mapper, # Adds absolute_path_mapper at the beginning.
    path_mappers
  )
  for (mapper in path_mappers) {
    path <- mapper(path)
  }
  return(path)
}

#' Append a timed directory with the format %Y%m%d_%H%M%S to the
#' target_directory path.
#'
#' @param target_directory The target directory.
#' @return Timed directory with the format %Y%m%d_%H%M%S appended to
#' target_directory path.
modifier_timed_directory <- function(target_directory) {
  return(
    file.path(target_directory, strftime(time, "%Y%m%d_%H%M%S"))
  )
}

with_scratch_directory <- function(
    source_directory,
    target_directory,
    modifier = modifier_timed_directory,
    callback_function,
    cleanup_after = FALSE,
    copy_pattern = NULL) {
  if (!dir.exists(source_directory)) {
    stop(paste("source_directory", source_directory, "doesn't exist."))
  }
  if (typeof(callback_function) != "closure") {
    stop("callback_function is not a closure.")
  }
  if (modifier != NULL) {
    if (typeof(modifier) != "closure") {
      stop("modifier is not a closure.")
    }
    target_directory <- modifier(target_directory)
  }
  # Create directory if doesn't exist.
  dir.create(target_directory, recursive = TRUE)

  # Calls the callback function.
  callback_function(target_directory)

  if (!cleanup_after) {
    return
  }
  # Copy back from target_directory to source_directory.
  file.copy(
    list.files(
      pattern = copy_pattern,
      path = target_directory,
      full.names = TRUE,
      all.files = TRUE,
      recursive = TRUE,
    ),
    source_directory
  )
  unlink(target_directory, recursive = TRUE)
}
