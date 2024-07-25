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

run_in_docker <- function(image_name,
                          volumes = list(),
                          additional_arguments = c()) {
  base_command <- "run --privileged=true --platform linux/amd64 --rm"
  for (volume in volumes) {
    base_command <- paste(base_command, "-v", volume)
    base_command <- paste(base_command, "-v", paste(
      volume[1],
      volume[2],
      sep = ":"
    ))
  }
  base_command <- paste(base_command, image_name)
  for (argument in additional_arguments) {
    base_command <- paste(base_command, argument)
  }
  system2("docker", args = base_command, stdout = TRUE)
}

#' Check if the script is running in a container.
#'
#' @returns A truthy value indicating the state.
is_running_in_docker <- function() {
  dockerenv_exists <- file.exists("/.dockerenv")
  cgroup_exists <- file.exists("/proc/1/cgroup")
  in_container_runtime <- FALSE
  if (cgroup_exists) {
    in_container_runtime <- any(
      grepl("docker", readLines("/proc/1/cgroup", warn = FALSE))
    )
  }
  return(dockerenv_exists || in_container_runtime)
}
