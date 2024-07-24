dimensions <- function(host_folder, input_file, separator) {
  # Type checking.
  if (typeof(host_folder) != "character") {
    stop(
      paste(
        "host_folder type is",
        paste0(typeof(host_folder), "."),
        "It should be \"character\""
      )
    )
  }
  if (typeof(input_file) != "character") {
    stop(
      paste(
        "input_file type is",
        paste0(typeof(input_file), "."),
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

  if (!file.exists(input_file)) {
    stop(paste("input_file:", input_file, "does not exist."))
  }
  if (!dir.exists(host_folder)) {
    stop(paste("host_folder:", host_folder, "does not exist."))
  }

  # Obtain parent folder from input_file.
  parent_folder <- dirname(input_file)

  # Obtain matrix name and file format from input_file.
  input_file_parts <- strsplit(basename(input_file), "\\.")[[1]]
  matrix_name <- input_file_parts[1]
  format <- input_file_parts[2] # Uses extension as an heuristic.

  run_in_docker(
    image_name = "docker.io/repbioinfo/r332.2017.01:latest",
    volumes = c(
      paste0(host_folder, ":/data"),
      paste0(parent_folder, ":/scratch")
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
                          volumes = c(),
                          additional_arguments = c()) {
  base_command <- "run --privileged=true --platform linux/amd64 --rm"
  for (volume in volumes) {
    base_command <- paste(base_command, "-v", volume)
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
