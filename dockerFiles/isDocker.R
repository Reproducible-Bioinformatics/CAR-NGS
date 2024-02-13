# Defining is_running_in_docker function
is_running_in_docker <- function() {
  # Check if running in a Docker container
  if (file.exists("/.dockerenv")) {
    return(TRUE)
  }

  if (file.exists("/proc/1/cgroup")) {
    content <- readLines("/proc/1/cgroup", warn = FALSE)
    return(any(grepl("docker", content)))
  }

  return(FALSE)
}
