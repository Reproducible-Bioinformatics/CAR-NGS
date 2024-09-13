#' @title Permutation Analysis
#' @description This function analyze the data that came up from permutationClustering script.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param range1, First number of cluster that has to be analyzed
#' @param range2, Last number of cluster that has to be analyzed
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param sp, minimun number of percentage of cells that has to be in common between two permutation to be the same cluster.
#' @param clusterPermErr, error that can be done by each permutation in cluster number depicting.Default = 0.05
#' @param maxDeltaConfidence, max value for Delta confidence for genes feature selection
#' @param minLogMean, min value for Log mean for genes feature selection
#' @author Luca Alessandri , alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return stability plot for each nCluster,two files with score information for each cell for each permutation.
#' @examples
#' \dontrun{
#' permAnalysis("docker", "path/to/scratch", "path/to/data/TOTAL", 3, 4, ",", 0.8) #
#' }
#' @export
permAnalysis <- function(group = c("sudo", "docker"), scratch.folderDOCKER, scratch.folderHOST, file, range1, range2, separator, sp, clusterPermErr = 0.05, maxDeltaConfidence = 0.01, minLogMean = 0.05) {
  data.folder <- dirname(file)
  positions <- length(strsplit(basename(file), "\\.")[[1]])
  matrixNameC <- strsplit(basename(file), "\\.")[[1]]
  matrixName <- paste(matrixNameC[seq(1, positions - 1)], collapse = "")
  format <- strsplit(basename(basename(file)), "\\.")[[1]][positions]

  # checking if the function is running from Docker
  isDocker <- is_running_in_docker()

  if (isDocker == TRUE) {
    scratch.folderHOST <- gsub("\\\\", "/", scratch.folderHOST)

    # creating a HOSTdocker variable for the datafolder
    host_parts <- unlist(strsplit(scratch.folderHOST, "/"))
    docker_parts <- unlist(strsplit(scratch.folderDOCKER, "/"))
    matches_path <- paste(matches, collapse = "/")
    HOSTpath <- gsub(matches_path, "", scratch.folderHOST)

    # checking if the datafolder is inside a shared folder
    wd_parts <- unlist(strsplit(data.folder, "/"))
    dmatches <- intersect(docker_parts, wd_parts)
    dmatches_path <- paste(dmatches, collapse = "/")
    d_path <- gsub(dmatches_path, "", data.folder)
    # creating the variable data.folderHOST
    data.folderHOST <- paste(HOSTpath, d_path, sep = "")
  }

  if (isDocker == FALSE) {
    scratch.folderDOCKER <- scratch.folderHOST
    data.folderHOST <- data.folder
  }

  # running time 1
  ptm <- proc.time()
  # setting the data.folder as working folder
  if (!file.exists(data.folder)) {
    cat(paste("\nIt seems that the ", data.folder, " folder does not exist\n"))
    return(2)
  }

  # storing the position of the home folder
  home <- getwd()
  setwd(data.folder)
  # initialize status
  exitStatus <- 0
  writeLines(as.character(exitStatus), "ExitStatusFile")

  # testing if docker is running
  test <- dockerTest()
  if (!test) {
    cat("\nERROR: Docker seems not to be installed in your system\n")
    exitStatus <- 10
    writeLines(as.character(exitStatus), "ExitStatusFile")
    setwd(home)
    return(10)
  }



  # check  if scratch folder exist
  if (!file.exists(scratch.folderDOCKER)) {
    cat(paste("\nIt seems that the ", scratch.folderDOCKER, " folder does not exist\n"))
    exitStatus <- 3
    writeLines(as.character(exitStatus), "ExitStatusFile")
    setwd(data.folder)
    return(3)
  }

  tmp.folder <- gsub(":", "-", gsub(" ", "-", date()))
  scrat_tmp.folderHOST <- file.path(scratch.folderHOST, tmp.folder)

  scrat_tmp.folderDOCKER <- file.path(scratch.folderDOCKER, tmp.folder)
  writeLines(scrat_tmp.folderDOCKER, paste(data.folder, "/tempFolderID", sep = ""))
  cat("\nCreating a folder in scratch folder\n")
  dir.create(file.path(scrat_tmp.folderDOCKER))

  # preprocess matrix and copying files
  if (separator == "\t") {
    separator <- "tab"
  }

  if (!file.exists(paste(data.folder, "/Results/", matrixName, "/", sep = ""))) {
    cat(paste("\nIt seems that some file are missing, check that your previously analysis results are still in the same folder,check Results folder!\n"))
    exitStatus <- 3
    writeLines(as.character(exitStatus), "ExitStatusFile")
    setwd(data.folder)
    return(3)
  }


  items1 <- list.files(paste(data.folder, "/Results"))
  for (item1 in items1) {
    file.copy(file.path(data.folder, "/Results/", item1), scrat_tmp.folderDOCKER, recursive = TRUE, overwrite = TRUE)
  }
  file.copy(paste(data.folder, "/", matrixName, ".", format, sep = ""), paste(scrat_tmp.folderDOCKER, sep = ""))

  dockerID_name <- "dockerID"
  nr_dockerID <- 0
  while (file.exists(dockerID_name)) {
    nr_dockerID <- nr_dockerID + 1
    dockerID_name <- paste0("dockerID", "_", nr_dockerID)
  }

  # executing the docker job
  params <- paste("--cidfile ", data.folder, "/", dockerID_name, " -v ", scrat_tmp.folderHOST, ":/scratch -v ", data.folderHOST, ":/data -d docker.io/repbioinfo/permutationanalysis Rscript /home/main.R ", matrixName, " ", range1, " ", range2, " ", format, " ", separator, " ", sp, " ", clusterPermErr, " ", maxDeltaConfidence, " ", minLogMean, sep = "")

  resultRun <- runDocker(group = group, params = params)

  # waiting for the end of the container work
  if (resultRun == 0) {
    #  system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
  }


  # saving log and removing docker container
  container.id <- readLines(paste(data.folder, "/", dockerID_name, sep = ""), warn = FALSE)
  file.copy(paste("docker logs ", substr(container.id, 1, 12), " &> "), paste(data.folder, "\\", substr(container.id, 1, 12), ".log", sep = ""))
  file.copy(paste("docker rm "), paste(container.id, sep = ""))



  # Copy result folder
  cat("Copying Result Folder")
  items <- list.files(scrat_tmp.folderDOCKER)
  for (item in items) {
    file.copy(file.path(scrat_tmp.folderDOCKER, item), data.folder, recursive = TRUE, overwrite = TRUE)
  }

  # removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  unlink(scrat_tmp.folderDOCKER, recursive = TRUE)
  file.remove(c("tempFolderID"))
  setwd(home)
}
