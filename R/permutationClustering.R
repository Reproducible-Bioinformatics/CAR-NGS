#' @title Permutations and Clustering
#' @description This function executes a ubuntu docker that produces a specific number of permutation to evaluate clustering.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param nPerm, number of permutations to perform the pValue to evaluate clustering
#' @param permAtTime, number of permutations that can be computes in parallel
#' @param percent, percentage of randomly selected cells removed in each permutation
#' @param range1, first number of cluster for k means algorithm
#' @param range2, last number of cluster for k means algorithm
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param logTen, 1 if the count matrix is already in log10, 0 otherwise
#' @param clustering, clustering method to use : "SIMLR" , "tSne", "griph"
#' @param perplexity, Number of close neighbors for each point. This parameter is specific for tSne. Default value is 10.Setting this parameter when use a clustering method different by tSne will be ignored.
#' @param seed, important value to reproduce the same results with same input
#' @param rK, 1 for rankGene algorithm 0 otherwise WARNING, very slow with this feature. This parameter is specific for SIMLR. Setting this parameter to 1 with other clustering methods will not give any different result to set the parameter to 0.

#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return VioPlot of silhouette cells value for each number of cluster used,clusterP file with clustering results for each permutation, killedCell file with removed cells in each permutation, clustering.output a sommarize file with general information for each cells.
#' @examples
#' \dontrun{
#'  permutationClustering("docker","/home/lucastormreig/CASC2.0/permutationClustering/scratch/","/home/lucastormreig/CASC2.0/permutationClustering/Data/TOTAL.csv",4,2,10,3,4,separator=",",logTen=0,clustering="SIMLR",perplexity=0)
#'}
#' @export
permutationClustering <- function(group=c("sudo","docker"), scratch.folderDOCKER, scratch.folderHOST, file, nPerm, permAtTime, percent, range1=3, range2=3, separator, logTen=0, clustering, perplexity=10 , seed=1111, rK=0){

  data.folder=dirname(file)
  positions=length(strsplit(basename(file),"\\.")[[1]])
  matrixNameC=strsplit(basename(file),"\\.")[[1]]
  matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
  format=strsplit(basename(basename(file)),"\\.")[[1]][positions]

  isDocker = is_running_in_docker()

    if (isDocker == TRUE){
      scratch.folderHOST = gsub("\\\\", "/", scratch.folderHOST)

      #creating a HOSTdocker variable for the datafolder
      host_parts = unlist(strsplit(scratch.folderHOST, "/"))
      docker_parts = unlist(strsplit(scratch.folderDOCKER, "/"))
      matches_path = paste(matches, collapse="/")
      HOSTpath = gsub(matches_path, "", scratch.folderHOST)

      #checking if the datafolder is inside a shared folder
      wd_parts = unlist(strsplit(data.folder, "/"))
      dmatches = intersect(docker_parts, wd_parts)
      dmatches_path = paste(dmatches, collapse="/")
      d_path = gsub(dmatches_path, "", data.folder)
      #creating the variable data.folderHOST
      data.folderHOST = paste(HOSTpath, d_path, sep="")
    }

    if(isDocker == FALSE){
      scratch.folderDOCKER = scratch.folderHOST
      data.folderHOST = data.folder
    }

  #running time 1
  ptm <- proc.time()
  #setting the data.folder as working folder
  if (!file.exists(data.folder)){
    cat(paste("\nIt seems that the ",data.folder, " folder does not exist\n"))
    return(2)
  }

  #storing the position of the home folder
  home <- getwd()
  setwd(data.folder)
  #initialize status
  exitStatus <- 0
  writeLines(as.character(exitStatus), "ExitStatusFile")

  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    exitStatus <- 10
    writeLines(as.character(exitStatus), "ExitStatusFile")
    setwd(home)
    return(10)
  }



  #check  if scratch folder exist
  if (!file.exists(scratch.folderDOCKER)){
    cat(paste("\nIt seems that the ",scratch.folderDOCKER, " folder does not exist\n"))
    exitStatus <- 3
    writeLines(as.character(exitStatus), "ExitStatusFile")
    setwd(data.folder)
    return(3)
  }

   tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folderHOST=file.path(scratch.folderHOST, tmp.folder)

  scrat_tmp.folderDOCKER=file.path(scratch.folderDOCKER, tmp.folder)
  writeLines(scrat_tmp.folderDOCKER,paste(data.folder,"/tempFolderID", sep=""))
  cat("\nCreating a folder in scratch folder\n")
  dir.create(file.path(scrat_tmp.folderDOCKER))

  #preprocess matrix and copying files

  if(separator=="\t"){
  separator="tab"
  }

  dir.create(paste(scrat_tmp.folderDOCKER,"/",matrixName,sep=""))
  dir.create(paste(data.folder,"/Results",sep=""))
  file.copy(paste(data.folder,"/",matrixName,".",format,sep=""), paste(scrat_tmp.folderDOCKER,"/",matrixName,sep=""))

  dockerID_name="dockerID"
  nr_dockerID = 0
  while (file.exists(dockerID_name)){
      nr_dockerID = nr_dockerID + 1
      dockerID_name =paste0("dockerID" ,"_", nr_dockerID)
  }


  #executing the docker job
  params <- paste("--cidfile ",data.folder,"/", dockerID_name," -v ",scrat_tmp.folderHOST,":/scratch -v ", data.folderHOST, ":/data -d docker.io/repbioinfo/permutationclustering Rscript /home/main.R ",matrixName," ",nPerm," ",permAtTime," ",percent," ",range1," ",range2," ",format," ",separator," ",logTen," ",clustering," ",rK," ",perplexity," ", seed, sep="")

  resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){
    #system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
  }

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/",dockerID_name, sep=""), warn = FALSE)
  file.copy(paste("docker logs ", substr(container.id,1,12), " &> "),paste(data.folder,"\\", substr(container.id,1,12),".log", sep=""))
  file.copy(paste("docker rm "),paste(container.id, sep=""))


  #Copy result folder
  cat("Copying Result Folder")
  items <- list.files(scrat_tmp.folderDOCKER)
          for (item in items) {
          file.copy(file.path(scrat_tmp.folderDOCKER, item), data.folder, recursive = TRUE, overwrite = TRUE)
        }

  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  unlink(scrat_tmp.folderDOCKER,recursive=TRUE)
  file.remove(c("tempFolderID"))
  setwd(home)
}
