#' @title Marker genes discovery with COMETSC
#' @description This function executes a ubuntu docker for cometsc (https://github.com/MSingerLab/COMETSC)
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param threads, integer refering to the max number of process run in parallel default 1 max the number of clusters under analysis, i.e. nCluster
#' @param X, from 0 to 1 argument for XL-mHG default 0.15, for more info see cometsc help.
#' @param K, the number of gene combinations to be considered., possible values 2, 3, 4, default 2. WARNING increasing the number of combinations makes the matrices very big
#' @param counts, if set to True it will graph the log(expression+1). To be used if unlogged data are provided
#' @param skipvis, set to True to skip visualizations
#' @param nCluster, number of interested cluster used for analysis
#' @param scratch.folder, temporary folder where calculation is made
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @return folders with prefix output. More info in output at https://hgmd.readthedocs.io/en/latest/Output.html
#' @author Raffaele Calogero,raffaele.calogero [at] unito [dot] it, University of Torino
#' 
#' @examples
#' \dontrun{
#'     #running cometsc
#'     cometsc(group="docker", file="/Users/raffaelecalogero/Desktop/AXLN1/data/topx_veanno.csv", 
#'            scratch.folder="/Users/raffaelecalogero/Desktop",
#'            threads=1, counts="True", skipvis="False", nCluster=8, separator=",") 
#' }
#'
#' @export
cometsc <- function(group=c("sudo","docker"), file, scratch.folderDOKER, scratch.folderHOST, threads=1,  X=0.15, K=2, counts=c("True", "False"), skipvis=c("True", "False"), nCluster, separator){

    id="results_cellranger"


  data.folder=dirname(file)
  positions=length(strsplit(basename(file),"\\.")[[1]])
  matrixNameC=strsplit(basename(file),"\\.")[[1]]
  matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
  format=strsplit(basename(basename(file)),"\\.")[[1]][positions]

  isDocker <- is_running_in_docker()

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
  #system("echo 0 > ExitStatusFile 2>&1")
  exitStatus <- 0
  writeLines(as.character(exitStatus), "ExitStatusFile")
  
  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    #system("echo 10 > ExitStatusFile 2>&1")
    exitStatus <- 10
    writeLines(as.character(exitStatus), "ExitStatusFile")
    setwd(home)
    return(10)
  }
  

  
  #check  if scratch folder exist
  if (!file.exists(scratch.folderDOCKER)){
    cat(paste("\nIt seems that the ",scratch.folderDOCKER, " folder does not exist\n"))
    #system("echo 3 > ExitStatusFile 2>&1")
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

  #system(paste("cp -r ",data.folder,"/Results/", matrixName,"/",nCluster, "/* ",scrat_tmp.folder,sep=""))
  file.copy(paste(data.folder,"/Results/", matrixName), paste (scrat_tmp.folder,sep=""))
  file.copy(paste(file), paste(scrat_tmp.folder,sep=""))
  
  dockerID_name="dockerID"
  nr_dockerID = 0
  while (file.exists(dockerID_name)){
      nr_dockerID = nr_dockerID + 1
      dockerID_name =paste0("dockerID" ,"_", nr_dockerID)
  }

  #executing the docker
  params <- paste("--cidfile ",data.folder,"/", dockerID_name," -v ",scrat_tmp.folderHOST,":/scratch -v ", data.folderHOST, ":/data -d docker.io/repbioinfo/cometsc.2020.01 bash /bin/cometsc.sh ", matrixName, " ", threads, " ", X, " ", K, " ", counts, " ", skipvis, " ", nCluster," ", separator, sep="")

  resultRun <- runDocker(group=group, params=params)
  
  #waiting for the end of the container work
  if(resultRun==0){
    source_path <- file.path(scrat_tmp.folder, matrixName, nCluster)
    files_to_copy <- list.files(source_path, pattern = "output.*")
    file.copy(paste(file.path(source_path, files_to_copy)), paste(data.folder,"/Results/",matrixName,"/",nCluster,"/", sep=""), recursive = TRUE)
  }

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/",dockerID_name, sep=""), warn = FALSE)
  file.copy(paste("docker logs ", substr(container.id,1,12), " &> "),paste(data.folder,"\\", substr(container.id,1,12),".log", sep=""))
  file.copy(paste("docker rm "),paste(container.id, sep=""))

  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  unlink(scrat_tmp.folderDOCKER,recursive=TRUE)
  file.remove(c("tempFolderID"))
  #system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
}
