permAnalysisSeurat <-
function(group=c("sudo","docker"), scratch.folderDOCKER, scratch.folderHOST, file,nCluster,separator,sp=0.8,sparse=FALSE,format="NULL"){

if(!sparse){
  data.folder=dirname(file)
positions=length(strsplit(basename(file),"\\.")[[1]])
matrixNameC=strsplit(basename(file),"\\.")[[1]]
matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
format=strsplit(basename(basename(file)),"\\.")[[1]][positions]
}else{
  matrixName=strsplit(dirname(file),"/")[[1]][length(strsplit(dirname(file),"/")[[1]])]
  data.folder=paste(strsplit(dirname(file),"/")[[1]][-length(strsplit(dirname(file),"/")[[1]])],collapse="/")
  if(format=="NULL"){
  stop("Format output cannot be NULL for sparse matrix")
  }
}
isDocker = is_running_in_docker()

    if (isDocker == TRUE){
      scratch.folderHOST = gsub("\\", "/", scratch.folderHOST)

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
    exitStatus <- 2
    writeLines(as.character(exitStatus), "ExitStatusFile")
    setwd(data.folder)
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
  if (!file.exists(scratch.folderHOST)){
    cat(paste("\nIt seems that the ",scratch.folder, " folder does not exist\n"))
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
 dockerID_name="dockerID"
  nr_dockerID = 0
  while (file.exists(dockerID_name)){
      nr_dockerID = nr_dockerID + 1
      dockerID_name =paste0("dockerID" ,"_", nr_dockerID)
  }

if(separator=="\t"){
separator="tab"
}
    if (!file.exists(paste(data.folder,"/Results/",matrixName,"/",sep=""))){
    cat(paste("\nIt seems that some file are missing, check that your previously analysis results are still in the same folder,check Results folder!\n"))
    exitStatus <- 3
    writeLines(as.character(exitStatus), "ExitStatusFile")
    setwd(data.folder)
    return(3)
  }
 file.copy(paste(data.folder,"/",matrixName,".",format,sep=""), paste(scrat_tmp.folderDOCKER,"/",sep=""))
 file.copy(paste(cl,sep=""), paste(scrat_tmp.folderDOCKER,"/",sep=""))


  #executing the docker job
    params <- paste("--cidfile ",data.folder,"/", dockerID_name," -v ",scrat_tmp.folderHOST,":/scratch -v ", data.folderHOST, ":/data -d docker.io/repbioinfo/seuratanalysis Rscript /home/main.R ",matrixName," ",nCluster," ",format," ",separator," ",sp," ",sparse,sep="")

resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){
  #  system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
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
  #system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  file.remove(c("tempFolderID"))
  setwd(home)
}
