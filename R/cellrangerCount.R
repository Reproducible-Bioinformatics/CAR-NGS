cellrangerCount <-
function(group=c("sudo","docker"),  transcriptome.folder,  fastq.folder,  sample=NULL, expect.cells=NULL, force.cells=NULL, nosecondary=TRUE, chemistry="auto", r1.length=NULL,  r2.length=NULL, lanes=NULL, localcores=NULL, localmem=NULL,  scratch.folderDOCKER, scratch.folderHOST, version="5"){

#checking if the function is running from Docker
isDocker <- is_running_in_docker()
    
if (isDocker == TRUE){
    scratch.folderHOST <- gsub("\\\\", "/", scratch.folderHOST)
    #creating a HOSTdocker variable for the transcriptome.folder
    host_parts = unlist(strsplit(scratch.folderHOST, "/"))
    docker_parts = unlist(strsplit(scratch.folderDOCKER, "/"))
    matches = intersect(host_parts, docker_parts)
    matches_path = paste(matches, collapse="/")
    HOSTpath = gsub(matches_path, "", scratch.folderHOST)
    #checking if the trascriptome.folderHOST is inside a shared folder
    tr_parts = unlist(strsplit(transcriptome.folder, "/"))
    tmatches = intersect(docker_parts, tr_parts)
    tmatches2 = intersect(host_parts, tr_parts)
    tmatches_path = paste(tmatches, collapse="/")
    tmatches_path2 = paste(tmatches2, collapse="/")
    t_path = gsub(tmatches_path, "", transcriptome.folder)
    #creating the variable trascriptome.folderHOST
    transcriptome.folderHOST = paste(HOSTpath, tmatches_path2, t_path, sep="")
}

if (isDocker == FALSE){
    scratch.folderDOCKER = scratch.folderHOST
    transcriptome.folderHOST = transcriptome.folder
}



  id="results_cellranger"
  #docker image
  if(version == "2"){
    dockerImage="docker.io/repbioinfo/cellranger"
  } else if(version == "3"){
    dockerImage="docker.io/repbioinfo/cellranger.2018.03"
  } else if(version == "5"){
    dockerImage="docker.io/repbioinfo/cellranger.2020.05"
  } else if(version == "7"){
    dockerImage="docker.io/repbioinfo/cellranger.2023.7.1.0"
  }
  

#storing the position of the home folder
  home <- getwd()


  #running time 1
  ptm <- proc.time()

  #setting the data.folder as working folder
  if (!file.exists(fastq.folder)){
    cat(paste("\nIt seems that the ",fastq.folder, " folder does not exist\n"))
    return(2)
  }
  setwd(fastq.folder)

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
      setwd(home)
    return(3)
  }

  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folderHOST=file.path(scratch.folderHOST, tmp.folder)

  scrat_tmp.folderDOCKER=file.path(scratch.folderDOCKER, tmp.folder)
  writeLines(scrat_tmp.folderDOCKER,paste(fastq.folder,"/tempFolderID", sep=""))
  cat("\nCreating a folder in scratch folder\n")
  dir.create(file.path(scrat_tmp.folderDOCKER))

  #cp fastq folder in the scrat_tmp.folder
  file.copy(list.files(fastq.folder, pattern = "*.gz", full.names = TRUE), scrat_tmp.folderDOCKER)

  #checking if a dockerID file exist and creating another one if necessary
  dockerID_name="dockerID"
  nr_dockerID = 0
  while (file.exists(dockerID_name)){
    nr_dockerID = nr_dockerID + 1
    dockerID_name =paste0("dockerID" ,"_", nr_dockerID)
    }

  #executing the docker job
  params <- paste("--cidfile ",fastq.folder,"/", dockerID_name," -v ",transcriptome.folderHOST,":/transcr -v ", scrat_tmp.folderHOST, ":/data -d ",dockerImage, " /bin/cellranger count  --id=",id," --transcriptome=/transcr --fastqs=/data", sep="")

  if(!is.null(sample)){
    params<-paste(params," --sample=",sample, sep="")

  }

  if (!is.null(expect.cells)){
   params<-paste(params," --expect-cells=",expect.cells, sep="")
  }

  if (!is.null(force.cells)){
   params<-paste(params," --force-cells=",force.cells, sep="")
  }

  if (nosecondary){
   params<-paste(params," --nosecondary", sep="")
  }

  if (!is.null(chemistry)){
   params<-paste(params," --chemistry=",chemistry, sep="")
  }

  if (!is.null(r1.length)){
   params<-paste(params," --r1-length=",r1.length, sep="")
  }

  if (!is.null(r2.length)){
   params<-paste(params," --r2-length=",r2.length, sep="")
  }

  if (!is.null(lanes)){
   params<-paste(params," --lanes=",lanes, sep="")
  }

  if (!is.null(localcores)){
   params<-paste(params," --localcores=",localcores, sep="")
  }

  if (!is.null(localmem)){
   params<-paste(params," --localmem=", localmem, sep="")
  }

  params.split <- unlist(strsplit(params, dockerImage))
  params0 <- paste(params.split[1], " ", dockerImage, " /bin/bash /data/script.sh", sep="")
  cat(params0,"\n")
  params1 <- NULL
  params1[1] <- "cd /data"
  params1[2] <- params.split[2]
  params1[3] <- paste("chmod 777 -R /data/", id, sep="")
  if(version=="2"){
    params1[4] <- paste("/bin/cellranger mat2csv /data/", id,"/outs/filtered_gene_bc_matrices ",id,".csv", sep="")
  }else if(version=="3"){
    params1[4] <- paste("/bin/cellranger mat2csv /data/", id,"/outs/filtered_feature_bc_matrix ",id,".csv", sep="")
  }else if(version=="5"){
    params1[4] <- paste("/bin/cellranger mat2csv /data/", id,"/outs/filtered_feature_bc_matrix ",id, ".csv", sep="")
  }else if(version=="7"){
    params1[4] <- paste("/bin/cellranger mat2csv /data/", id,"/outs/filtered_feature_bc_matrix ",id, ".csv", sep="")
  }
  



  fileConn<-file(paste(scrat_tmp.folderDOCKER,"/script.sh", sep=""), "w")
  writeLines(params1, fileConn)
  close(fileConn)
  #system(paste("chmod +x ", scrat_tmp.folder,"/script.sh", sep=""))
  script_path <- file.path(scrat_tmp.folderDOCKER, "script.sh")
  Sys.chmod(script_path, mode = "0777", use_umask = FALSE)

  #Run docker
  resultRun <- runDocker(group=group, params=params0, dockerID_name)
  
  #waiting for the end of the container work
  if(resultRun==0){
    file.copy(paste(scrat_tmp.folderDOCKER, "/", id, sep=""), paste(fastq.folder, sep=""), recursive = TRUE)
    file.copy(paste(scrat_tmp.folderDOCKER, "/results_cellranger.csv", sep=""), paste(fastq.folder, sep=""))
    #system(paste("sed \'s|,|\t|g\' ",fastq.folder,"/",id,".csv > ", fastq.folder,"/",id,".txt", sep=""))
    csv_file <- file.path(fastq.folder, paste0(id, ".csv"))
    data <- read.csv(csv_file)
    txt_file <- file.path(fastq.folder, paste0(id, ".txt"))
    write.table(data, txt_file, sep = "\t", row.names = FALSE, col.names = TRUE)

    cat("\nCellranger analysis is finished\n")
  }

  #saving log and removing docker container
  container.id <- readLines(paste(fastq.folder,"/",dockerID_name, sep=""), warn = FALSE)
  file.copy(paste("docker logs ", substr(container.id,1,12), " &> "),paste(fastq.folder,"\\", substr(container.id,1,12),".log", sep=""))
  file.copy(paste("docker rm "),paste(container.id, sep=""))

  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  unlink(scrat_tmp.folderDOCKER,recursive=TRUE)
  file.remove(c("tempFolderID"))

  setwd(home)

}
