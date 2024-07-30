#' @title Shifted logarithm transformation
#' @description This function executes a ubuntu docker that performs shifted logarithm transformation (Ahlmann-Eltze Nature Methods, vol. 20, pg. 665, 2023; transformGamPoi Bioconductor package)
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param file, a character string indicating the path of the count matrix
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param overdispersion, the simplest count model is the Poisson model. However, the Poisson model assumes that variance = mean. For many applications this is too rigid and the Gamma-Poisson allows a more flexible mean-variance relation (variance = mean + mean^2 * overdispersion).
#' @param pseudo_count, it is set to 4, which is use internally in pseudo-count = 1/( pseudo_count * overdispersion)).
#' @param size_factors,TRUE or False in large scale experiments, each sample is typically of different size (for example different sequencing depths). A size factor is an internal mechanism of GLMs to correct for this effect. normed_sum and poscounts are fairly simple methods and can lead to suboptimal results. For the best performance, I recommend to use deconvolution.
#' @param minimum_overdispersion, the 'acosh_transform' converges against 2 * sqrt(x) for 'overdispersion == 0'. However, the 'shifted_log_transform' would just become '0', thus here we apply the 'minimum_overdispersion' to avoid this behavior.
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param outputFolder, a character string indicating the path of the output folder
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return plot
#' @examples
#'\dontrun{
#'dir.create("scratch")
#'logTransform(group="docker",file="/home/gianluca/newFunction/setA.csv",separator=",",scratch.folder="/home/gianluca/newFunction/scratch/","/home/gianluca/newFunction/res/")
#' }
#' @export
logTransform <- function(group=c("sudo","docker"),file,separator, overdispersion=0.05,pseudo_count=4,size_factors=TRUE,minimum_overdispersion=0.001, scratch.folderDOCKER,scatch.folderHOST,outputFolder){

  data.folder=dirname(file)
  positions1=length(strsplit(basename(file),"\\.")[[1]])
  matrixNameC1=strsplit(basename(file),"\\.")[[1]]
  matrixName1=paste(matrixNameC1[seq(1,positions1-1)],collapse="")
  format1=strsplit(basename(basename(file)),"\\.")[[1]][positions1]
  data.folder=data.folder
    
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
  cat("\ncreating a folder in scratch folder\n")
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

  file.copy(paste(file," "),paste(scrat_tmp.folderDOCKER,"/",sep=""))
  dir.create(outputFolder)
  #executing the docker job
  params <- paste("--cidfile ",data.folder,"/", dockerID_name," -v ",scrat_tmp.folderHOST,":/scratch -v ", data.folder, ":/data -d repbioinfo/transformgampoi Rscript /home/main.R ",separator," ",format1," ",matrixName1," ",overdispersion," ",pseudo_count," ",minimum_overdispersion," ",size_factors,sep="")

  resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){
      items <- list.files(scrat_tmp.folderDOCKER)
          for (item in items) {
          file.copy(file.path(scrat_tmp.folderDOCKER, item),data.folder, recursive = TRUE, overwrite = TRUE)
        }
    #system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
  }
 

  #saving log and removing docker container
  container.id <- readLines(paste(data.folder,"/",dockerID_name, sep=""), warn = FALSE)
  file.copy(paste("docker logs ", substr(container.id,1,12), " >& "),paste(data.folder,"/", substr(container.id,1,12),".log", sep=""))
  file.copy(paste("docker rm "),paste(container.id, sep=""))


  #Copy result folder
  cat("Copying Result Folder")
  file.copy(paste(scrat_tmp.folder,"/* "),paste(outputFolder,sep=""))
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  file.remove(c("tempFolderID"))
  setwd(home)
}
