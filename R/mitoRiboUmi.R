#' @title MitoRiboUmi
#' @description This function executes a ubuntu docker showing ribosomal mitochondrial protein genes fraction in each cell
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folder, a character string indicating the path of the scratch folder
#' @param file, a character string indicating the path of the file, with file name and extension included
#' @param separator, separator used in count file, e.g. '\\t', ','
#' @param gtf.name, name for the gtf file tu be used.
#' @param bio.type, ENSEMBL biotype of interest, default for single-cell protein_coding
#' @param umiXgene, a integer defining how many UMI are required to call a gene as present. default: 3
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return a pdf file called Ribo_mito.pdf
#' @examples
#'\dontrun{
#'    system("wget ftp://ftp.ensembl.org/pub/release-94/gtf/mus_musculus/Mus_musculus.GRCm38.94.gtf.gz")
#'    system("gzip -d Mus_musculus.GRCm38.94.gtf.gz")
#'    system("wget http://130.192.119.59/public/testSCumi_mm10.csv.zip")
#'    unzip("testSCumi_mm10.csv.zip")
#'    mitoRiboUmi(group="docker", file=paste(getwd(), "testSCumi_mm10.csv", sep="/"), 
#'                 scratch.folder="/data/scratch", separator=",", umiXgene=3, 
#'                 gtf.name="Mus_musculus.GRCm38.94.gtf", bio.type="protein_coding")

#'}
#' @export
mitoRiboUmi <- function(group=c("sudo","docker"),scratch.folderDOCKER, scratch.folderHOST, file,separator, gtf.name, bio.type, umiXgene){


  data.folder=dirname(file)
positions=length(strsplit(basename(file),"\\.")[[1]])
matrixNameC=strsplit(basename(file),"\\.")[[1]]
matrixName=paste(matrixNameC[seq(1,positions-1)],collapse="")
format=strsplit(basename(basename(file)),"\\.")[[1]][positions]

 #checking if the function is running from Docker
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
  

  items1 <- list.files(paste(data.folder))
  for (item1 in items1) {
      file.copy(file.path(data.folder, item1), scrat_tmp.folderDOCKER, recursive = TRUE, overwrite = TRUE)
      }

  #executing the docker job
 params <- paste("--cidfile ",data.folder,"/", dockerID_name," -v ",scrat_tmp.folderHOST,":/scratch -v ", data.folderHOST, ":/data -d docker.io/repbioinfo/ribomitoumi Rscript /home/main.R ",matrixName,"  ",format," ",separator," ",gtf.name," ",bio.type," ",umiXgene,sep="")

resultRun <- runDocker(group=group, params=params)

  #waiting for the end of the container work
  if(resultRun==0){
  #  system(paste("cp ", scrat_tmp.folder, "/* ", data.folder, sep=""))
  }

  items <- list.files(paste(data.folder,"/Result"))
          for (item in items) {
          file.copy(file.path(paste(data.folder,"/Result"), item), scrat_tmp.folderDOCKER, recursive = TRUE, overwrite = TRUE)
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
