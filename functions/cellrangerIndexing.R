#' @title Cellranger indexing
#' @description This function creates the indexing for 10Xgenomics
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folderDOCKER, a character string indicating the path of the scratch folder inside the docker
#' @param scratch.folderHOST, a character string indicating the path of the scratch folder inside the host. If not running from docker, this is the character string that indicates the path of the scratch.folder
#' @param gtf.url, a character string indicating the URL from ENSEMBL ftp for the GTF for genome of interest
#' @param fasta.url, a character string indicating the URL from ENSEMBL ftp for the unmasked genome sequence of interest
#' @param genomeFolder,  path for the genome folder
#' @param bio.type, ENSEMBL biotype to filter the GTF
#' @param nThreads, number of cores for parallelization
#' @param version,  cellranger version: 2, 3 or 5.
#' @author Luca Alessandrì, Sebastian Bucatariu, Agata D'Onofrio
#'
#'
#' @return an indexed genome compliant with 10XGenomics cellranger
#' @examples
#' \dontrun{
#' library(rCASC)
#' setwd("/data/genomes/hg38refcellranger")
#'
#' cellrangerIndexing(group="docker", scratch.folderDOCKER="/data/scratch", scratch.folderHOST="/data/scratch",
#'             gtf.url="ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz",
#'             fasta.url="ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
#'             genomeFolder = getwd(), bio.type="protein_coding", nThreads = 8)
#' }
#'
#'
#' @export

cellrangerIndexing <- function(group=c("sudo","docker"),scratch.folderDOCKER, scratch.folderHOST,genomeFolder,gtf.url,fasta.url,bio.type=c("protein_coding","unitary_pseudogene",
                                                           "unprocessed_pseudogene","processed_pseudogene",
                                                           "transcribed_unprocessed_pseudogene","processed_transcript",
                                                           "antisense","transcribed_unitary_pseudogene",
                                                           "polymorphic_pseudogene","lincRNA",
                                                           "sense_intronic","transcribed_processed_pseudogene",
                                                           "sense_overlapping","IG_V_pseudogene",
                                                           "pseudogene","TR_V_gene",
                                                           "3prime_overlapping_ncRNA","IG_V_gene",
                                                           "bidirectional_promoter_lncRNA","snRNA",
                                                           "miRNA","misc_RNA",
                                                           "snoRNA","rRNA",
                                                           "IG_C_gene","IG_J_gene",
                                                           "TR_J_gene","TR_C_gene",
                                                           "TR_V_pseudogene","TR_J_pseudogene",
                                                           "IG_D_gene","ribozyme",
                                                           "IG_C_pseudogene","TR_D_gene",
                                                           "TEC","IG_J_pseudogene",
                                                           "scRNA","scaRNA",
                                                           "vaultRNA","sRNA",
                                                           "macro_lncRNA","non_coding","IG_pseudogene"),nThreads, version="5"){

  #checking if the function is running in docker
  isDocker <- is_running_in_docker()
    if (isDocker == TRUE){
    scratch.folderHOST <- gsub("\\\\", "/", scratch.folderHOST)
    }
    if(isDocker == FALSE){
   scratch.folderDOCKER = scratch.folderHOST
   }
    
  id="results_cellranger"

  #docker image
    if(version == "2"){
      dockerImage="docker.io/repbioinfo/cellranger"
    } else if(version == "3"){
      dockerImage="docker.io/repbioinfo/cellranger.2018.03"
    } else if(version == "5"){
      dockerImage="docker.io/repbioinfo/cellranger.2020.05"
    }



  #storing the position of the home folder
  home <- getwd()


  #running time 1
  ptm <- proc.time()

  #setting the data.folder as working folder
  if (!file.exists(genomeFolder)){
    cat(paste("\nIt seems that the ",genomeFolder, " folder does not exist\n"))
    return(2)
  }

  setwd(genomeFolder)


  #initialize status
  exitStatus <- 0
  writeLines(as.character(exitStatus), "ExitStatusFile")


  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    exitStatus <- 0
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
    setwd(home)
    return(3)
  }

  #creating a temporary folder
  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folderHOST=file.path(scratch.folderHOST, tmp.folder)
  
  scrat_tmp.folderDOCKER=file.path(scratch.folderDOCKER, tmp.folder)
  writeLines(scrat_tmp.folderDOCKER,paste(genomeFolder,"/tempFolderID", sep=""))
  cat("\nCreating a folder in scratch folder\n")
  dir.create(file.path(scrat_tmp.folderDOCKER))

  #checking if a dockerID file exist and creating another one if necessary
  dockerID_name="dockerID"
  nr_dockerID = 0
  while (file.exists(dockerID_name)){
      nr_dockerID = nr_dockerID + 1
      dockerID_name =paste0("dockerID" ,"_", nr_dockerID)
    } 
    
  #executing the docker job

  params <- paste("--cidfile ",genomeFolder,"/", dockerID_name," -v ", scrat_tmp.folderHOST, ":/data -d ",dockerImage, " /home/indexing.sh ",gtf.url," ",fasta.url," ",bio.type," ",nThreads,sep="")


  #Run docker
  resultRun <- runDocker(group=group, params=params, dockerID_name)

  #waiting for the end of the container work
     if(resultRun==0){
        items <- list.files(scrat_tmp.folderDOCKER)
          for (item in items) {
          file.copy(file.path(scrat_tmp.folderDOCKER, item), genomeFolder, recursive = TRUE, overwrite = TRUE)
        }
    }

  #saving log and removing docker container
  container.id <- readLines(paste(genomeFolder,"/dockerID", sep=""), warn = FALSE)
  file.copy(paste("docker logs ", substr(container.id,1,12), " &> "),paste(genomeFolder,"\\", substr(container.id,1,12),".log", sep=""))
  file.copy(paste("docker rm "),paste(container.id, sep=""))

  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  unlink(scrat_tmp.folderDOCKER,recursive=TRUE)
  file.remove(c("tempFolderID"))


  setwd(home)

}
