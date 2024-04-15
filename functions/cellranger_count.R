#' @title Cellranger count
#' @description This function takes FASTQ files from cellranger mkfastq and performs alignment, filtering, barcode counting, and UMI counting.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param scratch.folderDOCKER, a character string indicating the path of the scratch folder inside the docker
#' @param scratch.folderHOST, a character string indicating the path of the scratch folder inside the host. If not running from docker, this is the character string that indicates the path of the scratch.folder
#' @param transcriptome.folder,  path to the Cell Ranger compatible transcriptome reference e.g. for a human and mouse mixture sample, use refdata-cellranger-hg19-and-mm10-1.2.0
#' @param fastq.folder,  path of the fastq path folder in fastq folder the  fastq must have the format SAMPLENAME_S1_L001_R1_001.fastq.gz
#' @param sample, fastq name, if fastq name is subject1_S1_L001_R1_001.fastq.gz sample is subject1
#' @param expect.cells,  optional setting the number of recovered cells. Default: 3000 cells.
#' @param force.cells,  optional to force pipeline to use this number of cells, bypassing the cell detection algorithm. Use this if the number of cells estimated by Cell Ranger is not consistent with the barcode rank plot.
#' @param nosecondary,  optional flag, default TRUE, to skip secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Set this if you plan to use cellranger reanalyze or your own custom analysis.
#' @param chemistry,  optional assay configuration. One of: auto for autodetection (default), threeprime for Single Cell 3end, fiveprime for Single Cell 5end, SC3Pv1 for Single Cell 3end v1, SC3Pv2 for Single Cell 3end v2, SC5P-PE for Single Cell 5end paired-end (both R1 and R2 are used for alignment), SC5P-R2 for Single Cell 5end R2-only (where only R2 is used for alignment).
#' @param r1.length,  optional hard-trim the input R1 sequence to this length. Note that the length includes the Barcode and UMI sequences so do not set this below 26 for Single Cell 3end v2 or Single Cell 5end. This and --r2-length are useful for determining the optimal read length for sequencing.
#' @param r2.length,  optional hard-trim the input R2 sequence to this length.
#' @param lanes,  optional, lanes associated with this sample
#' @param localcores,  restricts cellranger to use specified number of cores to execute pipeline stages. By default, cellranger will use all of the cores available on your system.
#' @param localmem,  restricts cellranger to use specified amount of memory, in GB, to execute pipeline stages. By default, cellranger will use 90\% of the memory available on your system. Please note that cellranger requires at least 16 GB of memory to run all pipeline stages.
#' @param version,  cellranger version: 2, 3, 5 or 7. 
#' @author Greta Romano, romano [dot] greta [at] gmail [dot] com, University of Torino
#'
#'
#' @return a folder called results_cellranger, more info on the structure of this folder at https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview . In /somewhewre_in_your_computer/results_cellranger/outs/filtered_gene_bc_matrices the cells counts matrices results_cellranger.cvs and results_cellranger.txt are saved for further use.
#'
#' @examples
#' \dontrun{
#' home <- getwd()
#' library(rCASC)
#' setwd("/data/genomes/cellranger_hg19mm10")
#' #getting the human and mouse cellranger index
#' system("wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-hg19-and-mm10-2.1.0.tar.gz")
#' untar("refdata-cellranger-hg19-and-mm10-2.1.0.tar.gz") 
#' setwd(home)
#' # 100 cells 1:1 Mixture of Fresh Frozen Human (HEK293T) and Mouse (NIH3T3) Cells
#' system("wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/hgmm_100/hgmm_100_fastqs.tar")
#' untar("hgmm_100_fastqs.tar")
#' home=paste(home,"/fastqs",sep="")
#' cellrangerCount(group="docker",  transcriptome.folder="/data/genomes/cellranger_hg19mm10/refdata-cellranger-hg19_and_mm10-2.1.0",  fastq.folder=getwd(),  expect.cells=100, nosecondary=TRUE, scratch.folder="/data/scratch", version="2")
#' 
#' sraDownload(group = "docker", sra.name = "SRR7762358", data.folder = getwd(), scratch.folder = "/data/scratch", threads = 8)
#' system("mv ./SRR7762358/SRR7762358.fastq.gz ./SRR7762358/SRR7762358_S1_L001_R1_001.fastq.gz")
#' cellrangerCount(group="docker",  transcriptome.folder="/data/genomes/refdata-cellranger-GRCh38-3.0.0",  fastq.folder=getwd(), sample="SRR7762358",  nosecondary=TRUE, scratch.folder="/data/scratch", version="3")
#'
#' }
#'
#'
#' @export

cellrangerCount <- function(group=c("sudo","docker"),  transcriptome.folder,  fastq.folder,  sample=NULL, expect.cells=NULL, force.cells=NULL, nosecondary=TRUE, chemistry="auto", r1.length=NULL,  r2.length=NULL, lanes=NULL, localcores=NULL, localmem=NULL,  scratch.folderDOCKER, scratch.folderHOST, version="5"){

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
    exitStatus <- 0
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
  }
    #else if(version=="5"){
    #params1[4] <- paste("/bin/cellranger mat2csv /data/outs/filtered_feature_bc_matrix.csv", sep="")
  #}
  



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
    file.copy(paste(scrat_tmp.folderDOCKER, "/", id), paste(fastq.folder, sep=""), recursive = TRUE)
    file.copy(paste(scrat_tmp.folderDOCKER, "/results_cellranger.csv"), paste(fastq.folder, sep=""))
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
