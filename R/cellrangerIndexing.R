#' Cellranger indexing
#'
#' This function creates the indexing for 10Xgenomics
#'
#' @param result_dir_path, a character string indicating the path where the result files are saved
#' @param gtf.url, a character string indicating the URL from ENSEMBL ftp for the GTF for genome of interest
#' @param fasta.url, a character string indicating the URL from ENSEMBL ftp for the unmasked genome sequence of interest
#' @param bio.type, ENSEMBL biotype to filter the GTF
#' @param nThreads, number of cores for parallelization
#' @author Luca Alessandrì,  Sebastian Bucatariu, Agata D'Onofrio
#'
#'
#' @examples
#' \dontrun{
#' cellrangerIndexing(
#'   result_dir_path = "/the/result/dir",
#'   gtf.url = "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz",
#'   fasta.url = "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz",
#'   bio.type = "protein_coding",
#'   nThreads = 8
#' )
#' }
#'
#' @export

cellrangerIndexing <- function(result_dir_path, gtf.url, fasta.url, bio.type = c(
                                 "protein_coding", "unitary_pseudogene",
                                 "unprocessed_pseudogene", "processed_pseudogene",
                                 "transcribed_unprocessed_pseudogene", "processed_transcript",
                                 "antisense", "transcribed_unitary_pseudogene",
                                 "polymorphic_pseudogene", "lincRNA",
                                 "sense_intronic", "transcribed_processed_pseudogene",
                                 "sense_overlapping", "IG_V_pseudogene",
                                 "pseudogene", "TR_V_gene",
                                 "3prime_overlapping_ncRNA", "IG_V_gene",
                                 "bidirectional_promoter_lncRNA", "snRNA",
                                 "miRNA", "misc_RNA",
                                 "snoRNA", "rRNA",
                                 "IG_C_gene", "IG_J_gene",
                                 "TR_J_gene", "TR_C_gene",
                                 "TR_V_pseudogene", "TR_J_pseudogene",
                                 "IG_D_gene", "ribozyme",
                                 "IG_C_pseudogene", "TR_D_gene",
                                 "TEC", "IG_J_pseudogene",
                                 "scRNA", "scaRNA",
                                 "vaultRNA", "sRNA",
                                 "macro_lncRNA", "non_coding", "IG_pseudogene"
                               ), nThreads, version = "5") {
  id <- "results_cellranger"

  # docker image
  if (version == "2") {
    dockerImage <- "docker.io/repbioinfo/cellranger"
  } else if (version == "3") {
    dockerImage <- "docker.io/repbioinfo/cellranger.2018.03"
  } else if (version == "5") {
    dockerImage <- "docker.io/repbioinfo/cellranger.2020.05"
  }

  # Type checking.
  if (typeof(result_dir_path) != "character") {
    stop(paste("result_dir_path type is", paste0(typeof(result_dir_path), "."), "It should be \"character\""))
  }
  if (typeof(gtf.url) != "character") {
    stop(paste("gtf.url type is", paste0(typeof(gtf.url), "."), "It should be \"character\""))
  }
  if (typeof(fasta.url) != "character") {
    stop(paste("fasta.url type is", paste0(typeof(fasta.url), "."), "It should be \"character\""))
  }
  if (typeof(bio.type) != "character") {
    stop(paste("bio.type type is", paste0(typeof(bio.type), "."), "It should be \"character\""))
  }
  if (!is.numeric(nThreads)) {
    stop(paste("nThreads type is", paste0(typeof(nThreads), "."), "It should be \"double\" or \"integer\""))
  }
  if (typeof(version) != "character") {
    stop(paste("version type is", paste0(typeof(version), "."), "It should be \"character\""))
  }

  # Check if result_dir_path exists
  if (!rrundocker::is_running_in_docker()) {
    if (!dir.exists(result_dir_path)) {
      stop(paste("result_dir_path:", result_dir_path, "does not exist."))
    }
  }

  # Executing the docker job
  rrundocker::run_in_docker(
    image_name = paste0(dockerImage, ":latest"),
    volumes = list(
      c(result_dir_path, "/data")
    ),
    additional_arguments = c(
      "/home/indexing.sh",
      gtf.url,
      fasta.url,
      bio.type,
      nThreads
    )
  )
}
