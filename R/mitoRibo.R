#' Mitochondrial and Ribosomal Filtering Script
#'
#' This script is responsible for filtering cells based on their mitochondrial and ribosomal gene content, providing quality control in single-cell RNA-seq experiments,
#' It returns a filtered matrix excluding mitochondrial and ribosomal genes inside the parent_folder and a plot of mitochondrial vs. ribosomal gene expressionas a PNG file
#'
#' @param input_file_path, a character string indicating the path of the count matrix file, which can be both dense (.csv/.txt) or sparse (.mtx)
#' @param mitoMin, start range for mitochondrial percentage, cells within the range are retained
#' @param mitoMax, end range for mitochondrial percentage, cells within the range are retained
#' @param riboMin, start range for ribosomal percentage, cells within the range are retained
#' @param riboMax, end range for ribosomal percentage, cells within the range are retained
#' @param separator, separator used in the count table
#' @param genes_file, a character string indicating the path of the genes name files necessary for the analysis of a sparse matrix ("*genes.tsv")
#' @param barcodes_file, a character string indicating the path of the barcodes file necessary for the analysis of a sparse matrix ("*barcodes.tsv")
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @examples
#'\dontrun{
#' # Dense matrix analysis
#'   mitoRibo(
#'     input_file_path = "/the/input/file.csv",
#'     mitoMin=0,
#'     mitoMax=10,
#'     riboMin=0,
#'     riboMax=20,
#'     separator=","
#'   )
#'
#' # Dense matrix analysis
#'   mitoRibo(
#'     input_file_path = "/the/input/file.mtx",
#'     mitoMin=0,
#'     mitoMax=10,
#'     riboMin=0,
#'     riboMax=20,
#'     genes_file="combined_filtered_with_sample_genes.tsv",
#'     barcodes_file="combined_filtered_with_sample_barcodes.tsv"
#'  )
#'}
#' @export
mitoRibo <- function(input_file_path, mitoMin, mitoMax, riboMin, riboMax, separator=NULL, genes_file=NULL, barcodes_file=NULL){

  # Type checking.
  if (typeof(input_file_path) != "character") {
    stop(paste("input_file_path type is", paste0(typeof(input_file_path), "."), "It should be \"character\""))
  }
  if (!is.numeric(mitoMin)) {
    stop(paste("mitoMin type is", paste0(typeof(mitoMin), "."), "It should be \"double\" or \"integer\""))
  }
  if (!is.numeric(mitoMax)) {
    stop(paste("mitoMax type is", paste0(typeof(mitoMax), "."), "It should be \"double\" or \"integer\""))
  }
  if (!is.numeric(riboMin)) {
    stop(paste("riboMin type is", paste0(typeof(riboMin), "."), "It should be \"double\" or \"integer\""))
  }
  if (!is.numeric(riboMax)) {
    stop(paste("riboMax type is", paste0(typeof(riboMax), "."), "It should be \"double\" or \"integer\""))
  }
  if (!is.null(separator) && typeof(separator) != "character") {
    stop(paste("separator type is", paste0(typeof(separator), "."),"It should be \"character\""))
  }
  if (!is.null(genes_file) && typeof(genes_file) != "character") {
    stop(paste("genes_file type is", paste0(typeof(genes_file), "."),"It should be \"character\""))
  }
  if (!is.null(barcodes_file) && typeof(barcodes_file) != "character") {
    stop(paste("barcodes_file type is", paste0(typeof(barcodes_file), "."),"It should be \"character\""))
  }

  # Check if the given paths exist
  if (!is_running_in_docker()) {
    if (!file.exists(input_file_path)) {
      stop(paste("input_file_path:", input_file_path, "does not exist."))
    }
  }

  # Obtain parent folder from input_file_path.
  parent_folder <- dirname(input_file_path)

  # Obtain matrix name and file format from input_file_path to create the variable matrix_file
  input_file_path_parts <- strsplit(basename(input_file_path), "\\.")[[1]]
  matrix_name <- input_file_path_parts[1]
  format <- input_file_path_parts[2] # Uses extension as an heuristic.
  matrix_file <- paste0(matrix_name, ".", format)

  # Executing the docker job
  rrundocker::run_in_docker(
    image_name = paste0("repbioinfo/singlecelldownstream:latest"),
    volumes = list(
      c(parent_folder, "/scratch")
    ),
    additional_arguments = c(
      "Rscript /home/mitoRiboFilter.R",
      paste0("/scratch/",matrix_name),
      mitoMin,
      mitoMax,
      riboMin,
      riboMax,
      separator,
      paste0("/scratch/",genes_file),
      paste0("/scratch/",barcodes_file)
    )
  )
}
