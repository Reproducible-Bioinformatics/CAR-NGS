#' MitoRiboUmi
#'
#' This function executes a ubuntu docker showing ribosomal mitochondrial protein genes fraction in each cell,
#' providing a pdf file called Ribo_mito.pdf inside result_dir_path
#'
#' @param result_dir_path, a character string indicating the path where the result files are saved
#' @param input_file_path, a character string indicating the path of the file, with file name and extension included
#' @param separator, separator used in the count table
#' @param gtf.name, name for the gtf file to be used.
#' @param bio.type, ENSEMBL biotype of interest, default for single-cell protein_coding
#' @param umiXgene, a integer defining how many UMI are required to call a gene as present. default: 3
#' @author Luca Alessandri, alessandri [dot] luca1991 [at] gmail [dot] com, University of Torino
#'
#' @return a pdf file called Ribo_mito.pdf
#' @examples
#'\dontrun{
#' mitoRiboUmi(
#'   result_dir_path = "/the/result/dir",
#'   input_file_path = "/the/input/file.csv",
#'   separator = ",",
#'   gtf.name="Homo_sapiens.GRCh38.112.gtf",
#'   bio.type="protein_coding",
#'   umiXgene=3
#' )
#'}
#' @export
mitoRiboUmi <- function(result_dir_path, input_file_path, separator, gtf.name, bio.type, umiXgene){

  # Type checking.
  if (typeof(result_dir_path) != "character") {
    stop(paste("result_dir_path type is", paste0(typeof(result_dir_path), "."), "It should be \"character\""))
  }
  if (typeof(input_file_path) != "character") {
    stop(paste("input_file_path type is", paste0(typeof(input_file_path), "."), "It should be \"character\""))
  }
  if (typeof(separator) != "character") {
    stop(paste("separator type is", paste0(typeof(separator), "."),"It should be \"character\""))
  }
   if (typeof(gtf.name) != "character") {
    stop(paste("gtf.name type is", paste0(typeof(gtf.name), "."), "It should be \"character\""))
  }
  if (typeof(bio.type) != "character") {
    stop(paste("bio.type type is", paste0(typeof(bio.type), "."), "It should be \"character\""))
  }
  if (typeof(umiXgene) != "integer") {
    stop(paste("umiXgene type is", paste0(typeof(umiXgene), "."), "It should be \"integer\""))
  }

  # Check if the given paths exist
  if (!rrundocker::is_running_in_docker()) {
    if (!file.exists(input_file_path)) {
      stop(paste("input_file_path:", input_file_path, "does not exist."))
    }
    if (!dir.exists(result_dir_path)) {
      stop(paste("result_dir_path:", result_dir_path, "does not exist."))
    }
  }

  # Obtain parent folder from input_file_path.
  parent_folder <- dirname(input_file_path)

  # Obtain matrix name and file format from input_file_path.
  input_file_path_parts <- strsplit(basename(input_file_path), "\\.")[[1]]
  matrix_name <- input_file_path_parts[1]
  format <- input_file_path_parts[2] # Uses extension as an heuristic.

  # Executing the docker job
  rrundocker::run_in_docker(
    image_name = "docker.io/repbioinfo/ribomitoumi:latest",
    volumes = list(
      c(result_dir_path, "/data"),
      c(parent_folder, "/scratch")
    ),
    additional_arguments = c(
      "Rscript /home/main.R",
      matrix_name,
      format,
      paste0('"', separator, '"'),
      gtf.name,
      bio.type,
      umiXgene
    )
  )
}
