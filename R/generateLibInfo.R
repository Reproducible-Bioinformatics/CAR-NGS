#' generateLibInfo
#'
#' This function is used to generate sequencing library information from an XML file,
#' it returns two output files in the parent_folder, libseqInfo.txt and libseqInfo2.txt.
#'
#' @param xml_file_path, a character string indicating the path to the XML file containing library information.
#' @param configType, a character string indicating the configuration type, as defined in config_types.json.
#' Examples include "HTGTS_mouse", "HTGTS_human", "CELTICSseq", "polyA"
#'
#' @author Luca Alessandrì, Agata D'Onofrio
#'
#' @examples
#' \dontrun{
#'  generateLibInfo(
#'     xml_file_path = "/the/xml/file/path",
#'     configType = "HTGTS_mouse"
#'  )
#' }
#' @export
generateLibInfo <- function(xml_file_path, configType){

  # Type checking.
  if (typeof(xml_file_path) != "character") {
    stop(
      paste(
        "xml_file_path type is", paste0(typeof(xml_file_path), "."),
        "It should be \"character\""
      )
    )
  }
  if (typeof(configType) != "character") {
    stop(
      paste(
        "configType type is", paste0(typeof(configType), "."),
        "It should be \"character\""
      )
    )
  }

  # Check if the given path exists
  if (!is_running_in_docker()) {
    if (!file.exists(xml_file_path)) {
      stop(paste("xml_file_path:", xml_file_path, "does not exist."))
    }
  }

  # Obtain parent folder from xml_file_path
  parent_folder <- dirname(xml_file_path)

  # Obtain xml name and file format from xml_file_path to create the
  # variable xml_file
  xml_file_path_parts <- strsplit(basename(xml_file_path), "\\.")[[1]]
  xml_name <- xml_file_path_parts[1]
  format <- xml_file_path_parts[2] # Uses extension as an heuristic.
  xml_file <- paste0(xml_name, ".", format)

  # Executing the docker job
  run_in_docker(
    image_name = paste0("repbioinfo/htgts_pipeline_lts_v16:latest"),
    volumes = list(
      c(parent_folder, "/Data")
    ),
    additional_arguments = c(
      "python3 /Algorithm/sample_sheetTolibInfo.py",
      paste0("/Data/", xml_file),
      "/Data/libseqInfo.txt",
      "/Data/libseqInfo2.txt",
      configType
    )
  )
}
