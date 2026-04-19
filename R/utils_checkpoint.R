# ==============================================================================
# File: /R/utils_checkpoint.R
# Purpose: Safely manage local simulation state to prevent data loss.
#          Includes atomic saving, resume-checking, and final chunk compilation.
# ==============================================================================
#
#' Safely save an R object to disk (Atomic Save)
#'
#' Writes to a temporary file first, then renames it. This prevents file 
#' corruption if the R session crashes or the computer loses power exactly 
#' while the file is being written.
#'
#' @param object The R object (usually a data.frame) to save.
#' @param file_path The final destination path (e.g., "../output/temp/chunk_1.rds").
#' @export
safe_save_rds <- function(object, file_path) {
  # Ensure the directory exists
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create a temporary filename
  temp_file <- paste0(file_path, ".tmp")
  
  # Save to the temporary file
  saveRDS(object, file = temp_file)
  
  # Atomically rename to the final file path (overwrites if exists)
  file.rename(from = temp_file, to = file_path)
}

#' Check if a simulation chunk is already completed
#'
#' Used at the top of a loop to skip iterations that were already computed
#' in a previous session.
#'
#' @param file_path The path to the expected .rds file.
#' @return Logical TRUE if the file exists and is valid, FALSE otherwise.
#' @export
is_computed <- function(file_path) {
  # If the file doesn't exist, we definitely need to compute it
  if (!file.exists(file_path)) {
    return(FALSE)
  }
  
  # Double-check that it's a valid, readable RDS file and not a 0-byte corrupted file
  is_valid <- tryCatch({
    # Just reading the header to confirm it's structurally sound
    invisible(readRDS(file_path))
    TRUE
  }, error = function(e) {
    warning(sprintf("Found corrupted chunk, will overwrite: %s", file_path))
    FALSE
  })
  
  return(is_valid)
}

#' Merge all temporary chunks into the final dataset
#'
#' Reads all chunked .rds files matching a pattern in a directory, binds 
#' them into a single data.frame, and saves the final result.
#'
#' @param temp_dir Directory containing the temporary chunks.
#' @param pattern Regex pattern to match files (e.g., "^02_.*\\.rds$").
#' @param final_output_path Where to save the merged dataset.
#' @param clear_temp Logical: If TRUE, deletes the chunks after successful merge.
#' @export
compile_checkpoints <- function(temp_dir, pattern, final_output_path, clear_temp = FALSE) {
  
  files_to_merge <- list.files(path = temp_dir, pattern = pattern, full.names = TRUE)
  
  if (length(files_to_merge) == 0) {
    stop(sprintf("No files found in %s matching pattern '%s'", temp_dir, pattern))
  }
  
  cat(sprintf("Merging %d chunk files...\n", length(files_to_merge)))
  
  # Read and bind all rows safely
  merged_data <- purrr::map_dfr(files_to_merge, readRDS)
  
  # Safe save the final monolithic dataset
  safe_save_rds(merged_data, final_output_path)
  
  cat(sprintf("Successfully saved final dataset to: %s\n", final_output_path))
  
  # Clean up temp files to save disk space, if requested
  if (clear_temp) {
    file.remove(files_to_merge)
    cat("Temporary chunks cleared.\n")
  }
  
  return(invisible(merged_data))
}