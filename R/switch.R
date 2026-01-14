#' Fixes link-check problem
#' - Finds files that use the string "(.mailto:" and corrects to "(mailto:"
#' - Updates links to `.qmd` files to point to `.md` files.
#'
#' @param file_list vector of filenames
#'
#' @return messages. Side effect: modifies files.
modify_files <- function(file_list) {
  # Create an empty vector to store file names that need modifications
  matching_files <- c()

  # Iterate over each file to check for issues
  for (file in file_list) {
    # Read the contents of the file
    file_contents <- readLines(file)

    # Check if the file contains the string "(.mailto:" OR references to `.qmd` files
    if (any(grepl("\\(\\.mailto:", file_contents)) || any(grepl("\\.qmd\\)", file_contents))) {
      # Add the file name to the vector
      matching_files <- c(matching_files, file)
    }
  }

  # Iterate over the matching files to apply fixes
  for (file in matching_files) {
    # Read the contents of the file
    file_contents <- readLines(file)

    # Fix email links
    modified_contents <- gsub("\\(\\.mailto:", "(mailto:", file_contents)

    # Update links to `.qmd` files to point to `.md` files
    # Matches patterns like `[text](filename.qmd)` and replaces `.qmd` with `.md`
    modified_contents <- gsub("\\.qmd\\)", ".md)", modified_contents)

    # Write the modified contents back to the file
    writeLines(modified_contents, file)

    # Print a message indicating what modifications have been made
    message("Modified file:", file, "\n")
  }

  # Print a list of matching files that were modified
  message("Matching files:", matching_files, "\n")
}

# Get all `.qmd` files
all_qmd <- list.files(full.names = FALSE, all.files = FALSE, pattern = ".qmd$", recursive = TRUE)

# Modify files to fix email links and update `.qmd` references
modify_files(all_qmd)

# Generate a list of `.md` filenames that will replace `.qmd` files
all_md <- gsub(".qmd$", ".md", all_qmd)

# Rename all `.qmd` files to `.md`
file.rename(all_qmd, all_md)