#' fixes link-check problem
#' Finds whether any files use the string: "(.mailto:" for linking email adresses
#' then just replaces with: "(mailto:" as is markdown standard
#'
#' @param file_list vector of filenames
#'
#' @return messages only. side effect is: changes files.
modify_files <- function(file_list) {
  # Create an empty vector to store the file names that contain the string
  matching_files <- c()

  # Iterate over each file in the directory
  for (file in file_list) {
    # Read the contents of the file
    file_contents <- readLines(file)

    # Check if the file contains the string "(.mailto:"
    if (any(grepl("\\(\\.mailto:", file_contents))) {
      # Add the file name to the vector
      matching_files <- c(matching_files, file)
    }
  }

  # Iterate over the matching files
  for (file in matching_files) {
    # Read the contents of the file
    file_contents <- readLines(file)

    # Remove the "." from each line that contains the string
    modified_contents <- gsub("\\(\\.mailto:", "(mailto:", file_contents)

    # Write the modified contents back to the file
    writeLines(modified_contents, file)

    # Print a message indicating the modification has been made
    message("Modified file:", file, "\n")
  }

  # Print the list of matching files
  message("Matching files:", matching_files, "\n")
}

# get all qmd files
all_qmd <- list.files(full.names = FALSE, all.files = FALSE, pattern = ".qmd$", recursive = TRUE)

# modify if needed the link to email problem for link checker
modify_files(all_qmd)
# get filenames ending with .md
all_md <- gsub(".qmd$", ".md", all_qmd)
# rename all files
file.rename(all_qmd, all_md)
