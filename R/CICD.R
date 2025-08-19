# Load required packages
library(spelling)
library(readr)
library(stringr)
library(dplyr)

#-------------------------- Spell-check Without Modifying Files --------------------------

# Function to extract non-code content while keeping comments from code blocks
clean_text_for_spellcheck <- function(file) {
  text <- read_lines(file)

  text_clean <- text %>%
    # Remove Markdown links but keep link text
    str_replace_all("\\[([^\\]]+)\\]\\([^)]*\\)", "\\1") %>%
    # Remove standalone URLs
    str_replace_all("https?://[^\\s)\"']+", "")

  return(text_clean)
}

#-------------------------- Spell-check Without Modifying Files --------------------------

# Get all .qmd files
qmd_files <- list.files(pattern = ".*\\.qmd$", recursive = TRUE)
file <- qmd_files[63]
# Run spell check per file while ignoring code chunks

spell_check_results <- lapply(qmd_files, function(file) {
  cleaned_text <- clean_text_for_spellcheck(file)
  words <- spelling::spell_check_text(cleaned_text, ignore = read_lines("inst/WORDLIST"))

  if (nrow(words) > 0) {
    words$file <- file # Add filename column
  }

  return(words)
})

# Combine results into a single dataframe
all_typos <- bind_rows(spell_check_results)

# Print results
if (nrow(all_typos) > 0) {
  print(all_typos %>% select(file, word, found))
} else {
  message("No spelling errors found!")
}

#-------------------------- Add Words to Wordlist (If Needed) --------------------------

# Uncomment if you want to manually add words to the wordlist
# write(all_typos$word, file = "inst/WORDLIST", append = TRUE)
# sort_words <- sort(read_lines("inst/WORDLIST")) %>% unique()
# write_lines(sort_words, "inst/WORDLIST")

#-------------------------- Style-check ----------------------------------------

# This is what happens in CI/CD:

styler::style_dir(dry = "fail")

# Fix it locally:
styler::style_dir()

# Now it should pass:
styler::style_dir(dry = "fail")
