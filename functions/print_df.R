library(reactable)
library(reactablefmtr)

print_df <- function(dataset, n = 10) {
  out <- dataset

  reactable(
    head(out, n),
    compact = TRUE,
    wrap = FALSE,
    bordered = TRUE,
    striped = TRUE,
    highlight = TRUE,
    defaultPageSize = 5,
    theme = reactableTheme(
      color = "#333",
      backgroundColor = "white",
      borderColor = "#ddd",
      stripedColor = "#f8f9fa",
      highlightColor = "#dfe6f2",
      cellPadding = "8px"
    )
  ) %>% 
    add_title("Sample of Data")
}
