#| class: output-block
# file: script.R
library(dplyr)

file_path <- "script.log"

average_mpg <- mtcars %>%
  group_by(cyl) %>%
  summarise(average_mpg = mean(mpg, na.rm = TRUE))

model <- lm(mpg ~ wt, data = mtcars)
summary(model)

# log file path: script.log
