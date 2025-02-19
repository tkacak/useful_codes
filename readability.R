# Word Count:
library(devtools)
devtools::install_github("benmarwick/wordcountaddin", type = "source", dependencies = TRUE)

# How to use
# 
# Open a Rmd file in RStudio.
# Select some text, it can include YAML, code chunks and inline code
# Go to Tools > Addins in RStudio and click on Word count or Readability. Computing Readability may take a few moments on longer documents because it has to count syllables for some of the stats.
# Look in the console for the output
library(wordcountaddin)
