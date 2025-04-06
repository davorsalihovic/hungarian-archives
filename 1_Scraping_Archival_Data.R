#  Quantitative Analysis of the Pre-Mohacs Collection at the Hungarian National Archives ------------
#    Davor Salihovic, University of Antwerp                                                 
#                      davor.salihovic@uantwerpen.be

# Clean the R environment
rm(list = ls())

library(rvest); library(tidyr)
library(dplyr); library(purrr)

# 1. SCRAPING DATA FROM THE HUNGARIAN NATIONAL ARCHIVES ----

scrape_document <- function(doc_id) {
  url <- paste0("https://archives.hungaricana.hu/en/charters/", doc_id, "/")
  page <- tryCatch(read_html(url), error = function(e) return(NULL))
  if (is.null(page)) return(NULL)
  data <- page %>% html_nodes("table") %>% html_table(fill = TRUE)
  if (length(data) == 0) return(NULL)
  doc_info <- as.data.frame(data[[1]], stringsAsFactors = FALSE)
  relevant_categories <- c("DL-DF", "Date", "Place of dating", "Document type", "Issuer of charter", "Survival form", "Alternative date",
                           "Language", "Index", "Subject", "Abstract", "The old reference of the records")
  extracted <- doc_info %>%
    filter(X1 %in% relevant_categories) %>%
    setNames(c("Category", "Value")) %>%
    spread(key = Category, value = Value)
  
  extracted$Doc_ID <- doc_id
  return(extracted)
}

start.time <- Sys.time()
doc_ids <- seq(0, 318752, 1)
documents <- doc_ids %>%
  map_dfr(scrape_document)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)

save.image("data/documents.RData")