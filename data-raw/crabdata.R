crabdata <- read.csv("data-raw/crabdata.csv", stringsAsFactors = FALSE)
save(crabdata, file = "data/crabdata.RData")