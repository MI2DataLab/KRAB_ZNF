library(knitr)
library(XLConnect)
library(xlsx)
write.csv(illumina_humanmethylation_27_data[1:100,], file = "illumina_humanmethylation_27.csv")

raw_data <- read.xlsx2("Raw_data.xlsx", sheetIndex = 1)
wb = loadWorkbook("Raw_data.xlsx")
df = readWorksheet(wb, sheet = "Sheet1", header = TRUE)