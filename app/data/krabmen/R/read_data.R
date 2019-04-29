# from boxplots
load.datasets <- function() {

  load('data/krabmen/boxplots/resources/all.Rdata')
  assign('all.tumors.long2', all.tumors.long, inherits = T)
  
  load('data/krabmen/boxplots/resources/pvalues.Rdata')
  assign('pvalues2', pvalues, inherits = T)
}

load.tumor.names <- function() {
  read.csv('data/krabmen/boxplots/resources/tumor-names.csv', stringsAsFactors = F)$name
}

load.gene.list <- function() {
  read.csv('data/krabmen/boxplots/resources/gene-names.csv', stringsAsFactors = F)$name
}

# from boxplot-subtypes
.tumor.data <- function() {
  
  names <- read.csv('data/krabmen/boxplots-subtypes/resources/tumor.names.csv', stringsAsFactors = F)$name
  
  f <- sapply(names, function(name) {
    
    file.name <- file.path('data/krabmen/boxplots-subtypes/resources', paste0(name, '.RData'))
    
    local(get(load(file.name)))
  })
  
  #print('tumor data done')
  
  return(f)
}