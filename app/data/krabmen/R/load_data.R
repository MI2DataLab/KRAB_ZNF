.tumor.names <- function() {
  as.character(read.csv('resources/tumor.names.csv')$name)
}

.gene.names <- function() {
  as.character(read.csv('resources/gene-names.csv')$name)
}

.tumor.data <- function() {
  
  names <- .tumor.names()
  
  f <- sapply(names, function(name) {
    
    file.name <- file.path('resources', paste0(name, '.RData'))
    
    local(get(load(file.name)))
  })
  
  print('tumor data done')
  
  return(f)
}
