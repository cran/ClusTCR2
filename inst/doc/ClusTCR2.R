## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=TRUE---------------------------------------------------------
library(ClusTCR2)

## ----echo=TRUE----------------------------------------------------------------
# Example usage of mcl_cluster function with a stored file
example_file <- read.csv(system.file("extdata", "my_data.csv",package = "ClusTCR2"))
# Perform clustering using mcl_cluster function
step1 <- ClusTCR(example_file,allele = FALSE)
# Print the result
step1[1:6,1:6]

## ----echo=TRUE----------------------------------------------------------------
# Example usage of mcl_cluster function with a stored file
step2 <- mcl_cluster(step1)
step2[[1]][1:3,1:3]
head(step2[[2]][1:6,1:6])

## ----echo=TRUE----------------------------------------------------------------
# Visualization of the network plot
netplot_ClusTCR2(step2, label = "Name_selected", Clust_selected = 1)

## ----echo=TRUE----------------------------------------------------------------
# Visualization of the network plot
# step2[[1]]
subset(step2[[1]],step2[[1]]$Clust_size_order == 1)
motif_plot(step2,Clust_selected = 1)

