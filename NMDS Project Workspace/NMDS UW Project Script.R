# install packages
install.packages("vegan")
install.packages("permute")
install.packages("tidyverse")
install.packages("lattice")
install.packages("ggplot")


# load libraries
library(vegan)
library(permute)
library(tidyverse)
library(lattice)
library(ggplot)

# analyze species proportions
setwd("C:/Users/nickd/OneDrive/Desktop/NMDS Project Workspace")
species_count <- table(species_data$Species)

species_prop <- species_count / sum(species_count)

species_df <- as.data.frame(species_prop)
colnames(species_df) <- c("Species", "Proportion")



