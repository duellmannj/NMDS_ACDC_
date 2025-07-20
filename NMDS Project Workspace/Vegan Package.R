#This code was created by Nicholas Duellman to explore the ACDC project data for community shifts and plot NMDS, ANOSIM and PERMANOVA tests as well
#as diversity. For any questions about the code, please reach out to nickduellman@gmail.com. 
#This code largely follows vegan package. Data format for species abundance can be found in vegan package for NMDS
# ordination and to create a species matrix (this is critical!). 

#I want to thank everybody in the AEFL at UW-GB for the oppurtunity this summer and to have access to this data




# install packages


install.packages("vegan")
install.packages("permute")
install.packages("tidyverse")
install.packages("lattice")
install.packages("ggplot2")
install.packages("dplyr")

# load libraries
library(vegan)
library(permute)
library(tidyverse)
library(lattice)
library(ggplot2)
library(dplyr)
library(ecodist)

set.seed(1108)

Nick Duellman 


"Evidence to community shifts across years"

setwd("C:/Users/nickd/OneDrive/Desktop/NMDS Project Workspace/data") #set a working directory for all the files (species abundance and environmental data can be seperate files, make sure the rows match up)

#Import data as .csv format 

env_data <- read.csv("C:/Users/nickd/OneDrive/Desktop/NMDS Project Workspace/data/2020-2022 Ash. Creek Environmental and Pop. Data/Ashw. Creek Sites and Fish Data 2020 - 2022 - Environmental Data and Temporal Data.csv")

# edit and trim data #any N/A becomes a 0
data_fish <- abund_fish_data
data_fish[is.na(data_fish)] <- 0

env_data <- env_data_ash[-(43:70), -(9:12)]
view(env_data)

# extract species presence from fish dataset
fish_community <- data_fish[,4:25]

# data frame to matrix for bray-curtis dissimilarity
fish_matrix <- as.matrix(fish_community)

# distance matrix for fish community
fish.dist <- vegdist(fish_matrix, method = "bray")

# ANOSIM - test to see if there are statistical differences between groups are significant 

fish.anosim <- with(env_data, anosim(fish.dist, env_data$Site), distance = "bray", strata = env_data$Year)
summary(fish.anosim)


fish.anosim <- anosim(fish.dist, env_data$Site, distance = "bray", permutations = 999)
  
fish.anosim <- anosim(fish.dist, env_data$Year, distance = "bray", permutations = 999)

fish.anosim <- anosim(fish.dist, env_data$Season, distance = "bray", permutations = 999)

fish.anosim

"R stat of 0.02648 and P > 0.05 (P = 0.178) indicate that there is no statistical significance to suggest that fish community changes across years"

"R stat of 0.4574 and P = 0.001 indicate that there is statistical significance to suggest that fish community changes across sites (permutations = 999).




plot(fish.anosim)

Site
P = 0.001 R stat = 0.4574

Year
P = 0.175 R stat = 0.02648


--------

# NMDS Ordination

fish.NMDS <- metaMDS(fish_matrix, 
                    distance =  
                      "bray", k = 3, autotransform = 
                      TRUE, trymax = 
                      100)

# check stress value
fish.NMDS$stress

Value less than 0.20 [Value = 0.1202]

# data visualization
plot(fish.NMDS)
ordiplot(fish.NMDS, type = "n")
orditorp(fish.NMDS, display = "species", col = "blue", air = 0.10)
orditorp(fish.NMDS, display = "sites", cex = 1.25, air = 0.10)

# ellipse plot (site)
ordiplot(fish.NMDS)
orditorp(fish.NMDS, display = "species", col = "blue", air = 0.10)
ordiellipse(fish.NMDS, env_data$Site, label = FALSE,
            col=c("darkorchid1",
                  "darkslategray1",
                  "bisque",
                  "brown3", 
                  "yellow"
                  ), draw = "polygon", alpha=120)
legend("topright", title = "Site",
       c("2", "3", "4", "6", "8"), fill =c("darkorchid1",
                                       "darkslategray1",
                                       "bisque",
                                       "brown3",
                                       "yellow"), 
       horiz=FALSE, cex = 0.9)


# ellipse plot (year)
ordiplot(fish.NMDS)
orditorp(fish.NMDS, display = "species", col = "blue", air = 0.10)
ordiellipse(fish.NMDS, env_data$Year, label = FALSE,
            col=c("darkorchid1",
                  "darkslategray1",
                  "bisque"
            ), draw = "polygon", alpha=120)

legend("topright", title = "Year",
       c("2020", "2021", "2022"), fill =c("darkorchid1",
                                           "darkslategray1",
                                           "bisque"
                                          ), 
       horiz=FALSE, cex = 0.9)


# season plot

ordiplot(fish.NMDS)
orditorp(fish.NMDS, display = "species", col = "blue", air = 0.10)
ordiellipse(fish.NMDS, env_data$Season, label = FALSE,
            col=c("green",
                  "lightblue1",
                  "red"
            ), draw = "polygon", alpha=120)
legend("topright", title = "Season",
       c("Spring", "Summer", "Fall"), fill =c("green",
                                           "lightblue1",
                                           "red"
                                          ), 
       horiz=FALSE, cex = 0.9)

#ordicluster(fish.NMDS,hclust(vegdist(fish_matrix,"bray"))) 


stressplot(fish.NMDS)

# betadisper and envfit to test ellipse significance

# betadisper looks at significant differences in groups 
bd <- betadisper(fish.dist, env_data$Season)
bd <- betadisper(fish.dist, env_data$Year)
bd <- betadisper(fish.dist, env_data$Site)
anova(bd)
permutest(bd, permutations = 999, pairwise = TRUE) #compares all variables with another

"betadisper function looks at if there is significant differences among grouping types"

"we can even tukey test post hoc the results to gain confirm significance"
"Anova table of env_data$Site indicated that there is significant differences in sites 
with the fish distance matrix p = 0.0352. We can confirm this with running a permutest
this looks at post hoc Tukey Test and confirms that the differences are significant"

"We can even run a pairwise permutest to see the differences across variables"

"permutest with permutations and pairwise with env_data$Sites can show us that the 
differences in sites are both significant and which sites are most significant from 
another"

"Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
         2        3        4        6     8
2          0.699000 0.611000 0.813000 0.022
3 0.698892          0.851000 0.554000 0.007
4 0.617391 0.868440          0.506000 0.021
6 0.797024 0.560698 0.513633          0.042
8 0.012408 0.010189 0.018288 0.035914  "





# envfit looks at if the groups have a significant impact on the community composition
ef <- envfit(fish.NMDS ~ env_data$Site, permutations = 999)
ef <- envfit(fish.NMDS ~ env_data$Season, permutations = 999)
ef <- envfit(fish.NMDS ~ env_data$Year, permutations = 999)


(print(ef)


"Results from envfit show that Site has a significant impact on community composition, 
R2 = .6745 with P < 0.001 and at 999 permutations. Season is not significant on community composition, 
R2 = 0.0245 and P = 0.75. Year is not significant, R2 = 0.0758 and P = 0.213, Year
explains 7.5% of the variation in community composition"

"'Envfit' function can compare the significance of centroid groupings on the species composition matrix"

en = envfit(fish.NMDS, env_data, permutations = 999, na.rm = TRUE)
en

"Will the null (p greater than 0.05) change if we add more years of site data? 
"Will year or season become significant factors that influence community composition?"
"We permutate 999 times to get 999 random reorgnizations of data to confirm our significance"


fish.NMDS$points






"data analysis of NMDS"

#fish_permanova <- adonis2(as.matrix(fish_community) ~ env_data$Year + env_data$Season,
                          data = env_data,
                          permutations = 999, 
                          method = "bray")

#fish_permanova

#fish.inv.dist <- vegdist(fish_matrix, method = "bray")
#fish.inv.dispersion <- betadisper(fish.inv.dist, group = env_data$Year)
#permutest(fish.inv.dispersion)

-------------------------------------------------------------------------------------
set.seed(1108)

# this looks across all sites (strata)
fish_permanova <- adonis2(as.matrix(fish_community) ~ Year, 
                          data = env_data,permutations = 999, strata = env_data$Site,
                          method = "bray",by = "margin")

#fish_permanova

"Year is significant across sites but with low effect although significant (P < 0.045)"


fish_permanova <- adonis2(as.matrix(fish_community) ~ Site, 
                          data = env_data,permutations = 999, strata = env_data$Year,
                          method = "bray",by = "margin")

fish_permanova

"Species abundance differences in sites are significant across years, R2 = 0.18698 P = 0.001"




fish_permanova <- adonis2(fish.dist ~ Site, 
                          data = env_data,permutations = 999, strata = env_data$Year,
                          method = "bray", by = 'margin')


with(env_data, adonis2(fish.dist ~ env_data$Oxygen..mg.L + env_data$Temperature + env_data$Site, permutations = 999,
                       strata = env_data$Year, method = "bray"))


with(env_data, adonis2(fish.dist ~ env_data$Oxygen..mg.L*env_data$Temperature:env_data$Site, permutations = 999,
                       strata = env_data$Year, method = "bray"))


fish_adonis <- with(env_data, adonis2(fish.dist ~ 
                                        env_data$Oxygen..mg.L*env_data$Temperature + env_data$Site,
                                      data = env_data, 
                                      strata = env_data$Site,
                                      method = "bray", permutations = 999))

fish_adonis <- with(env_data, adonis2(fish.dist ~ 
                                        env_data$Oxygen..mg.L*env_data$Temperature + env_data$Site,
                                      data = env_data, 
                                      strata = env_data$Site,
                                      method = "bray", permutations = 999))

fish_adonis
                    
# look at the interaction of site and temperature as predictors and site as strata, we got P = 0.59 and R2 = 0.24, 
# not good or significant predictors even with season and year as strata

fish_adonis <- with(env_data, adonis2(fish.dist ~ 
                                        env_data$Oxygen..mg.L * env_data$Temperature + env_data$Site,
                                      data = env_data, 
                                      strata = env_data$Year,
                                      method = "bray", permutations = 999))

fish_adonis

"The interaction of Oxygen and Temperature with the constant of Site with year as strata is significant and R2 = 0.24 and P = 0.001"

"adonis2(formula = fish.dist ~ env_data$Oxygen..mg.L * env_data$Temperature + env_data$Site, data = env_data, permutations = 999, method = "bray", strata = env_data$Year)
Df SumOfSqs      R2      F Pr(>F)    
Model     4   2.7251 0.24239 2.9595  0.001 ***
  Residual 37   8.5173 0.75761                  
Total    41  11.2424 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
  


"A*B*C (test for significance then pairwise, then A*B*C - A:B:C and then A+B+C then look for pairwise in all)
-------------------------------------------------------------------------------------

fish_adonis %>% as.data.frame() %>% 
  rownames_to_column("term") %>% 
  mutate_if(is.numeric, round, digits = 3) %>% 
  knitr::kable()



library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
?pairwise.adonis2

"Pairwise Adonis Test"


--------------------------------------------------


# Betadisper

fish.betadisper <- betadisper(fish.dist, group = env_data$Oxygen..mg.L) 
fish.betadisper <- betadisper(fish.dist, group = env_data$Temperature) 
fish.betadisper <- betadisper(fish.dist, group = env_data$Site) 
fish.betadisper <- betadisper(fish.dist, group = env_data$Year) 

permutest(fish.betadisper)

"P = 0.65, non-significant for Temperature
P = 0.742, non-significant for Oxygen
P = 0.036, significant for Site
P = 0.661, non-significant for Year"

--------------------------------------------------------------

"diversity function in vegan"

H <- diversity(fish_community)
simp <- diversity(fish_community, "simpson")

S <- specnumber(fish_community)
J <- H/log(S)

(alpha <- with(env_data, tapply(specnumber(fish_community), Site, mean)))
(gamma <- with(env_data, specnumber(fish_community, Site)))
gamma/alpha - 1


fish_2020 <- fish_community[ 2:14, , drop = FALSE ]
fish_2021 <- fish_community[15:28, , drop = FALSE ]
fish_2022 <- fish_community[29:42, , drop = FALSE ]

data_2020 <- env_data[ 2:14, , drop = FALSE]
data_2021 <- env_data[ 15:28, , drop = FALSE]
data_2022 <- env_data[ 29:43, , drop = FALSE]

(alpha <- with(data_2020, tapply(specnumber(fish_2020), Site, mean)))
(gamma <- with(data_2020, specnumber(fish_2020, Site)))
gamma/alpha - 1

(alpha <- with(data_2021, tapply(specnumber(fish_2021), Site, mean)))
(gamma <- with(data_2021, specnumber(fish_2021, Site)))
gamma/alpha - 1

(alpha <- with(data_2022, tapply(specnumber(fish_2022), Site, mean)))
(gamma <- with(data_2022, specnumber(fish_2022, Site)))
gamma/alpha - 1

# Shannon is Alpha





---------------------------------------------------------------------------------

# Species Proportion

library(dplyr)
library(tidyverse)

#df <- data.frame(value = c(5,3,4,6,5))

head(fish_community)

rownum(fish_community)

# 22 distinct species

# 2,163 total speciments

'df_spec' <- data.frame(colSums(fish_community)) # data frame that gets sums of the colums of species abundance
'total_species_df' <- data.frame(sum(df_spec$colSums.fish_community.)) # total species collected


# get the proportion of each species from the total 
prop.df <- data.frame(df_spec$colSums.fish_community./total_species_df$sum.df_spec.colSums.fish_community.)
View(prop.df) 


# rename the data for simplicity
prop.df <- prop.df %>% rename( fish_prop = 	
df_spec.colSums.fish_community..total_species_df.sum.df_spec.colSums.fish_community. 
)

# combine the columns into one data frame
df_fish_new <- cbind(df_spec, Proportion = prop.df$fish_prop * 100) # get percentages
View(df_fish_new)

knitr::kable(df_fish_new) # view data frame


-----------------------------------------------------------------------------------------

# individual examination of groupings with envfit

(This was already demonstrated)
envfit(fish.NMDS ~ env_data$Season, data = env_data, permutations = 999)

envfit(fish.NMDS ~ Site, data = env_data, permutations = 999)

envfit(fish.NMDS ~ env_data$Year, data = env_data, permutations = 999)



-------------------------------------------------------------------------------------

# Species diversity and abundance (see line 365 to 358 for splitting up matrix into years)

fish_2020 <- fish_2020[ , colSums(fish_2020 == 0) < nrow(fish_2020) ]
fish_2020

rownum(fish_2020)

head(fish_2020)

ncol(fish_2020)


fish_2021 <- fish_2021[ , colSums(fish_2021 == 0) < nrow(fish_2021) ]
fish_2021



head(fish_2021)

ncol(fish_2021)

fish_2022 <- fish_2022[ , colSums(fish_2022 == 0) < nrow(fish_2022) ]

fish_2022


head(fish_2022)

ncol(fish_2022)


species_abundance_2020 <- sum(as.matrix(fish_2020), na.rm = TRUE)
species_abundance_2020

species_abundance_2021 <- sum(as.matrix(fish_2021), na.rm = TRUE)
species_abundance_2021

species_abundance_2022 <- sum(as.matrix(fish_2022), na.rm = TRUE)
species_abundance_2022




