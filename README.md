# GWAS_project

**All data there are in folders, but they're in  in zip format**

## First step: prepare the data for use by the model.

Conversion of genotype data from vcf format to hapmap was performed in the [TASSEL 5](https://tassel.bitbucket.io/) program. 
The conversion from hapmap to numerical was done with the R package [GAPIT3](https://zzlab.net/GAPIT/gapit_help_document.pdf).
```ruby
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
```
Setting the working directory:
```ruby
setwd('Your/working/directory')
```
Conversion data:
Sesame data
```ruby
myG_sesame <- read.table("sesame_complicated_genotypes.hmp.txt", fill = TRUE)
myGAPIT_sesame <- GAPIT(G=myG_sesame, output.numerical=TRUE)
myGD_sesame= myGAPIT_sesame$GD
myGM_sesame= myGAPIT_sesame$GM
```
Soybean data
```ruby
myG_soybean <- read.table("soybean_simple_genotypes.hmp.txt", fill = TRUE)
myGAPIT_soybean <- GAPIT(G=myG_soybean, output.numerical=TRUE)
myGD_soybean= myGAPIT_soybean$GD
myGM_soybean= myGAPIT_soybean$GM
```
## Step two: Preparing phenotype data.
This step was performed in R.

**Soybean phenotype data**:
```ruby 
soy_phenotypes <- read.table('soybean_simple_phenotypes.tsv', sep = '\t') 
sum(is.na(soy_phenotypes)) # 0
soy_phenotypes[2:250, 3:17] <- lapply(soy_phenotypes[2:250, 3:17], as.numeric) # convert all values into a numeric format
soy_phenotypes_leu <- cbind(soy_phenotypes$V1, soy_phenotypes$V10) # column 1 with the taxon name and column 10 with the amino acid leucine are selected
colnames(soy_phenotypes_leu) <- c ("Taxa", "Leu") 
```
For the **sesame data**, a selection was made according to the Plant height characteristic for 2018, after which the average value for all groups was taken.
```ruby
sesame_phenotypes <- read.xlsx('sesame_complicated_phenotypes.xlsx', na.strings = TRUE)
sesame_phenotypes_2018 <- subset(sesame_phenotypes, Year=='2018')
sesame_phenotypes_2018_Plant_height <- data.frame(cbind(sesame_phenotypes_2018$Genotype, sesame_phenotypes_2018$Plant.height))
colnames ( sesame_phenotypes_2018_Plant_height ) <- c ("Genotype", "Plant.height")
sesame_phenotypes_2018_Plant_height$Genotype <- as.numeric(gsub("S-","", sesame_phenotypes_2018_Plant_height$Genotype))
sesame_phenotypes_2018_Plant_height$Plant.height <-  as.numeric(sesame_phenotypes_2018_Plant_height$Plant.height)
sesame_phenotypes_2018_Plant_height_means <- sesame_phenotypes_2018_Plant_height %>% group_by(Genotype) %>% summarise(mean.Plant.height = sum(Plant.height, na.rm=TRUE)/7)
sesame_phenotypes_2018_Plant_height_means$Genotype <- paste("S", sesame_phenotypes_2018_Plant_height_means$Genotype, sep="-")
sesame_phenotypes_2018_Plant_height$Genotype <- as.numeric(sesame_phenotypes_2018_Plant_height$Genotype)
```

