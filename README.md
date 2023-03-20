# GWAS_project

Time for all model:

|                   | FarmCPU| BLINK  | SUPER  | Plink2 glm |
|-------------------|--------|--------|------------|--------|
| Sesame            | 39.36s | 33.59s | 126.01s    | 0.277s |
| Sesame generated  | 51.53s | 53.19s | 515.22s    | 1.578s |
| Soybean           | 26.54s | 22.34s | 100.06s    | 2.689s |
| Soybean generated | 22.56s | 22.21s | 100.75s    | 0.952s |

# GAPIT (FarmCPU, MLM, Blink)

**All the data are in the folders, but in zip format**

## First step: prepare the data for use by the model.

Conversion of genotype data from vcf format to hapmap was performed in the [TASSEL 5](https://tassel.bitbucket.io/) program. 
The conversion from hapmap to numerical was done with the R package [GAPIT3](https://zzlab.net/GAPIT/gapit_help_document.pdf).
```ruby
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
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

#Process phenotype file
myY <- read.table("processed_soy_phenotypes_leu.tsv", sep = '\t', header = TRUE)
write.table(myY, file='processed_soy_phenotypes_leu.txt', quote=FALSE, sep='\t', row.names = F)

#Process hapmap file: manually deleted two "#" symbols in the header
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

## Soybean FarmCPU model 

```ruby
setwd('C:/../../../FarmCPU_soybean_2')
myY  <- read.table("C://..//..//..//..//processed_soy_phenotypes_leu.txt", sep = '\t', head = TRUE)
myG <- read.delim("C://..//..//..//..//soybean_simple_genotypes.hmp.txt", head = FALSE)
myG[,1]=c("rs",1:(nrow(myG)-1))
system.time(myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=3, 
  model="FarmCPU",
  SNP.FDR = 0.05
))
```
![image](https://user-images.githubusercontent.com/109213422/224039902-3e6f7eba-a96f-44d8-880e-c54c52ea2f15.png)

| SNP                   | Chr         | Pos      | P.value              | MAF               | nobs | H&B.P.Value        | Effect            |
|-----------------------|-------------|----------|----------------------|-------------------|------|--------------------|-------------------|
| SGLYMACHR_18_61818402 | GLYMACHR_18 | 61818402 | 1.88033730003591e-06 | 0.5               | 97   | 0.0069391406918542 | 1.83250994205456  |
| SGLYMACHR_18_61819070 | GLYMACHR_18 | 61819070 | 2.38478930899706e-06 | 0.494845360824742 | 97   | 0.0069391406918542 | -1.8435566021963  |
| SGLYMACHR_18_61846089 | GLYMACHR_18 | 61846089 | 2.38478930899706e-06 | 0.494845360824742 | 97   | 0.0069391406918542 | -1.8435566021963  |
| SGLYMACHR_18_61846097 | GLYMACHR_18 | 61846097 | 2.38478930899706e-06 | 0.494845360824742 | 97   | 0.0069391406918542 | -1.8435566021963  |
| SGLYMACHR_18_61846199 | GLYMACHR_18 | 61846199 | 2.38478930899706e-06 | 0.494845360824742 | 97   | 0.0069391406918542 | -1.8435566021963  |
| SGLYMACHR_18_61846240 | GLYMACHR_18 | 61846240 | 2.38478930899509e-06 | 0.494845360824742 | 97   | 0.0069391406918542 | 1.84355660219631  |
| SGLYMACHR_18_61846255 | GLYMACHR_18 | 61846255 | 2.38478930899509e-06 | 0.494845360824742 | 97   | 0.0069391406918542 | 1.84355660219631  |
| SGLYMACHR_18_61846357 | GLYMACHR_18 | 61846357 | 2.38478930899509e-06 | 0.494845360824742 | 97   | 0.0069391406918542 | 1.84355660219631  |
| SGLYMACHR_18_61832336 | GLYMACHR_18 | 61832336 | 9.06381141379862e-06 | 0.494845360824742 | 97   | 0.0234430446767116 | 1.74061809163047  |
| SGLYMACHR_20_18390575 | GLYMACHR_20 | 18390575 | 2.07679479332065e-05 | 0.123711340206186 | 97   | 0.0483436291989182 | -2.47779977697192 |
| SGLYMACHR_18_61915849 | GLYMACHR_18 | 61915849 | 3.52642338719254e-05 | 0.438144329896907 | 97   | 0.0746255305518799 | -1.48761165405494 |

## Soybean SUPER model 

```ruby
setwd('C:/../../../SUPER_soybean_2')
myY  <- read.table("C://..//..//..//..//processed_soy_phenotypes_leu.txt", sep = '\t', head = TRUE)
myG <- read.delim("C://..//..//..//..//soybean_simple_genotypes.hmp.txt", head = FALSE)
system.time(myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=3, 
  model="SUPER",
  SNP.FDR = 0.05
))
```
![image](https://user-images.githubusercontent.com/109213422/224086028-9331de96-40bd-436f-ac26-5a124a0665e6.png)

| SNP                   | Chr         | Pos      | P.value              | MAF               | nobs | H&B.P.Value        | Effect |
|-----------------------|-------------|----------|----------------------|-------------------|------|--------------------|--------|
| SGLYMACHR_18_61818402 | GLYMACHR_18 | 61818402 | 2.71734020996048e-05 | 0.5               | 92   | 0.0899196719844409 | NA     |
| SGLYMACHR_18_61819070 | GLYMACHR_18 | 61819070 | 3.09028858095853e-05 | 0.494845360824742 | 92   | 0.0899196719844409 | NA     |
| SGLYMACHR_18_61846089 | GLYMACHR_18 | 61846089 | 3.09028858095853e-05 | 0.494845360824742 | 92   | 0.0899196719844409 | NA     |
| SGLYMACHR_18_61846097 | GLYMACHR_18 | 61846097 | 3.09028858095853e-05 | 0.494845360824742 | 92   | 0.0899196719844409 | NA     |
| SGLYMACHR_18_61846199 | GLYMACHR_18 | 61846199 | 3.09028858095853e-05 | 0.494845360824742 | 92   | 0.0899196719844409 | NA     |
| SGLYMACHR_18_61846240 | GLYMACHR_18 | 61846240 | 3.09028858095796e-05 | 0.494845360824742 | 92   | 0.0899196719844409 | NA     |
| SGLYMACHR_18_61846255 | GLYMACHR_18 | 61846255 | 3.09028858095796e-05 | 0.494845360824742 | 92   | 0.0899196719844409 | NA     |
| SGLYMACHR_18_61846357 | GLYMACHR_18 | 61846357 | 3.09028858095796e-05 | 0.494845360824742 | 92   | 0.0899196719844409 | NA     |

## Soybean BLINK model 

```ruby
setwd('C:/../../../BLINK_soybean_2')
myY  <- read.table("C://..//..//..//..//processed_soy_phenotypes_leu.txt", sep = '\t', head = TRUE)
myG <- read.delim("C://..//..//..//..//soybean_simple_genotypes.hmp.txt", head = FALSE)
system.time(myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=3, 
  model="BLINK",
  SNP.FDR = 0.05
))
```
![image](https://user-images.githubusercontent.com/109213422/224045617-98c7c3d6-11da-4b5f-909e-cb281e018007.png)

| SNP                   | Chr         | Pos      | P.value              | MAF               | nobs | H&B.P.Value         | Effect            |
|-----------------------|-------------|----------|----------------------|-------------------|------|---------------------|-------------------|
| SGLYMACHR_18_61818402 | GLYMACHR_18 | 61818402 | 1.87561573471076e-06 | 0.5               | 97   | 0.00692203007624291 | 1.83250994205456  |
| SGLYMACHR_18_61819070 | GLYMACHR_18 | 61819070 | 2.37890886716828e-06 | 0.494845360824742 | 97   | 0.00692203007624291 | -1.8435566021963  |
| SGLYMACHR_18_61846089 | GLYMACHR_18 | 61846089 | 2.37890886716828e-06 | 0.494845360824742 | 97   | 0.00692203007624291 | -1.8435566021963  |
| SGLYMACHR_18_61846097 | GLYMACHR_18 | 61846097 | 2.37890886716828e-06 | 0.494845360824742 | 97   | 0.00692203007624291 | -1.8435566021963  |
| SGLYMACHR_18_61846199 | GLYMACHR_18 | 61846199 | 2.37890886716828e-06 | 0.494845360824742 | 97   | 0.00692203007624291 | -1.8435566021963  |
| SGLYMACHR_18_61846240 | GLYMACHR_18 | 61846240 | 2.37890886716632e-06 | 0.494845360824742 | 97   | 0.00692203007624291 | 1.84355660219631  |
| SGLYMACHR_18_61846255 | GLYMACHR_18 | 61846255 | 2.37890886716632e-06 | 0.494845360824742 | 97   | 0.00692203007624291 | 1.84355660219631  |
| SGLYMACHR_18_61846357 | GLYMACHR_18 | 61846357 | 2.37890886716632e-06 | 0.494845360824742 | 97   | 0.00692203007624291 | 1.84355660219631  |
| SGLYMACHR_18_61832336 | GLYMACHR_18 | 61832336 | 9.04376415752922e-06 | 0.494845360824742 | 97   | 0.0233911935621072  | 1.74061809163047  |
| SGLYMACHR_20_18390575 | GLYMACHR_20 | 18390575 | 2.07252902319826e-05 | 0.123711340206186 | 97   | 0.0482443306020092  | -2.47779977697192 |
| SGLYMACHR_18_61915849 | GLYMACHR_18 | 61915849 | 3.5195353949114e-05  | 0.438144329896907 | 97   | 0.0744797681115887  | -1.48761165405494 |

## Sesame FarmCPU model 

```ruby
myY  <- read.table("C://..//..//..//..//processed_sesame_phenotypes.txt", sep = '\t', head = TRUE)
myG <- read.delim("C://..//..//..//..//sesame_complicated_genotypes.hmp.txt", head = FALSE)
```
```ruby
setwd('C:/../../../FarmCPU_sesame_2')
system.time(myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=3, 
  model="FarmCPU",
  SNP.FDR = 0.05
))
```
![image](https://user-images.githubusercontent.com/109213422/226430317-1d86c0fe-dfb7-4306-9d70-b36546f47514.png)

| SNP   | Chr| Pos      | P.value              | MAF               |nobs| H&B.P.Value          | Effect            |
|-------|----|----------|----------------------|-------------------|----|----------------------|-------------------|
| 4414  | 3  | 45355483 | 2.58248046574916e-08 | 0.139175257731959 | 97 | 0.000409013256165352 | 103.164182890069  |
| 17486 | 12 | 39996708 | 4.03720395661647e-08 | 0.164948453608247 | 97 | 0.000426274908432611 | -72.7242386983184 |
| 30544 | 20 | 36271980 | 1.50995533501268e-07 | 0.474226804123711 | 97 | 0.00119573362979654  | 37.6661411116701  |
| 18909 | 13 | 39682867 | 2.30843202806337e-07 | 0.139175257731959 | 97 | 0.0014624378584187   | -72.2674495288569 |
| 2255  | 2  | 36639483 | 9.49534236032572e-06 | 0.185567010309278 | 97 | 0.0420920837856796   | 71.4424426947608  |
| 6435  | 4  | 49624157 | 9.51982159034987e-06 | 0.139175257731959 | 97 | 0.0420920837856796   | 56.9102746665879  |
| 29920 | 20 | 2039705  | 1.06306563418815e-05 | 0.293814432989691 | 97 | 0.0420920837856796   | 38.4965768677324  |

## Sesame SUPER model 

```ruby
setwd('C:/../../../SUPER_sesame_2')
system.time(myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=3, 
  model="SUPER",
  SNP.FDR = 0.05
))
```
![image](https://user-images.githubusercontent.com/109213422/226431303-7cf5ecad-796b-48de-95a5-904643c4e75a.png)

| SNP   | Chr| Pos      | P.value              | MAF                |nobs| H&B.P.Value          | 
|-------|----|----------|----------------------|--------------------|----|----------------------|
| 6779  | 5  | 4355701  | 4.48713308168275e-08 | 0.128865979381443  | 92 | 0.000448121012182729 |
| 10209 | 7  | 17902694 | 1.1317616168272e-07  | 0.108247422680412  | 92 | 0.000448121012182729 |
| 10210 | 7  | 17902749 | 1.1317616168272e-07  | 0.108247422680412  | 92 | 0.000448121012182729 |
| 10211 | 7  | 17902783 | 1.13176161682527e-07 | 0.108247422680412  | 92 | 0.000448121012182729 |
| 15396 | 10 | 40282492 | 1.0414651577763e-07  | 0.0824742268041238 | 92 | 0.000448121012182729 |
| 15403 | 10 | 40347104 | 7.61246075250568e-08 | 0.0773195876288659 | 92 | 0.000448121012182729 |
| 15404 | 10 | 40378675 | 1.1269186323058e-07  | 0.0721649484536082 | 92 | 0.000448121012182729 |
| 16483 | 11 | 33243330 | 2.18711202178033e-08 | 0.11340206185567   | 92 | 0.000448121012182729 |
| 6803  | 5  | 4684569  | 2.14014740636511e-07 | 0.154639175257732  | 92 | 0.000677913092440211 |
| 29653 | 19 | 47795676 | 2.09679193515062e-07 | 0.103092783505155  | 92 | 0.000677913092440211 |
| 15394 | 10 | 40256295 | 3.00669925736352e-07 | 0.0979381443298969 | 92 | 0.000865820051602244 |
| 29654 | 19 | 47796136 | 3.40319629354327e-07 | 0.108247422680412  | 92 | 0.000898330381618971 |
| 3026  | 3  | 4891054  | 4.33599195226515e-07 | 0.0979381443298969 | 92 | 0.000900575291925032 |
| 3027  | 3  | 4891056  | 4.33599195226515e-07 | 0.0979381443298969 | 92 | 0.000900575291925032 |
| 4414  | 3  | 45355483 | 5.11755122321334e-07 | 0.139175257731959  | 92 | 0.000900575291925032 |
| 7398  | 5  | 35824337 | 5.06838199221747e-07 | 0.108247422680412  | 92 | 0.000900575291925032 |
| 10304 | 7  | 23283363 | 5.04056014627574e-07 | 0.170103092783505  | 92 | 0.000900575291925032 |
| 15410 | 10 | 40511307 | 4.11417040454048e-07 | 0.139175257731959  | 92 | 0.000900575291925032 |
| 10396 | 7  | 27170909 | 5.60100721874187e-07 | 0.123711340206186  | 92 | 0.000933776340320356 |
| 10194 | 7  | 17555663 | 8.85125435460452e-07 | 0.180412371134021  | 92 | 0.00127441969516569  |
| 10195 | 7  | 17555788 | 8.85125435460359e-07 | 0.180412371134021  | 92 | 0.00127441969516569  |
| 26075 | 18 | 3591942  | 8.48421482431062e-07 | 0.154639175257732  | 92 | 0.00127441969516569  |
| 29324 | 19 | 39915150 | 1.09260814291269e-06 | 0.123711340206186  | 92 | 0.00150475893630011  |
| 13928 | 10 | 571336   | 1.51980805988163e-06 | 0.118556701030928  | 92 | 0.00200589333770044  |
| 13890 | 9  | 49775192 | 2.71094726686298e-06 | 0.123711340206186  | 92 | 0.00343487862500607  |
| 29213 | 19 | 37737773 | 3.64608318660933e-06 | 0.118556701030928  | 92 | 0.00444205119303989  |
| 6801  | 5  | 4684401  | 6.36045254128675e-06 | 0.144329896907217  | 92 | 0.00719548909635032  |
| 6802  | 5  | 4684477  | 6.36045254128706e-06 | 0.144329896907216  | 92 | 0.00719548909635032  |
| 6739  | 5  | 3679594  | 9.60777758309513e-06 | 0.0979381443298969 | 92 | 0.00980921091173332  |
| 10192 | 7  | 17555583 | 9.44905577149218e-06 | 0.154639175257732  | 92 | 0.00980921091173332  |
| 10193 | 7  | 17555591 | 9.44905577149218e-06 | 0.154639175257732  | 92 | 0.00980921091173332  |
| 13926 | 10 | 567658   | 9.90954505541945e-06 | 0.103092783505155  | 92 | 0.00980921091173332  |
| 29337 | 19 | 39978117 | 1.05706754820938e-05 | 0.0618556701030928 | 92 | 0.0101465671688123   |
| 6432  | 4  | 49620427 | 1.17942159545304e-05 | 0.11340206185567   | 92 | 0.0106741024164487   |
| 18909 | 13 | 39682867 | 1.17370229836713e-05 | 0.139175257731959  | 92 | 0.0106741024164487   |
| 20628 | 14 | 47961602 | 1.36012883269268e-05 | 0.134020618556701  | 92 | 0.011967622473437    |
| 6800  | 5  | 4684378  | 1.4503353964599e-05  | 0.139175257731959  | 92 | 0.0124164389238551   |
| 6433  | 4  | 49620465 | 1.60893149787009e-05 | 0.108247422680412  | 92 | 0.0134117142438245   |
| 8096  | 6  | 10080487 | 1.72402456493851e-05 | 0.237113402061856  | 92 | 0.0136525505297481   |
| 8097  | 6  | 10080514 | 1.72402456493851e-05 | 0.237113402061856  | 92 | 0.0136525505297481   |
| 14897 | 10 | 23267088 | 2.05228792007225e-05 | 0.11340206185567   | 92 | 0.0158556761356607   |
| 10596 | 7  | 36429122 | 2.61933999528336e-05 | 0.0979381443298969 | 92 | 0.0192953985326967   |
| 10597 | 7  | 36429129 | 2.61933999528336e-05 | 0.0979381443298969 | 92 | 0.0192953985326967   |
| 17486 | 12 | 39996708 | 2.82681800546619e-05 | 0.164948453608247  | 92 | 0.0203505198048062   |
| 5993  | 4  | 38758952 | 3.09011786743085e-05 | 0.15979381443299   | 92 | 0.0203922028268208   |
| 5994  | 4  | 38758986 | 3.09011786743085e-05 | 0.15979381443299   | 92 | 0.0203922028268208   |
| 5995  | 4  | 38759005 | 3.09011786742914e-05 | 0.15979381443299   | 92 | 0.0203922028268208   |
| 5997  | 4  | 38759009 | 3.09011786742914e-05 | 0.15979381443299   | 92 | 0.0203922028268208   |
| 2255  | 2  | 36639483 | 5.75373819352543e-05 | 0.185567010309278  | 92 | 0.0371949818404309   |
| 6505  | 4  | 50849376 | 7.45048428062922e-05 | 0.45360824742268   | 92 | 0.0472003080146422   |
| 10629 | 7  | 37267885 | 7.78602713306452e-05 | 0.206185567010309  | 92 | 0.0483588618562651   |

## Sesame BLINK model 

```ruby
setwd('C:/../../../BLINK_sesame_2')
system.time(myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=3, 
  model="BLINK",
  SNP.FDR = 0.05
))
```
![image](https://user-images.githubusercontent.com/109213422/226431869-75e7f824-c363-4963-a511-86f12dcf789e.png)

| SNP   |Chr| Pos      | P.value             | MAF               |nobs| H&B.P.Value          | Effect           |
|-------|---|----------|---------------------|-------------------|----|----------------------|------------------|
| 13039 | 9 | 28490006 | 4.3864778038606e-08 | 0.221649484536082 | 97 | 0.000463153569716961 | 105.772035506385 |

# FaST-LMM
## First step: prepare the data (on the example of sesame)
### Sesame phenotypes preprocessing (R)
```ruby
sesame_phenotypes <- read.table('C:/Users/Roman/Desktop/IB/GWAS_project/data/sesame_phenotypes_2018_Plant_height_means.tsv ', header = TRUE)
sesame_phenotypes$V1 <- sesame_phenotypes$Genotype
sesame_phenotypes <- sesame_phenotypes[c("Genotype", "V1", "mean.Plant.height")]

library("dplyr")

sesame_phenotypes <- sesame_phenotypes %>% 
        rename("FID" = "Genotype",
               "IID" = "V1",
               "Height" = "mean.Plant.height")

           
write.table(sesame_phenotypes, file='C:/Users/Roman/Desktop/IB/GWAS_project/Fastlmm/sesame/raw_data/sesame_phenotypes_fast.txt', quote=FALSE,
            sep='\t', row.names = F)
```
### Sesame genotypes preprocessing (command line)
```ruby
sed 's/LG//g' sesame_complicated_genotypes_copy.vcf > sesame_complicated_genotypes_no_chrom.vcf
```

Create env (command line):
```ruby
python3.11 -m venv FaST-LMM
source FaST-LMM/bin/activate
```
Installing (command line):
```ruby
pip install fastlmm
```

Installing and using plink (to convert vcf fotmat in bed, bim, fam) (command line):
```ruby
conda install -c bioconda plink2
plink2 --vcf sesame_complicated_genotypes.vcf.gz --make-bed --allow-extra-chr --out sesame_complicated_genotypes
plink2 --vcf soybean_simple_genotypes.vcf.gz --make-bed --allow-extra-chr --out sesame_complicated_genotypes
```

# BOLT-LMM
**Prepare the data for use by the model (command line)**
```ruby
bcftools view --header-only soybean_simple_genotypes.vcf | sed 's/##contig=<ID=GLYMAchr_/##contig=<ID=/' | sed 's/##contig=<ID=scaffold_/##contig=<ID=/' > soybean_rename
awk '{gsub(/^GLYMAchr_/,""); print}' soybean_simple_genotypes.vcf > soybean_rename_chr.vcf
# Replace the description in soybean_rename_chr.vcf with soybean_rename
awk '{gsub(/^LG/,""); print}' sesame_complicated_genotypes.vcf > sesame_rename_chr.vcf
```

**Make bed, bim and fam file (command line)**

```ruby
plink2 --vcf ~/IB/GWAS_progect/Raw_data/soybean_rename_chr.vcf --make-bed --out soybean_rename_chs
./gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile soybean_rename_chs --ld-wind 5000 --ld-sig 0.05 --out soybean_ld
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' ../Raw_data/genotypes_sesame.vcf > genotypes_sesame_named.vcf
plink2 --vcf genotypes_sesame_named.vcf -pheno ~/IB/GWAS_progect/Raw_data/phenotypes_sesame.tsv --max-alleles 2 --make-bed --allow-extra-chr --out sesame_2
~/IB/GWAS_progect/BOLT_LMM/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile sesame_2 --ld-wind 5000 --ld-sig 0.05 --out sesame_ld
```
# PLINK2 glm

**For soybean (command line):**
```ruby
plink2 --glm allow-no-covars --allow-extra-chr --bed Raw_data/Bed_bim_fam/soybean_simpe_genotypes.bed --bim Raw_data/Bed_bim_fam/soybean_simpe_genotypes.bim --fam Raw_data/Bed_bim_fam/soybean_simpe_genotypes.fam --pheno Raw_data/soy_phenotypes_leu_plink_2col.txt --adjust cols='chrom','pos','alt','a1','ref','gc','fdrbh' --out plink_glm_2/soy_plink_result --covar-variance-standardize --freq --threads 32 --memory 100000
cd plink_glm_2
awk '{gsub(/^GLYMAchr_/,""); print}' soy_plink_result.Leu.glm.linear.adjusted  > soy_plink_result_chr.Leu.glm.linear.adjusted
awk '{gsub(/^GLYMAchr_/,""); print}' soy_plink_result.Leu.glm.linear > soy_plink_result_chr.Leu.glm.linear
```
*Python code:*
```ruby
# this section is enough to execute once
!pip3 install qmplot
import pandas as pd
import matplotlib.pyplot as plt
from qmplot import manhattanplot
```
```ruby
df_soy = pd.read_table('plink_glm_2/soy_plink_result_chr.Leu.glm.linear', sep="\t")
df_soy_ad = pd.read_table('plink_glm_2/soy_plink_result_chr.Leu.glm.linear.adjusted', sep="\t")
df_soy_merge = df_soy.merge(df_soy_ad, on=['#CHROM', 'POS', 'REF', 'ALT', 'A1'], how='left')
df_soy_merge = df_soy_merge.drop(['ID_y'], axis='columns')
df_soy_merge = df_soy_merge.rename({'ID_x':'ID'}, axis='columns')
df_soy_merge
```
Your table must look something like this

| #CHROM | POS | ID       | REF | ALT | A1  | TEST | OBS_CT | BETA | SE        | T_STAT   | P         | ERRCODE  | GC  | FDR_BH   |
|--------|-----|----------|-----|-----|-----|------|--------|------|-----------|----------|-----------|----------|-----|----------|
| 0      | 1   | 41579    | .   | T   | A   | A    | ADD    | 234  | 0.147398  | 0.219730 | 0.670813  | 0.503007 | .   | 0.693188 | 0.716952 |
| 1      | 1   | 69607    | .   | G   | A   | A    | ADD    | 237  | 0.867344  | 0.242347 | 3.578930  | 0.000419 | .   | 0.037701 | 0.025466 |
| 2      | 1   | 93302    | .   | A   | G   | G    | ADD    | 237  | 0.848699  | 0.246706 | 3.440120  | 0.000688 | .   | 0.045557 | 0.032772 |
| 3      | 1   | 123428   | .   | A   | G   | G    | ADD    | 238  | 0.659877  | 0.226817 | 2.909290  | 0.003968 | .   | 0.089719 | 0.076670 |
| 4      | 1   | 123445   | .   | A   | G   | G    | ADD    | 238  | 0.640807  | 0.223982 | 2.860970  | 0.004602 | .   | 0.095083 | 0.081543 |
| ...    | ... | ...      | ... | ... | ... | ...  | ...    | ...  | ...       | ...      | ...       | ...      | ... | ...      | ...      |
| 23273  | 20  | 46740364 | .   | C   | T   | T    | ADD    | 233  | 0.305600  | 0.212748 | 1.436440  | 0.152229 | .   | 0.399029 | 0.388883 |
| 23274  | 20  | 46740390 | .   | A   | G   | G    | ADD    | 234  | 0.332001  | 0.212400 | 1.563090  | 0.119394 | .   | 0.358956 | 0.341908 |
| 23275  | 20  | 46751278 | .   | A   | G   | G    | ADD    | 238  | 0.191220  | 0.229718 | 0.832409  | 0.406020 | .   | 0.624516 | 0.644343 |
| 23276  | 20  | 46751325 | .   | T   | C   | T    | ADD    | 238  | -0.627696 | 0.213608 | -2.938540 | 0.003625 | .   | 0.086596 | 0.073956 |
| 23277  | 20  | 46755301 | .   | T   | C   | C    | ADD    | 238  | -0.067125 | 0.239267 | -0.280544 | 0.779306 | .   | 0.868890 | 0.886079 |

**Plotting**

```ruby
if __name__ == "__main__":

    df_soy_merge = df_soy_merge.dropna(how="any", axis=0)  # clean data
    ax = manhattanplot(data=df_soy_merge, pv='FDR_BH', suggestiveline=1e-2)
    plt.savefig("output_manhattan_plot_soy.png", dpi=300)
```
**Here, the red line separates the SNP c p-value with the Benjamini-Hochberg FDR correction $-log_{10}{(5×10^{-2})}$**

![output_manhattan_plot_soy](https://user-images.githubusercontent.com/109213422/225274969-cb8a5c0d-757f-480f-9fdd-9f36baf76b67.png)

**For sesame:**

(command line)
```ruby
plink2 --glm allow-no-covars --allow-extra-chr --bed sesame_2.bed --bim sesame_2.bim --fam sesame_2.fam --pheno ~/IB/GWAS_progect/Raw_data/phenotypes_sesame.tsv --adjust --out plink_glm_2/sesame_plink_result --covar-variance-standardize --freq --threads 32 --memory 100000
```
**Python code:**

```ruby
df_sesame = pd.read_table('plink_glm_2/sesame_plink_result_chr.PH.glm.linear', sep="\t")
df_sesame_ad = pd.read_table('plink_glm_2/sesame_plink_result_chr.PH.glm.linear.adjusted', sep="\t")
df_sesame_merge = df_sesame.merge(df_sesame_ad, on='ID', how='left')
df_sesame_merge = df_sesame_merge.drop(['#CHROM_y', 'A1_y'], axis='columns')
df_sesame_merge = df_sesame_merge.rename({'#CHROM_x':'#CHROM', 'A1_x':'A1'}, axis='columns')
df_sesame_merge
```
Your table must look something like this

| #CHROM | POS      | ID     | REF               | ALT | A1_x | TEST | OBS_CT | BETA | SE        | ...     | P   | ERRCODE  | UNADJ | GC       | BONF     | HOLM | SIDAK_SS | SIDAK_SD | FDR_BH | FDR_BY   |
|--------|----------|--------|-------------------|-----|------|------|--------|------|-----------|---------|-----|----------|-------|----------|----------|------|----------|----------|--------|----------|
| 0      | 1        | 40522  | 1_40522_C_T       | C   | T    | T    | ADD    | 83   | -7.28654  | 39.2311 | ... | 0.853118 | .     | 0.853118 | 0.876121 | 1.0  | 1.0      | 1.0      | 1.0    | 0.969641 | 1.0 |
| 1      | 1        | 41592  | 1_41592_T_A       | T   | A    | A    | ADD    | 94   | -44.78850 | 38.9581 | ... | 0.253264 | .     | 0.253264 | 0.336078 | 1.0  | 1.0      | 1.0      | 1.0    | 0.777471 | 1.0 |
| 2      | 1        | 111641 | 1_111641_G_C      | G   | C    | C    | ADD    | 87   | 53.16680  | 27.0588 | ... | 0.052697 | .     | 0.052697 | 0.102832 | 1.0  | 1.0      | 1.0      | 1.0    | 0.593615 | 1.0 |
| 3      | 1        | 111675 | 1_111675_G_C      | G   | C    | C    | ADD    | 87   | 53.16680  | 27.0588 | ... | 0.052697 | .     | 0.052697 | 0.102832 | 1.0  | 1.0      | 1.0      | 1.0    | 0.593615 | 1.0 |
| 4      | 1        | 111688 | 1_111688_T_C      | T   | C    | C    | ADD    | 87   | 49.46560  | 26.7467 | ... | 0.067874 | .     | 0.067874 | 0.124208 | 1.0  | 1.0      | 1.0      | 1.0    | 0.619146 | 1.0 |
| ...    | ...      | ...    | ...               | ... | ...  | ...  | ...    | ...  | ...       | ...     | ... | ...      | ...   | ...      | ...      | ...  | ...      | ...      | ...    | ...      | ... |
| 31361  | KZ848124 | 2376   | KZ848124_2376_A_G | A   | G    | A    | ADD    | 97   | -25.54030 | 20.8585 | ... | 0.223810 | .     | 0.223810 | 0.305716 | 1.0  | 1.0      | 1.0      | 1.0    | 0.759581 | 1.0 |
| 31362  | KZ848124 | 2400   | KZ848124_2400_G_A | G   | A    | G    | ADD    | 97   | 2.97387   | 25.5964 | ... | 0.907753 | .     | 0.907753 | 0.922278 | 1.0  | 1.0      | 1.0      | 1.0    | 0.981726 | 1.0 |
| 31363  | KZ848124 | 2416   | KZ848124_2416_T_A | T   | A    | T    | ADD    | 90   | -26.24970 | 25.7602 | ... | 0.310995 | .     | 0.310995 | 0.393629 | 1.0  | 1.0      | 1.0      | 1.0    | 0.804824 | 1.0 |
| 31364  | KZ848124 | 2425   | KZ848124_2425_T_A | T   | A    | T    | ADD    | 97   | -5.23986  | 26.1919 | ... | 0.841864 | .     | 0.841864 | 0.866594 | 1.0  | 1.0      | 1.0      | 1.0    | 0.966697 | 1.0 |
| 31365  | KZ848124 | 2648   | KZ848124_2648_A_T | A   | T    | A    | ADD    | 76   | -26.17170 | 33.8941 | ... | 0.442477 | .     | 0.442477 | 0.517848 | 1.0  | 1.0      | 1.0      | 1.0    | 0.849949 | 1.0 |

```ruby
if __name__ == "__main__":

    df_sesame_merge = df_sesame_merge.dropna(how="any", axis=0)  # clean data
    ax = manhattanplot(data=df_sesame_merge, pv='FDR_BH')
    plt.savefig("output_manhattan_plot_sesame.png", dpi=300, suggestiveline=1e-2)
```

**Here, the red line separates the SNP c p-value with the Benjamini-Hochberg FDR correction $-log_{10}{(5×10^{-2})}$**

![image](https://user-images.githubusercontent.com/109213422/226433922-3430a3f5-4c43-4b88-aefc-60ffa73492a1.png)

**For soybean processed (command line):**
```ruby
plink2 --glm allow-no-covars --allow-extra-chr --bed Raw_data/Bed_bim_fam/soybean_simpe_genotypes.bed --bim Raw_data/Bed_bim_fam/soybean_simpe_genotypes.bim --fam Raw_data/Bed_bim_fam/soybean_simpe_genotypes.fam --pheno processed_gen_soy_phenotypes_leu.txt --adjust cols='chrom','pos','alt','a1','ref','gc','fdrbh' --out plink_glm_2/soy_processed_plink_result --covar-variance-standardize --freq --threads 32 --memory 100000
cd plink_glm_2
awk '{gsub(/^GLYMAchr_/,""); print}' soy_processed_plink_result.Leu.glm.linear.adjusted  > soy_processed_plink_result_chr.Leu.glm.linear.adjusted
awk '{gsub(/^GLYMAchr_/,""); print}' soy_processed_plink_result.Leu.glm.linear > soy_processed_plink_result_chr.Leu.glm.linear
```

*Python code:*

```ruby
#soy processed
df_soy_processed = pd.read_table('plink_glm_2/soy_processed_plink_result_chr.Leu.glm.linear', sep="\t")
df_soy_processed_ad = pd.read_table('plink_glm_2/soy_processed_plink_result_chr.Leu.glm.linear.adjusted', sep="\t")
df_soy_processed_merge = df_soy_processed.merge(df_soy_processed_ad, on=['#CHROM', 'POS', 'REF', 'ALT', 'A1'], how='left')
df_soy_processed_merge = df_soy_processed_merge.drop(['ID_y'], axis='columns')
df_soy_processed_merge = df_soy_processed_merge.rename({'ID_x':'ID'}, axis='columns')
df_soy_processed_merge
```

**Plotting**

```ruby
if __name__ == "__main__":

    df_soy_processed_merge = df_soy_processed_merge.dropna(how="any", axis=0)  # clean data
    ax = manhattanplot(data=df_soy_processed_merge, pv='FDR_BH', suggestiveline=1e-2)
    plt.savefig("output_manhattan_plot_soy_processed.png", dpi=300)
```
**Here, the red line separates the SNP c p-value with the Benjamini-Hochberg FDR correction $-log_{10}{(5×10^{-2})}$**

![output_manhattan_plot_soy_processed](https://user-images.githubusercontent.com/109213422/225308786-1af5f9b1-f994-4a04-a4ae-c370ce63d521.png)

**For sesame processed:**
(command line)
```ruby
plink2 --glm allow-no-covars --allow-extra-chr --bed Raw_data/Bed_bim_fam/sesame_complicated_genotypes.bed --bim Raw_data/Bed_bim_fam/sesame_complicated_genotypes.bim --fam Raw_data/Bed_bim_fam/sesame_complicated_genotypes.fam --pheno processed_gen_sesame_phenotypes.txt --adjust --out plink_glm_2/sesame_processed_plink_result --covar-variance-standardize --freq --threads 32 --memory 100000
cd plink_glm_2
awk '{gsub(/^LG/,""); print}' sesame_processed_plink_result.PH.glm.linear.adjusted > sesame_processed_plink_result_chr.PH.glm.linear.adjusted
awk '{gsub(/^LG/,""); print}' sesame_processed_plink_result.PH.glm.linear > sesame_processed_plink_result_chr.PH.glm.linear
```
**Python code:**

```ruby
#sesame processed
df_sesame_processed = pd.read_table('plink_glm_2/sesame_processed_plink_result_chr.PH.glm.linear', sep="\t")
df_sesame_processed_ad = pd.read_table('plink_glm_2/sesame_processed_plink_result_chr.PH.glm.linear.adjusted', sep="\t")
df_sesame_processed_merge = df_sesame_processed.merge(df_sesame_processed_ad, on='ID', how='left')
df_sesame_processed_merge = df_sesame_processed_merge.drop(['#CHROM_y', 'A1_y'], axis='columns')
df_sesame_processed_merge = df_sesame_processed_merge.rename({'#CHROM_x':'#CHROM', 'A1_x':'A1'}, axis='columns')
df_sesame_processed_merge
```

```ruby
if __name__ == "__main__":

    df_sesame_processed_merge = df_sesame_processed_merge.dropna(how="any", axis=0)  # clean data
    ax = manhattanplot(data=df_sesame_processed_merge, pv='FDR_BH', suggestiveline=1e-2)
    plt.savefig("output_manhattan_plot_sesame_processed.png", dpi=300)
```
![output_manhattan_plot_sesame_processed](https://user-images.githubusercontent.com/109213422/225309195-fbc6c982-2936-4b71-ac5c-be049dfbc382.png)

