# GWAS_project

# GAPIT (FarmCPU, MLM, Blink)

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
![image](https://user-images.githubusercontent.com/109213422/224057912-16b45025-3a11-4a07-982f-4cbaaed3880f.png)

## Sesame SUPER model 

```ruby
myY  <- read.table("C://..//..//..//..//processed_sesame_phenotypes.txt", sep = '\t', head = TRUE)
myG <- read.delim("C://..//..//..//..//sesame_complicated_genotypes.hmp.txt", head = FALSE)
```
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
![image](https://user-images.githubusercontent.com/109213422/224059719-b53bef48-6978-49ab-a104-d4c7cbff72a5.png)

## Sesame BLINK model 

```ruby
myY  <- read.table("C://..//..//..//..//processed_sesame_phenotypes.txt", sep = '\t', head = TRUE)
myG <- read.delim("C://..//..//..//..//sesame_complicated_genotypes.hmp.txt", head = FALSE)
```
```ruby
setwd('C:/../../../BLINK_sesame_2')
system.time(myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=3, 
  model="BLONK",
  SNP.FDR = 0.05
))
```
![image](https://user-images.githubusercontent.com/109213422/224060054-f872017b-e449-469d-8a7d-01936889e1bc.png)

Time for all model:

|                   | FarmCPU | BLINK  | SUPER    |
|-------------------|---------|--------|----------|
| Sesame            | 22.48   | 22.78  | 523.51   |
| Sesame generated  | 51.53   | 53.19  | 515.22   |
| Soybean           | 26.54   | 22.34  |  100.06  |
| Soybean generated | 22.56   | 22.21  | 100.75   |

# FaST-LMM
Create env:
```ruby
python3.11 -m venv FaST-LMM
source FaST-LMM/bin/activate
```
Installing:
```ruby
pip install fastlmm
```

Installing and using plink (to convert vcf fotmat in bed, bim, fam):
```ruby
conda install -c bioconda plink2
plink2 --vcf sesame_complicated_genotypes.vcf.gz --make-bed --allow-extra-chr --out sesame_complicated_genotypes
plink2 --vcf soybean_simple_genotypes.vcf.gz --make-bed --allow-extra-chr --out sesame_complicated_genotypes
```

# BOLT-LMM
Prepare the data for use by the model
```ruby
bcftools view --header-only soybean_simple_genotypes.vcf | sed 's/##contig=<ID=GLYMAchr_/##contig=<ID=/' | sed 's/##contig=<ID=scaffold_/##contig=<ID=/' > soybean_rename
awk '{gsub(/^GLYMAchr_/,""); print}' soybean_simple_genotypes.vcf > soybean_rename_chr.vcf
# Replace the description in soybean_rename_chr.vcf with soybean_rename
awk '{gsub(/^LG/,""); print}' sesame_complicated_genotypes.vcf > sesame_rename_chr.vcf
```

Make bed, bim and fam file

```ruby
plink2 --vcf ~/IB/GWAS_progect/Raw_data/sesame_rename_chr.vcf --make-bed --out sesame_rename_chs
plink2 --vcf ~/IB/GWAS_progect/Raw_data/soybean_rename_chr.vcf --make-bed --out soybean_rename_chs

./gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile sesame_rename_chs --ld-wind 5000 --ld-sig 0.05 --out sesame_ld
./gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile soybean_rename_chs --ld-wind 5000 --ld-sig 0.05 --out soybean_ld
```
# PLINK2 glm

for soybean:
```ruby
plink2 --glm allow-no-covars --allow-extra-chr --bed Raw_data/Bed_bim_fam/soybean_simpe_genotypes.bed --bim Raw_data/Bed_bim_fam/soybean_simpe_genotypes.bim --fam Raw_data/Bed_bim_fam/soybean_simpe_genotypes.fam --pheno Raw_data/soy_phenotypes_leu_plink_2col.txt --out plink_result --mac 20 --covar-variance-standardize --freq --threads 32 --memory 100000
 ```
 
for sesame:
```ruby
plink2 --glm allow-no-covars --allow-extra-chr --bed Raw_data/Bed_bim_fam/sesame_complicated_genotypes.bed --bim Raw_data/Bed_bim_fam/sesame_complicated_genotypes.bim --fam Raw_data/Bed_bim_fam/sesame_complicated_genotypes.fam --pheno Raw_data/sesame_phenotypes_leu_BOLT.txt --out plink2_glm/sesame_plink_result --mac 20 --covar-variance-standardize --freq --threads 32 --memory 100000
```
