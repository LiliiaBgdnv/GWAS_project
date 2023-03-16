# GWAS_project

Time for all model:

|                   | FarmCPU| BLINK  | SUPER  | Plink2 glm |
|-------------------|--------|--------|------------|--------|
| Sesame            | 22.48s | 22.78s | 523.51s    | 0.277s |
| Sesame generated  | 51.53s | 53.19s | 515.22s    | 1.578s |
| Soybean           | 26.54s | 22.34s | 100.06s    | 2.689s |
| Soybean generated | 22.56s | 22.21s | 100.75s    | 0.952s |

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
  model="BLINK",
  SNP.FDR = 0.05
))
```
![image](https://user-images.githubusercontent.com/109213422/224060054-f872017b-e449-469d-8a7d-01936889e1bc.png)

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
plink2 --vcf ~/IB/GWAS_progect/Raw_data/sesame_rename_chr.vcf --make-bed --out sesame_rename_chs
plink2 --vcf ~/IB/GWAS_progect/Raw_data/soybean_rename_chr.vcf --make-bed --out soybean_rename_chs

./gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile sesame_rename_chs --ld-wind 5000 --ld-sig 0.05 --out sesame_ld
./gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --bfile soybean_rename_chs --ld-wind 5000 --ld-sig 0.05 --out soybean_ld
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
plink2 --glm allow-no-covars --allow-extra-chr --bed Raw_data/Bed_bim_fam/sesame_complicated_genotypes.bed --bim Raw_data/Bed_bim_fam/sesame_complicated_genotypes.bim --fam Raw_data/Bed_bim_fam/sesame_complicated_genotypes.fam --pheno Raw_data/sesame_phenotypes_leu_BOLT.txt --adjust --out plink_glm_2/sesame_plink_result --covar-variance-standardize --freq --threads 32 --memory 100000
cd plink_glm_2
awk '{gsub(/^LG/,""); print}' sesame_plink_result.PH.glm.linear.adjusted > sesame_plink_result_chr.PH.glm.linear.adjusted
awk '{gsub(/^LG/,""); print}' sesame_plink_result.PH.glm.linear > sesame_plink_result_chr.PH.glm.linear
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

| #CHROM | POS | ID       | REF          | ALT | A1_x | TEST | OBS_CT | BETA | SE        | ...      | P   | ERRCODE  | UNADJ | GC       | BONF     | HOLM | SIDAK_SS | SIDAK_SD | FDR_BH | FDR_BY   |
|--------|-----|----------|--------------|-----|------|------|--------|------|-----------|----------|-----|----------|-------|----------|----------|------|----------|----------|--------|----------|
| 0      | 1   | 385      | 1:138:+      | T   | C    | C    | ADD    | 56   | -1.342850 | 1.299340 | ... | 0.305984 | .     | 0.305984 | 0.663770 | 1.0  | 1.0      | 1.0      | 1.0    | 0.471993 | 1.0 |
| 1      | 1   | 6852     | 28:36:-      | G   | A    | A    | ADD    | 71   | -0.854213 | 1.581950 | ... | 0.590954 | .     | 0.590954 | 0.819465 | 1.0  | 1.0      | 1.0      | 1.0    | 0.734041 | 1.0 |
| 2      | 1   | 8182     | 40:101:+     | G   | A    | A    | ADD    | 62   | -0.736148 | 1.074830 | ... | 0.496045 | .     | 0.496045 | 0.772524 | 1.0  | 1.0      | 1.0      | 1.0    | 0.655826 | 1.0 |
| 3      | 1   | 8237     | 45:36:+      | C   | A    | A    | ADD    | 77   | -0.737742 | 1.231620 | ... | 0.550978 | .     | 0.550978 | 0.800098 | 1.0  | 1.0      | 1.0      | 1.0    | 0.702149 | 1.0 |
| 4      | 1   | 8259     | 45:58:+      | C   | T    | T    | ADD    | 77   | -0.737742 | 1.231620 | ... | 0.550978 | .     | 0.550978 | 0.800098 | 1.0  | 1.0      | 1.0      | 1.0    | 0.702149 | 1.0 |
| ...    | ... | ...      | ...          | ... | ...  | ...  | ...    | ...  | ...       | ...      | ... | ...      | ...   | ...      | ...      | ...  | ...      | ...      | ...    | ...      | ... |
| 84142  | 9   | 12298737 | 529391:145:+ | C   | A    | A    | ADD    | 119  | -0.690236 | 0.996285 | ... | 0.489802 | .     | 0.489802 | 0.769310 | 1.0  | 1.0      | 1.0      | 1.0    | 0.650541 | 1.0 |
| 84143  | 9   | 12376846 | 529498:146:- | A   | G    | G    | ADD    | 42   | 3.198200  | 2.104780 | ... | 0.136505 | .     | 0.136505 | 0.527200 | 1.0  | 1.0      | 1.0      | 1.0    | 0.260864 | 1.0 |
| 84144  | 9   | 12377410 | 529502:145:- | A   | G    | G    | ADD    | 21   | 1.960410  | 5.943250 | ... | 0.745119 | .     | 0.745119 | 0.890202 | 1.0  | 1.0      | 1.0      | 1.0    | 0.847391 | 1.0 |
| 84145  | 9   | 12377781 | 529515:127:- | T   | C    | C    | ADD    | 32   | -0.522238 | 2.890830 | ... | 0.857854 | .     | 0.857854 | 0.939373 | 1.0  | 1.0      | 1.0      | 1.0    | 0.919861 | 1.0 |
| 84146  | 9   | 12377983 | 529512:79:+  | G   | T    | T    | ADD    | 21   | 2.864840  | 3.547890 | ... | 0.429385 | .     | 0.429385 | 0.737186 | 1.0  | 1.0      | 1.0      | 1.0    | 0.596036 | 1.0 |

```ruby
if __name__ == "__main__":

    df_sesame_merge = df_sesame_merge.dropna(how="any", axis=0)  # clean data
    ax = manhattanplot(data=df_sesame_merge, pv='FDR_BH')
    plt.savefig("output_manhattan_plot_sesame.png", dpi=300)
```

**The red line at the $-log_{10}{(10^{-5})}$ level for "putative" association and the green line $-log_{10}{(5×10^{-8})}$ for the "genome-wide significance" threshold (with the Benjamini-Hochberg FDR).**

![output_manhattan_plot_sesame](https://user-images.githubusercontent.com/109213422/225277372-b0a9b543-3cbc-4f17-ab57-b31b015ddd80.png)

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

