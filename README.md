# Optimizing the search for significant SNPs <img width="458" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/4df1ba44-5a7d-46e0-bb75-4eb892320d8c">

**Students:**
>Liliia Bogdanova([github](https://github.com/LiliiaBgdnv), [telegram](http://t.me/bt_despair_and_hope))
>
>Maria Molodova (github, telegram)

**Supervisors:**
> Elena Grigoreva, NOVA PLANT
> 
> Lavrentiy Danilov, NOVA PLANT

## Introduction

The identification of the genetic loci associated with agronomically important phenotypes is a crucial task in plant selection, since it helps to optimize breeding strategies and accelerate crop breeding [1]. Genome-wide association study (GWAS) is a powerful tool commonly used to identify genetic variants responsible for the phenotype of interest. Most GWAS models are developed and optimized using human datasets; however, they have limited application in plant studies, especially when mapping complex traits such as stress tolerance and yield [2]. GWAS in plant species should account for such specific properties of plants as polyploidy, a complex genome, and a complex population structure. 

In this work, we aim to analyze the existing approaches for GWAS analysis in plants, to test them on real data of plant genotypes and phenotypes, and to select the optimal approaches.

1.  Brachi B. Genome-wide association studies in plants: the missing heritability is in the field / B. Brachi, G.P. Morris, J.O. Borevitz // Genome Biology. – 2011. – Vol. 12. – Genome-wide association studies in plants. – № 10. – P. 232.
2.  Prioritized candidate causal haplotype blocks in plant genome-wide association studies / X. Wu [et al.] // PLOS Genetics. – 2022. – Vol. 18. – № 10. – P. e1010437.

## Running tools

To get the data and code, clone the git repository:

```ruby
git clone https://github.com/LiliiaBgdnv/GWAS_project.git
cd GWAS_project
```

All **.bed, .bim, .fam** files are made using the program [plink2](https://www.cog-genomics.org/plink/2.0/). You can install it on the website.

### [GAPIT](https://zzlab.net/GAPIT/) (FarmCPU, MLM, Blink)
All data for the methods are in the "[GAPIT](https://github.com/LiliiaBgdnv/GWAS_project/tree/main/GAPIT)" folder, and there is also an [**.rmd** file](https://github.com/LiliiaBgdnv/GWAS_project/blob/main/GAPIT/GAPIT.Rmd) with all the work of the three models.

### PLINK2 glm
All data is presented in the appropriate [folder](https://github.com/LiliiaBgdnv/GWAS_project/tree/main/plink_glm)

**For simple trait (command line):**
```ruby
plink2 --glm allow-no-covars --allow-extra-chr --bed /plink_glm/inputs/simple_trait.bed --bim /plink_glm/inputs/simple_trait.bim --fam /plink_glm/inputs/simple_trait.fam --pheno /plink_glm/inputs/phenotype_simple.txt --adjust cols='chrom','pos','alt','a1','ref','gc','fdrbh' --out /plink_glm/results/simple_glm_result_chr --covar-variance-standardize --freq --threads 32 --memory 100000
awk '{gsub(/^GLYMAchr_/,""); print}' plink_glm/results/simple_glm_result_chr.Leu.glm.linear.adjusted  > plink_glm/results/simple_glm_result_chr.Leu.glm.linear.adjusted
awk '{gsub(/^GLYMAchr_/,""); print}' plink_glm/results/simple_glm_result_chr.Leu.glm.linear > plink_glm/results/simple_glm_result_chr.Leu.glm.linear
```

**For complex trait (command line):**

```ruby
plink2 --glm allow-no-covars --allow-extra-chr --bed /plink_glm/inputs/complex_trait.bed --bim /plink_glm/inputs/complex_trait.bim --fam /plink_glm/inputs/complex_trait.fam --pheno /plink_glm/inputs/phenotypes_complex.tsv --adjust --out /plink_glm/results/complex_glm_result --covar-variance-standardize --freq --threads 32 --memory 100000
```

**For simple trait generated data (command line):**
```ruby
plink2 --glm allow-no-covars --allow-extra-chr --bed /plink_glm/inputs/simple_trait.bed --bim /plink_glm/inputs/simple_trait.bim --fam /plink_glm/inputs/simple_trait.fam --pheno /plink_glm/inputs/gen_phenotype_simple.txt --adjust cols='chrom','pos','alt','a1','ref','gc','fdrbh' --out /plink_glm/results/simple_glm_gen_result_chr --covar-variance-standardize --freq --threads 32 --memory 100000
awk '{gsub(/^GLYMAchr_/,""); print}' plink_glm/results/simple_glm_gen_result_chr.Leu.glm.linear.adjusted  > plink_glm/results/simple_glm_gen_result_chr.Leu.glm.linear.adjusted
awk '{gsub(/^GLYMAchr_/,""); print}' plink_glm/results/simple_glm_gen_result_chr.Leu.glm.linear > plink_glm/results/simple_glm_gen_result_chr.Leu.glm.linear
```

**For complex trait generated data (command line):**

```ruby
plink2 --glm allow-no-covars --allow-extra-chr --bed /plink_glm/inputs/complex_trait.bed --bim /plink_glm/inputs/complex_trait.bim --fam /plink_glm/inputs/complex_trait.fam --pheno /plink_glm/inputs/gen_phenotypes_complex.tsv --adjust --out /plink_glm/results/complex_glm_gen_result --covar-variance-standardize --freq --threads 32 --memory 100000
```

### [GEMMA](https://github.com/genetics-statistics/GEMMA)

Installation:
```ruby
wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz
``
### BSLMM

**For simple trait (command line):**
```ruby
./gemma -bfile inplut/simple_trait -bslmm 1 -o gemma_bslmm_output_simple
```

**For complex trait (command line):**

```ruby
./gemma -bfile inplut/complex_trait -bslmm 1 -o gemma_bslmm_output_complex
```

**For simple trait generated data (command line):**
```ruby
./gemma -bfile inplut/simple_generated -bslmm 1 -o gemma_bslmm_output_simple_gen
```

**For complex trait generated data (command line):**

```ruby
./gemma -bfile inplut/complex_generated -bslmm 1 -o gemma_bslmm_output_complex_gen
```
## Analysis

The [Jupiter notebook](https://github.com/LiliiaBgdnv/GWAS_project/blob/main/Benjamini_Yekutieli.ipynb) contains the code for calculating the Benjamini-Yekutieli correction and plotting the Manhattan and QQ plots for GAPIT tools and GLM.
For the generated data, the graphics can be found in the Jupiter notebook.
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/bc15a182-bf22-4bb1-82fd-bd864ec3c068)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/9220e670-5565-4468-b0fc-3da522892f02)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/42db7742-a88d-4f0b-9d88-129143b2b8b3)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/81c953f6-eeb2-41a8-a2a4-d43070776873)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/5cdbbf64-f639-4093-b518-19e3575ecb71)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/1afa407c-768d-49f6-b393-43844a587fcd)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/2c21ea70-23ab-43ee-9f0b-c9fe378fe463)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/60192c43-3ffa-45a4-b4db-fb7c17618897)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/e3a2fc4d-1d14-4865-8c95-85d53f2a94d0)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/44103b0a-b004-4ed6-a0ae-06698f5680ab)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/248cf1c7-1b5c-42eb-8faf-5e46ca27c8d5)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/5736ece9-5d33-4181-b44c-4c290650419e)


The processing of the results of the BSLMM model differed from other models, because the result of the Bayesian statistics does not contain the p-value calculated for each SNP, but rather contains posterior samples of parameters: random effects (alpha), fixed effects (beta), and sparse effects (gamma). Based on the gamma values for each file, the effect sizes for each SNP were calculated and the posterior probability was calculated. The posterior inclusion probability (PIP) was calculated as a measure of the strength of the association between the SNP and the phenotype. Based on the PIP, the distribution of linked SNPs across the genome was visualized. The size of the dots reflects the magnitude of the effect.

The results of the processing are in the [folder](https://github.com/LiliiaBgdnv/GWAS_project/tree/main/GEMMA_BSLMM/results), [script](https://github.com/LiliiaBgdnv/GWAS_project/blob/main/GEMMA_BSLMM/visualization.R) by processing.

**Simple trait**
<img width="1388" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/23a2986c-5f5e-4d15-b1a5-744dcac58b03">

**Complex trait**
<img width="1388" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/b4020276-a790-4433-9291-8776beef578b">

As you can see, with a threshold of **0.05 only 1 SNPs were found for the complex trait, and no SNPs were found at all for the simple trait**, so we had to lower the threshold to 0.02 for the rest of the analysis in order to compare the SNPs found by each model.

**Simple trait**
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/b19e31bc-b854-4569-a280-54f7971d8bb3)

**Complex trait**
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/ea6cfe92-ba6d-4597-a73f-f55d16a023b8)


### Time for all model in seconds:
| Model      | Dataset 1 (simple trait) | Dataset 2 (complex trait) |
|------------|--------------------------|---------------------------|
| FarmCPU    | 26.54                    | 39.36                     |
| BLINK      | 22.34                    | 33.59                     |
| SUPER      | 100.06                   | 126.01                    |
| Plink2 glm | 2.689                    | 0.277                     |
| BSLMM      | 189261                   | 250576.8                  |
| MLMM       | 1.68                     | 6.875                     |


### Number of intersections of the significant SNPs identified by the tools.
<img width="615" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/ff5b4c62-0783-4c15-bb95-fe8faab181e6">
