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



## Running tools

All **.bed, .bim, .fam** files are made using the program [plink2](https://www.cog-genomics.org/plink/2.0/). You can install it on the website.

### [GAPIT](https://zzlab.net/GAPIT/) (FarmCPU, MLM, Blink)
All data for the methods are in the "[GAPIT](https://github.com/LiliiaBgdnv/GWAS_project/tree/main/GAPIT)" folder, and there is also an **.rmd** file with all the work of the three models.

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

**BSLMM**

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

Time for all model:

|                   | FarmCPU| BLINK  | SUPER  | Plink2 glm |
|-------------------|--------|--------|------------|--------|
| Sesame            | 39.36s | 33.59s | 126.01s    | 0.277s |
| Sesame generated  | 24.53 s | 21.90s | 126.31s    | 1.578s |
| Soybean           | 26.54s | 22.34s | 100.06s    | 2.689s |
| Soybean generated | 22.56s | 22.21s | 100.75s    | 0.952s |
