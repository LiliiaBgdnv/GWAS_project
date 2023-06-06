# Optimizing the search for significant SNPs <img width="458" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/4df1ba44-5a7d-46e0-bb75-4eb892320d8c">

**Students:**
>Liliia Bogdanova([github](https://github.com/LiliiaBgdnv), [telegram](http://t.me/bt_despair_and_hope))
>
>Maria Molodova ([github](https://github.com/MaryM12), [telegram](http://t.me/maria_molodova))

**Supervisors:**
> Elena Grigoreva, NOVA PLANT
> 
> Lavrentiy Danilov, NOVA PLANT

## Introduction

The identification of the genetic loci associated with agronomically important phenotypes is a crucial task in plant selection, since it helps to optimize breeding strategies and accelerate crop breeding [1]. Genome-wide association study (GWAS) is a powerful tool commonly used to identify genetic variants responsible for the phenotype of interest. Most GWAS models are developed and optimized using human datasets; however, they have limited application in plant studies, especially when mapping complex traits such as stress tolerance and yield [2]. GWAS in plant species should account for such specific properties of plants as polyploidy, a complex genome, and a complex population structure. 

In this work, we aim to analyze the existing approaches for GWAS analysis in plants, to test them on real data of plant genotypes and phenotypes, and to select the optimal approaches.

The 6 methods are tested on **2 datasets for soybeans (Glycine max)**, on a simple and a complex trait. Simple traits develop as a result of some type of interaction between a limited number of genes, as we assume. Complex traits develop as a result of the interaction of a large number of genes involved in many pathways.

**A simple trait is the content of the amino acid leucine in beans,the data for working with a simple sign was taken from the [article](https://www.frontiersin.org/articles/10.3389/fpls.2019.01445/full). The complex trait is a trade secret :)**

Preparation of the data for the simple trait: we selected the leucine content for the work, so for further work it is necessary to select the column 'Leu' from the phenotype data for the simple trait.
```ruby
cut -f 1,2,10 ./Raw_data/phenotypes_simple_trait_full.tsv > ./Raw_data/phenotypes_simple_trait.tsv
```

Also to check for false positives, files were generated for the models in which these phenotypes were represented as a normal distribution. The script for the generation is presented in this [file](https://github.com/LiliiaBgdnv/GWAS_project/blob/main/GAPIT/GAPIT.Rmd).

**We have analyzed the literature and selected six models for comparison:**

**BLINK** (Bayesian Linkage and Association with IBD-sharing) is a method that uses Bayesian models to determine the relationship between genotype and phenotype, considering shared inheritance. It can be used for both association and linkage analysis. uses shared inheritance (IBD-sharing) information between individuals for linkage analysis. This can improve the power of analysis, especially in cases where the population structure is complex or unknown [3].

**FarmCPU** (Fixed and random model Circulating Probability Unification) is an MLM method that accounts for differences in the distribution of effects between genes and also uses random effects to control for population structure. uses mixed models that account for both fixed and random effects for linkage analysis. It can also use a pre-trained model for faster analysis [4].

**SUPER** (Set-based Unified P-value for Association and Linkage Using Related Samples) is a method that uses MLM with a reduced covariance matrix to account for correlation between traits. It also accounts for linkage information in the sample. uses linkage information and MLM with a reduced covariance matrix to account for correlation between traits. It can also use gene group information for more accurate analysis [5].

**GLM** is based on generalizing linear regression to cases where the dependent variable (phenotype) does not have a normal distribution. Instead, the dependent variable can have any distribution from an exponential family of distributions, such as a binomial, Poisson, or Gamma distribution. The GLM assumes that the relationship between genotype and phenotype is described not only by a linear relationship, but also by other factors such as age, sex, population structure, etc. In the GLM model, these factors can be accounted for by adding the relevant variables as factors in the regression equation [6].

**GEMMA (BSLMM)** - uses a Bayesian approach to estimate gene effects. BSLMM incorporates the LASSO model to reduce the number of markers analyzed, which can improve the accuracy of the analysis [7].

**GEMMA (MVLMM)** - uses multilayer models to account for heterogeneity of gene effects. This can increase the power of the analysis, especially in cases where there are differences in gene effects in different subsets of the sample [8].

*Each method was run with default parameters.*


## Getting Started

1. Fork the repository.
2. If you have a personal token that allows you to work with private repositories through the command line, fine, if not, then generate it `Settings -> Developer settings -> Personal access tokens -> Generate new token
3. 
```ruby
git clone https://<PERSONAL-ACCESS Token>/<user_name>/GWAS_project.git
```
4.
```ruby
cd GWAS_project
```

## Running the tools

### [GAPIT (Version 3)](https://zzlab.net/GAPIT/) (FarmCPU, SUPER, Blink)
All the data for the methods is in the "[GAPIT](https://github.com/LiliiaBgdnv/GWAS_project/tree/main/GAPIT)" folder, and there is also an [**.rmd** file](https://github.com/LiliiaBgdnv/GWAS_project/blob/main/GAPIT/GAPIT.Rmd) with all the work of the three models. Here we also use files with the extension ".hmp.txt". Conversion of genotype data from vcf format to hapmap was performed in the [TASSEL 5](https://tassel.bitbucket.io/) program. 
 
### PLINK2 glm
**All the data is presented in the appropriate [folder](https://github.com/LiliiaBgdnv/GWAS_project/tree/main/plink_glm)**

Here we use the `allow-no-covars` flag since we do not have a file with principal component covariates and `--allow-extra-chr` in order to allow work with non-standard chromosome names. The `--adjust flag is used to calculate the p-value with the Bonferroni one-step correction, the Sidak one-step correction, the Bonferroni step-by-step method, the Benjamini/Hochberg FDR correction, the Benjamini/Yekutieli FDR correction

Here **.bed, .bim, .fam** files are made using the program [plink2](https://www.cog-genomics.org/plink/2.0/). You can install it on the website.

```ruby
cd plink_glm
```
Run [run_plink_glm.sh](https://github.com/LiliiaBgdnv/GWAS_project/blob/add-bash_scripts/plink_glm/run_plink_glm.sh)  script to preprocess the phenotype data and run plink glm tool.


### [GEMMA  v0.98.6](https://github.com/genetics-statistics/GEMMA)

For this model, we prepared other genotype and phenotype files, changing the chromosome names so that they contain only numbers, the processed files are in the [folder](https://github.com/LiliiaBgdnv/GWAS_project/tree/main/GEMMA_BSLMM/data_using_for_create_bim).

Installation:
```ruby
wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz
gunzip gemma-0.98.5-linux-static-AMD64.gz
mv gemma-0.98.5-linux-static-AMD64 gemma
chmod 777 gemma

```

#### BSLMM
```ruby
cd ../GEMMA_BSLMM
```
To preprocess the data and run BSLMM, run [run_gemma_bslmm.sh](https://github.com/LiliiaBgdnv/GWAS_project/blob/main/GEMMA_BSLMM/run_gemma_bslmm.sh) script

#### MVLMM
```ruby
cd ../GEMMA_MVLMM
```
Here we used [plink1.9](https://www.cog-genomics.org/plink/) to generate .bed files. Install it from the website.
To preprocess the data and run BSLMM, run [run_gemma_mvlmm.sh](https://github.com/MaryM12/GWAS_project/blob/main/GEMMA_MVLMM/run_gemma_mvlmm.sh) script


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
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/037af1bd-6fa0-4b55-ad05-faabb22f871c)
![image](https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/8f341bb3-af1b-402f-9f73-32f91ff1455c)


The processing of the results of the BSLMM model differed from other models, because the result of the Bayesian statistics does not contain the p-value calculated for each SNP, but rather contains posterior samples of parameters: random effects (alpha), fixed effects (beta), and sparse effects (gamma). Based on the gamma values for each file, the effect sizes for each SNP were calculated and the posterior probability was calculated. The posterior inclusion probability (PIP) was calculated as a measure of the strength of the association between the SNP and the phenotype. Based on the PIP, the distribution of linked SNPs across the genome was visualized. The size of the dots reflects the magnitude of the effect.

The results of the processing are in the [folder](https://github.com/LiliiaBgdnv/GWAS_project/tree/main/GEMMA_BSLMM/results), [script](https://github.com/LiliiaBgdnv/GWAS_project/blob/main/GEMMA_BSLMM/visualization.R) for the processing.

**Simple trait**
<img width="1388" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/23a2986c-5f5e-4d15-b1a5-744dcac58b03">

**Complex trait**
<img width="1388" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/b4020276-a790-4433-9291-8776beef578b">

As evident from the plot, with a threshold of **0.05 only 1 SNPs were found for the complex trait, and no SNPs were found at all for the simple trait**, so we had to lower the threshold to 0.02 for the rest of the analysis in order to compare the SNPs found by each model.

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
| MVLMM      | 1.68                     | 6.875                     |

<img width="1057" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/96f35edf-fe84-4127-97a8-95e66a0405ee">

## Annotation of the identified SNPs
The SNPs found by the program were annotated in the following steps:
1) Selection of found SNPs with p-value < 0.05 for the complex trait or p-value < 0.1 for the simple trait (for the BSLMM model PIP > 0.02). 
2) For each SNP, its chromosome number, its coordinate - 500bp and its coordinate + 500bp are taken and written into a .bed file. We investigate exactly the gaps around the SNPs found, because the model can predict not only SNPs in the protein-coding gene, but also the associated genes if the SNPs found are in the decoding gene. This and the previous step are done in [Jupiter notebook](https://github.com/LiliiaBgdnv/GWAS_project/blob/main/Benjamini_Yekutieli.ipynb).
3) Download the annotation of Glycine max:
```ruby
cd ../for_annotation
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/gff3/glycine_max/Glycine_max.Glycine_max_v2.1.56.chr.gff3.gz
gzip -d Glycine_max.Glycine_max_v2.1.56.chr.gff3.gz
```
4) Intersecting our .bed files and annotation using [bedtools v2.31.0](https://bedtools.readthedocs.io/en/latest/index.html):
Optional:
```ruby
#Install bedtools
#Linux
apt-get install bedtools
#macOS
brew tap homebrew/science
brew install bedtools
```
The intersecting complex trait
```ruby
bedtools intersect -wao -nonamecheck -a blink_complex.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/blink_complex_output.txt
bedtools intersect -wao -nonamecheck -a farmcpu_complex.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/farmcpu_complex_output.txt
bedtools intersect -wao -nonamecheck -a super_complex.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/super_complex_output.txt
bedtools intersect -wao -nonamecheck -a glm_complex.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/glm_complex_output.txt
bedtools intersect -wao -nonamecheck -a bslmm_complex.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/bslmm_complex_output.txt
```
The intersection of a simple feature is somewhat more complicated:
5) Upload the file to the [site](https://www.ncbi.nlm.nih.gov/genome/tools/remap#tab=asm&src_org=Glycine%20max&src_asm=GCF_000004515.3&tgt_asm=GCF_000004515.6&min_ratio=0.5&max_ratio=2.0&allow_locations=true&merge_fragments=false&in_fmt=guess&out_fmt=guess&genome_workbench=true) with the parameters 
> Source Organism : Glycine max 
> 
> Source Assembly: V1.1 
> 
> Target Assembly: Glycine_max_v4.0
> 
<img width="890" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/b32f63b2-8e00-401a-b0a8-5f1a1cfa5615">

upload our **.bed** file

Press `submit`.

6) Download the resulting file `Download Full Mapping Report`
7) In the resulting file, leave the columns **source_id, mapped_start, mapped_stop**, delete their name and the name of the table (if you have one).
8) The intersecting simple trait
```ruby
bedtools intersect -wao -nonamecheck -a r_blink_simple.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/blink_simple_output.txt
bedtools intersect -wao -nonamecheck -a r_farmcpu_simple.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/farmcpu_simple_output.txt
bedtools intersect -wao -nonamecheck -a r_bslmm_simple.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/bslmm_simple_output.txt
bedtools intersect -wao -nonamecheck -a r_glm_simple.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/glm_simple_output.txt
```

9) Extract protein IDs from these files:
```ruby
grep "CDS:" output/blink_complex_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_blink_complex.txt
grep "CDS:" output/farmcpu_complex_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_farmcpu_complex.txt     
grep "CDS:" output/super_complex_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_super_complex.txt
grep "CDS:" output/glm_complex_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_glm_complex.txt
grep "CDS:" output/bslmm_complex_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_bslmm_complex.txt
grep "CDS:" output/blink_simple_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_blink_simple.txt     
grep "CDS:" output/farmcpu_simple_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_farmcpu_simple.txt
grep "CDS:" output/bslmm_simple_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_bslmm_simple.txt
grep "CDS:" output/glm_simple_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_glm_simple.txt
```

10) For the last stage of the annotation, you need to see which protein is encoded under that ID at that [site](http://go.pantherdb.org) choosing the Glycine max organism.

On this [site](http://go.pantherdb.org/tools/gxIdsList.do?list=upload_1&organism=Glycine%20max) we visualize biological processes in which found proteins are used. In 1 point, the names from the abstract are given, and the next two points are indicated by the picture.
<img width="637" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/09873b74-c99f-4fce-9786-61f7eab3abec">

## Simple trait:
<img width="1174" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/9a1b974d-d403-4c42-8687-762dfb3445c6">

## Complex trait
<img width="1382" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/122a704e-9bc6-439f-b7fe-e2e06dd4cea8">


### Number of intersections of the significant SNPs identified by the tools.

For visualization, we built a heatmap. Preparation of **.bed, .bim, .fam** fills is done as in the *annotation part*, but with the difference that in the 2nd point they take coordinates with +- 1000bp.

```ruby
cd GWAS_ptoject/heatmap
```
We used tool [intervene](https://intervene.readthedocs.io/en/latest/). You can download it with one of these options:
```ruby
conda install -c bioconda intervene
# or
pip install intervene
```

```ruby
# simple trait
./intervene pairwise -i ./simple/*.bed
# complex trait
./intervene pairwise -i ./complex/*.bed
```

The heatmaps show the total number number of the SNPs, identified by each model (on the diagonal), as well the number of intersections between the models.

<img width="615" alt="image" src="https://github.com/LiliiaBgdnv/GWAS_project/assets/109213422/ff5b4c62-0783-4c15-bb95-fe8faab181e6">


## CONCLUSIONS:
1. **FarmCPU**, implemented in the GAPIT package, performed best for the **complex trait** dataset. It identified 8 SNPs with a relevant annotation and demonstrated the greatest statistical power and the highest degree of consistency with other models.
2. For the **simple trait** dataset, **FarmCPU** method also performed best. It identified 8 SNPs, with the highest statistical power and the highest number of overlaps with other models. But in general, few SNPs were identified, because amino acid content is not a very contrasting trait.

 To improve the quality of the analysis, we plan to further test the models with the fine-tuned parameters.

## Software Requirements

* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/python.svg height=20> Python 3.10
* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/ubuntu.svg height = 20> Ubuntu 21.04
* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/macos.svg height = 20> Ventura 13.2.1
* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/gnubash.svg height=20> Bash
* <img src=https://github.com/simple-icons/simple-icons/blob/develop/icons/r.svg height=20> R 4.2.3

| R using packages | Python using packages |
|------------------|-----------------------|
| .GlobalEnv       | statsmodels           |
| lme4             | pandas                |
| scatterplot3d    | matplotlib            |
| EMMREML          | qmplot                |
| Matrix           | scipy                 |
| bigmemory        | numpy                 |
| grid             | typing                |
| compiler         |                       |
| ape              |                       |
| genetics         |                       |
| mvtnorm          |                       |
| MASS             |                       |
| gtools           |                       |
| gdata            |                       |
| combinat         |                       |
| gplots           |                       |
| tools:rstudio    |                       |
| stats            |                       |
| graphics         |                       |
| grDevices         |                       |

## References

1.  Brachi B. Genome-wide association studies in plants: the missing heritability is in the field / B. Brachi, G.P. Morris, J.O. Borevitz // Genome Biology. – 2011. – Vol. 12. – Genome-wide association studies in plants. – № 10. – P. 232.
2.  Prioritized candidate causal haplotype blocks in plant genome-wide association studies / X. Wu [et al.] // PLOS Genetics. – 2022. – Vol. 18. – № 10. – P. e1010437.
3.  Meng Huang and others, BLINK: a package for the next level of genome-wide association studies with both individuals and markers in the millions, GigaScience, Volume 8, Issue 2, February 2019, giy154, https://doi.org/10.1093/gigascience/giy154
4.  Iterative Usage of Fixed and Random Effect Models for Powerful and Efficient Genome-Wide Association Studies, Liu X, Huang M, Fan B, Buckler ES, Zhang Z (2016) Iterative Usage of Fixed and Random Effect Models for Powerful and Efficient Genome-Wide Association Studies. PLOS Genetics 12(2): e1005767. https://doi.org/10.1371/journal.pgen.1005767
5.  A SUPER Powerful Method for Genome Wide Association Study, Wang Q, Tian F, Pan Y, Buckler ES, Zhang Z (2014) A SUPER Powerful Method for Genome Wide Association Study. PLOS ONE 9(9): e107684. https://doi.org/10.1371/journal.pone.0107684
6.  Chu BB, Keys KL, German CA, Zhou H, Zhou JJ, Sobel EM, Sinsheimer JS, Lange K. Iterative hard thresholding in genome-wide association studies: Generalized linear models, prior weights, and double sparsity. Gigascience. 2020 Jun 1;9(6):giaa044. doi: 10.1093/gigascience/giaa044. PMID: 32491161; PMCID: PMC7268817.
7.  Xiang Zhou, Peter Carbonetto and Matthew Stephens (2013). Polygenic modeling with bayesian sparse linear mixed models. PLoS Genetics 9, e1003264.
8.  Xiang Zhou and Matthew Stephens (2014). Efficient multivariate linear mixed model algorithms for genome-wide association studies. Nature Methods 11, 407–409.
