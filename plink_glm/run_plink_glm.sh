# prepare phenotype files 
# select 2 columns needed for analysis
cut -f2- ../Raw_data/phenotypes_simple_trait.tsv > ../Raw_data/phenotypes_simple_trait_2col.tsv
cut -f2- ../Raw_data/phenotypes_simple_trait_gen.tsv > ../Raw_data/phenotypes_simple_trait_2col_gen.tsv
../plink2 --vcf ../Raw_data/genotypes_complex_trait.vcf -pheno ../Raw_data/phenotypes_complex_trait.tsv --make-bed --allow-extra-chr --max-alleles 2 --out ./inputs/complex_trait
../plink2 2 --vcf ../Raw_data/genotypes_simple_trait.vcf -pheno ../Raw_data/phenotypes_simple_trait_2col.tsv --make-bed --allow-extra-chr --out ./inputs/simple_trait
../plink2 --vcf ../Raw_data/genotypes_complex_trait.vcf -pheno ../Raw_data/phenotypes_complex_trait_gen.tsv --make-bed --allow-extra-chr --max-alleles 2 --out ./inputs/complex_trait_gen
../plink2 --vcf ../Raw_data/genotypes_simple_trait.vcf -pheno ../Raw_data/phenotypes_simple_trait_gen_2col_gen.tsv --make-bed --allow-extra-chr --out ./inputs/simple_trait_gen

# Run plink2
# simple trait
../plink2 --glm allow-no-covars --allow-extra-chr --bed ./inputs/simple_trait.bed --bim ./inputs/simple_trait.bim --fam ./inputs/simple_trait.fam --pheno ../Raw_data/phenotypes_simple_trait_2col.tsv --adjust cols='chrom','pos','alt','a1','ref','gc','fdrbh' --out ./results/simple_glm_result_chr --covar-variance-standardize --freq --threads 32 --memory 100000
awk '{gsub(/^GLYMAchr_/,""); print}' ./results/simple_glm_result_chr.Leu.glm.linear.adjusted  > ./results/simple_glm_result_chr.Leu.glm.linear.adjusted
awk '{gsub(/^GLYMAchr_/,""); print}' ./results/simple_glm_result_chr.Leu.glm.linear > plink_glm/results/simple_glm_result_chr.Leu.glm.linear

# complex trait
../plink2 --glm allow-no-covars --allow-extra-chr --bed ./inputs/complex_trait.bed --bim ./inputs/complex_trait.bim --fam ./inputs/complex_trait.fam --pheno ../Raw_data/phenotypes_complex_trait.tsv --adjust --out ./results/complex_glm_result --covar-variance-standardize --freq --threads 32 --memory 100000

# simple trait generated
../plink2 --glm allow-no-covars --allow-extra-chr --bed ./inputs/simple_trait.bed --bim ./inputs/simple_trait.bim --fam ./inputs/simple_trait.fam --pheno ../Raw_data/phenotypes_simple_trait_gen_2col_gen.tsv --adjust cols='chrom','pos','alt','a1','ref','gc','fdrbh' --out ./results/simple_glm_gen_result_chr --covar-variance-standardize --freq --threads 32 --memory 100000
awk '{gsub(/^GLYMAchr_/,""); print}' ./results/simple_glm_gen_result_chr.Leu.glm.linear.adjusted  > ./results/simple_glm_gen_result_chr.Leu.glm.linear.adjusted
awk '{gsub(/^GLYMAchr_/,""); print}' ./results/simple_glm_gen_result_chr.Leu.glm.linear > ./results/simple_glm_gen_result_chr.Leu.glm.linear

# complex trait generated
plink2 --glm allow-no-covars --allow-extra-chr --bed ./inputs/complex_trait.bed --bim ./inputs/complex_trait.bim --fam ./inputs/complex_trait.fam --pheno ../Raw_data/phenotypes_complex_trait_gen.tsv --adjust --out ./results/complex_glm_gen_result --covar-variance-standardize --freq --threads 32 --memory 100000
