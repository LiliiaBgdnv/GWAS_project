# preprocess the data
sed 's/KZ//g' genotypes_complex_trait.vcf > genotypes_complex_trait_processed.vcf
# create bed bim fam files
./plink1.9 --vcf ../Raw_data/genotypes_complex_trait_processed.vcf -pheno ../Raw_data/soy2_phenotypes_fast.txt --make-bed --allow-extra-chr --out ./input/soy2
./plink1.9 --vcf ../Raw_data/genotypes_simple_trait.vcf -pheno ../Raw_data/soybean_simple_genotypes_phen_recode1.9.txt --make-bed --allow-extra-chr --out ./input/soybean_simple_genotypes_phen_recode1.9
./plink1.9 --vcf ../Raw_data/genotypes_complex_trait_processed.vcf -pheno ../Raw_data/soy2_gen_phenotypes_fast.txt --make-bed --allow-extra-chr --out ./input/soy2gen
./plink1.9 --vcf ../Raw_data/genotypes_simple_trait.vcf -pheno ../Raw_data/gen_soy_phenotypes_leu.tsv --make-bed --allow-extra-chr --out ./input/soybean_generated
# simple trait
## generating the kinship matrix
../gemma -bfile ./input/soybean_simple_genotypes_phen_recode1.9 -gk 1 -o kinship_soybean_matrix_center
## running mvlmm
../gemma -bfile ./input/soybean_simple_genotypes_phen_recode1.9 -k ./output/kinship_soybean_matrix_center.cXX.txt -lmm 4 -n 1 -o simple_soy_mvlmm

# complex trait
## generating the kinship matrix
../gemma -bfile ./input/soy2 -gk 1 -o kinship_soy2_matrix_center
## running mvlmm
../gemma -bfile ./input/soy2 -k ./output/kinship_soy2_matrix_center.cXX.txt -lmm 4 -n 1 -o complex_soy2_mvlmm

# simple trait generated data
## generating the kinship matrix
../gemma -bfile ./input/soy2gen -gk 1 -o kinship_soybean_matrix_center
## running mvlmm
../gemma -bfile ./input/soy2gen -k ./output/kinship_soybean_matrix_center.cXX.txt -lmm 4 -n 1 -o simple_soy_gen_mvlmm

# complex trait generated data
## generating the kinship matrix
../gemma -bfile ./input/soybean_generated -gk 1 -o kinship_soy2_matrix_center 
## running mvlmm
../gemma -bfile ./input/soybean_generated -k ./output/kinship_soy2_matrix_center.cXX.txt -lmm 4 -n 1 -o complex_soy2_gen_mvlmm.assoc