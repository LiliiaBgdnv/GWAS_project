# preprocess the data
../plink2 --vcf ./data_using_for_create_bim/genotypes_complex_3_1_3.vcf -pheno ./data_using_for_create_bim/phenotypes_soy_complex_3_1_3.tsv --make-bed --allow-extra-chr --max-alleles 2 --out ./input/complex_trait
../plink2 --vcf ./data_using_for_create_bim/soybean_rename_chr.vcf -pheno ../Raw_data/phenotypes_simple_trait_2col.tsv --make-bed --allow-extra-chr --out ./input/simple_trait
../plink2 --vcf ./data_using_for_create_bim/genotypes_complex_3_1_3.vcf -pheno ./data_using_for_create_bim/phenotypes_soy_complex_gen.tsv --make-bed --allow-extra-chr --max-alleles 2 --out ./input/complex_trait_gen
../plink2 --vcf ./data_using_for_create_bim/soybean_rename_chr.vcf -pheno ../Raw_data/phenotypes_simple_trait_2col_gen.tsv --make-bed --allow-extra-chr --out ./input/simple_trait_gen

# run BSLMM
# simple trait 
../gemma -bfile ./input/simple_trait -bslmm 1 -o gemma_bslmm_output_simple
# complex trait 
../gemma -bfile ./input/complex_trait -bslmm 1 -o gemma_bslmm_output_complex
# simple trait generated
../gemma -bfile ./input/simple_generated -bslmm 1 -o gemma_bslmm_output_simple_gen
# complex trait generated
../gemma -bfile ./input/complex_generated -bslmm 1 -o gemma_bslmm_output_complex_gen
