# Intersect SNPs for the complex trait
bedtools intersect -wao -nonamecheck -a blink_complex.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/blink_complex_output.txt
bedtools intersect -wao -nonamecheck -a farmcpu_complex.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/farmcpu_complex_output.txt
bedtools intersect -wao -nonamecheck -a super_complex.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/super_complex_output.txt
bedtools intersect -wao -nonamecheck -a glm_complex.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/glm_complex_output.txt
bedtools intersect -wao -nonamecheck -a bslmm_complex.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/bslmm_complex_output.txt

# Intersect SNPs for the simple trait
bedtools intersect -wao -nonamecheck -a r_blink_simple.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/blink_simple_output.txt
bedtools intersect -wao -nonamecheck -a r_farmcpu_simple.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/farmcpu_simple_output.txt
bedtools intersect -wao -nonamecheck -a r_bslmm_simple.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/bslmm_simple_output.txt
bedtools intersect -wao -nonamecheck -a r_glm_simple.bed -b Glycine_max.Glycine_max_v2.1.56.chr.gff3 > output/glm_simple_output.txt

# Extract protein IDs from the resulting files
grep "CDS:" output/blink_complex_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_blink_complex.txt
grep "CDS:" output/farmcpu_complex_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_farmcpu_complex.txt     
grep "CDS:" output/super_complex_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_super_complex.txt
grep "CDS:" output/glm_complex_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_glm_complex.txt
grep "CDS:" output/bslmm_complex_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_bslmm_complex.txt
grep "CDS:" output/blink_simple_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_blink_simple.txt     
grep "CDS:" output/farmcpu_simple_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_farmcpu_simple.txt
grep "CDS:" output/bslmm_simple_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_bslmm_simple.txt
grep "CDS:" output/glm_simple_output.txt | awk -F'CDS:|;' '{print $2}' > output/protein/protein_glm_simple.txt