### Calculate HET for each individual for all SNPs in the putative inversion

module load bcftools 
module load vcftools

file="/.../pa_vo_sf5_final.vcf.gz"
output="/.../het/${1}"

# subset the vcf file to the inversion region
bcftools filter ${file} --threads 12 --regions ${1} -o ${output}

# calcuate heterozygosity per individual 
vcftools \
    --vcf ${output} \
    --het \
    --out /.../het/${1}