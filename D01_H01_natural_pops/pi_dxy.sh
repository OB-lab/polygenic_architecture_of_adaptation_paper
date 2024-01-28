### Calculate Pi and Dxy for each putative inversion

module load bcftools 

DIR="/.../pi_dxy"
FILE="/.../pa_as_sf5_final.recode.vcf.gz"

# filter to the chromosome
bcftools filter ${FILE} --regions ${1} --output ${DIR}/${2}.vcf 
bgzip -@ 12 ${DIR}/${2}.vcf
tabix ${DIR}/${2}.vcf.gz

# for each inversion check the genotype, and assign the genotype value as the 'population' for that individual
POPS="${DIR}/${2}_POPS.txt"
genotypes="/home/564/km6006/Scripts/avneet_paper/inversion_genotypes.txt"

while IFS=$'\t' read -r inversion name genotype population; do
    echo "current $inversion has "$genotype" genotype"
    if [ "$inversion" == "$2" ] && [ "$genotype" -eq 0 ]; then
        echo -e "$name\t$genotype" >> "$POPS"
    elif [ "$inversion" == "$2" ] && [ "$genotype" -eq 1 ]; then
        echo -e "$name\t$genotype" >> "$POPS"
    elif [ "$inversion" == "$2" ] && [ "$genotype" -eq 2 ]; then
        echo -e "$name\t$genotype" >> "$POPS"        
    else
        echo "no match"
    fi
done < "$genotypes"

# Run pixy - calculate pi and dxy values
pixy --stats pi dxy --vcf ${DIR}/${2}.vcf.gz --populations ${POPS} --window_size 10000 --output_prefix ${2} --output_folder ${DIR} --n_cores 12