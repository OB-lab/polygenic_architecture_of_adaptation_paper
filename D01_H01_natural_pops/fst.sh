### calculate Fst for each putative inversion

module load bcftools
module load vcftools 

FILE="/.../pa_vo_sf5_final.recode_RS.vcf.gz"
DIR="/.../fst"
OUTPUT="/.../${2}_fst"

# filter to the chromosome
bcftools filter ${FILE} --regions ${1} --output ${DIR}/${2}.vcf 

# create list variables for each homozygous genotype group
G0="${DIR}/${2}_G0.txt"
G2="${DIR}/${2}_G2.txt"

# read in the sample genotype data for all inversions
genotypes="/.../inversion_genotypes.txt"

# for the current inversion add each homozygous sample to the correct list, based on their genotype
while IFS=$'\t' read -r inversion name genotype population; do
    echo "current $inversion has $genotype genotype"
    if [ "$inversion" == "$2" ] && [ "$genotype" -eq 0 ]; then
        echo -e "$name\t$population" >> "$G0"
    elif [ "$inversion" == "$2" ] && [ "$genotype" -eq 2 ]; then
        echo -e "$name\t$population" >> "$G2"
    else
        echo "het sample"
    fi
done < "$genotypes"

wait

# calculate Weir and Cockerham’s Fst between the two homozygous genotypes
vcftools \
    --vcf ${DIR}/${2}.vcf \
    --weir-fst-pop ${G0} \
    --weir-fst-pop ${G2} \
    --fst-window-size 10000 \
    --fst-window-step 10000 \
    --out ${OUTPUT}