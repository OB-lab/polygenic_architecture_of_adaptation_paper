### Calculate Pi, Dxy and Fst between populations, for all chromosomes

source ~/miniconda3/etc/profile.d/conda.sh
conda activate bio

module load bcftools 

FST_IN="/.../pa_vo_sf5_final.vcf.gz"
D01="/.../D01.txt"
H01="/.../H01.txt"

PIXY_IN="/.../pa_as_sf5_final.vcf.gz"
PIXY_DIR="/.../pvp_popgen"
POPS="/.../D01H01.txt"

# Calculate Pi and Dxy between populations 
tabix ${PIXY_IN}
pixy --stats pi dxy --vcf ${PIXY_IN} --populations ${POPS} --window_size 10000 --output_prefix D01H01 --output_folder ${PIXY_DIR} --n_cores 12

# Calculate Fst between populations 
/home/564/km6006/bin/vcftools \
    --gzvcf ${FST_IN} \
    --weir-fst-pop ${D01} \
    --weir-fst-pop ${H01} \
    --fst-window-size 10000 \
    --fst-window-step 10000 \
    --out /.../pvp_popgen/D01H01
