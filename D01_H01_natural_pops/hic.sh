dir="/.../..."

# trim 
homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 -q 15 ${dir}/${1}_HiC_R1.fastq.gz
pigz ${dir}/${1}_HiC_R1_trimmed.fastq
homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 -q 15 ${dir}/${1}_HiC_R2.fastq.gz
pigz ${dir}/1_${1}_HiC_R2_trimmed.fastq

# align 
ngm -r /.../SLv141Asm_Ch20RN.fasta -t 24 -q ${dir}/${1}_HiC_R1_trimmed.fastq.gz -o ${dir}/${1}_HiC_R1.sam; 
ngm -r /.../SLv141Asm_Ch20RN.fasta -t 24 -q ${dir}/${1}_HiC_R2_trimmed.fastq.gz -o ${dir}/${1}_HiC_R2.sam;

# make tag directories
makeTagDirectory ${dir}/${1}_HiC-TAG ${dir}/${1}_HiC_R1.sam,${dir}/${1}_HiC_R2.sam -tbp 1 -mapq 10’ 

# create Hi-C interaction matricies and convert to a table 
declare -a arr=("SLv141Ch1" "SLv141Ch2" "SLv141Ch3" "SLv141Ch4" "SLv141Ch5" "SLv141Ch6" "SLv141Ch7" "SLv141Ch8" "SLv141Ch9" "SLv141Ch10" "SLv141Ch11" "SLv141Ch12" "SLv141Ch13" "SLv141Ch14" "SLv141Ch15" "SLv141Ch16" "SLv141Ch17" "SLv141Ch18" "SLv141Ch19" "SLv141Ch20")

for i in "${arr[@]}";
do
    analyzeHiC ${dir}/1_tags/${1}_HiC-TAG -chr ${i} -res 100000 -coverageNorm > ${dir}/${1}_${i}_hic_matrix.txt
    cat ${dir}/${1}_${i}_hic_matrix.txt | perl /.../hicmatrix2table.pl ${1} > ${dir}/${1}_${i}_hic_table.txt
done

wait