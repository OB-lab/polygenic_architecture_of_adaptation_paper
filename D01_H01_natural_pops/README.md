### Analysis of WGS and HiC data from D01 and HO1 natural populations 


- heterozygosity.sh, het_plot.R 

   Calculate heterozygosity per sample for all SNPs from each putative inversion

- pca.R

   Create PCA with kmeans clustering of all SNPs from each putative inversion 

- fst.sh 

   Calculate fst between homozygous genotype clusters for each putative inversion identified by local PCA. 

- pi_dxy.sh 

    Calculte pi and dxy between homozygous genotype clusters for each putative inversion identified by local PCA.

- pi_dxy_pop.sh 

    Calculte pi, dxy and fst between populations for every chromosome.

- hic.sh 

    Create hic interaction matricies for each chromosome at a resolution of 100kb windows.