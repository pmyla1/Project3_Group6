####This script was written by Luke Archer (2024) to select only biallelic variants from a 
##filtered VCF file to remove indels/large mutations and subsequently run the 
##poly_sfs.c script sourced from Tuomas Hamala (2023) Github 
##https://github.com/thamala/polySV/blob/main/poly_sfs.c 

#########################
##CONFIGURATION
##source bash profile
source $HOME/.bash_profile

##cd to filtered_VCFs_for_faststructure
cd /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/filtered_VCFs_for_faststructure/
##########################

###########################
##GATK SELECT VARIANTS TO INCLUDE BIALLELIC SITES ONLY
##activate GATK environment
conda activate /shared/apps/conda/bio2/
########################

#########GATK SELECT VARIANTS
gatk SelectVariants -R ../lyrata.fasta -V 150324_filtered_pops.vcf.gz --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O 160324_SFS_filtered_pops.vcf.gz
##########################

#############################
##POLY_SFS CONFIGURATION AND RUNNING
##########compile and run poly_sfs.c by Tuomas Hamala
#compiling
gcc ../scripts/poly_sfs.c -o poly_sfs -lm

##running the poly_sfs script on the VCF file using a text file containing the individuals in the VCF
./poly_sfs -vcf 160324_SFS_filtered_pops.vcf.gz -inds individuals.txt -mis 0.9 > SFS_output.sfs
############################

conda deactivate

echo "DONE!!"

##lastline
