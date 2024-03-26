###################
#This script is written by Luke Archer (2024) and uses the poly_freq.c script from Tuomas Hämälä (2023) 
##which can be accessed on: https://github.com/thamala/polySV/blob/main/poly_sfs.c
###############

################
##CONFIGURATION
source $HOME/.bash_profile

cd /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/filtered_VCFs_for_faststructure/

conda activate /shared/conda/shared/
################

#################
##compiling the poly_sfs.c script
#gcc ./poly_sfs.c -o poly_sfs -lm
################

################
##use the script on the 200324_tets_only.vcf
./poly_sfs -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -inds ./individuals_200324.txt -mis 0.8 > ./200324_poly_sfs_out.tsv
################

##################
##compile and use poly_freq.c script
gcc ./poly_freq.c -o poly_freq -lm

###specify the population file, proportion of missing data as 0.8 (only 20% missing) and minor allele freq of 0.05
##-out 0 specifies allele frequencies rather than counts
./poly_freq -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pops ./new_pops_200324.txt -mis 0.8 -maf 0.05 -out 0 > ./200324_poly_freq_out
##################

conda deactivate

echo "DONE!!"

##lastline
