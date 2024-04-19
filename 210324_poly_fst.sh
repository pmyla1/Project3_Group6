########################
##This script is written by Luke Archer (2024) and uses Tuomas Hämälä's (2023) script to compile poly_fst.c and calculate Fst on a mixed ploidy vcf
##source of the poly_fst script https://github.com/thamala/polySV/blob/main/poly_fst.c
##############

#####################
##CONFIGURATION
source $HOME/.bash_profile

cd /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/filtered_VCFs_for_faststructure/

conda activate /shared/conda/shared/
############################

#############################
##uses grep and bcftools to obtain the individuals from each population for the fst scans 
bcftools query -l ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf.gz | grep "BZD" > ./BZD_pop.txt 

bcftools query -l ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf.gz | grep "OCH" > ./OCH_pop.txt

bcftools query -l ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf.gz | grep "MOD" > ./MOD_pop.txt

bcftools query -l ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf.gz | grep "KEH" > ./KEH_pop.txt

bcftools query -l ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf.gz | grep "LIC" > ./LIC_pop.txt

bcftools query -l ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf.gz | grep "FRE" > ./FRE_pop.txt

bcftools query -l ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf.gz | grep "HAB" > ./HAB_pop.txt
#######################


########################
##compile the poly_fst.c script
gcc ./poly_fst.c -o poly_fst -lm
######################

#######################
#run poly_fst on the 2002324_tets_only_filtered_vcf and BZD vs OCH population contrast
./poly_fst -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pop1 ./BZD_pop.txt -pop2 ./OCH_pop.txt -mis 0.8 > ./BZD_OCH_Fst_output.fst 
##now run poly_fst on BZD vs KEH
./poly_fst -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pop1 ./BZD_pop.txt -pop2 ./KEH_pop.txt -mis 0.8 > ./BZD_KEH_Fst_output.fst
##OCH vs KEH
./poly_fst -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pop1 ./OCH_pop.txt -pop2 ./KEH_pop.txt -mis 0.8 > ./OCH_KEH_Fst_output.fst
##OCH vs LIC
./poly_fst -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pop1 ./OCH_pop.txt -pop2 ./LIC_pop.txt -mis 0.8 > ./OCH_LIC_Fst_output.fst
##OCH vs MOD
./poly_fst -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pop1 ./OCH_pop.txt -pop2 ./MOD_pop.txt -mis 0.8 > ./OCH_MOD_Fst_output.fst
##OCH vs HAB
./poly_fst -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pop1 ./OCH_pop.txt -pop2 ./HAB_pop.txt -mis 0.8 > ./OCH_HAB_Fst_output.fst
##OCH vs FRE
./poly_fst -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pop1 ./OCH_pop.txt -pop2 ./FRE_pop.txt -mis 0.8 > ./OCH_FRE_Fst_output.fst
##KEH vs LIC
./poly_fst -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pop1 ./KEH_pop.txt -pop2 ./LIC_pop.txt -mis 0.8 > ./KEH_LIC_Fst_output.fst
##KEH vs MOD 
./poly_fst -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pop1 ./KEH_pop.txt -pop2 ./MOD_pop.txt -mis 0.8 > ./KEH_MOD_Fst_output.fst
##KEH vs HAB
./poly_fst -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pop1 ./KEH_pop.txt -pop2 ./HAB_pop.txt -mis 0.8 > ./KEH_HAB_Fst_output.fst
##KEH vs FRE
./poly_fst -vcf ./200324_tets_only_VCFs/200324_cleaned_lyrata_tets_only.vcf -pop1 ./KEH_pop.txt -pop2 ./FRE_pop.txt -mis 0.8 > ./KEH_FRE_Fst_output.fst
######################

conda deactivate

echo "DONE!!"

##lastline
