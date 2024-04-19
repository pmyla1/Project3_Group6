#################
##This script is for the compilation of poly_fst.c to calculate Fst and Dxy on mixed
##ploidy vcfs, written by Tuomas Hamala (2023), source https://github.com/thamala/polySV/blob/main/poly_fst.c
##############

#############
##Configuration
source $HOME/.bash_profile

cd /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/290324_whole_pipeline_VCFs/

conda activate /shared/conda/shared/
############

##############
##Make a new directory for the populations to go into
mkdir 190424_Fst_populations

##make new directory for fst output
mkdir 190424_Fst_output
#############

############
#use grep and bcftools to obtain the individuals from each population for the fst scans 
bcftools query -l ./290324_tetraploids_only.vcf.gz | grep "BZD" > ./190424_Fst_populations/BZD_pop.txt 

bcftools query -l ./290324_tetraploids_only.vcf.gz | grep "OCH" > ./190424_Fst_populations/OCH_pop.txt

bcftools query -l ./290324_tetraploids_only.vcf.gz | grep "MOD" > ./190424_Fst_populations/MOD_pop.txt

bcftools query -l ./290324_tetraploids_only.vcf.gz | grep "KAG" > ./190424_Fst_populations/KAG_pop.txt

bcftools query -l ./290324_tetraploids_only.vcf.gz | grep "KEH" > ./190424_Fst_populations/KEH_pop.txt

bcftools query -l ./290324_tetraploids_only.vcf.gz | grep "LOI" > ./190424_Fst_populations/LOI_pop.txt

bcftools query -l ./290324_tetraploids_only.vcf.gz | grep "LIC" > ./190424_Fst_populations/LIC_pop.txt

bcftools query -l ./290324_tetraploids_only.vcf.gz | grep "FRE" > ./190424_Fst_populations/FRE_pop.txt

bcftools query -l ./290324_tetraploids_only.vcf.gz | grep "HAB" > ./190424_Fst_populations/HAB_pop.txt
############


############
##compile the poly_fst.c script
gcc ../scripts/poly_fst.c -o ./poly_fst -lm
###########

##########
#run poly_fst on the 2002324_tets_only_filtered_vcf and BZD vs OCH population contrast
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/BZD_pop.txt -pop2 ./190424_Fst_populations/OCH_pop.txt -mis 0.8 > ./190424_Fst_output/BZD_OCH_Fst_output.fst 
##now run poly_fst on BZD vs KEH
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/BZD_pop.txt -pop2 ./190424_Fst_populations/KEH_pop.txt -mis 0.8 > ./190424_Fst_output/BZD_KEH_Fst_output.fst
##OCH vs KEH
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/OCH_pop.txt -pop2 ./190424_Fst_populations/KEH_pop.txt -mis 0.8 > ./190424_Fst_output/OCH_KEH_Fst_output.fst
##OCH vs LIC
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/OCH_pop.txt -pop2 ./190424_Fst_populations/LIC_pop.txt -mis 0.8 > ./190424_Fst_output/OCH_LIC_Fst_output.fst
##OCH vs MOD
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/OCH_pop.txt -pop2 ./190424_Fst_populations/MOD_pop.txt -mis 0.8 > ./190424_Fst_output/OCH_MOD_Fst_output.fst
##OCH vs HAB
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/OCH_pop.txt -pop2 ./190424_Fst_populations/HAB_pop.txt -mis 0.8 > ./190424_Fst_output/OCH_HAB_Fst_output.fst
##OCH vs FRE
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/OCH_pop.txt -pop2 ./190424_Fst_populations/FRE_pop.txt -mis 0.8 > ./190424_Fst_output/OCH_FRE_Fst_output.fst
##KEH vs LIC
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/KEH_pop.txt -pop2 ./190424_Fst_populations/LIC_pop.txt -mis 0.8 > ./190424_Fst_output/KEH_LIC_Fst_output.fst
##KEH vs MOD 
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/KEH_pop.txt -pop2 ./190424_Fst_populations/MOD_pop.txt -mis 0.8 > ./190424_Fst_output/KEH_MOD_Fst_output.fst
##KEH vs HAB
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/KEH_pop.txt -pop2 ./190424_Fst_populations/HAB_pop.txt -mis 0.8 > ./190424_Fst_output/KEH_HAB_Fst_output.fst
##KEH vs FRE
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/KEH_pop.txt -pop2 ./190424_Fst_populations/FRE_pop.txt -mis 0.8 > ./190424_Fst_output/KEH_FRE_Fst_output.fst
##KEH vs KAG
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/KEH_pop.txt -pop2 ./190424_Fst_populations/KAG_pop.txt -mis 0.8 > ./190424_Fst_output/KEH_KAG_Fst_output.fst
##KEH vs LOI
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/KEH_pop.txt -pop2 ./190424_Fst_populations/LOI_pop.txt -mis 0.8 > ./190424_Fst_output/KEH_LOI_Fst_output.fst
##BZD vs KAG
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/BZD_pop.txt -pop2 ./190424_Fst_populations/KAG_pop.txt -mis 0.8 > ./190424_Fst_output/BZD_KAG_Fst_output.fst
##BZD vs LOI
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/BZD_pop.txt -pop2 ./190424_Fst_populations/LOI_pop.txt -mis 0.8 > ./190424_Fst_output/BZD_LOI_Fst_output.fst
##########


conda deactivate

echo "DONE!!"

##lastline
