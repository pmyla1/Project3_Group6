##############
##This script was written by Luke Archer (2024) and uses GATK to select individuals from the original VCF file and goes
##through the  whole pipeline from selecting variants to faststructure

##########
##CONFIGURATION
source $HOME/.bash_profile

cd /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/

##activate GATK env
conda activate /shared/apps/conda/bio2/
###############

##########
##GATK select variants to include one pure arenosa (KEH), one pure lyrata (MOD) and hybrid tetraploids only 
gatk SelectVariants -V Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz -sn FRE-03tl -sn FRE-04tl -sn FRE-05tl -sn FRE-06tl -sn FRE-07tl -sn FRE-08tl -sn HAB-01tl -sn HAB-02tl -sn HAB-03tl -sn KEH-05tl -sn KEH-06tl -sn KEH-07tl -sn KEH-08tl -sn KEH-09tl -sn KEH-10tl -sn LOI-01tl -sn LOI-02tl -sn LOI-03tl -sn MOD-01tl -sn MOD-02tl -sn MOD-03tl -sn MOD-04tl -sn OCH-03tl -sn OCH-04tl -sn OCH-05tl -sn OCH-06tl -sn OCH-07tl -sn OCH-08tl -sn PIL-01tl -sn PIL-02tl -sn PIL-03tl -O ./290324_whole_pipeline_VCFs/290324_tetraploids_only.vcf.gz 

conda deactivate
#################

#################
conda activate /shared/conda/shared/

##copy the vcf file then unzip for further analyses with poly_sfs, poly_freq, Cochlearia_create_structure_file.py
cp ./290324_whole_pipeline_VCFs/290324_tetraploids_only.vcf.gz ./290324_whole_pipeline_VCFs/290324_tetraploids_only_copy.vcf.gz

gunzip ./290324_whole_pipeline_VCFs/290324_tetraploids_only_copy.vcf.gz > ./290324_whole_pipeline_VCFs/290324_tetraploids_only_copy.vcf

conda deactivate
#################

###################
##activate faststructure env
conda activate /shared/conda/faststructure/

##################
##poly_sfs script for site-frequency spectra estimates
./filtered_VCFs_for_faststructure/poly_sfs -vcf ./290324_whole_pipeline_VCFs/290324_tetraploids_only_copy.vcf -inds ./290324_whole_pipeline_VCFs/individuals_final.txt -mis 0.8 > ./290324_whole_pipeline_VCFs/290324_SFS_output.sfs

#########
##poly_freq script for estimating allele frequencies from mixed ploidy vcf files
./filtered_VCFs_for_faststructure/poly_freq -vcf ./290324_whole_pipeline_VCFs/290324_tetraploids_only_copy.vcf -pops ./290324_whole_pipeline_VCFs/populations_final.txt -mis 0.8 -maf 0.02 -out 0 > ./290324_whole_pipeline_VCFs/290324_poly_freq_output

######
##poly_fst for estimating pairwise Fst and Dxy from mixed ploidy VCFs


###################
##CREATE STRUCTURE APPROPRIATE FILES 
##Cochlearia_create_structure.py to create faststructure appropriate files .str
python3 ./scripts/Cochlearia_create_structure_file.py -v ./290324_whole_pipeline_VCFs/ -o 290324_structure_files -s true 

##remove the first and last lines from the diploidized structure file
sed "1d" ./290324_whole_pipeline_VCFs/vcf_to_str/290324_structure_files.StructureInputDiploidized.str | sed "$d" ./290324_whole_pipeline_VCFs/vcf_to_str/290324_structure_files.StructureInputDiploidized.str > ./290324_whole_pipeline_VCFs/290324_first_last_removed.StructureInputDiploidized.str

###################
##reorder the structure files into alphabetical order

grep "FRE" ./290324_whole_pipeline_VCFs/290324_first_last_removed.StructureInputDiploidized.str > ./290324_whole_pipeline_VCFs/FRE.str

grep "HAB" ./290324_whole_pipeline_VCFs/290324_first_last_removed.StructureInputDiploidized.str > ./290324_whole_pipeline_VCFs/HAB.str

grep "KEH" ./290324_whole_pipeline_VCFs/290324_first_last_removed.StructureInputDiploidized.str > ./290324_whole_pipeline_VCFs/KEH.str

grep "LOI" ./290324_whole_pipeline_VCFs/290324_first_last_removed.StructureInputDiploidized.str > ./290324_whole_pipeline_VCFs/LOI.str

grep "MOD" ./290324_whole_pipeline_VCFs/290324_first_last_removed.StructureInputDiploidized.str > ./290324_whole_pipeline_VCFs/MOD.str

grep "OCH" ./290324_whole_pipeline_VCFs/290324_first_last_removed.StructureInputDiploidized.str > ./290324_whole_pipeline_VCFs/OCH.str

grep "PIL" ./290324_whole_pipeline_VCFs/290324_first_last_removed.StructureInputDiploidized.str > ./290324_whole_pipeline_VCFs/PIL.str


##concatenate into a reordered .str file
cd ./290324_whole_pipeline_VCFs/

cat FRE.str HAB.str KEH.str LOI.str MOD.str OCH.str PIL.str > ./290324_reordered_structure.str 

cd ../
######################

############
#run structure with different K values
python /shared/conda/faststructure/bin/structure.py -K 2 --input=./290324_whole_pipeline_VCFs/290324_reordered_structure --output=./290324_whole_pipeline_VCFs/290324_K2_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 3 --input=./290324_whole_pipeline_VCFs/290324_reordered_structure --output=./290324_whole_pipeline_VCFs/290324_K3_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 4 --input=./290324_whole_pipeline_VCFs/290324_reordered_structure --output=./290324_whole_pipeline_VCFs/290324_K4_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 5 --input=./290324_whole_pipeline_VCFs/290324_reordered_structure --output=./290324_whole_pipeline_VCFs/290324_K5_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 6 --input=./290324_whole_pipeline_VCFs/290324_reordered_structure --output=./290324_whole_pipeline_VCFs/290324_K6_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 7 --input=./290324_whole_pipeline_VCFs/290324_reordered_structure --output=./290324_whole_pipeline_VCFs/290324_K7_out --format=str --full
#############

#############
##run new_distruct.py
python ./scripts/new_distruct.py -K 2 --input=./290324_whole_pipeline_VCFs/290324_K2_out --output=./290324_whole_pipeline_VCFs/290324_K2_plot --popfile=./290324_whole_pipeline_VCFs/populations_final.txt --title="Arabidopsis lyrata admixture: K=2"

python ./scripts/new_distruct.py -K 3 --input=./290324_whole_pipeline_VCFs/290324_K3_out --output=./290324_whole_pipeline_VCFs/290324_K3_plot --popfile=./290324_whole_pipeline_VCFs/populations_final.txt --title="Arabidopsis lyrata admixture: K=3"

python ./scripts/new_distruct.py -K 4 --input=./290324_whole_pipeline_VCFs/290324_K4_out --output=./290324_whole_pipeline_VCFs/290324_K4_plot --popfile=./290324_whole_pipeline_VCFs/populations_final.txt --title="Arabidopsis lyrata admixture: K=4"

python ./scripts/new_distruct.py -K 5 --input=./290324_whole_pipeline_VCFs/290324_K5_out --output=./290324_whole_pipeline_VCFs/290324_K5_plot --popfile=./290324_whole_pipeline_VCFs/populations_final.txt --title="Arabidopsis lyrata admixture: K=5"

python ./scripts/new_distruct.py -K 6 --input=./290324_whole_pipeline_VCFs/290324_K6_out --output=./290324_whole_pipeline_VCFs/290324_K6_plot --popfile=./290324_whole_pipeline_VCFs/populations_final.txt --title="Arabidopsis lyrata admixture: K=6"

python ./scripts/new_distruct.py -K 7 --input=./290324_whole_pipeline_VCFs/290324_K7_out --output=./290324_whole_pipeline_VCFs/290324_K7_plot --popfile=./290324_whole_pipeline_VCFs/populations_final.txt --title="Arabidopsis lyrata admixture: K=7"
###############

conda deactivate 

echo "DONE!!"

##lastline

