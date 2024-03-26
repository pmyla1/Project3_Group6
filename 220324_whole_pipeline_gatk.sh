####################################
##This script was written by Luke Archer (2024) and can be used to select individuals from the original VCF
##and subsequently produce site-frequency spectra, allele frequencies, and potentially Fst.
##The script also creates fastSTRUCTURE appropriate files and runs structure.py on these files to produce structure output,
##and finally runs a slightly modified distruct.py script called new_distruct.py to create admixture plots for varying K values.

####################
##CONFIGURATION
source $HOME/.bash_profile

cd /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/

##activate GATK env
conda activate /shared/apps/conda/bio2/
########################

#######################
#GATK SELECT VARIANTS TO INCLUDE ONLY TETRAPLOID POPS IN THE WACHAU HYBRID ZONE
gatk SelectVariants -V Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz -sn BZD-01tl -sn BZD-02tl -sn BZD-03tl -sn BZD-04tl -sn BZD-05tl -sn BZD-06tl -sn BZD-07tl -sn BZD-08tl -sn FRE-03tl -sn FRE-04tl -sn FRE-05tl -sn FRE-06tl -sn FRE-07tl -sn FRE-08tl -sn HAB-01tl -sn HAB-02tl -sn HAB-03tl -sn KAG-01tl -sn KAG-02tl -sn KAG-03tl -sn KAG-04tl -sn KAG-05tl -sn KEH-05tl -sn KEH-06tl -sn KEH-07tl -sn KEH-08tl -sn KEH-09tl -sn KEH-10tl -sn LIC-01tl -sn LIC-02tl -sn LIC-03tl -sn LOI-01tl -sn LOI-02tl -sn LOI-03tl -sn MAU-01tl -sn MAU-02tl -sn MAU-03tl -sn MOD-01tl -sn MOD-02tl -sn MOD-03tl -sn MOD-04tl -sn OCH-03tl -sn OCH-04tl -sn OCH-05tl -sn OCH-06tl -sn OCH-07tl -sn OCH-08tl -sn PIL-01tl -sn PIL-02tl -sn PIL-03tl -O ./220324_whole_pipeline_VCFs/220324_tetraploids_only.vcf.gz 
######################

conda deactivate
###########################
##activate faststructure env
conda activate /shared/conda/faststructure/

##########################
##RUN POLY_SFS TO GET SITE FREQUENCY SPECTRA FOR ALL TETRAPLOIDS
./filtered_VCFs_for_faststructure/poly_sfs -vcf ./220324_whole_pipeline_VCFs/220324_tetraploids_only_copy.vcf -inds ./220324_whole_pipeline_VCFs/individuals.txt -mis 0.8 > ./220324_whole_pipeline_VCFs/220324_SFS_output.sfs

##########################
##RUN POLY_FREQ SCRIPT FROM TUOMAS HÄMÄLÄ (2023) TO ESTIMATE ALLELE FREQUENCIES FROM MIXED PLOIDY VCF FILE
./filtered_VCFs_for_faststructure/poly_freq -vcf ./220324_whole_pipeline_VCFs/220324_tetraploids_only_copy.vcf -pops ./220324_whole_pipeline_VCFs/populations.txt -mis 0.8 -maf 0.02 -out 0 > ./220324_whole_pipeline_VCFs/220324_poly_freq_output

######
##poly_fst for estimating pairwise Fst and Dxy from mixed ploidy VCFs

###################
##RUN COCHLEARIA_CREATE_STRUCTURE.PY TO create faststructure appropriate files .str
python3 ./scripts/Cochlearia_create_structure_file.py -v ./220324_whole_pipeline_VCFs/ -o 220324_structure_files -s true 

##remove the first and last lines from the diploidized structure file
sed "1d" ./220324_whole_pipeline_VCFs/vcf_to_str/220324_structure_files.StructureInputDiploidized.str | sed "$d" ./220324_whole_pipeline_VCFs/vcf_to_str/220324_structure_files.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str

#########
##reorder the structure files into alphabetical order
grep "BZD" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/BZD.str

grep "FRE" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/FRE.str

grep "HAB" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/HAB.str

grep "KAG" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/KAG.str

grep "KEH" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/KEH.str

grep "LIC" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/LIC.str

#grep "LOI" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/LOI.str

#grep "MAU" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/MAU.str

grep "MOD" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/MOD.str

grep "OCH" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/OCH.str

grep "PIL" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/PIL.str
##########################

##concatenate into a reordered .str file
cd ./220324_whole_pipeline_VCFs/

cat BZD.str FRE.str HAB.str KAG.str KEH.str LIC.str LOI.str MAU.str MOD.str OCH.str PIL.str > ./220324_reordered_structure.str 

cd ../
######################

########################
##RUN STRUCTURE.PY WITH DIFFERENT K VALUES (2-9)
python /shared/conda/faststructure/bin/structure.py -K 2 --input=./220324_whole_pipeline_VCFs/220324_reordered_structure --output=./220324_whole_pipeline_VCFs/220324_K2_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 3 --input=./220324_whole_pipeline_VCFs/220324_reordered_structure --output=./220324_whole_pipeline_VCFs/220324_K3_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 4 --input=./220324_whole_pipeline_VCFs/220324_reordered_structure --output=./220324_whole_pipeline_VCFs/220324_K4_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 5 --input=./220324_whole_pipeline_VCFs/220324_reordered_structure --output=./220324_whole_pipeline_VCFs/220324_K5_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 6 --input=./220324_whole_pipeline_VCFs/220324_reordered_structure --output=./220324_whole_pipeline_VCFs/220324_K6_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 7 --input=./220324_whole_pipeline_VCFs/220324_reordered_structure --output=./220324_whole_pipeline_VCFs/220324_K7_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 8 --input=./220324_whole_pipeline_VCFs/220324_reordered_structure --output=./220324_whole_pipeline_VCFs/220324_K8_out --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 9 --input=./220324_whole_pipeline_VCFs/220324_reordered_structure --output=./220324_whole_pipeline_VCFs/220324_K9_out --format=str --full
###########################

##########################
##RUN NEW_DISTRUCT.PY TO CREATE STRUCTURE/ADMIXTURE PLOTS
python ./scripts/new_distruct.py -K 2 --input=./220324_whole_pipeline_VCFs/220324_K2_out --output=./220324_whole_pipeline_VCFs/220324_K2_plot --popfile=./220324_whole_pipeline_VCFs/populations.txt --title="Arabidopsis lyrata admixture: K=2"

python ./scripts/new_distruct.py -K 3 --input=./220324_whole_pipeline_VCFs/220324_K3_out --output=./220324_whole_pipeline_VCFs/220324_K3_plot --popfile=./220324_whole_pipeline_VCFs/populations.txt --title="Arabidopsis lyrata admixture: K=3"

python ./scripts/new_distruct.py -K 4 --input=./220324_whole_pipeline_VCFs/220324_K4_out --output=./220324_whole_pipeline_VCFs/220324_K4_plot --popfile=./220324_whole_pipeline_VCFs/populations.txt --title="Arabidopsis lyrata admixture: K=4"

python ./scripts/new_distruct.py -K 5 --input=./220324_whole_pipeline_VCFs/220324_K5_out --output=./220324_whole_pipeline_VCFs/220324_K5_plot --popfile=./220324_whole_pipeline_VCFs/populations.txt --title="Arabidopsis lyrata admixture: K=5"

python ./scripts/new_distruct.py -K 6 --input=./220324_whole_pipeline_VCFs/220324_K6_out --output=./220324_whole_pipeline_VCFs/220324_K6_plot --popfile=./220324_whole_pipeline_VCFs/populations.txt --title="Arabidopsis lyrata admixture: K=6"

python ./scripts/new_distruct.py -K 7 --input=./220324_whole_pipeline_VCFs/220324_K7_out --output=./220324_whole_pipeline_VCFs/220324_K7_plot --popfile=./220324_whole_pipeline_VCFs/populations.txt --title="Arabidopsis lyrata admixture: K=7"

python ./scripts/new_distruct.py -K 8 --input=./220324_whole_pipeline_VCFs/220324_K8_out --output=./220324_whole_pipeline_VCFs/220324_K8_plot --popfile=./220324_whole_pipeline_VCFs/populations.txt --title="Arabidopsis lyrata admixture: K=8"

python ./scripts/new_distruct.py -K 9 --input=./220324_whole_pipeline_VCFs/220324_K9_out --output=./220324_whole_pipeline_VCFs/220324_K9_plot --popfile=./220324_whole_pipeline_VCFs/populations.txt --title="Arabidopsis lyrata admixture: K9"
##########################

#conda deactivate 

echo "DONE!!"

##lastlinee
