##########CONFIG
source $HOME/.bash_profile

cd /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/filtered_VCFs_for_faststructure/190324_filtered_VCFs/

conda activate /shared/conda/faststructure/
##############

##################
##convert VCFs to structure format with Cochlearia_create_structure_file.py
#python3 /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/scripts/Cochlearia_create_structure_file.py -v /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/filtered_VCFs_for_faststructure/190324_filtered_VCFs/ -o structure_modified_VCF -s true
###############

##############
##run structure.py with varying K values
#python /shared/conda/faststructure/bin/structure.py -K 10 --input=./vcf_to_str/structure_modified_VCF_no_heeaders_lastline_removed.StructureInputDiploidized --output=./faststructure_K10_190324 --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 9 --input=./vcf_to_str/structure_modified_VCF_no_heeaders_lastline_removed.StructureInputDiploidized --output=./faststructure_K9_200324 
#############

############
##run new_distruct.py to create admixture plots
python /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/scripts/new_distruct.py -K 10 --input=./faststructure_K10_190324 --output=./faststructure_K10_190324_plot --popfile=/workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/filtered_VCFs_for_faststructure/populations_individuals.txt --title="Arabidopsis lyrata admixture: K=10"
#############

conda deactivate

##########
echo "DONE!!"

##lastline

