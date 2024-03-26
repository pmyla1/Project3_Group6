########################
###This script was written by Luke Archer (2024) and can be used to reorder the structure input for easier plotting

##################
##CONFIGURATION
source $HOME/.bash_profile

cd /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/filtered_VCFs_for_faststructure/200324_tets_only_VCFs/

conda activate /shared/conda/faststructure/
###################

####################
##REORDERING THE POPULATIONS ACCORDING TO THE DESIRED PLOTTING LAYOUT

grep "BZD" ./200324_tets_first_last_removed.str > ./BZD.str

grep "FRE" ./200324_tets_first_last_removed.str > ./FRE.str

grep "HAB" ./200324_tets_first_last_removed.str > ./HAB.str

grep "KEH" ./200324_tets_first_last_removed.str > ./KEH.str

grep "LIC" ./200324_tets_first_last_removed.str > ./LIC.str

grep "MOD" ./200324_tets_first_last_removed.str > ./MOD.str

grep "OCH" ./200324_tets_first_last_removed.str > ./OCH.str
###################

####################
##concatenate into a reordered structure file
cat ./BZD.str ./FRE.str ./HAB.str ./KEH.str ./LIC.str ./MOD.str ./OCH.str > ./220324_reordered_structure.str
#####################

#####################
##run structure.py with different K values using the reordered .str file
python /shared/conda/faststructure/bin/structure.py -K 2 --input=./220324_reordered_structure --output=./220324_structure/K2_admixture --format=str --full 

python /shared/conda/faststructure/bin/structure.py -K 3 --input=./220324_reordered_structure --output=./220324_structure/K3_admixture --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 4 --input=./220324_reordered_structure --output=./220324_structure/K4_admixture --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 5 --input=./220324_reordered_structure --output=./220324_structure/K5_admixture --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 6 --input=./220324_reordered_structure --output=./220324_structure/K6_admixture --format=str --full

python /shared/conda/faststructure/bin/structure.py -K 7 --input=./220324_reordered_structure --output=./220324_structure/K7_admixture --format=str --full
###################

###################
##run new_distruct.py with different K values
python ../../scripts/new_distruct.py -K 2 --input=./220324_structure/K2_admixture --output=./220324_structure/K2_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=2"

python ../../scripts/new_distruct.py -K 3 --input=./220324_structure/K3_admixture --output=./220324_structure/K3_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=3"

python ../../scripts/new_distruct.py -K 4 --input=./220324_structure/K4_admixture --output=./220324_structure/K4_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=4"

python ../../scripts/new_distruct.py -K 5 --input=./220324_structure/K5_admixture --output=./220324_structure/K5_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=5"

python ../../scripts/new_distruct.py -K 6 --input=./220324_structure/K6_admixture --output=./220324_structure/K6_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=6"

python ../../scripts/new_distruct.py -K 7 --input=./220324_structure/K7_admixture --output=./220324_structure/K7_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=7"
#######################

conda deactivate

echo "DONE!!"

##lastline



