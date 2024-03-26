################
###This script can be used to reorder the structure input for easier plotting

#######
##configuration
source $HOME/.bash_profile

cd /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/filtered_VCFs_for_faststructure/200324_tets_only_VCFs/

conda activate /shared/conda/faststructure/
########

########
##reorder the populations according to desired plotting layout

#grep "BZD" ./200324_tets_first_last_removed.str > ./BZD.str

#grep "FRE" ./200324_tets_first_last_removed.str > ./FRE.str

#grep "HAB" ./200324_tets_first_last_removed.str > ./HAB.str

#grep "KEH" ./200324_tets_first_last_removed.str > ./KEH.str

#grep "LIC" ./200324_tets_first_last_removed.str > ./LIC.str

#grep "MOD" ./200324_tets_first_last_removed.str > ./MOD.str

#grep "OCH" ./200324_tets_first_last_removed.str > ./OCH.str
########

#########
##concatenate into a reordered structure file
#cat ./BZD.str ./FRE.str ./HAB.str ./KEH.str ./LIC.str ./MOD.str ./OCH.str > ./220324_reordered_structure.str
########

############
##run structure.py again with different K values using the reordered .str file
#python /shared/conda/faststructure/bin/structure.py -K 2 --input=./220324_reordered_structure --output=./220324_structure/K2_admixture --format=str --full 

#python /shared/conda/faststructure/bin/structure.py -K 3 --input=./220324_reordered_structure --output=./220324_structure/K3_admixture --format=str --full

#python /shared/conda/faststructure/bin/structure.py -K 4 --input=./220324_reordered_structure --output=./220324_structure/K4_admixture --format=str --full

#python /shared/conda/faststructure/bin/structure.py -K 5 --input=./220324_reordered_structure --output=./220324_structure/K5_admixture --format=str --full

#python /shared/conda/faststructure/bin/structure.py -K 6 --input=./220324_reordered_structure --output=./220324_structure/K6_admixture --format=str --full

#python /shared/conda/faststructure/bin/structure.py -K 7 --input=./220324_reordered_structure --output=./220324_structure/K7_admixture --format=str --full
#########

###########
##run distruct2.3.py on these to produce admixture plots
#python ../../scripts/distruct2.3.py -K 2 --input=./220324_structure/K2_admixture.meanQ --output=./220324_structure/K2_plot --poporder=../poporder.txt --popcolors=../popcolors.txt --title="A. lyrata admixture: K=2"
######################
##try new_distruct.py 
python ../../scripts/new_distruct.py -K 2 --input=./220324_structure/K2_admixture --output=./220324_structure/K2_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=2"

python ../../scripts/new_distruct.py -K 3 --input=./220324_structure/K3_admixture --output=./220324_structure/K3_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=3"

python ../../scripts/new_distruct.py -K 4 --input=./220324_structure/K4_admixture --output=./220324_structure/K4_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=4"

python ../../scripts/new_distruct.py -K 5 --input=./220324_structure/K5_admixture --output=./220324_structure/K5_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=5"

python ../../scripts/new_distruct.py -K 6 --input=./220324_structure/K6_admixture --output=./220324_structure/K6_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=6"

python ../../scripts/new_distruct.py -K 7 --input=./220324_structure/K7_admixture --output=./220324_structure/K7_plot --popfile=../new_pops_200324.txt --title="A. lyrata admixture: K=7"
###################

###################
##distruct2.3.py - CANNOT GET THIS WORKING
python ../../scripts/distruct2.3.py -K 2 --input=./220324_structure/K2_admixture.meanQ --output=./220324_structure/K2_plot --popfile=../new_pops.txt --title="A. lyrata	admixture: K=2"

#python ../../scripts/distruct2.3.py -K 3 --input=./220324_structure/K3_admixture.meanQ --output=./220324_structure/K3_plot.png --poporder=../poporder.txt --popcolors=../popcolors.txt --title="A. lyrata admixture: K=3"

#python ../../scripts/distruct2.3.py -K 4 --input=./220324_structure/K4_admixture.meanQ --output=./220324_structure/K4_plot.png --popfile=../distruct2_pops.txt --poporder=../poporder.txt --popcolors=../popcolors.txt --title="A. lyrata admixture: K=4"

#python ../../scripts/distruct2.3.py -K 5 --input=./220324_structure/K5_admixture.meanQ --output=./220324_structure/K5_plot.png --popfile=../distruct2_pops.txt --poporder=../poporder.txt --popcolors=../popcolors.txt --title="A. lyrata admixture: K=5"

#python ../../scripts/distruct2.3.py -K 6 --input=./220324_structure/K6_admixture.meanQ --output=./220324_structure/K6_plot.png --popfile=../distruct2_pops.txt --poporder=../poporder.txt --popcolors=../popcolors.txt --title="A. lyrata admixture: K=6"

#python ../../scripts/distruct2.3.py -K 7 --input=./220324_structure/K7_admixture.meanQ --output=./220324_structure/K7_plot.png --popfile=../distruct2_pops.txt --poporder=../poporder.txt --popcolors=../popcolors.txt --title="A. lyrata admixture: K=7"
###########

conda deactivate

echo "DONE!!"

##lastline



