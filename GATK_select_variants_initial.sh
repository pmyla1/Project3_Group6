#######################
##This script was written by Luke Archer (2024) to index a mixed ploidy VCF file and subsequently 
##subsequently select only the populations that were suspected to be 50/50 hybrids, 100% lyrata, or
##100% arenosa. THIS SCRIPT WAS RUN ON THE HPC.

#######################
##CONFIGURATION
##source bash profile
source $HOME/.bash_profile

##cd to lyrata_VCF
cd /workhere/students_2023/Rot2_DYL/LEVI_project/lyrata_VCF/

#activate gatk conda environment
conda activate /shared/apps/conda/bio2/
#########################

########################
##INDEX VCF AND CREATE SEQUENCE DICTIONARY
##make an indexed VCF file using gatk IndexFeatureFile
gatk IndexFeatureFile -I Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz 

##create a sequence dictionary file for the lyrata.fasta reference 
gatk CreateSequenceDictionary -R lyrata.fasta

#DEACTIVATE GATK ENVIRONMENT 
conda deactivate
########################

########################
##INDEX THE REFERENCE FASTA USING SAMTOOLS
##now create the fasta index file using samtools faidx
conda deactivate

#first activate the samtools environment
conda activate samtools

#use samtools faidx 
samtools faidx lyrata.fasta

#now deactivate samtools environment and activate the shared gatk environment once again
conda deactivate
#######################

########################
##ACTIVATE GATK ENVIRONMENT AND SELECT POPULATIONS OF INTEREST
conda activate /shared/apps/conda/bio2/

#use gatk SelectVariants to select samples we want
gatk SelectVariants -R lyrata.fasta -V Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz -sn HAB-01tl -sn HAB-02tl -sn HAB-03tl -sn FRE-03tl -sn FRE-04tl -sn FRE-05tl -sn FRE-06tl -sn FRE-07tl -sn FRE-08tl -sn FRE-01tl -sn FRE-02tl -sn OCH-01tl -sn OCH-02tl -sn OCH-03tl -sn OCH-04tl -sn OCH-05tl -sn OCH-06tl -sn OCH-07tl -sn OCH-08tl -sn KEH-01tl -sn KEH-04tl -sn KEH-05tl -sn KEH-06tl -sn KEH-07tl -sn KEH-08tl -sn KEH-09tl -sn KEH-10tl -sn LIC-01tl -sn LIC-02tl -sn LIC-03tl -sn MOD-01tl -sn MOD-02tl -sn MOD-03tl -sn MOD-04tl -sn KAG-01tl -sn KAG-02tl -sn KAG-03tl -sn KAG-04tl -sn KAG-05tl -sn BZD-01tl -sn BZD-02tl -sn BZD-03tl -sn BZD-04tl -sn BZD-05tl -sn BZD-06tl -sn BZD-07tl -sn BZD-08tl -sn PIZ-06dl -sn PIZ-08dl -sn PIZ-09dl -sn PIZ-11dl -sn PIZ-02dl -sn PIZ-03dl -sn PIL-01tl -sn PIL-02tl -sn PIL-03tl -sn LOI-01tl -sn LOI-02tl -sn LOI-03tl -sn MAU-01tl -sn MAU-02tl -sn MAU-03tl -sn ANI-01dl -sn ANI-02dl -sn BRA-01dl -sn BRA-02dl -sn BRA-03dl -sn BRA-04dl -sn BRA-05dl -sn BRA-06dl -sn BRA-07dl -sn BRA-07dl -sn BRA-08dl -sn BRA-09dl -sn BRA-10dl -O new_filtered_pops.vcf.gz

echo "DONE!!!!!"

##lastline
