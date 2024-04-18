###This python script was written by Luke Archer (2024) to trim the arenosa_632.txt
##to include only SNPs found in both arenosa_632.txtand in lyrata_272_with_some_hybrids.txt
##Please change directory to the same folder as the text files to be merged

###############
##IMPORT MODULES 
import pandas as pd

##############
##STEP 1 - DEFINE THE INPUT AND OUTPUT FILES
#define input and output files
in1=input("Please enter the filename of your A. arenosa file: ")

in2=input("Please enter the filename of your A. lyrata file: ")

#create an empty output file 
output=[]
##############

##############
##STEP 2 - USE PANDAS MODULE TO READ THE INPUT FILES 
##read the input files with pd.read_csv()
input1=pd.read_csv(in1,sep="\t")

input2=pd.read_csv(in2,sep="\t")

##convert input files to a pandas dataframe for the pd.merge() to work
input1=pd.DataFrame(input1)

input2=pd.DataFrame(input2)

##use pd.merge with an inner join (for the intersection between 2 files) i.e. common values of CHROM and POS
output=pd.merge(input1,input2,how='inner',on=['CHROM','POS'],sort=True)

##rename the AF columns to AF_arenosa and AF_lyrata
output=output.rename(columns={'AF_x':'AF_arenosa','AF_y':'AF_lyrata'})

##drop allele count (AC_x, AC_y) allele number (AN_x, AC_y), REF, and ALT columns from output to retain allele frequencies only
output=output.drop(['AC_x','AC_y','AN_x','AN_y','REF_x','REF_y','ALT_x','ALT_y'],axis=1)

##write the output to_csv
output.to_csv('CommonSNPs_lyrata_arenosa.tsv',sep='\t',index=True)

##check by printing the output to screen
print(output)

print("Thank you for using CommonSNPsIdentifier.v2, have a good day!!")
##################

##lastline
