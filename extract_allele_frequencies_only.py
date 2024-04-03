##This script was written by Luke Archer (2024) and takes the output tsv file from 250324_combined_lyrata_arenosa.py
##and creates a new pandas dataframe including the allele frequency columns only.  

##import pandas module 
import pandas as pd 

##create an empty output file
output=[]

in3=input("Please enter your input tsv file for allele frequency extraction: ")

##use pd.read_csv on the input file
input=pd.read_csv(in3,sep="\t")

##drop the alternate and reference allele columns from the dataframe
output=input.drop(['ALT_x','ALT_y','REF_x','REF_y'],axis=1)

##write the output file to a new tsv file
output.to_csv('AFs_only_lyrata_arenosa.tsv',sep="\t",index=True)

##print to screen
print(output)

print("Thank you for using extract_allele_frequencies_only.py :)")

######lastline
