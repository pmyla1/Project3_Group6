######This script will take the output tsv file from 250324_combined_lyrata_arenosa.py
######and extract only the allele frequency columns 

##import argparse module to define input and output file
import pandas as pd 

##set up the argument  parser
#parser=argparse.ArgumentParser(description='Retains only allele frequency data from a data frame')

##add arguments
#parser.add_argument('-i',metavar='--input',type=list,nargs='+',required=True,help='Input filename')

#parser.add_argument('-o',metavar='--output',type=str,required=False,help='Output file prefix')

#args=parser.parse_args()

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