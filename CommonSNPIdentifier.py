###python script to trim the arenosa_632.txt file to include only SNPs found in both arenosa_632.txt
##and in lyrata_272_with_some_hybrids.txt
##Please change directory to the same folder as the text files to be merged

##import pandas module and argparse 
import pandas as pd
#import argparse

##set up the argument parser
#parser=argparse.ArgumentParser(prog="CommonSNPIdentifier",
#                               description="Compares SNPs between 2 files and retains only common SNPs")

##add arguments
#parser.add_argument("-i1","--input1",help="Path to input file 1",type=str,required=True)

#parser.add_argument("-i2","--input2",help="Path to input file 2",type=str,required=True)

#parser.add_argument("-o","--output",help="Output file prefix",type=str,required=True)

#args=parser.parse_args()
#######
##set up the pandas arguments
#in1=open(args.i1,"r")
#in2=open(args.i2,"r")
#output=open(args.o,"w")

##############
##STEP 1 - DEFINE THE INPUT AND OUTPUT FILES
#define input and output files
in1=input("Please enter the filename of your first input file: ")

in2=input("Please enter the filename of your second input file: ")

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
##use pd.merge with an inner join (for the intersection between 2 files)
output=pd.merge(input1,input2,how='outer',on=['CHROM','POS','REF','ALT','AF','AC','AN'],sort=True)

output.to_csv('CommonSNPs_lyrata_arenosa.tsv',sep='\t',index=True)

##check by printing the output to screen
print(output)

print("Thank you for using CommonSNPsIdentifier, have a good day!!")
##################

##lastline