###########################
##This script was provided by Levi Yant (2024) and can be used to create a fastSTRUCTURE appropriate
##structure file (ending in .str) from an input VCF file
##usage: python Cochlearia_create_structure.py -v <path/to/VCF> -o <output/prefix> -s true 
##the -s true flag creates pseudo-diploid data from polyploid data

import os, sys, subprocess, argparse, random, numpy, csv, scipy, gzip

def transpose(i, o=None, d=','):
    f = open (i, 'r')
    file_contents = f.readlines ()
    f.close ()

    out_data = map((lambda x: d.join([y for y in x])),
                   zip(* [x.strip().split(d) for x in file_contents if x]))
    if o:
        f = open (o,'w')

        # here we map a lambda, that joins the first element of a column, the 
        # header, to the rest of the members joined by a comma and a space. 
        # the lambda is mapped against a zipped comprehension on the 
        # original lines of the csv file. This groups members vertically
        # down the columns into rows. 
        f.write ('\n'.join (out_data))
        f.close ()

    return out_data

# create variables that can be entered in the command line
parser = argparse.ArgumentParser(description = 'This program selects one random SNP per defined window, and outputs structure-like format')
parser.add_argument('-v', type = str, metavar = 'vcf_path', required = True, help = 'path to directory with vcfs')
parser.add_argument('-o', type = str, metavar = 'Output_Prefix', required = True, help = 'Vcfs retain original scaffold name but with this prefix')
parser.add_argument('-s', type = str, metavar = 'Subset?', required = False, default = 'false', help = 'if true, this will subsample polyploid data to create psuedo-diploid data')
args = parser.parse_args()

#Create folder for output if it doesn't already exist
if os.path.exists(args.v + '/vcf_to_str/') == False:
    os.mkdir(args.v + '/vcf_to_str/')
vcf_list = []
print(vcf_list)

#get names of vcf files in args.v directory
for file in os.listdir(args.v):
    if file[-4:] == ".vcf":
        vcf_list.append(file)

markernames=[]
tot_sites=0

#create output and temp files
structtempfile= open(args.v+ '/vcf_to_str/'+args.o+"TransposedStruct.str",'w')
subtempfile= open(args.v+ '/vcf_to_str/'+args.o+"TransposedStructSubSample.str",'w')
structfile= open(args.v+ '/vcf_to_str/'+args.o+".StructureInpu.str",'w')
structfile.write("0\t0\t0\t0\t\t\t")

#create aditional file if diploidizing
if args.s=='true':
    subfile= open(args.v+ '/vcf_to_str/'+args.o+".StructureInputDiploidized.str",'w')
    subfile.write("0\t0\t0\t0\t\t\t")

#cycle over vcf's in args.v directory
for iii,vcf in enumerate(vcf_list):
    src = open(args.v + vcf)

    line_num = 0
    chosen_lines = []
    site_num = 0
    names=[]
    ploidy=[]
    first_site=True
    
    # evaluate contents of each line of input file
    for line_idx, line in enumerate(src): #Cycle over lines in the VCF file
        if line_idx % 10000 == 0:
            print(line_idx)
        cols = line.replace('\n', '').split('\t') #Split each line of vcf
        if len(cols) < 2: # This should be info just before header
            pretend_function=1
        elif cols[0] == "#CHROM": #This should be header
            for j in cols[9:]: #get names of individuals in vcf
                names.append(j)
        else: #parse important info from each line    
            genos = []
            dgenos=[]
            scaff = cols[0]
            position = int(cols[1])
            ref_base = cols[3]
            alt_base = cols[4]
            info = cols[7].split(";")
            AN = float(info[2].split("=")[1])
            AC = float(info[0].split("=")[1])
            newsite=[]
            #Convert alleles to Structure input
            if ref_base == "A":
                ref_base = 1
            elif ref_base == "T":
                ref_base = 2
            elif ref_base == "G":
                ref_base = 3
            elif ref_base == "C":
                ref_base = 4
    
            if alt_base == "A":
                alt_base = 1
            elif alt_base == "T":
                alt_base = 2
            elif alt_base == "G":
                alt_base = 3
            elif alt_base == "C":
                alt_base = 4
            #initial setup of infor for output files
            if first_site == True and iii == 0:
            	#Determine ploidy of each individual
                for i,geno in enumerate(cols[9:]):
                    geno=geno.split(":")[0]
                    geno=geno.split("/")
                    ploidy.append(len(geno))
                #Write individual name information for the temporary file that is to be transposed.
                for j,item in enumerate(names):
                    for jj in range(0,ploidy[j]):
                        structtempfile.write("%s\t" % item)
                structtempfile.write("\n")
                #Write PopFlag for each individual which is a unique integer for each population.
                for j,item in enumerate(names):
                    if j==0: 
                            oldname=item[:3]
                            popcount=1
                            for jj in range(0,ploidy[j]):
                                structtempfile.write("%s\t" % str(popcount))
                    else:
                        if oldname!=item[:3]:
                            popcount+=1
                            oldname=item[:3]
                            for jj in range(0,ploidy[j]):
                                structtempfile.write("%s\t" % str(popcount))
                        else:
                            for jj in range(0,ploidy[j]):
                                structtempfile.write("%s\t" % str(popcount))
                structtempfile.write("\n")
                #Same task as steps immediately above, but adjusted to accomodate subsetting.
                if args.s=='true':
                    for j,item in enumerate(names):
                        for jj in range(0,2):
                            subtempfile.write("%s\t" % item)
                    subtempfile.write("\n")
                    #Write PopFlag for each individual which is a unique integer for each population.
                    for j,item in enumerate(names):
                        if j==0: 
                                oldname=item[:3]
                                popcount=1
                                for jj in range(0,2):
                                    subtempfile.write("%s\t" % str(popcount))
                        else:
                            if oldname!=item[:3]:
                                popcount+=1
                                oldname=item[:3]
                                for jj in range(0,2):
                                    subtempfile.write("%s\t" % str(popcount))
                            else:
                                for jj in range(0,2):
                                    subtempfile.write("%s\t" % str(popcount))
                    subtempfile.write("\n")

                num_alleles=float(sum(ploidy))
                first_site=False
                    
            #to each line of data
            if cols[0] != "#CHROM" and len(cols) > 2:
                line_num += 1
                site_num = 1
                markernames.append(str(cols[0])+"_"+str(cols[1]))
                #select sites here
                site = cols
                for geno in site[9:]:
                    geno=geno.replace("|", "/")
                    geno=geno.split(":")
                    geno=geno[0].split("/")
                    for allele in geno:
                        if allele == '.':
                            genos.append(-9)
                        elif allele == '0':
                            genos.append(ref_base)
                        elif allele == '1':
                            genos.append(alt_base)
                        else:
                            print("allele not matched")
                    for allele in numpy.random.choice(geno,2,replace = False):
                        if allele == '.':
                            dgenos.append(-9)
                        elif allele == '0':
                            dgenos.append(ref_base)
                        elif allele == '1':
                            dgenos.append(alt_base)
                        else:
                            print("allele not matched")
                    # filling chosen_sites array with random choices
                for item in genos:
                    structtempfile.write("%s\t" % item)       
                structtempfile.write("\n")
                if args.s == 'true':
                    for item in dgenos:
                        subtempfile.write("%s\t" % item)       
                    subtempfile.write("\n")
                tot_sites+=1


#This section transposes the temporary files created during this replicate and adds headers, so that they are formatted according to structure
structtempfile.close()
subtempfile.close()

#Transposes file
jj=transpose(i=args.v+ '/vcf_to_str/'+args.o+"TransposedStruct.str",d="\t",)
kk=transpose(i=args.v+ '/vcf_to_str/'+args.o+"TransposedStructSubSample.str",d="\t",)

#Write header names for each marker
for marker in markernames:
    structfile.write("%s\t" % marker)
    if args.s == 'true':
        subfile.write("%s\t" % marker)
structfile.write("\n0\t0\t0\t0\t")
if args.s == 'true':
    subfile.write("\n0\t0\t0\t0\t")
for item in jj:
    structfile.write("%s" % item)       
    structfile.write("\n0\t0\t0\t0\t")
if args.s=='true':
    for item in kk:
        subfile.write("%s" % item)       
        subfile.write("\n0\t0\t0\t0\t")

#remove the temporary files that contained the info that needed to be transposed
os.remove(args.v+ '/vcf_to_str/'+args.o+"TransposedStruct.str")
os.remove(args.v+ '/vcf_to_str/'+args.o+"TransposedStructSubSample.str")

