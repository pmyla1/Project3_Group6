# Gene Flow & Hybridisation Between *Arabidopsis arenosa* & *Arabidopsis lyrata* in the Austrian Forealps

This Github page offers an introduction into manipulating **mixed ploidy VCF files** and how to calculate various population genetics metrics from polyploid VCFs such as **site frequency spectra**, and **Fst**.

# Introduction

## Polyploidy & Whole Genome Duplication 
**Polyploidy** and **whole genome duplications** (WGD) occur throughout **all kingdoms** of life and are especially **ubiquitous in plants**. WGD is a **major mutation** that disrupts **ionomic**, **cellular**, and **meiotic** processes, and **neo-polyploids** must overcome various challenges including **genomic instability** and **chromosomal mis-segregation** during meiosis (Margburger et al., 2019 [Nature](https://www.nature.com/articles/s41467-019-13159-5)). One of the immediate challenges in neo-polyploids is the **formation of multivalent crossovers** between homologous chromosomes during meiosis, which can lead to **entanglement** and **chromosomal breakage** at anaphase I (Bray et al., 2023). If neo-polyploids **overcome** the initial challenges related to meiosis and genome instability, they can become **established** as a **polyploid lineage**. 


Some of the genetic adaptations to WGD in polyploids have been characterised, however, despite **process-level convergence**, there appears to be **low convergence** at the **gene/orthologue level** (Bray et al, 2023 [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.09.27.559727v2)). For example, many of the genes **under selection in polyploid** lineages of *Cochlearia*, *Cardamine amara*, and *Arabidopsis arenosa* belong to common functional categories involving **DNA repair**, **cell division**, and **ion homeostasis**, but there are **no orthologous genes** in common between lineages (Bray et al., 2023).

There are **two main types** of polyploids: **autopolyploids** (where all subgenomes originate from the **same species** without hybridisation) and **allopolyploids** (where the subgenomes are inherited from **different species** through **hybridisation**). Autopolyploids and allopolyploids display distinct characteristic site frequency spectra (SFS). **Autopolyploid** SFS follow a **Poisson distribution** with a **high proportion** of **low frequency** variants, whereas **allopolyploid** SFS follow a characteristic **trimodal distribution** with **high proportions** of **low**, **intermediate**, and **high frequency** variants.

## Arabidopsis arenosa & Arabidopsis lyrata
***Arabdopsis arenosa*** is a biennial **outcrossing plant** which produces distinct **diploid** and **tetraploid** lineages throughout **Central Europe** and is closely related to the **model species** ***Arabidopsis thaliana*** (Margburger et al., 2019). Similarly, the sister species ***Arabidopsis lyrata*** also forms distinct **diploid** and **tetraploid** lineages across its **distribution** range, and there have been reports of **gene flow** between the **tetraploid** lineages of ***A. arenosa*** and ***A. lyrata*** where these species have **overlapping ranges** (Jørgensen et al., 2011 [BMC Ecology & Evolution](https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-11-346). 

The **rationale** for this project was the discovery of **bidirectional gene flow** between **autotetraploid** lineages of *A. arenosa* and *A. lyrata* which may have **facilitated** the **stabilisation of polyploidy** post WGD (Marburger et al., 2019). In this project, we selected specific **tetraploid populations** of *A. arenosa* and *A. lyrata* from the **Austrian Forealps** and the **Wachau hybrid zone** suspected to be hybrids and assessed the degree of **genetic admixture** using **fastSTRUCTURE**. We aimed to discover whether **hybridisation** between *A. arenosa* and *A. lyrata* has created any **allopolyploid populations** in the **Wachau hybrid zone** consisting of approximately **50:50 admixture**. A **graphical summary** of the process of **speciation**, followed by subsequent **hybridisation** and **adaptive introgression** between *A. arenosa* & *A. lyrata* can be seen **below**.

<img width="946" alt="Screenshot 2024-04-22 at 13 35 32" src="https://github.com/pmyla1/Project3_Group6/assets/151543531/1de2717b-bf10-4f01-aa79-4655cd3c401e">


Source: [Schmickl, R. & Yant, L. (2021). Adaptive introgression: how polyploidy reshapes gene flow landscapes. *New Phytologist* 230(2), 457-461.](https://doi.org/10.1111/nph.17204)


The original VCF file **Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz** contains samples from **diploid** and **tetraploid** *A. arenosa* and *A. lyrata* populations sampled **throughout Europe** including lineages from the **Austrian Forealps** and others from the **Wachau hybrid zone**. The **Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz** was filtered by Levi Yant (2024) and includes **57,613** single nucleotide polymorphisms **(SNPs)** on **Chromosome 1 only**. The **population map** can be accessed via the following link [Population Map](https://www.google.com/maps/d/viewer?mid=1HAhM5y-bYMJbXCtMSZaubk1qe0wX6JI&ll=48.09350708234271%2C15.809612499999968&z=9).

## Some of the Populations in the VCF

 | Species        | Ploidy           | 3-letter pop code(s) |
 | ------------- | :-------------: | -----:|
 | pure *A. arenosa* | 4x | KEH, BZD |
 | pure *A. lyrata* | 4x | KAG, LIC, MOD, MAU |
 | pure *A. lyrata* | 2x | PIZ, PLE |
 | hybrid *A. arenosa* x *A. lyrata* | 4x | FRE, HAB, OCH |  


# **Software Used**

## **1) fastSTRUCTURE**.
An algorithm to infer population structure: sourced from [fastSTRUCTURE](https://rajanil.github.io/fastStructure/). 

Citation: Anil Raj, Matthew Stephens, and Jonathan K. Pritchard. fastSTRUCTURE: Variational Inference of Population Structure in Large SNP Data Sets , (Genetics) June 2014 197:573-589.

## Installation and dependencies

### fastSTRUCTURE requires:
[Numpy](https://numpy.org/)   

[SciPy](https://scipy.org/) 

[Cython](https://cython.org/) 

[GNU Scientific Library](https://www.gnu.org/software/gsl/)  

## Obtaining the Source code from GitHub

You can obtain the source code from the [fastSTRUCTURE](https://rajanil.github.io/fastStructure/) github page by cloning the repository into a new directory called fastSTRUCTURE:
```
## first make a new directory called fastSTRUCTURE (or whatever you want to call the folder)
mkdir ~/fastSTRUCTURE

## change into the new directory
cd ~/fastSTRUCTURE

## clone the github repository for fastSTRUCTURE
git clone https://github.com/rajanil/fastStructure

```
Furthermore, the code can also be obtained with the following wget command

```
wget --no-check-certificate https://github.com/rajanil/fastStructure/archive/master.tar.gz
```

## Building Python extensions

**Before** building the Python extensions for fastSTRUCTURE, identify your **path** to the **library files** ***libgsl.so*** and ***libgslcblas.so***, in addition to the **header file** ***gsl/gsl_sf_psi.h*** that are part of your **GNU Scientific Library installation**. 

If you have used **default options** to donwload and install GSL, the **libraries** (*.so* files) should be located in **/usr/local/lib**, and the **header files** (*.h* files) should be located in **/usr/local/include**. If you have successfully located your GSL library and header files, and they are in the correct location, **add** the following lines of code to your ***.bashrc*** file in your home directory **(~)**.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export CFLAGS="-I/usr/local/include"
export LDFLAGS="-L/usr/local/lib"
```

Then execute the following command to **set** these **environmental variables**
```
source ~/.bashrc
```
Next, to **build** the library **extensions**, you can use the following commands
```
## assuming you created a directory called fastSTRUCTURE
cd ~/fastSTRUCTURE/fastStructure/vars
## run the setup.py script to build the library extensions
python setup.py build_ext --inplace
```

Now, to **compile** the main **`cython`** scripts, run the following code
```
## change directory into your newly created fastStructure folder
cd ~/fastSTRUCTURE/fastStructure
## compile the cython scripts using setup.py
python setup.py build_ext --inplace
```
**Errors** suggesting the **build failed** may be due to using the **wrong compiler** or **incorrectly set environmental variables**.
If you would like to use a **specific gcc compiler**, you can do the following:
```
CC=/path/to/compiler python setup.py build_ext --inplace 
```

## **2) The Genome Analysis Toolkit (GATK) `v4.2.2.0`**. 
For support and documentation click on the following link [GATK](https://software.broadinstitute.org/gatk/) 

## Dependencies for GATK v4
To configure, build, and run GATK:

You require **`JAVA 17`** to run or build GATK.
**`JAVA 17`** can be installed using the following instructions 
1) on **macOS**, you can install the [Homebrew package manager](https://docs.brew.sh/Installation) and then run the following command
   ```
   brew tap homebrew/cask-versions
   ## then to install the Eclipse Foundation's OpenJDK 17 run 
   brew install --cask temurin17
   ```
Furthermore, you must have a **Python** version greater than or equal to **Python 2.6**, and an **R version of at least 3.2.5**

To **build GATK** you must have [git-lfs](https://git-lfs.com/) of **1.1.0 or greater**. This is required to download the large files that are used to build GATK, including the test files. 
After downloading git-lfs, run the following command 
```
git lfs install
```
You also **require Gradle 5.6** - you can use the `./gradlew` script which will download and use an **appropriate Gradle version automatically**.

## Installation of `GATK v4.2.2.0`

## Python dependencies
GATK uses **Conda** to establish and manage the Python environment and dependencies required by the GATK tools that use Python.
Firstly, please **uninstall any previous GATK versions** installed on your device using the following commands:
```
source deactivate
## assuming you called your GATK environment gatk
conda env remove -n gatk
```

Next, **update conda** to the **latest version** for your base environment
```
conda update -n base conda
```

To install the **GATK v4** package

First **create** and then **activate** a Conda environment for GATK using either **Miniconda** or **Conda**.

If you are running from a **zip** or **tar** distribution, then navigate to the **directory** where you have **stored your GATK jars** and the **gatk wrapper** script, ensure the **gatkcondaenv.yml** is present, and run the following command:
```
conda env create -n gatk -f gatkcondaenv.yml

source activate gatk
```
Or if you are running from a **cloned repository**, run the following command:
```
./gradlew localDevCondaEnv

source activate gatk
```

To check if your GATK environment is **properly installed** run: 
```
conda list
```
Which should return a list of packages you have installed, **gatkpythonpackages** should be one of these.

## **3) R and RStudio - R version 4.3.1 (2023-06-16)**

# Installation 

## To install R-4.3.1 on Windows PC:

1) **Uninstall previous versions** of R. 

2) Click on the following link [R installation for Windows](https://cran.r-project.org/bin/windows/base/) 

3) Click **Download R-4.3.1 for Windows**. 

4) Open the installer and follow the instructions using **default options**. 


## To install Rtools-4.3 on Windows PC:

1) **Uninstall previous versions** of Rtools. 

2) Click on the following link [Rtools installation for Windows](https://cran.r-project.org/bin/windows/Rtools/)  

3) Click **RTools 4.3**. 

4) Click **Rtools43 installer**. 

5) Open the installer and follow the instructions using **default options**.


## To install RStudio-4.3.1 on Windows PC:

1) **Uninstall previous versions** of RStudio. 

2) Click on the following link [RStudio installation for Windows](https://posit.co/download/rstudio-desktop/) 

3) Follow the downloads instructions from **Step 2 onwards**. 

4) Open the installer and then follow the instructions using **default options**. 


## To install R version 4.3.1. on Apple macOS

1) **Uninstall previous installations** of R by navigating to your **applications folder** and moving **R** and **XQuartz** to the **Bin**. 

2) To install R, **check which Apple macOS version** you have - briefly, click the **Apple Logo** in the **top left corner** of your screen and then **choose** ***About this mac*** > note what the **Processor line** says. 

3) Click on the following link [R installation for Apple macOS](https://cran.r-project.org/bin/macosx/)  

  
   If your Processor line has **"Intel"** then follow the download **instructions** for **`R-4.3.1-x86_64.pkg`**. 
   
   **Otherwise**, follow the download instructions for **`R-4.3.1-arm64.pkg`**. 
  

4) Open the installer and follow the instructions using **default options**. 

5) To install XQuartz - **download XQuartz** from the following link [XQuartz download for macOS](https://www.xquartz.org/), open the installer using **default options**. 


## To install RStudio-4.3.1 on Apple macOS

1) **Uninstall previous versions** of RStudio. 
2) Click on the following link [RStudio Version 4.3.1 Installation for macOS](https://posit.co/download/rstudio-desktop/) and follow the download instructions from **step 2 onwards**. 
3) Open the installer and follow the instructions using **default options**.

## To install SplitsTree for Phylogenetic Analyses

1) To install **SplitsTree**, go to the **University of Tübingen** Website [SplitsTree](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html) and click on the **`Download`** link.

2) Select the **correct program installer** for you specific device and follow the instructions for installation/download using **default options**.

## To install samtools using Conda

To **install samtools** using **conda** you can use the **following commands**:

   ```
   ## create the samtools environment
   conda create -n samtools
   ## activate the environment
   conda activate samtools
   ## configure the bioconda channel
   conda config --add channels bioconda
   ## configure conda-forge
   conda config --add channels conda-forge
   ## install samtools
   conda install -c bioconda samtools
   ```

## To install bcftools using conda 

To install bcftools using conda, use the following command
```
conda install -c bioconda bcftools
```


# Scripts

## GATK_select_variants_initial.sh

**GATK_select_variants_initial.sh** was executed on the HPC using a **shared GATK environment**, and a custom **samtools environment** created using the install steps above.

The `gatk IndexFeatureFile -I`  command is used to create an **indexed fasta** file from the **Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz**.

```
##make an indexed VCF file using gatk IndexFeatureFile
gatk IndexFeatureFile -I Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz
```

Subsequently, a **sequence dictionary** is created using `gatk CreateSequenceDictionary`, and then `samtools faidx` is used to index the previously created fasta file. 

```
##create a sequence dictionary
gatk CreateSequenceDictionary -R lyrata.fasta

##index lyrata.fasta with samtools faidx
samtools faidx lyrata.fasta
```

Finally, `gatk SelectVariants` with the `-sn` flag is used to **select specific individuals** from **Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz** and to make a new .vcf.gz with only these individuals called **new_pops_filtered.vcf.gz**.

## 220324_whole_pipeline_gatk.sh 

**220324_whole_pipeline_gatk.sh** was executed on the **HPC** using a **shared conda environment** `/shared/conda/shared/` and a **shared gatk environment** `/shared/apps/conda/bio2/`.

`gatk SelectVariants` with the `-sn` flag is used to select only the tetraploid populations of interest (either pure *A. lyrata*, pure *A. arenosa*, or expected to be hybrids) and creates a new filtered vcf called 220324_filtered_pops.vcf.gz.

```
##Use gatk SelectVariants  to select specific populations 
gatk SelectVariants -V Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz -sn BZD-01tl -sn ...
```

The **site-frequency spectra** and the **allele frequencies** are calculated using the unzipped 220324_filtered_pops.vcf and Tuomas Hämälä's (2023) [**poly_sfs.c**](https://github.com/thamala/polySV/blob/main/poly_sfs.c) and [**poly_freq.c**](https://github.com/thamala/polySV/blob/main/poly_freq.c) scripts, respectively.

```
##Run poly_sfs to get site frequency spectra for all tetraploids
./filtered_VCFs_for_faststructure/poly_sfs -vcf ./220324_whole_pipeline_VCFs/220324_tetraploids_only_copy.vcf -inds ./220324_whole_pipeline_VCFs/individuals.txt -mis 0.8 > ./220324_whole_pipeline_VCFs/220324_SFS_output.sfs

##Run poly_freq to estimate allele frequencies from a tetraploid VCF
./filtered_VCFs_for_faststructure/poly_freq -vcf ./220324_whole_pipeline_VCFs/220324_tetraploids_only_copy.vcf -pops ./220324_whole_pipeline_VCFs/populations.txt -mis 0.8 -maf 0.02 -out 0 > ./220324_whole_pipeline_VCFs/220324_poly_freq_output
```

Next, the polyploid VCF is prepared for **fastSTRUCTURE** using Yant et al (2023) [**Cochlearia_create_structure_file.py**](https://github.com/pmyla1/Project3_Group6/blob/main/Cochlearia_create_structure_file.py) script for polyploid data, using the **'-s true'** flag to subsample the data to make a pseudo-diploid structure output. 

```
##run cochlearia_create_structure_file.py to create faststructure appropriate files
python3 ./scripts/Cochlearia_create_structure_file.py -v ./220324_whole_pipeline_VCFs/ -o 220324_structure_files -s true
```

The populations are then **rearranged** into **alphabetical order** for plotting purposes. 

```
##example code to move all BZD individuals into a single .str file
grep "BZD" ./220324_whole_pipeline_VCFs/220324_first_last_removed.StructureInputDiploidized.str > ./220324_whole_pipeline_VCFs/BZD.str
```

Next, the [**structure.py**](https://github.com/rajanil/fastStructure/blob/master/structure.py) script from **fastSTRUCTURE** is used to assess population genetic structure with varying K-values (k=2-9). **Only K=2 command shown**.

```
##Run structure.py with K=2 
python /shared/conda/faststructure/bin/structure.py -K 2 --input=./220324_whole_pipeline_VCFs/220324_reordered_structure --output=./220324_whole_pipeline_VCFs/220324_K2_out --format=str --full
```

Finally, a slightly **modified distruct.py** script called [**new_distruct.py**](https://github.com/pmyla1/Project3_Group6/blob/main/new_distruct.py) is used to produce fastSTRUCTURE **admixture plots** for varying K-values (k=2-9). **Only K=2 shown**.

```
##Run new_distruct.py to create admixture plots, using the output from structure.py as input
python ./scripts/new_distruct.py -K 2 --input=./220324_whole_pipeline_VCFs/220324_K2_out --output=./220324_whole_pipeline_VCFs/220324_K2_plot --popfile=./220324_whole_pipeline_VCFs/populations.txt --title="Arabidopsis lyrata admixture: K=2"
```
## 290324_whole_pipe.sh

The **290324_whole_pipe.sh** was executed on the **HPC** using a **shared conda environment** `/shared/conda/shared/` and a **shared GATK environment** `/shared/apps/conda/bio2/`.

Similarly to the **220324_whole_pipeline_gatk.sh**, `gatk SelectVariants -sn <sample in vcf>` was used to select tetraploid only individuals from **Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz**. 

```
##GATK select variants to include tetraploids only 
gatk SelectVariants -V Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz -sn BZD-01tl ...
```

Site-frequency spectra and allele frequencies are calculated using Tuomas Hämälä's (2023) [**poly_sfs.c**](https://github.com/thamala/polySV/blob/main/poly_sfs.c) and [**poly_freq.c**](https://github.com/thamala/polySV/blob/main/poly_freq.c) scripts, respectively.

```
##################
##poly_sfs for site-frequency spectra calculation
./filtered_VCFs_for_faststructure/poly_sfs -vcf ./290324_whole_pipeline_VCFs/290324_tetraploids_only_copy.vcf -inds ./290324_whole_pipeline_VCFs/individuals_final.txt -mis 0.8 > ./290324_whole_pipeline_VCFs/290324_SFS_output.sfs

#########
##poly_freq for estimating allele frequencies from mixed ploidy vcf files
./filtered_VCFs_for_faststructure/poly_freq -vcf ./290324_whole_pipeline_VCFs/290324_tetraploids_only_copy.vcf -pops ./290324_whole_pipeline_VCFs/populations_final.txt -mis 0.8 -maf 0.02 -out 0 > ./290324_whole_pipeline_VCFs/290324_poly_freq_output
```

The VCF is subsequently prepared for **fastSTRUCTURE** using Yant et al (2023) [**Cochlearia_create_structure_file.py**](https://github.com/pmyla1/Project3_Group6/blob/main/Cochlearia_create_structure_file.py) python script for polyploids, using the **`-s true`** flag to subsample the data to make a pseudo-diploid output structure file. The structure files are arranged according to populations. 

```
##Cochlearia_create_structure.py to create faststructure appropriate files .str
python3 ./scripts/Cochlearia_create_structure_file.py -v ./290324_whole_pipeline_VCFs/ -o 290324_structure_files -s true 
```

The fastSTRUCTURE script [**structure.py**](https://github.com/rajanil/fastStructure/blob/master/structure.py) is used to infer the admixture proportions of the individual samples in **290324_filtered_pops.vcf.gz** with K-values ranging from 2-7. Subsequently, a modified **distruct.py** script called [**new_distruct.py**](https://github.com/pmyla1/Project3_Group6/blob/main/new_distruct.py) is used to produce fastSTRUCTURE plots for K-values ranging from 2-7. **Only K=2 shown**.

```
##Run structure.py when K=2
python /shared/conda/faststructure/bin/structure.py -K 2 --input=./290324_whole_pipeline_VCFs/290324_reordered_structure --output=./290324_whole_pipeline_VCFs/290324_K2_out --format=str --full

##Run new_distruct.py when K=2
python ./scripts/new_distruct.py -K 2 --input=./290324_whole_pipeline_VCFs/290324_K2_out --output=./290324_whole_pipeline_VCFs/290324_K2_plot --popfile=./290324_whole_pipeline_VCFs/populations_final.txt --title="Arabidopsis lyrata admixture: K=2"
```

## 250324_combined_lyrata_arenosa.py 

To determine whether **hybridisation** betweeen *A. arenosa* and *A. lyrata* had produced an **allotetraploid** lineage, we leveraged allele frequency information from **text files** containing **4-fold degenerate SNP** data from a **larger number of samples** from both species.

**250324_combined_lyrata_arenosa.py** is a python script which utilises the [pandas](https://pandas.pydata.org/docs/user_guide/index.html) python package to convert the input **arenosa_632.txt** and **lyrata_272_with_some_hybrids.txt** files into pandas dataframes with the `pd.dataframe()` command. The script requires **user input** and asks for the *A. arenosa* and *A. lyrata* text files in that order.

## Input Text File Structure

 | **CHROM** | **POS** | **REF** | **ALT** | **AF** | **AC** | **AN** |
 | :-------: | :-----: | :-----: | :-----: | :----: | :----: | :----: |
 | scaffold_1 | 32 | C | A | 0.00053 | 1 | 1886 | 
 | scaffold_1 | 38 | A | T | 0.558 | 977 | 1750 |
 | scaffold_1 | 160 | C | A | 0.00974 | 18 | 1848 |

*Key: CHROM, chromosome; POS, position; REF, reference allele; ALT, alternative allele; AF, allele frequency; AC, allele count; AN, allele  number*

The structure of the **final output file** used to calculate the **allele frequency differences** between species can be visualised **below**.

The converted input files are then merged using an **inner join** for the **intersection** between 2 files, based on the **CHROM** and **POS** columns using **`pd.merge(how='inner',on=['CHROM','POS'])`** to only include sites that are **shared**. 

The script then **renames** the **AF_x** & **AF_y** columns as **AF_arenosa** & **AF_lyrata**, respectively, and drops the allele count, allele number, reference, and  alternative columns from the output file to **retain the allele frequencies only**. 

```
import pandas as pd

#define input and output files
in1=input("Please enter the filename of your A. arenosa file: ")
in2=input("Please enter the filename of your A. lyrata file: ")

#create an empty output file 
output=[]

##read the input files with pd.read_csv()
input1=pd.read_csv(in1,sep="\t")
input2=pd.read_csv(in2,sep="\t")

##convert input files to a pandas dataframe for the pd.merge() to work
input1=pd.DataFrame(input1)
input2=pd.DataFrame(input2)

##use pd.merge with an inner join for the intersection/common values of CHROM and POS
output=pd.merge(input1,input2,how='inner',on=['CHROM','POS'],sort=True)

##rename the AF columns to AF_arenosa and AF_lyrata
output=output.rename(columns={'AF_x':'AF_arenosa','AF_y':'AF_lyrata'})

##drop AC, AN, REF, and ALT columns from output to retain allele frequencies only
output=output.drop(['AC_x','AC_y','AN_x','AN_y','REF_x','REF_y','ALT_x','ALT_y'],axis=1)

##write the output to_csv
output.to_csv('Common_SNPs.tsv',sep='\t',index=True)
```
## Final Output File Structure

 | **CHROM** | **POS** | **AF_arenosa** | **AF_lyrata** | 
 | :-------: | :-----: | :------------: | :-----------: |
 | scaffold_1 | 32 | 0.154 | 0.00053 |
 | scaffold_1 | 160 | 0.232 | 0.558 |

*Key: CHROM, chromosome; POS, position; AF_arenosa, allele frequency in A. arenosa; AF_lyrata, allele frequency in A. lyrata* 

## 250324_common_SNPs.R

**250324_common_SNPs.R** was **executed locally** and is an R script that uses ggplot2 to produce a range of different plots from the output of **250324_combined_lyrata_arenosa.py**. After setting your **working directory** to the filepath where your **Common_SNPs.tsv output** file is and **loading** the appropriate **libraries** specified in the script,the data is **read** into R and **converted** into a **tibble**. 

```
##Read Common_SNPs.tsv into R
arenosa_lyrata_AFs<-read_tsv('Common_SNPs.tsv')
##convert to a tibble for easier data extraction
arenosa_lyrata_AFs<-as_tibble(arenosa_lyrata_AFs)
##drop ...1 and. `Unnamed: 0` columns from the df using dplyr::select 
cleaned_arenosa_lyrata_AFs<-dplyr::select(arenosa_lyrata_AFs,-c(...1,`Unnamed: 0`))
```

Firstly, the data is subsetted into each chromosome scaffold using `subset()`. **Chromosome 1 subsetting shown only**.

```
##subset the data per chromosome/scaffold 
chrom1<-subset(cleaned_arenosa_lyrata_AFs,CHROM=='scaffold_1')
#compute mean AFs per site
arenosa_mean1<-mean(chrom1$AF_arenosa)
lyrata_mean1<-mean(chrom1$AF_lyrata)
```

Next, the **genome-wide** allele frequency distributions are plotted for *A. arenosa* and *A. lyrata*, and subsequently, the allele frequency distributions **per scaffold** are plotted using ggplot2. **Genome-wide plots for *A. arenosa* shown only.**

```
###################
##PLOTS OF THE GENOME WIDE INFORMATION
mean_arenosa<-mean(cleaned_arenosa_lyrata_AFs$AF_arenosa)

##GENOME WIDE PLOTS
arenosa<-ggplot(cleaned_arenosa_lyrata_AFs,aes(AF_arenosa))+
  geom_histogram(fill='orange',colour='black',bins=100)+
  geom_vline(xintercept=mean_arenosa,linetype=2,colour=2)+
  theme_bw()+
  labs(title="Genome wide Allele frequency\ndistribution for A. arenosa SNPs",
       x="A. arenosa AF",y='Count')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor.x=element_blank())
```

Next, the scaffolds are filtered for minor allele frequencies >0.025, and the **allele frequency differences** between *A. arenosa* and *A. lyrata* are calculated, and then plotted per scaffold (**CHROM**) as a Manhattan plot using **ggplot2**. **Chromosome 1 plotting shown only**.

```
##This command calculates allele frequency differences between arenosa/lyrata
arenosa_lyrata_MAF2.5$AF_difference<-arenosa_lyrata_MAF2.5$AF_arenosa-arenosa_lyrata_MAF2.5$AF_lyrata

##choose a threshold for plotting outliers/extreme allele frequency differences
chrom1$threshold<-chrom1$AF_difference>0.85

##Plot allele frequency differences along chromosome 1
chrom1_AF_diff<-ggplot(chrom1,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 1 AF differences\nat 4-fold neutral sites",
       x='Chromosome 1 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

```

Subsequently, the **top 1% outlier** allele frequency differences are calculated by using **`dplyr::arrange(desc(AF_difference))`** and then taking the **top 1% rows**. This was performed with an aim of visualizing the **"fixed" allele frequency differences** between *A. arenosa* and *A. lyrata* at common SNPs. **Calculations and plotting for Chromosome 1 shown only**.

```
##CALCULATING SITES WITH THE TOP 1% OUTLIERS IN TERMS OF AF DIFFERENCES
top_AF_diff_chrom1<-chrom1%>%arrange(desc(chrom1$AF_difference))
##calculate the top 1% of rows 
top_1PCT_rows_chrom1<-round(0.01*nrow(top_AF_diff_chrom1)) ##289
##use this to calculate the top outliers
top_1PCT_AF_outliers_chrom1<-top_AF_diff_chrom1%>%top_n(289,AF_difference)

##PLOTTING THE TOP 1% ALLELE FREQUENCY DIFFERENCES
chrom1_1PCT_diff<-ggplot(top_1PCT_AF_outliers_chrom1,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.8)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.7,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 1 top 1% AF\ndifferences at 4-fold sites",
       x='Chromosome 1 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')
```

## 170424_pop_genetics.R

This script was executed locally and can be used to perfom exploratory genetic analysis. Firstly, the **290324_tetraploids_only.vcf.gz** created by **290324_whole_pipe.sh** is loaded into RStudio and converted into a genlight object by the **vcf2genlight.tetra** function. Subsequently, principal component analysis **(PCA)** is performed on the genlight object utilising the **glPcaFast** function. Next, discriminant analysis of principal components **(DAPC)** is performed on the genlight object. Finally, **Nei's genetic distances** are calculated for both the **individuals** and **populations** in the genlight object. 

Subsequently, the **Nei's genetic distance** files can be loaded into **SplitsTree** and used to **create phylogenetic networks**, both of the **individual** samples and the **populations**.

The following code block was provided by Ana (2023)

```
## Calculate Nei's distances between individuals/pops
aa.D.ind<-stamppNeisD(aa.genlight,pop=FALSE) # Nei's 1972 distance between indivs
# export matrix - for SplitsTree
stamppPhylip(aa.D.ind,file="290324_individuals.phy.dst")

aa.D.pop<-stamppNeisD(aa.genlight,pop=TRUE)   # Nei's 1972 distance between pops
# export matrix - for SplitsTree
stamppPhylip(aa.D.pop,file="290324_populations.phy.dst") 

### create the dist objects
colnames(aa.D.ind)<-rownames(aa.D.ind)
aa.D.ind.dist<-as.dist(aa.D.ind,diag=T)
attr(aa.D.ind.dist,"Labels")<-rownames(aa.D.ind)   # name the rows of a matrix  

colnames(aa.D.pop)<-rownames(aa.D.pop) 
aa.D.pop.dist<-as.dist(aa.D.pop,diag=T)
attr(aa.D.pop.dist,"Labels")<-rownames(aa.D.pop)   # name the rows of a matrix  

aa.D.ind # individuals
aa.D.pop # populations
Dgen_ind<-aa.D.ind.dist
Dgen_pop<-aa.D.pop.dist
```

## PCA_alternative_script.R 

This script was used to peform an **alternative PCA** utilising Tuomas Hämälä's (2023) [est_adapt_pca.R](https://github.com/thamala/polySV/blob/main/est_adapt_dist.r) adapted PCA script. 

The input file is the **290324_tetraploids_only.vcf.gz** created by the **290324_whole_pipe.sh** script.

## 190424_poly_fst.sh

**190324_poly_fst.sh** can be used to calculate **pairwise Fst** values between individuals from different populations. Firstly, **grep & bcftools** are used to **create text files** for each population. Subsequently, the **[poly_fst.c](https://github.com/thamala/polySV/blob/main/poly_fst.c)** script from Tuomas Hämälä (2023) is **compiled** and **executed** to calculate pairwise polyploid  Fst scores. **Fst calculations for the BZD-OCH pairwise population contrast are shown below**.

```
#use grep and bcftools to obtain the individuals from each population for the fst scans 
bcftools query -l ./290324_tetraploids_only.vcf.gz | grep "BZD" > ./190424_Fst_populations/BZD_pop.txt

bcftools query -l ./290324_tetraploids_only.vcf.gz | grep "OCH" > ./190424_Fst_populations/OCH_pop.txt

##compile the poly_fst.c script
gcc ../scripts/poly_fst.c -o ./poly_fst -lm

##run poly_fst on the 290324_tets_only_filtered_vcf and BZD vs OCH population contrast
./poly_fst -vcf ./290324_tetraploids_only_copy.vcf -pop1 ./190424_Fst_populations/BZD_pop.txt -pop2 ./190424_Fst_populations/OCH_pop.txt -mis 0.8 > ./190424_Fst_output/BZD_OCH_Fst_output.fst

```

## 190424_polyploid_Fst.R

This script takes the output .fst files from **190424_poly_fst.sh** as input and **creates Manhattan plots** of the Fst scores along **Chromosome 1** for various pairwise population contrasts. **Fst scores** along Chromosome 1 for the **BZD-OCH** population contrast **shown only.**

```
####BZD vs OCH contrast
BZD_OCH_Fst<-read_tsv("BZD_OCH_Fst_output.fst")
##rename columns for easier plotting
BZD_OCH_Fst<-BZD_OCH_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.260870`) 
##make a threshold for colouring points
BZD_OCH_Fst$threshold<-BZD_OCH_Fst$Fst>=0.6
##calculate the mean Fst value
BZD_OCH_Fst$mean_Fst<-mean(BZD_OCH_Fst$Fst)
##plot
BZD_OCH_plot<-ggplot(BZD_OCH_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="BZD vs OCH Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
```
Links to the multipanel **Fst Manhattan plots** can be found **below**: 

[BZD_Fst_scans](https://github.com/pmyla1/Project3_Group6/blob/main/BZD_Fst_scans.png) 

[KEH_Fst_scans](https://github.com/pmyla1/Project3_Group6/blob/main/KEH_Fst_scans.png)

[OCH_Fst_scans](https://github.com/pmyla1/Project3_Group6/blob/main/OCH_Fst_contrasts.png)


# Phylogenetic analyses with SplitsTree

Next, we looked at the **relationships** between the individuals and populations with **SplitsTree**, following the correct **download instructions** for your **specific device** on the **University of Tübingen** website in the following link [SplitsTree](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html). 

Briefly, after converting **290324_tetraploids_only.vcf.gz** into a **genlight** object using the **vcf2genlight** function from the **170424_pop_genetics.R script**, the genlight object can be converted into **Nei's genetic distance** data and subsequently converted into a **phylogentic distance file** (.phy.dst) before being loaded into **SplitsTree**. Both the **individual** and the **population** data can be used to produce **phylogenetic networks** in SplitsTree.

## How to produce phylogenetic networks in SplitsTree

1) After successfully downloading SplitsTree, **open the application** and select **`File` >> `Open`** on the **toolbar**.
2) Navigate to your **.phy.dst** files created in the **170424_pop_genetics.R** script, and select **Open**.
3) Your phylogenetic network should now be loaded into **SplitsTree** and can be saved by selecting **`File` >> `Export Image`**, then selecting **`Format: PNG(*.png)`**, and finally ticking the **`Save visible region`** option.

# Population structure analysis

Alternative population structure figures were created using the GUI **[Omicsspeaks Structure Plot V2.0](http://omicsspeaks.com/strplot2/).** 

The **input file** is a comma separated value **(CSV)** file containing the **individual name** from the VCF (e.g BZD-01tl), followed by the **population** in the VCF (e.g. BZD), followed by the **fastSTRUCTURE output** (two floating numbers) for the genetic admixture proportions when **K=2** (e.g 0.95,0.05).

## Omics Speaks Structure Input File

 |           |         |                |               |
 | :-------: | :-----: | :------------: | :-----------: |
 | BZD-01tl | BZD | 0.95 | 0.05 |
 | BZD-02tl | BZD | 0.97 | 0.03 |
 | BZD-03tl | BZD | 0.98 | 0.02 |
 | BZD-04tl | BZD | 0.96 | 0.04 |

The link to the input file when K=2 can be accessed here [OmicsSpeaks input](https://github.com/pmyla1/Project3_Group6/files/14791197/K2_omics_speaks_input.csv)

## Steps taken to produce OmicsSpeaks Structure Plot

1) Click on the following link for [OmicsSpeaks](http://omicsspeaks.com/) and select the **Structure Plot V2.0** button.
2) **Scroll down** and select the green **Take me to the application** button.
3) Click on the blue **Select input file** near the **top left** of the page.
4) Select your **input .csv** file with the correct structure as shown above, keep the **default options** for now, and click the green **Submit** button.
5) You can **customise the colour** palette if you wish.
6) Once you are happy with the colour palette, select the green **Preview & Download** button, select **png**, and then click on the green **Download** button to download the figure as a png.


