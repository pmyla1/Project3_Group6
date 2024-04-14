# Project3_Group6
This Github page offers an introduction into manipulating **mixed ploidy VCF files** and how to calculate various population genetics metrics from polyploid VCFs such as **site frequency spectra**, and **Fst**.

## **Introduction**
**Polyploidy** and **whole genome duplications** (WGD) occur throughout **all kingdoms** of life and are especially **ubiquitous in plants**. WGD is a **major mutation** that disrupts **ionomic**, **cellular**, and **meiotic** processes, and **neo-polyploids** must overcome various challenges including **genomic instability** and **chromosomal mis-segregation** during meiosis (Margburger et al., 2019 [Nature](https://www.nature.com/articles/s41467-019-13159-5)). One of the immediate challenges in neo-polyploids is the **formation of multivalent crossovers** between homologous chromosomes during meiosis, which can lead to **entanglement** and **chromosomal breakage** at anaphase I (Bray et al., 2023). If neo-polyploids **overcome** the initial challenges related to meiosis and genome instability, they can become **established** as a **polyploid lineage**. 

Some of the genetic adaptations to WGD in polyploids have been characterised, however, despite **process-level convergence**, there appears to be **low convergence** at the **gene/orthologue level** (Bray et al, 2023 [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.09.27.559727v2)). For example, many of the genes **under selection in polyploid** lineages of *Cochlearia*, *Cardamine amara*, and *Arabidopsis arenosa* belong to functional categories involving **DNA repair**, **cell division**, and **ion homeostasis** (Bray et al., 2023).

There are **two main types** of polyploids: **autopolyploids** (all subgenomes originate from the **same species** without hybridisation) and **allopolyploids** (subgenomes are inherited from **different species** through **hybridisation**). Autopolyploids and allopolyploids display distinct characteristic site frequency spectra (SFS). **Autopolyploid** SFS have a **Poisson distribution** with a **high proportion** of **low frequency** variants, whereas **allopolyploid** SFS have a characteristic **trimodal distribution** with **high proportions** of **low**, **intermediate**, and **high frequency** variants.


***Arabdopsis arenosa*** is a biennial **outcrossing plant** which produces distinct **diploid** and **tetraploid** lineages throughout **Central Europe** and is closely related to the **model species** ***Arabidopsis thaliana*** (Margburger et al., 2019). Similarly, the sister species ***Arabidopsis lyrata*** also forms distinct **diploid** and **tetraploid** lineages across its **distribution** range, and there have been reports of **gene flow** between the **tetraploid** lineages of ***A. arenosa*** and ***A. lyrata*** where these species have **overlapping ranges** (Jørgensen et al., 2011 [BMC Ecology & Evolution](https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-11-346])). 

The **rationale** for this project was the discovery of **bidirectional gene flow** between **autotetraploid** lineages of *A. arenosa* and *A. lyrata* which may have **facilitated** the **stabilisation of polyploidy** post WGD (Marburger et al., 2019)[Nature Communications](https://www.nature.com/articles/s41467-019-13159-5). In this project, we selected specific **tetraploid populations** of *A. arenosa* and *A. lyrata* from the **Austrian Forealps** and the **Wachau hybrid zone** suspected to be 50:50 hybrids and to assess the degree of **genetic admixture** using **fastSTRUCTURE**. We aimed to discover whether **hybridisation** between *A. arenosa* and *A. lyrata* has created any **allopolyploid lineages** in the **Wachau hybrid zone** consisting of approximately **50:50 admixture**. If there were 50:50 hybrids in the VCF, we expected the **site-frequency spectra** of the admixed *A. arenosa* and *A. lyrata* **populations** to show a **characteristic allopolyploid SFS** with peaks at low, intermediate, and high allele frequencies. 


The original VCF file **Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz** contains samples from **diploid** and **tetraploid** *A. arenosa* and *A. lyrata* populations sampled **throughout Europe** including lineages from the **Austrian Forealps** and others from the **Wachau hybrid zone**. The population map can be accessed via the following link [Population Map](https://www.google.com/maps/d/viewer?mid=1HAhM5y-bYMJbXCtMSZaubk1qe0wX6JI&ll=48.09350708234271%2C15.809612499999968&z=9).



### **Software Used**
#### **1) fastSTRUCTURE**.
An algorithm to infer population structure: sourced from [fastSTRUCTURE](https://rajanil.github.io/fastStructure/). 

Citation: Anil Raj, Matthew Stephens, and Jonathan K. Pritchard. fastSTRUCTURE: Variational Inference of Population Structure in Large SNP Data Sets , (Genetics) June 2014 197:573-589.

#### Installation and dependencies

#### fastSTRUCTURE depends on the following:
[Numpy](https://numpy.org/)   

[SciPy](https://scipy.org/) 

[Cython](https://cython.org/) 

[GNU Scientific Library](https://www.gnu.org/software/gsl/)  

#### Obtaining the Source code from GitHub

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

### Building Python extensions

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

#### **2) The Genome Analysis Toolkit (GATK) `v4.2.2.0`**. 
For support and documentation click on the following link [GATK](https://software.broadinstitute.org/gatk/) 

### Dependencies for GATK v4
To configure, build, and run GATK:

You require **`JAVA 17`** to run or build GATK.
**`JAVA 17`** can be installed using the following instructions 
1) on **macOS**, you can install the [Homebrew package manager](https://docs.brew.sh/Installation) and then run the following command
   ```
   brew tap homebrew/cask-versions
   ## then to install the Eclipse Foundation's OpenJDK 17 run 
   brew install --cask temurin17
   ```
Furthermore, you must have a **`Python`** version greater than or equal to **`Python 2.6`**, and an **`R version of at least 3.2.5`**

To **build GATK** you must have [git-lfs](https://git-lfs.com/) of **1.1.0 or greater**. This is required to download the large files that are used to build GATK, including the test files. 
After downloading git-lfs, run the following command 
```
git lfs install
```
You also **require Gradle 5.6** - you can use the `./gradlew` script which will download and use an **appropriate Gradle version automatically**.

### Python dependencies
GATK uses **Conda** as a **package manager** to establish and manage the Python environment and dependencies required by GATK tools that require Python.

You must first **create** and then **activate** a `Conda` environment for `GATK` using either **`Miniconda`** or **`Conda`**.
To create the GATK Conda environment, there are **two main methods** that can be used:
If running from a **zip or tar** distribution, run the following commands to create the `gatk` environment:
```
conda env create -f gatkcondanev.yml
```
Or if you are running from a **cloned repository**, run the following command:
```
./gradlew localDevCondaEnv
```
You can **activate** the conda environment by executing the following **shell command**.  
```
source activate gatk
```

### Installation of `GATK v4.2.2.0`

Firstly, if you have any **previous GATK versions** installed on your device, please **uninstall** using the following commands:
```
source deactivate
## assuming you called your GATK environment gatk
conda env remove -n gatk
```

Next, **update your conda** package manager to the **latest version** for your base environment
```
conda update -n base conda
```

To install the **GATK v4** package

Navigate to the **directory** where you have **stored your GATK jars** and the **gatk wrapper** script, ensure the **`gatkcondaenv.yml`** is present, then run the following command:
```
conda env create -n gatk -f gatkcondaenv.yml

source activate gatk
```
To check if your GATK environment is **properly installed** run: 
```
conda list
```
Which should return a list of packages you have installed, **`gatkpythonpackages`** should be one of these.

## **3) R and RStudio - R version 4.3.1 (2023-06-16)**

### Installation 

### To install R-4.3.1 on Windows PC:

1) **Uninstall previous versions** of R. 

2) Click on the following link [R installation for Windows](https://cran.r-project.org/bin/windows/base/) 

3) Click **Download R-4.3.1 for Windows**. 

4) Open the installer and follow the instructions using **default options**. 


### To install Rtools-4.3 on Windows PC:

1) **Uninstall previous versions** of Rtools. 

2) Click on the following link [Rtools installation for Windows](https://cran.r-project.org/bin/windows/Rtools/)  

3) Click **RTools 4.3**. 

4) Click **Rtools43 installer**. 

5) Open the installer and follow the instructions using **default options**.


### To install RStudio-4.3.1 on Windows PC:

1) **Uninstall previous versions** of RStudio. 

2) Click on the following link [RStudio installation for Windows](https://posit.co/download/rstudio-desktop/) 

3) Follow the downloads instructions from **Step 2 onwards**. 

4) Open the installer and then follow the instructions using **default options**. 


### To install R version 4.3.1. on Apple macOS

1) **Uninstall previous installations** of R by navigating to your **applications folder** and moving **R** and **XQuartz** to the **Bin**. 

2) To install R, **check which Apple macOS version** you have - briefly, click the **Apple Logo** in the **top left corner** of your screen and then **choose** ***About this mac*** > note what the **Processor line** says. 

3) Click on the following link [R installation for Apple macOS](https://cran.r-project.org/bin/macosx/)  

  
   If your Processor line has **"Intel"** then follow the download **instructions** for **`R-4.3.1-x86_64.pkg`**. 
   
   **Otherwise**, follow the download instructions for **`R-4.3.1-arm64.pkg`**. 
  

4) Open the installer and follow the instructions using **default options**. 

5) To install XQuartz - **download XQuartz** from the following link [XQuartz download for macOS](https://www.xquartz.org/), open the installer using **default options**. 


### To install RStudio-4.3.1 on Apple macOS

1) **Uninstall previous versions** of RStudio. 
2) Click on the following link [RStudio Version 4.3.1 Installation for macOS](https://posit.co/download/rstudio-desktop/) and follow the download instructions from **step 2 onwards**. 
3) Open the installer and follow the instructions using **default options**.

## To install SplitsTree

1) To install SplitsTree, go to the University of Tübingen Website [SplitsTree](https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/) and click on the `Download' link.

2) Select the correct program installer for you specific device and follow the instructions for installation/download using default options. 

## Scripts

1) **GATK_select_variants_initial.sh** was executed on the HPC using a **shared GATK environment**, and a custom **samtools environment**.

The script uses the `gatk IndexFeatureFile -I`  command to create an **indexed fasta** file from the original **Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz**.

Subsequently, the script creates a **sequence dictionary** using `gatk CreateSequenceDictionary`, and then uses `samtools faidx` to index the fasta file created in a previous step of the script. 

Finally, the script uses `gatk SelectVariants` with the `-sn` flag to **select specific individuals** from **Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz** and to make a new .vcf.gz with only these individuals called **new_pops_filtered.vcf.gz**.

2) **220324_whole_pipeline_gatk.sh** was executed on the HPC using a shared conda environment `/shared/conda/shared/ and a shared gatk environment `/shared/apps/conda/bio2/`.

This script uses `gatk SelectVariants` with the `-sn` flag to select only the tetraploid populations of interest (either pure *A. lyrata*, pure *A. arenosa*, or suspected to be hybrids) and creates a new filtered vcf called 220324_filtered_pops.vcf.gz.

The script then calculates the **site-frequency spectra** on 220324_filtered_pops.vcf.gz, and **allele frequencies** using Tuomas Hämälä's (2023) **`poly_sfs.c`** and **`poly_freq.c`** scripts, respectively.

Next, the VCF file is prepared for **`Cochlearia_create_structure_file.py`** using the **'-s true'** flag for pseudo-diploidization, and the populations are rearranged into alphabetical order for plotting purposes. 


### Some of the populations in Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf.gz 


### *A. arenosa* and *A. lyrata* designated populations and ploidy levels

 | Species        | Ploidy           | 3-letter pop code(s) |
 | ------------- | :-------------: | -----:|
 | pure *A. arenosa* | 4x | KEH, BZD |
 | pure *A. lyrata* | 4x | KAG, LIC, MOD, MAU |
 | pure *A. lyrata* | 2x | PIZ, PLE |
 | hybrid *A. arenosa* x *A. lyrata* | 4x | FRE, HAB, OCH |  


## Exploratory genetic analyses with PCA 

We performed **exploratory population genetic analyses** using two different PCA techniques, (1) **Adegenet**, and (2) using Tuomas Hämälä's (2023) [est_adapt_pca.R](https://github.com/thamala/polySV/blob/main/est_adapt_dist.r) adapted PCA script and using our filtered vcf.gz as the input file for the PCA. 

### Without BZD population

<img width="468" alt="Alt_PCA_tets_only_no_BZD" src="https://github.com/pmyla1/Project3_Group6/assets/151543531/c0c26237-cabe-4620-8b3e-a095750fddab">


*This figure shows the tetraploid only populations retained in the last VCF file used to investigate population structure and admixture between A. arenosa and A. lyrata. KEH samples ***KEH-06 and KEH-08*** form a cluster with ***OCH-05 and FRE-06*** along ***PC1***, which explains ***33% of the variance*** in the dataset. Conversely, ***KEH-07 and KEH-09*** form a separate cluster with ***FRE-08*** scoring positively along PC1 with the same magnitude as the previously mentioned cluster, but being ***differentiated along PC2**. Lastly, ***KEH-05 and KEH-10*** form a separate cluster with ***FRE-05***, with extremely ***negative scores along PC1***. The other samples appear to form clusters with individuals from the same population, e.g. MOD with MOD, etc.*

### With BZD plus additional populations (MAU, LIC, KAG)

<img width="475" alt="Alt_PCA_tets_only_with_BZD" src="https://github.com/pmyla1/Project3_Group6/assets/151543531/c2ae0068-0773-4919-a283-7bc49135b286">

## Exploratory genetic analyses with discriminant analysis of principal components (DAPC)

We subsequently performed a discrimininant analysis of principal components (DAPC) on individuals in the filtered VCF in order to determine the number of population clusters and to discern if there were any sample mix ups. Our cluster analysis suggested that there were only 3 population clusters (K = 3). (**ADD MORE INFORMATION HERE FROM YOUR R SCRIPTS AND OUTPUT**).


## Phylogenetic analyses with SplitsTree

Next, we looked at the **relationships** between the individuals and populations with SplitsTree, following the download instructions from the **University of Tübingen** website in the link [SplitsTree](https://software-ab.cs.uni-tuebingen.de/download/splitstree4/welcome.html). 

Briefly, after converting the latest VCF file (**290324_tetraploids_only.vcf.gz**) into a **genlight** object using the **vcf2genlight** function in the **290324_populations.R script**, the genlight object can be converted into **Nei's genetic distance** data and subsequently converted into a **phylogentic distance file** (.phy.dst) before being loaded into SplitsTree. Both the **individual** and the **population** data were used to produce **phylogenetic networks** in SplitsTree, and can be visualized below.

### Phylogenetic network - Relationship between Individuals
![290324_individuals](https://github.com/pmyla1/Project3_Group6/assets/151543531/b07c9e8a-4758-469f-b15e-36e3f2c15c07)


***Figure 3a*** *Phylogenetic network showing the relationship between individuals. Most individuals appear to group together with their respective populations, (e.g. MOD appears to form a single cluster) however, there are some sample mixups with KEH individuals forming three separate clusters. Furthermore, some OCH individuals cluster with individuals from different populations.*


<img width="475" alt="Phylo_tree_inds_with_BZD" src="https://github.com/pmyla1/Project3_Group6/assets/151543531/76d26c5d-2d06-49e5-b9da-1e492e65974f">

***Figure 3b*** *Phylogenetic network showing relationships between individuals, including BZD*. 

### Phylogenetic network - Relationship between Populations
![290324_populations](https://github.com/pmyla1/Project3_Group6/assets/151543531/b7d838f1-e7e6-4165-a55e-f3142c03a75d)

***Figure 3c*** *Phylogenetic network showing the relationship between populations.*


## Population structure analysis

### Structure plot of the tetraploid lineages with K = 2 

![K2_structure_plot](https://github.com/pmyla1/Project3_Group6/assets/151543531/cc49ac45-9aaa-494e-a258-691b162e312e)

This image was produced using the GUI [Omicsspeaks Structure Plot V2.0](http://omicsspeaks.com/strplot2/). 

The input file is a comma separated value (CSV) file containing the individual name from the VCF (e.g BZD-01tl), followed by the population in the VCF (e.g. BZD), followed by the fastSTRUCTURE output for the genetic admixture proportions when K=2 (e.g 0.95,0.05).

The link to the input file can be accessed here [OmicsSpeaks input]
[K2_omics_speaks_input.csv](https://github.com/pmyla1/Project3_Group6/files/14791197/K2_omics_speaks_input.csv)

**Orange** bars represent ***A. arenosa*** whereas **green** bars represent ***A. lyrata***.

Contrary to our expectations, when K = 2, **FRE** was estimated to be **pure *A. lyrata*** as opposed to a 50/50 hybrid. 

### Allele frequency differences on a larger cohort of *A. arenosa* and *A. lyrata* samples

In order to determine whether there was **hybridisation** betweeen *A. arenosa* and *A. lyrata* and the subsequent **formation** of an **allotetraploid** lineage (2 subgenomes: one from *A. arenosa*, the other from *A. lyrata*), we were given **text files** containing **4-fold degenerate** single nucleotide polymorphism (SNP) data from a **larger number of samples** from both species. The **structure** of the input files can be seen **below**.

 | **CHROM** | **POS** | **REF** | **ALT** | **AF** | **AC** | **AN** |
 | :-------: | :-----: | :-----: | :-----: | :----: | :----: | :----: |
 | scaffold_1 | 32 | C | A | 0.00053 | 1 | 1886 | 
 | scaffold_1 | 38 | A | T | 0.558 | 977 | 1750 |
 | scaffold_1 | 160 | C | A | 0.00974 | 18 | 1848 |

***Key: CHROM, chromosome; POS, position; REF, reference allele; ALT, alternative allele; AF, allele frequency; AC, allele count; AN, allele  number***

The structure of the **final output file** used to calculate **allele frequency differences** between species can be visualised **below**.

 | **CHROM** | **POS** | **AF_arenosa** | **AF_lyrata** | 
 | :-------: | :-----: | :------------: | :-----------: |
 | scaffold_1 | 32 | 0.154 | 0.00053 |

***Key: CHROM, chromosome; POS, position; AF_arenosa, allele frequency in A. arenosa; AF_lyrata, allele frequency in A. lyrata*** 

We obtained only the **common/shared SNPs** between **both species** by running the **250324_combined_lyrata_arenosa.py** script, and subsequently **compared the allele frequencies** between *A. arenosa* and *A. lyrata* at these common sites by calculating the **allele frequency difference**. Plots of the site frequency spectra per species can be seen below. 

<img width="475" alt="AF_spectrum_arenosa_lyrata" src="https://github.com/pmyla1/Project3_Group6/assets/151543531/b85b800f-6df0-47b7-ac73-ab0e6d46c6a0">


Next, we plotted the **allele frequency differences** per chromosome scaffold as a **Manhattan** plot using the **250324_common_SNPs.R script**, using an arbritary **threshold of 0.85** for SNPs expected to be **"fixed"** in one species relative to the others. The **orange dots** above the dashed red line represent **"fixed" allele frequency differences**, suggesting that these SNPs are **private to one species**. These plots can be seen **below**.

<img width="475" alt="Chrom_1_4_AF_differences" src="https://github.com/pmyla1/Project3_Group6/assets/151543531/d47d2218-269c-4948-b2fd-4a03af4d8f4e">


<img width="475" alt="Chrom_5_8_AF_differences" src="https://github.com/pmyla1/Project3_Group6/assets/151543531/b39061d5-245c-460d-b4ff-fbfc7ca751dd">


