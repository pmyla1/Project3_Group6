# Project3_Group6
An introduction into manipulating **mixed ploidy VCF files** and obtaining various population genetics metrics from these files including **site frequency spectra**, **Fst**, etc.

## **Overview**
**Polyploidy** and **whole genome duplications** (WGD) occur throughout **all kingdoms** of life, including in animals such as *Xenopus laevis*, and are **especially ubiquitous** in **plants**. WGD is a **major mutational** process that disrupts **ionomic**, **cellular**, and **meiotic** processes, and **neo-polyploids** must overcome various challenges including **genomic instability** and **chromosomal mis-segregation** during meiosis (Bray et al., 2023). One of the immediate novel challenges related to WGD is the **formation of multivalent crossovers** between homologous chromosomes during **metaphase I** of meiosis, which can result in **entanglement** and **chromosomal breakage** at anaphase I (Bray et al., 2023). If neo-polyploids are **able to overcome** these initial challenges with meiosis and genome instability, they can become **established** as a **stable polyploid lineage**. 

Some of the genetic adaptations to WGD in polyploid lineages have been characterised, however, despite **convergence at the process level**, there is **low/no convergence** at the level of the **gene/orthologue** (Bray et al, 2023 [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.09.27.559727v2)). For example, some genes found to be **under selection** in **polyploid** lineages of *Cochlearia*, *Cardamine amara*, and *Arabidopsis arenosa* were found to be enriched in functional categories involving **DNA repair**, **cell division**, and **ion homeostasis**, to name a few (Bray et al., 2023).

There are **two main types** of polyploidy: **autopolyploids** (where all subgenomes are from the **same species** without hybridisation) and **allopolyploids** (where the subgenomes are inherited from **different species** through **hybridisation**). Autopolyploids and allopolyploids display distinct characteristic site frequency spectra (SFS). **Autopolyploid** SFS show a **Poisson distribution** with a **high proportion** of **low frequency** variants, whereas **allopolyploid** SFS have a characteristic **trimodal distribution** with **high proportions** of **low**, **intermediate**, and **high frequency** variants.

### Example of an autopolyploid SFS

![lyrata_AF_spectrum](https://github.com/pmyla1/Project3_Group6/assets/151543531/54ede91e-def1-44f3-b91e-f4078e570b37)



The VCF file included **diploid** and **tetraploid** *Arabidopsis arenosa* and *Arabidopsis lyrata* populations sampled from **throughout Europe** including lineages from the **Austrian Forealps** and others from a well-established **hybrid zone** in **Wachau**. 

The **rationale** for this project was the discovery of **bidirectional gene flow** between **autotetraploid** lineages of *A. arenosa* and *A. lyrata* which may have **facilitated** the **stabilisation of polyploidy** post WGD (Marburger et al., 2019)[Nature Communications](https://www.nature.com/articles/s41467-019-13159-5). Our aim was to **select specific populations** from the **Austrian Forealps** and the **Wachau hybrid zone** and to look for **genetic admixture/hybridisation** between *A. arenosa* and *A. lyrata* **autotetraploid lineages**. We wanted to discover **whether hybridisation** had resulted in the formation of an **allopolyploid** lineage consisting of a **50/50 admixture** between *A. arenosa* and *A. lyrata*. If there were 50/50 hybrids in the VCF, we expected the **site-frequency spectra** of the admixed *A. arenosa* and *A. lyrata* **populations/hybrids** to show a **characteristic allopolyploid SFS** with peaks at low (~0.0-0.1), intermediate (~0.4-0.6), and high (~0.9-1.0) allele frequencies. 

### Example of an allopolyploid SFS


### Information on the *A. arenosa* and *A. lyrata* Populations in the VCF


| Species        | Ploidy           | 3-letter pop code(s) |
| ------------- |:-------------:| -----:|
| pure *A. arenosa* | 4x | KEH, BZD |
| pure *A. lyrata* | 4x | KAG, LIC, MOD, MAU |
| pure *A. lyrata* | 2x | PIZ, PLE |
| hybrid *A. arenosa* x *A. lyrata* | 4x | FRE, HAB, OCH |  



### **Software/Programs Used**
#### **1) fastSTRUCTURE**.
An algorithm to infer population structure: sourced from [fastSTRUCTURE](https://rajanil.github.io/fastStructure/). 

Citation: Anil Raj, Matthew Stephens, and Jonathan K. Pritchard. fastSTRUCTURE: Variational Inference of Population Structure in Large SNP Data Sets , (Genetics) June 2014 197:573-589.

#### **2) The Genome Analysis Toolkit (GATK) v4.2.2.0**. 
For support and documentation click on the following link [GATK](https://software.broadinstitute.org/gatk/) 

#### **3) R and RStudio - R version 4.3.1 (2023-06-16)**


## Installation and Dependencies 
#### 1) fastSTRUCTURE.
#### fastSTRUCTURE depends on the following:
[Numpy](https://numpy.org/)   

[SciPy](https://scipy.org/) 

[Cython](https://cython.org/) 

[GNU Scientific Library](https://www.gnu.org/software/gsl/)  

#### Obtaining the Source code from GitHub

In order to successfully install fastSTRUCTURE you can follow the guidelines given on the fastSTRUCTURE GitHub page in the following link [fastSTRUCTURE](https://rajanil.github.io/fastStructure/).


## Scripts

Some of the linked scripts were written on the HPC using various software packages on a shared conda environment, whereas others were written and executed on a local machine using other packages in a local miniconda environment. 


## Phylogenetic analyses with SplitsTree


## Population structure analysis

### Structure plot of the tetraploid lineages with K = 2 

![K2_structure_plot](https://github.com/pmyla1/Project3_Group6/assets/151543531/cc49ac45-9aaa-494e-a258-691b162e312e)

**Orange** bars represent ***A. arenosa*** whereas *green* bars represent ***A. lyrata***.
Contrary to our expectations, when K =2, **FRE** was estimated to be **pure *A. lyrata*** as opposed to a 50/50 hybrid. 





