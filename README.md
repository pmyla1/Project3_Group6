# Project3_Group6
An introduction into manipulating **mixed ploidy VCF files** and obtaining various population genetics metrics from these files including **site frequency spectra**, **Fst**, etc.

## **Overview**
**Polyploidy** and **whole genome duplications** (WGD) occur throughout **all kingdoms** of life, including in animals such as *Xenopus laevis*, and are **especially ubiquitous** in **plants**. WGD is a **major mutational** process that disrupts **ionomic**, **cellular**, and **meiotic** processes, and **neo-polyploids** must overcome various challenges including **genomic instability** and **chromosomal mis-segregation** during meiosis (Bray et al., 2023). One of the immediate novel challenges related to WGD is the **formation of multivalent crossovers** between homologous chromosomes during **metaphase I** of meiosis, which can result in **entanglement** and **chromosomal breakage** at anaphase I (Bray et al., 2023). If neo-polyploids are **able to overcome** these initial challenges with meiosis and genome instability, they can become **established** as a **stable polyploid lineage**. 

Some of the genetic adaptations to WGD in polyploid lineages have been characterised, however, despite **convergence at the process level**, there is **low/no convergence** at the level of the **gene/orthologue** (Bray et al, 2023 [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.09.27.559727v2)). For example, some genes found to be **under selection** in **polyploid** lineages of *Cochlearia*, *Cardamine amara*, and *Arabidopsis arenosa* were found to be enriched in functional categories involving **DNA repair**, **cell division**, and **ion homeostasis**, to name a few (Bray et al., 2023).

There are **two main types** of polyploidy: **autopolyploids** (where all subgenomes are from the **same species** without hybridisation) and **allopolyploids** (where the subgenomes are inherited from **different species** through **hybridisation**). Autopolyploids and allopolyploids display distinct characteristic site frequency spectra (SFS). **Autopolyploid** SFS show a **Poisson distribution** with a **high proportion** of **low frequency** variants, whereas **allopolyploid** SFS have a characteristic **trimodal distribution** with **high proportions** of **low**, **intermediate**, and **high frequency** variants.

## Example of an autopolyploid SFS
![plot](/Users/lukearcher/Desktop/230324_lyrata_arenosa/lyrata_AF_spectrum.png)  


The VCF file included **diploid** and **tetraploid** *Arabidopsis arenosa* and *Arabidopsis lyrata* populations sampled from **throughout Europe** (insert more geographical context here) including from a well-established **hybrid zone** in **Wachau** in the **Austrian Forealps**. 

The aim of the project was to **select specific populations** from this **hybrid zone** and **control** populations from distinct locations and to look for **genetic admixture/hybridisation** between *A. arenosa* and *A. lyrata* **autotetraploid lineages**. One of the expectations was that the **site-frequency spectra** of the mixed ploidy VCF files containing admixed *A. arenosa* and *A. lyrata* **populations/hybrids** would show a **characteristic bimodal distribution** with peaks at low (~0.0-0.1) and high (~0.9-1.0) allele frequencies, and a smaller peak at intermediate (~0.4-0.6) allele frequencies. 

### Information on the *A. arenosa* and *A. lyrata* Populations in the VCF

| Species        | Ploidy           | 3-letter pop code(s) |
| ------------- |:-------------:| -----:|
| pure *A. arenosa* | 4x | KEH, BZD |
| pure *A. lyrata*  | 4x | KAG, LIC, MOD, MAU |
| *A. arenosa* x *A. lyrata* hybrids | 4x | FRE, HAB, OCH |  


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



