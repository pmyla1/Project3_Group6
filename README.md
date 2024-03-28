# Project3_Group6
An introduction into manipulating **mixed ploidy VCF files** and obtaining various population genetics metrics from these files including **site frequency spectra**, **Fst**, etc.

## **Overview**
The VCF files in question include **diploid** and **tetraploid** *Arabidopsis arenosa* and *Arabidopsis lyrata* populations sampled from **throughout Europe** (insert more geographical context here) including from a well-established **hybrid zone** in **Wachau** in the **Austrian Forealps**. 

The aim of the project was to **select specific populations** from this **hybrid zone** and **control** populations from distinct locations and to look for **genetic admixture/hybridisation** between *A. arenosa* and *A. lyrata* **autotetraploid lineages**. One of the expectations was that the **site-frequency spectra** of the mixed ploidy VCF files containing admixed *A. arenosa* and *A. lyrata* **populations/hybrids** would show a **characteristic bimodal distribution** with peaks at low (~0.0-0.1) and high (~0.9-1.0) allele frequencies, and a smaller peak at intermediate (~0.4-0.6) allele frequencies. 

## **Software/Programs Used**
### **1) fastSTRUCTURE**.
An algorithm to infer population structure: sourced from [fastSTRUCTURE](https://rajanil.github.io/fastStructure/). 
Citation: Anil Raj, Matthew Stephens, and Jonathan K. Pritchard. fastSTRUCTURE: Variational Inference of Population Structure in Large SNP Data Sets , (Genetics) June 2014 197:573-589.
### **2) The Genome Analysis Toolkit (GATK) v4.2.2.0**. 
For support and documentation go to [GATK](https://software.broadinstitute.org/gatk/) 
### **3) R and RStudio - R version 4.3.1 (2023-06-16)** - Nickname Beagle Scouts

# Installation and Dependencies 
### 1) fastSTRUCTURE.
#### fastSTRUCTURE depends on the following:
Numpy [Numpy](https://numpy.org/)
SciPy [SciPy](https://scipy.org/)
Cython [Cython](https://cython.org/)
GNU Scientific Library [GNU Scientific Library](https://www.gnu.org/software/gsl/) 
### Obtaining the Source code from GitHub
In order to successfully install fastSTRUCTURE you can follow the guidelines given on the fastSTRUCTURE GitHub page in the following link [fastSTRUCTURE](https://rajanil.github.io/fastStructure/).



