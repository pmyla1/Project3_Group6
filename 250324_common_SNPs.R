#############################
##This R script was written by Luke Archer (2024) and can be used to compare the allele frequencies of common
##SNPs between Arabidopsis lyrata and Arabidopsis arenosa tetraploids.
##The input file is a .tsv file created using pandas pd.to_csv() with 
##columns containing allele frequencies in both A. lyrata and A. arenosa 
##at the same sites.

##set working directory and load appropriate libraries
library(tidyverse)
library(dplyr)
library(adegenet)
library(ade4)
library(MASS)
library(vcfR)
library(gridExtra)

setwd('<path/to/Common_SNPs.tsv/file>')

##define input files
arenosa_lyrata_AFs<-read_tsv('Common_SNPs.tsv')
##convert to a tibble for easier data extraction
arenosa_lyrata_AFs<-as_tibble(arenosa_lyrata_AFs)
##drop ...1 and. `Unnamed: 0` columns from the df using dplyr::select 
cleaned_arenosa_lyrata_AFs<-dplyr::select(arenosa_lyrata_AFs,-c(...1,`Unnamed: 0`))
######################

#################
##subset the data per chromosome/scaffold 
chrom1<-subset(cleaned_arenosa_lyrata_AFs,CHROM=='scaffold_1')
#compute mean AFs per site
arenosa_mean1<-mean(chrom1$AF_arenosa)
lyrata_mean1<-mean(chrom1$AF_lyrata)

chrom2<-subset(cleaned_arenosa_lyrata_AFs,CHROM=='scaffold_2')
arenosa_mean2<-mean(chrom2$AF_arenosa)
lyrata_mean2<-mean(chrom2$AF_lyrata)

chrom3<-subset(cleaned_arenosa_lyrata_AFs,CHROM=='scaffold_3')
arenosa_mean3<-mean(chrom3$AF_arenosa)
lyrata_mean3<-mean(chrom3$AF_lyrata)

chrom4<-subset(cleaned_arenosa_lyrata_AFs,CHROM=='scaffold_4')
arenosa_mean4<-mean(chrom4$AF_arenosa)
lyrata_mean4<-mean(chrom4$AF_lyrata)

chrom5<-subset(cleaned_arenosa_lyrata_AFs,CHROM=='scaffold_5')
arenosa_mean5<-mean(chrom5$AF_arenosa)
lyrata_mean5<-mean(chrom5$AF_lyrata)

chrom6<-subset(cleaned_arenosa_lyrata_AFs,CHROM=='scaffold_6')
arenosa_mean6<-mean(chrom6$AF_arenosa)
lyrata_mean6<-mean(chrom6$AF_lyrata)

chrom7<-subset(cleaned_arenosa_lyrata_AFs,CHROM=='scaffold_7')
arenosa_mean7<-mean(chrom7$AF_arenosa)
lyrata_mean7<-mean(chrom7$AF_lyrata)

chrom8<-subset(cleaned_arenosa_lyrata_AFs,CHROM=='scaffold_8')
arenosa_mean8<-mean(chrom8$AF_arenosa)
lyrata_mean8<-mean(chrom8$AF_lyrata)
###################

###################
##PLOTS OF THE GENOME WIDE INFORMATION
mean_arenosa<-mean(cleaned_arenosa_lyrata_AFs$AF_arenosa)
mean_lyrata<-mean(cleaned_arenosa_lyrata_AFs$AF_lyrata)

##GENOME WIDE PLOTS
arenosa<-ggplot(cleaned_arenosa_lyrata_AFs,aes(AF_arenosa))+
  geom_histogram(fill='orange',colour='black',bins=100)+
  geom_vline(xintercept=mean_arenosa,linetype=2,colour=2)+
  theme_bw()+
  labs(title="Genome wide Allele frequency\ndistribution for A. arenosa SNPs",
       x="A. arenosa AF",y='Count')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor.x=element_blank())

lyrata<-ggplot(cleaned_arenosa_lyrata_AFs,aes(AF_lyrata))+
  geom_histogram(fill='violet',colour='black',bins=100)+
  geom_vline(xintercept=mean_lyrata,linetype=2,colour=2)+
  theme_bw()+
  labs(title="Genome wide Allele frequency\ndistribution for A. lyrata SNPs",
       x="A. lyrata AF",y='Count')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor.x=element_blank())
##arrange into a multipanel figure
grid.arrange(arenosa,lyrata,nrow=2)
################

################
###filter for minor allele frequencies >0.025
arenosa_lyrata_MAF2.5<-subset(cleaned_arenosa_lyrata_AFs, AF_arenosa>0.025&AF_lyrata>0.025)
arenosa_mean<-mean(arenosa_lyrata_MAF2.5$AF_arenosa)
lyrata_mean<-mean(arenosa_lyrata_MAF2.5$AF_lyrata)

######plot the arenosa_MAF2.5 data as a frequency histogram
arenosa_MAF2.5<-ggplot(arenosa_lyrata_MAF2.5,aes(AF_arenosa))+
  geom_histogram(colour='black',fill='orange',bins=100)+
  geom_vline(xintercept=mean(arenosa_lyrata_MAF2.5$AF_arenosa),
             linetype=2,colour=2)+
  theme_bw()+
  labs(title="Allele frequency spectrum of\ncommon SNPs in A. arenosa",
       x='A. arenosa AF',y='Count')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank())

lyrata_MAF2.5<-ggplot(arenosa_lyrata_MAF2.5,aes(AF_lyrata))+
  geom_histogram(colour='black',fill='violet',bins=100)+
  geom_vline(xintercept=mean(arenosa_lyrata_MAF2.5$AF_lyrata),
             linetype=2,colour=2)+
  theme_bw()+
  labs(title="Allele frequency spectrum of\ncommon SNPs in A. lyrata",
       x='A. lyrata AF',y='Count')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank())

##make a multipanel plot
grid.arrange(arenosa_MAF2.5,lyrata_MAF2.5,nrow=2)
##############

##############
##plot position vs AF_difference
##CALCULATE AF DIFFERENCES BETWEEN A. ARENOSA AND A. LYRATA
##POSITIVE DIFFERENCES = ALLELE MORE COMMON IN ARENOSA AND VICE VERSA
arenosa_lyrata_MAF2.5$AF_difference<-arenosa_lyrata_MAF2.5$AF_arenosa-arenosa_lyrata_MAF2.5$AF_lyrata

####make thresholds for all scaffolds as 0.85
chrom1$threshold<-chrom1$AF_difference>0.85
chrom2$threshold<-chrom2$AF_difference>0.85
chrom3$threshold<-chrom3$AF_difference>0.85
chrom4$threshold<-chrom4$AF_difference>0.85
chrom5$threshold<-chrom5$AF_difference>0.85
chrom6$threshold<-chrom6$AF_difference>0.85
chrom7$threshold<-chrom7$AF_difference>0.85
chrom8$threshold<-chrom8$AF_difference>0.85
##############
##PLOTTING OF AF DIFFERENCES BETWEEN SPECIES AND SCAFFOLDS
chrom1_AF_diff<-ggplot(chrom1,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
   scale_y_continuous(limits=c(0.5,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 1 AF differences\nat 4-fold neutral sites",
       x='Chromosome 1 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom2_AF_diff<-ggplot(chrom2,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.5,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 2 AF differences\nat 4-fold neutral sites",
       x='Chromosome 2 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom3_AF_diff<-ggplot(chrom3,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.5,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 3 AF differences\nat 4-fold neutral sites",
       x='Chromosome 3 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom4_AF_diff<-ggplot(chrom4,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.5,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 4 AF differences\nat 4-fold neutral sites",
       x='Chromosome 4 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom5_AF_diff<-ggplot(chrom5,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.5,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 5 AF differences\nat 4-fold neutral sites",
       x='Chromosome 5 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom6_AF_diff<-ggplot(chrom6,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.5,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 6 AF differences\nat 4-fold neutral sites",
       x='Chromosome 6 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom7_AF_diff<-ggplot(chrom7,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.5,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 7 AF differences\nat 4-fold neutral sites",
       x='Chromosome 7 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom8_AF_diff<-ggplot(chrom8,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.4)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.5,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 8 AF differences\nat 4-fold neutral sites",
       x='Chromosome 8 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

##make 2 multi-panel plots from SCAFFOLD1-4 AND 5-8
multipanel1<-grid.arrange(chrom1_AF_diff,chrom2_AF_diff,chrom3_AF_diff,chrom4_AF_diff,nrow=2,ncol=2)

multipanel2<-grid.arrange(chrom5_AF_diff,chrom6_AF_diff,chrom7_AF_diff,chrom8_AF_diff,nrow=2,ncol=2)
#############


################
##CALCULATING THE SITES WITH THE TOP 1% OUTLIERS IN TERMS OF AF DIFFERENCES
##find the sites with the highest maximum AF differences per scaffold
top_AF_diff_chrom1<-chrom1%>%arrange(desc(chrom1$AF_difference))
##calculate the top 1% of rows 
top_1PCT_rows_chrom1<-round(0.01*nrow(top_AF_diff_chrom1)) ##289
##use this to calculate the top outliers
top_1PCT_AF_outliers_chrom1<-top_AF_diff_chrom1%>%top_n(289,AF_difference)


top_AF_diff_chrom2<-chrom2%>%arrange(desc(chrom2$AF_difference))
##calculate the top 1% of rows 
top_1PCT_rows_chrom2<-round(0.01*nrow(top_AF_diff_chrom2)) #146
##use this to calculate the top outliers
top_1PCT_AF_outliers_chrom2<-top_AF_diff_chrom2%>%top_n(146,AF_difference)

top_AF_diff_chrom3<-chrom3%>%arrange(desc(chrom3$AF_difference))
##calculate the top 1% of rows 
top_1PCT_rows_chrom3<-round(0.01*nrow(top_AF_diff_chrom3)) ##223
##use this to calculate the top outliers
top_1PCT_AF_outliers_chrom3<-top_AF_diff_chrom3%>%top_n(223,AF_difference)

top_AF_diff_chrom4<-chrom4%>%arrange(desc(chrom4$AF_difference))
##calculate the top 1% of rows 
top_1PCT_rows_chrom4<-round(0.01*nrow(top_AF_diff_chrom4)) #183
##use this to calculate the top outliers
top_1PCT_AF_outliers_chrom4<-top_AF_diff_chrom4%>%top_n(183,AF_difference)

top_AF_diff_chrom5<-chrom5%>%arrange(desc(chrom5$AF_difference))
##calculate the top 1% of rows 
top_1PCT_rows_chrom5<-round(0.01*nrow(top_AF_diff_chrom5)) ##170
##use this to calculate the top outliers
top_1PCT_AF_outliers_chrom5<-top_AF_diff_chrom5%>%top_n(170,AF_difference)

top_AF_diff_chrom6<-chrom6%>%arrange(desc(chrom6$AF_difference))
##calculate the top 1% of rows 
top_1PCT_rows_chrom6<-round(0.01*nrow(top_AF_diff_chrom6)) ##228
##use this to calculate the top outliers
top_1PCT_AF_outliers_chrom6<-top_AF_diff_chrom6%>%top_n(228,AF_difference)

top_AF_diff_chrom7<-chrom7%>%arrange(desc(chrom7$AF_difference))
##calculate the top 1% of rows 
top_1PCT_rows_chrom7<-round(0.01*nrow(top_AF_diff_chrom7)) ##225
##use this to calculate the top outliers
top_1PCT_AF_outliers_chrom7<-top_AF_diff_chrom7%>%top_n(225,AF_difference)

##calculate the sites with the greatest AF differences
top_AF_diff_chrom8<-chrom8%>%arrange(desc(chrom8$AF_difference))
##calculate the top 1% of rows 
top_1PCT_rows_chrom8<-round(0.01*nrow(top_AF_diff_chrom8))
##use this to calculate the top outliers
top_1PCT_AF_outliers_chrom8<-top_AF_diff_chrom8%>%top_n(173,AF_difference)
##################

###################
##PLOTTING THE TOP 1% ALLELE FREQUENCY DIFFERENCES
##now plot the top 1% AF differences between species
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

chrom2_1PCT_diff<-ggplot(top_1PCT_AF_outliers_chrom2,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.8)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.7,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 2 top 1% AF\ndifferences at 4-fold sites",
       x='Chromosome 2 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom3_1PCT_diff<-ggplot(top_1PCT_AF_outliers_chrom3,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.8)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.7,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 3 top 1% AF\ndifferences at 4-fold sites",
       x='Chromosome 3 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom4_1PCT_diff<-ggplot(top_1PCT_AF_outliers_chrom4,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.8)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.7,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 4 top 1% AF\ndifferences at 4-fold sites",
       x='Chromosome 4 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom5_1PCT_diff<-ggplot(top_1PCT_AF_outliers_chrom5,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.8)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.7,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 5 top 1% AF\ndifferences at 4-fold sites",
       x='Chromosome 5 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom6_1PCT_diff<-ggplot(top_1PCT_AF_outliers_chrom6,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.8)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.7,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 6 top 1% AF\ndifferences at 4-fold sites",
       x='Chromosome 6 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom7_1PCT_diff<-ggplot(top_1PCT_AF_outliers_chrom7,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.8)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.7,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 7 top 1% AF\ndifferences at 4-fold sites",
       x='Chromosome 7 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')

chrom8_1PCT_diff<-ggplot(top_1PCT_AF_outliers_chrom8,aes(x=POS,y=AF_difference,colour=threshold))+
  geom_point(alpha=0.8)+
  geom_hline(yintercept=c(0,0.85),linetype=2,colour=2)+
  scale_y_continuous(limits=c(0.7,1))+
  scale_colour_brewer(palette="Set2")+
  theme_bw()+
  labs(title="Chromosome 8 top 1% AF\ndifferences at 4-fold sites",
       x='Chromosome 8 position',y='AF difference')+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor=element_blank(),
        legend.position='none')
chrom8_1PCT_diff
###multi-panel graphs
multipanel1PCT<-grid.arrange(chrom1_1PCT_diff,chrom2_1PCT_diff,chrom3_1PCT_diff,chrom4_1PCT_diff,nrow=2,ncol=2)

multipanel2PCT<-grid.arrange(chrom5_1PCT_diff,chrom6_1PCT_diff,chrom7_1PCT_diff,chrom8_1PCT_diff,nrow=2,ncol=2)
###################
