#190424 Polyploid Fst output script written by Luke Archer (2024). 
##Utilises the output from 190424_poly_fst.sh to produce Fst Manhattan plots

##set working directory 

setwd("/Users/lukearcher/Desktop/LEVI_project/190424_poly_fst_output/")

###############
##load appropriate libraries
options(warn=1)

library(adegenet)
library(adegraphics) #not strictly necessary for all of this (hombrew r installs will interfere)
library(vcfR)
library(pegas)
library(StAMPP)
library(ade4)
library(MASS)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)
###############

###############
##BZD vs KAG
BZD_KAG_Fst<-read_tsv("BZD_KAG_Fst_output.fst")
head(BZD_KAG_Fst)
##rename columns for easier plotting purposes
BZD_KAG_Fst<-BZD_KAG_Fst %>% rename(Chromosome=NW_003302555.1) %>% rename(Position=`716`) %>% rename(Fst=`0.000000`) 
head(BZD_KAG_Fst)
##make a threshold for colouring purposes
BZD_KAG_Fst$threshold<-BZD_KAG_Fst$Fst>=0.6
head(BZD_KAG_Fst)
##calculate the mean Fst value
BZD_KAG_Fst$mean_Fst<-mean(BZD_KAG_Fst$Fst)
head(BZD_KAG_Fst)
##plot
BZD_KAG_plot<-ggplot(BZD_KAG_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="BZD vs KAG Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
##############

##############
####BZD vs KEH
BZD_KEH_Fst<-read_tsv("BZD_KEH_Fst_output.fst")
head(BZD_KEH_Fst)
##rename columns for easier plotting purposes
BZD_KEH_Fst<-BZD_KEH_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.000000`) 
head(BZD_KEH_Fst)
##make a threshold for colouring purposes
BZD_KEH_Fst$threshold<-BZD_KEH_Fst$Fst>=0.6
#head(BZD_KEH_Fst)
##calculate the mean Fst value
BZD_KEH_Fst$mean_Fst<-mean(BZD_KEH_Fst$Fst)
head(BZD_KEH_Fst)
##plot
BZD_KEH_plot<-ggplot(BZD_KEH_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="BZD vs KEH Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
#############

#############
####BZD vs OCH
BZD_OCH_Fst<-read_tsv("BZD_OCH_Fst_output.fst")
head(BZD_OCH_Fst)
##rename columns for easier plotting purposes
BZD_OCH_Fst<-BZD_OCH_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.260870`) 
head(BZD_OCH_Fst)
##make a threshold for colouring purposes
BZD_OCH_Fst$threshold<-BZD_OCH_Fst$Fst>=0.6
#head(BZD_OCH_Fst)
##calculate the mean Fst value
BZD_OCH_Fst$mean_Fst<-mean(BZD_OCH_Fst$Fst)
head(BZD_OCH_Fst)
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
###############

##############
####BZD vs LOI
BZD_LOI_Fst<-read_tsv("BZD_LOI_Fst_output.fst")
head(BZD_LOI_Fst)
##rename columns for easier plotting purposes
BZD_LOI_Fst<-BZD_LOI_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.260870`) 
head(BZD_LOI_Fst)
##make a threshold for colouring purposes
BZD_LOI_Fst$threshold<-BZD_LOI_Fst$Fst>=0.6
#head(BZD_LOI_Fst)
##calculate the mean Fst value
BZD_LOI_Fst$mean_Fst<-mean(BZD_LOI_Fst$Fst)
head(BZD_LOI_Fst)
##plot
BZD_LOI_plot<-ggplot(BZD_LOI_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="BZD vs LOI Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
##############

###############
####KEH vs FRE
KEH_FRE_Fst<-read_tsv("KEH_FRE_Fst_output.fst")
head(KEH_FRE_Fst)
##rename columns for easier plotting purposes
KEH_FRE_Fst<-KEH_FRE_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.123746`) 
head(KEH_FRE_Fst)
##make a threshold for colouring purposes
KEH_FRE_Fst$threshold<-KEH_FRE_Fst$Fst>=0.6
#head(KEH_FRE_Fst)
##calculate the mean Fst value
KEH_FRE_Fst$mean_Fst<-mean(KEH_FRE_Fst$Fst)
head(KEH_FRE_Fst)
##plot
KEH_FRE_plot<-ggplot(KEH_FRE_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="KEH vs FRE Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
##############

##############
####KEH vs FRE
KEH_HAB_Fst<-read_tsv("KEH_HAB_Fst_output.fst")
head(KEH_HAB_Fst)
##rename columns for easier plotting purposes
KEH_HAB_Fst<-KEH_HAB_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`716`) %>% 
  rename(Fst=`0.000000`) 
head(KEH_HAB_Fst)
##make a threshold for colouring purposes
KEH_HAB_Fst$threshold<-KEH_HAB_Fst$Fst>=0.6
#head(KEH_FRE_Fst)
##calculate the mean Fst value
KEH_HAB_Fst$mean_Fst<-mean(KEH_HAB_Fst$Fst)
head(KEH_HAB_Fst)
##plot
KEH_HAB_plot<-ggplot(KEH_HAB_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="KEH vs HAB Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
###############


##############
####KEH vs KAG
KEH_KAG_Fst<-read_tsv("KEH_KAG_Fst_output.fst")
head(KEH_KAG_Fst)
##rename columns for easier plotting purposes
KEH_KAG_Fst<-KEH_KAG_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`716`) %>% 
  rename(Fst=`0.000000`) 
head(KEH_KAG_Fst)
##make a threshold for colouring purposes
KEH_KAG_Fst$threshold<-KEH_KAG_Fst$Fst>=0.6
#head(KEH_KAG_Fst)
##calculate the mean Fst value
KEH_KAG_Fst$mean_Fst<-mean(KEH_KAG_Fst$Fst)
head(KEH_KAG_Fst)
##plot
KEH_KAG_plot<-ggplot(KEH_KAG_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="KEH vs KAG Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
###############


###############
####KEH vs LIC
KEH_LIC_Fst<-read_tsv("KEH_LIC_Fst_output.fst")
head(KEH_LIC_Fst)
##rename columns for easier plotting purposes
KEH_LIC_Fst<-KEH_LIC_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.000000`) 
head(KEH_LIC_Fst)
##make a threshold for colouring purposes
KEH_LIC_Fst$threshold<-KEH_LIC_Fst$Fst>=0.6
#head(KEH_FRE_Fst)
##calculate the mean Fst value
KEH_LIC_Fst$mean_Fst<-mean(KEH_LIC_Fst$Fst)
head(KEH_LIC_Fst)
##plot
KEH_LIC_plot<-ggplot(KEH_LIC_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="KEH vs LIC Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
###############

###############
####KEH vs FRE
KEH_LOI_Fst<-read_tsv("KEH_LOI_Fst_output.fst")
head(KEH_LOI_Fst)
##rename columns for easier plotting purposes
KEH_LOI_Fst<-KEH_LOI_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`716`) %>% 
  rename(Fst=`0.090909`) 
head(KEH_LOI_Fst)
##make a threshold for colouring purposes
KEH_LOI_Fst$threshold<-KEH_LOI_Fst$Fst>=0.6
#head(KEH_FRE_Fst)
##calculate the mean Fst value
KEH_LOI_Fst$mean_Fst<-mean(KEH_LOI_Fst$Fst)
head(KEH_LOI_Fst)
##plot
KEH_LOI_plot<-ggplot(KEH_LOI_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="KEH vs LOI Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
##############

##############
####KEH vs MOD
KEH_MOD_Fst<-read_tsv("KEH_MOD_Fst_output.fst")
head(KEH_MOD_Fst)
##rename columns for easier plotting purposes
KEH_MOD_Fst<-KEH_MOD_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.000000`) 
head(KEH_MOD_Fst)
##make a threshold for colouring purposes
KEH_MOD_Fst$threshold<-KEH_MOD_Fst$Fst>=0.6
#head(KEH_MOD_Fst)
##calculate the mean Fst value
KEH_MOD_Fst$mean_Fst<-mean(KEH_MOD_Fst$Fst)
head(KEH_MOD_Fst)
##plot
KEH_MOD_plot<-ggplot(KEH_MOD_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="KEH vs MOD Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
##############

##############
####OCH vs FRE
OCH_FRE_Fst<-read_tsv("OCH_FRE_Fst_output.fst")
head(OCH_FRE_Fst)
##rename columns for easier plotting purposes
OCH_FRE_Fst<-OCH_FRE_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`-0.038902`) 
head(OCH_FRE_Fst)
##make a threshold for colouring purposes
OCH_FRE_Fst$threshold<-OCH_FRE_Fst$Fst>=0.6
#head(OCH_FRE_Fst)
##calculate the mean Fst value
OCH_FRE_Fst$mean_Fst<-mean(OCH_FRE_Fst$Fst)
head(OCH_FRE_Fst)
##plot
OCH_FRE_plot<-ggplot(OCH_FRE_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="OCH vs FRE Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
##############

##############
####OCH vs FRE
OCH_HAB_Fst<-read_tsv("OCH_HAB_Fst_output.fst")
head(OCH_HAB_Fst)
##rename columns for easier plotting purposes
OCH_HAB_Fst<-OCH_HAB_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`716`) %>% 
  rename(Fst=`0.000000`) 
head(OCH_HAB_Fst)
##make a threshold for colouring purposes
OCH_HAB_Fst$threshold<-OCH_HAB_Fst$Fst>=0.6
#head(OCH_HAB_Fst)
##calculate the mean Fst value
OCH_HAB_Fst$mean_Fst<-mean(OCH_HAB_Fst$Fst)
head(OCH_HAB_Fst)
##plot
OCH_HAB_plot<-ggplot(OCH_HAB_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="OCH vs HAB Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
###############

###############
####OCH vs KEH
OCH_KEH_Fst<-read_tsv("OCH_KEH_Fst_output.fst")
head(OCH_KEH_Fst)
##rename columns for easier plotting purposes
OCH_KEH_Fst<-OCH_KEH_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.167562`) 
head(OCH_KEH_Fst)
##make a threshold for colouring purposes
OCH_KEH_Fst$threshold<-OCH_KEH_Fst$Fst>=0.6
#head(OCH_KEH_Fst)
##calculate the mean Fst value
OCH_KEH_Fst$mean_Fst<-mean(OCH_KEH_Fst$Fst)
head(OCH_KEH_Fst)
##plot
OCH_KEH_plot<-ggplot(OCH_KEH_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="OCH vs KEH Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
##############

##############
####OCH vs FRE
OCH_LIC_Fst<-read_tsv("OCH_LIC_Fst_output.fst")
head(OCH_LIC_Fst)
##rename columns for easier plotting purposes
OCH_LIC_Fst<-OCH_LIC_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.260870`) 
head(OCH_LIC_Fst)
##make a threshold for colouring purposes
OCH_LIC_Fst$threshold<-OCH_LIC_Fst$Fst>=0.6
#head(OCH_LIC_Fst)
##calculate the mean Fst value
OCH_LIC_Fst$mean_Fst<-mean(OCH_LIC_Fst$Fst)
head(OCH_LIC_Fst)
##plot
OCH_LIC_plot<-ggplot(OCH_LIC_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="OCH vs LIC Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
##############

##############
####OCH vs MOD
OCH_MOD_Fst<-read_tsv("OCH_MOD_Fst_output.fst")
head(OCH_MOD_Fst)
##rename columns for easier plotting purposes
OCH_MOD_Fst<-OCH_MOD_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.260870`) 
head(OCH_MOD_Fst)
##make a threshold for colouring purposes
OCH_MOD_Fst$threshold<-OCH_MOD_Fst$Fst>=0.6
#head(OCH_MOD_Fst)
##calculate the mean Fst value
OCH_MOD_Fst$mean_Fst<-mean(OCH_MOD_Fst$Fst)
head(OCH_MOD_Fst)
##plot
OCH_MOD_plot<-ggplot(OCH_MOD_Fst,aes(x=Position,y=Fst,colour=threshold))+
  geom_point(alpha=0.5)+
  scale_colour_brewer(palette='Set2')+
  theme_bw()+
  labs(title="OCH vs MOD Fst scan",x="Position",y="Fst",
       colour=element_blank())+
  theme(title=element_text(face='bold',size=11),
        panel.grid.minor = element_blank(),
        legend.position='none')
#############

###############
##This code block produces multipanel figures of the Fst scans 
##firstly, the BZD population contrasts
grid.arrange(BZD_KAG_plot,BZD_KEH_plot,BZD_OCH_plot,BZD_LOI_plot,nrow=2)
##next, the KEH population contrasts
grid.arrange(KEH_FRE_plot,KEH_HAB_plot,KEH_KAG_plot,KEH_LIC_plot,KEH_LOI_plot,KEH_MOD_plot,nrow=3)
##finally, the OCH population contrasts
grid.arrange(OCH_FRE_plot,OCH_HAB_plot,OCH_KEH_plot,OCH_LIC_plot,OCH_MOD_plot,nrow=3)
################