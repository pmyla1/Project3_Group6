#190424 Polyploid Fst output script written by Luke Archer (2024). 
##Utilises the output from 190424_poly_fst.sh to produce Fst Manhattan plots

##set working directory 

setwd("<path/to/.fst/files/produced/by/190324_poly_fst.sh>")

###############
##load appropriate libraries
library(tidyverse)
library(ggplot2)
library(dplyr) ##for manipulation of dataframes
library(gridExtra) ##for multipanel plotting 
###############

###############
##BZD vs KAG
BZD_KAG_Fst<-read_tsv("BZD_KAG_Fst_output.fst")
##rename columns for easier plotting
BZD_KAG_Fst<-BZD_KAG_Fst %>% rename(Chromosome=NW_003302555.1) %>% rename(Position=`716`) %>% rename(Fst=`0.000000`) 
##make a threshold for colouring
BZD_KAG_Fst$threshold<-BZD_KAG_Fst$Fst>=0.6
##calculate the mean Fst value
BZD_KAG_Fst$mean_Fst<-mean(BZD_KAG_Fst$Fst)
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
##rename columns for easier plotting
BZD_KEH_Fst<-BZD_KEH_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.000000`) 
##make a threshold for colouring
BZD_KEH_Fst$threshold<-BZD_KEH_Fst$Fst>=0.6
##calculate the mean Fst value
BZD_KEH_Fst$mean_Fst<-mean(BZD_KEH_Fst$Fst)
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
##rename columns for easier plotting
BZD_OCH_Fst<-BZD_OCH_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.260870`) 
##make a threshold for colouring
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
###############

##############
####BZD vs LOI
BZD_LOI_Fst<-read_tsv("BZD_LOI_Fst_output.fst")
head(BZD_LOI_Fst)
##rename columns for easier plotting
BZD_LOI_Fst<-BZD_LOI_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.260870`) 
##make a threshold for colouring
BZD_LOI_Fst$threshold<-BZD_LOI_Fst$Fst>=0.6
##calculate the mean Fst value
BZD_LOI_Fst$mean_Fst<-mean(BZD_LOI_Fst$Fst)
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
##rename columns for easier plotting
KEH_FRE_Fst<-KEH_FRE_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.123746`) 
##make a threshold for colouring
KEH_FRE_Fst$threshold<-KEH_FRE_Fst$Fst>=0.6
##calculate the mean Fst value
KEH_FRE_Fst$mean_Fst<-mean(KEH_FRE_Fst$Fst)
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
##rename columns for easier plotting
KEH_HAB_Fst<-KEH_HAB_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`716`) %>% 
  rename(Fst=`0.000000`) 
##make a threshold for colouring
KEH_HAB_Fst$threshold<-KEH_HAB_Fst$Fst>=0.6
##calculate the mean Fst value
KEH_HAB_Fst$mean_Fst<-mean(KEH_HAB_Fst$Fst)
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
##############

##############
####KEH vs KAG
KEH_KAG_Fst<-read_tsv("KEH_KAG_Fst_output.fst")
##rename columns for easier plotting
KEH_KAG_Fst<-KEH_KAG_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`716`) %>% 
  rename(Fst=`0.000000`) 
##make a threshold for colouring
KEH_KAG_Fst$threshold<-KEH_KAG_Fst$Fst>=0.6
##calculate the mean Fst value
KEH_KAG_Fst$mean_Fst<-mean(KEH_KAG_Fst$Fst)
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
##rename columns for easier plotting
KEH_LIC_Fst<-KEH_LIC_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.000000`) 
##make a threshold for colouring
KEH_LIC_Fst$threshold<-KEH_LIC_Fst$Fst>=0.6
##calculate the mean Fst value
KEH_LIC_Fst$mean_Fst<-mean(KEH_LIC_Fst$Fst)
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
##rename columns for easier plotting
KEH_LOI_Fst<-KEH_LOI_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`716`) %>% 
  rename(Fst=`0.090909`) 
##make a threshold for colouring
KEH_LOI_Fst$threshold<-KEH_LOI_Fst$Fst>=0.6
##calculate the mean Fst value
KEH_LOI_Fst$mean_Fst<-mean(KEH_LOI_Fst$Fst)
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
##rename columns for easier plotting
KEH_MOD_Fst<-KEH_MOD_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.000000`) 
##make a threshold for colouring
KEH_MOD_Fst$threshold<-KEH_MOD_Fst$Fst>=0.6
##calculate the mean Fst value
KEH_MOD_Fst$mean_Fst<-mean(KEH_MOD_Fst$Fst)
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
##rename columns for easier plotting
OCH_FRE_Fst<-OCH_FRE_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`-0.038902`) 
##make a threshold for colouring
OCH_FRE_Fst$threshold<-OCH_FRE_Fst$Fst>=0.6
##calculate the mean Fst value
OCH_FRE_Fst$mean_Fst<-mean(OCH_FRE_Fst$Fst)
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
##rename columns for easier plotting
OCH_HAB_Fst<-OCH_HAB_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`716`) %>% 
  rename(Fst=`0.000000`) 
##make a threshold for colouring
OCH_HAB_Fst$threshold<-OCH_HAB_Fst$Fst>=0.6
##calculate the mean Fst value
OCH_HAB_Fst$mean_Fst<-mean(OCH_HAB_Fst$Fst)
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
##rename columns for easier plotting
OCH_KEH_Fst<-OCH_KEH_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.167562`) 
##make a threshold for colouring
OCH_KEH_Fst$threshold<-OCH_KEH_Fst$Fst>=0.6
##calculate the mean Fst value
OCH_KEH_Fst$mean_Fst<-mean(OCH_KEH_Fst$Fst)
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
##rename columns for easier plotting
OCH_LIC_Fst<-OCH_LIC_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.260870`) 
##make a threshold for colouring
OCH_LIC_Fst$threshold<-OCH_LIC_Fst$Fst>=0.6
##calculate the mean Fst value
OCH_LIC_Fst$mean_Fst<-mean(OCH_LIC_Fst$Fst)
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
##rename columns for easier plotting
OCH_MOD_Fst<-OCH_MOD_Fst %>% rename(Chromosome=NW_003302555.1) %>% 
  rename(Position=`1317`) %>% 
  rename(Fst=`0.260870`) 
##make a threshold for colouring
OCH_MOD_Fst$threshold<-OCH_MOD_Fst$Fst>=0.6
##calculate the mean Fst value
OCH_MOD_Fst$mean_Fst<-mean(OCH_MOD_Fst$Fst)
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
