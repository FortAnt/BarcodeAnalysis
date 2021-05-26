library(rncl)
library(splits)
library(dplyr)
library(ggplot2)
setwd("/path/to/directory/")
#load tree from BEAST+treannotator output   
yule_tr<-read.nexus("./barcode_BEASTout.trees.nex")

#run species delimitation, with single threshold
yule_gmyc<-gmyc(yule_tr)

#summary of GMYC run
summary(yule_gmyc)

#create species groups
GMYC_species<-spec.list(yule_gmyc)
#write GMYC species table for export
write.table(GMYC_species,file="GMYC_species.txt",quote=F,row.names = F)
#Our fasta headers files are in this format: NCBIaccession_Genus_species
#We only want to keep the Genus and Species info, by removing the NCBI accession number, using the sub command (remove everything before the first "_"):
GMYC_species$species<-sub("^.*?_","",GMYC_species$sample_name)
#transform GMYC_spec into factor
GMYC_species$GMYC_spec<-as.factor(GMYC_species$GMYC_spec)
#Count number of sequences in each GMYC species, and each species names within:
tally<-GMYC_species %>% group_by(GMYC_spec,species) %>% tally()

#create new tibble containing each GMYC cluster, and count number of different sequences names within:
tally_perc<-GMYC_species %>% group_by(GMYC_spec) %>% tally()


#loop through GMYC clusters, extract the species names with the highest number of sequences, and divide by the total number of sequences within GMYC cluster:
for (i in levels(tally_perc$GMYC_spec)){
  #add max percentage agreement within GMYC cluster
  tally_perc[tally_perc$GMYC_spec==i,"perc_agreement"]<-max(tally[tally$GMYC_spec==i,"n"])/sum(tally[tally$GMYC_spec==i,"n"])*100
  #add list of all species names found within GMYC cluster
  tally_perc[tally_perc$GMYC_spec==i,"species_names"]<-paste(tally[tally$GMYC_spec==i,"species"],collapse=",")
}

#exclude all GMYC clusters with less than two sequences
tally_perc<-subset(tally_perc,n>1)
#remove unused levels
tally_perc$GMYC_spec<-droplevels(tally_perc$GMYC_spec)
#clean string from species names, by removing "_"
tally_perc$species_names<-gsub("_"," ",tally_perc$species_names)


#Plot species names agreement within GMYC cluster with list of species names within GMYC cluster
pdf("Species_names_agreement_BARCODE.pdf")
ggplot()+
  geom_bar(data=tally_perc,aes(perc_agreement,GMYC_spec),
           stat="identity",col="grey",fill="#e0e0eb")+
  annotate("text",y=1:nlevels(tally_perc$GMYC_spec),
           x=140,
           label=tally_perc$species_names,
           fontface = 'italic',
           hjust=1,
           size=9)+
  ylab("GMYC cluster")+
  xlab("Species names agreement within GMYC cluster [%]")+
  coord_cartesian(expand=c(0,0))
dev.off()

