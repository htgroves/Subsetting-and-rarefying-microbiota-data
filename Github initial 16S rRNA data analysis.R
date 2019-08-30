setwd("M_final_files_for_r_HG0115")

library(Biostrings) 
library(ggplot2)
library(plyr)
library(stringr)
library(gridExtra)
library(reshape2)
library(knitr)
library(vegan)
library(DESeq2)
library(phyloseq)

meta_data <- read.csv("HG0115_qiime_mapping_file.csv", header = T, row.names = 1)
names(meta_data)
all_data <-import_biom(BIOMfilename = "otu_table_json.biom",refseqfilename = "rep_set_for_r.fasta",treefilename = "HG0115.tre")
sample_data(all_data)=sample_data(meta_data)
sample_sums(all_data)
all_data
all_data = subset_taxa(all_data, Rank1 != "Unclassified")


makeTaxlabel <- function(OTU, mydata){ # Function from Aaron Saunders which does the OTU relabelling. Makes a string label using the lowest information tax level. First arguemtn is the OTU number, second is mydata which indicates the phyloseq object with the tax table
  OTU <- as.character(OTU) # To make sure the OTU numbers are stored as a character rather than an integer
  taxstrings <- as.character(tax_table(mydata)[OTU])
  empty_strings <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  tax_name <- NA 
  tax_level <- length(taxstrings)
  
  while(is.na(tax_name)|
        (tax_name %in% empty_strings)){
    tax_name <- taxstrings[tax_level] # start at the lowest level
    tax_level <- tax_level -1
  }
  tax_name
}
tax_table(all_data) =gsub("s__uncultured_bacterium", as.character(NA), tax_table(all_data)) # These all adjust some of the useless names from the taxonomy so that the function ignores them
tax_table(all_data) =gsub("s__uncultured_organism", as.character(NA), tax_table(all_data))
tax_table(all_data) =gsub("g__uncultured", as.character(NA), tax_table(all_data))
tax_table(all_data) =gsub("g__uncultured_bacterium", as.character(NA), tax_table(all_data))
tax_table(all_data) =gsub("f__uncultured_bacterium", as.character(NA), tax_table(all_data))
mynames = NULL # This runs the function on each row of the tax able and out the mynames object, which is the list of new names
for (i in 1:length(taxa_names(all_data))){
  mynames <- rbind(mynames, c(makeTaxlabel(taxa_names(all_data)[i],all_data)))
}
mynames =gsub("s__", "", mynames) #editing my names to get rid of some guff
mynames =gsub("g__", "", mynames)
mynames =gsub("f__", "", mynames)
mynames =gsub("o__", "", mynames)
mynames =gsub("c__", "", mynames)
mynames =gsub("p__", "", mynames)
length(taxa_names(all_data))
OTUID = str_c(mynames[,1], "_", seq(1,15512,1))
tax_table(all_data) <- cbind(tax_table(all_data), mynames=OTUID)
head(tax_table(all_data)) 
colnames(tax_table(all_data)) = c("Kingdon", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTUID") #renames the tax table headings to something real
taxa_names(all_data) <- tax_table(all_data)[,8]

## These extract individual columns from your dataset
king <- tax_table(all_data)[,1]
phyl <- tax_table(all_data)[,2]
cla <- tax_table(all_data)[,3]
ord <- tax_table(all_data)[,4]
fam <- tax_table(all_data)[,5]
gen <- tax_table(all_data)[,6]
spe <- tax_table(all_data)[,7]
id <- tax_table(all_data)[,8]

## These delete the unnecessary prefix
king <- gsub("k__", "", king)
phyl <- gsub("p__", "", phyl)
cla <- gsub("c__", "", cla)
ord <- gsub("o__", "", ord)
fam <- gsub("f__", "", fam)
gen <- gsub("g__", "", gen)
spe <- gsub("s__", "", spe)

## This substitudes your dataset with hte 'fixed' dataset
tax_table(all_data) <- cbind(king, phyl, cla, ord, fam, gen, spe, id)


save.image("HG0115_basic")

######### Subsetting

TN = subset_samples(all_data, Dosed =="TN") # subsetting out the naives
RSV = subset_samples(all_data, Dosed =="RSV") #subsetting out the RSV infected
PBS = subset_samples(all_data, Dosed =="PBS") # subsetting out the PBS dosed


RSV_wt = subset_samples(RSV, Genetics =="wt") # Just RSV infected wt mice

PBS_wt = subset_samples(PBS, Genetics =="wt") # Just RSV infected wt mice


wt = subset_samples(all_data, Genetics =="wt")

RSV_a = subset_samples(RSV_wt, Group=="A")

RSV_e = subset_samples(RSV_wt, Group=="E")


save.image("HG0115_subset_by_gen")

#### rarefying data

min(sample_sums(TN))
sample_sums(TN) #Finds the minium number of reads in the TN experiment
set.seed(111) # set the random number generator to 111.  Always type this before doing rarefication or anything that invovles a random number generation
TN_rare <- rarefy_even_depth(TN, 29000)

min(sample_sums(RSV))
set.seed(111)
RSV_rare = rarefy_even_depth(RSV, 29000)

min(sample_sums(PBS))
sample_sums(PBS)
set.seed(111)
PBS_rare = rarefy_even_depth(PBS, 29000)

min(sample_sums(RSV_wt))
set.seed(111)
RSV_wt_rare = rarefy_even_depth(RSV_wt, 29000)

min(sample_sums(PBS_wt))
set.seed(111)
PBS_wt_rare = rarefy_even_depth(PBS_wt, 29000)


min(sample_sums(RSV_a))
set.seed(111)
RSV_a_rare = rarefy_even_depth(RSV_a, 29000)

min(sample_sums(RSV_e))
set.seed(111)
RSV_e_rare = rarefy_even_depth(RSV_e, 29000)


save.image("HG0115_rarefied")
