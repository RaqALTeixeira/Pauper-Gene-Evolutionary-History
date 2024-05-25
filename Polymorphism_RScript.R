# We have diploid sequences in our fasta file, but for the analysis in pegas
# we need to have haploid sequences by converting IUPAC codes to ACGT. 
#This involves creating a new fasta file
#where each (diploid) individual is represented by 2 (haploid) sequences. The 
#haplotypes will be named LupXXXX_1 and LupXXXX_2.
#using the script diploid_to_randomphased.pl 

#!/bin/bash

./diploid_to_randomphased.pl -in allExons.fas -out allExons_d2rp.fas


##In R:

#If the package is not installed, this will install
#pegas and load the library
install.packages("pegas")
library(pegas)
library(openxlsx)
library(readxl)

#read the fasta file (with the IUPAC converted do acgt)
pgene <- read.dna(file="allExons_d2rp",format="fas")
head(pgene)

#Create a loop that reads the
#table with the species assignment of each
#individual, then extracts just those sequences to a new object

###
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
library(readxl)

#Read the excel with the sample info
LupSamples <- readxl::read_excel("Lupinus_reseq_SAMPLES.xlsx", sheet = "updateJan24")

#Extract all species corresponding to the North American Perinials clade
#and check for occurences of each species
NAP_SPECIES <- LupSamples[grep("NAP", LupSamples$clade), ]
table(NAP_SPECIES$Species) -> species_count

#Repeat for North American Annuals
NAA_SPECIES <- LupSamples[grep("NAA", LupSamples$clade), ]
table(NAA_SPECIES$Species) -> Aspecies_count

#check species with occurence higher than 3
speciesall <- names(species_count[species_count >= 3])

Aspeciesall <- names(Aspecies_count[Aspecies_count >= 3])

#Create a loop in which, for each species
#in speciesall (25 species), checks the Bioinf code
#and searches the pgene dna list for the
#corresponding sequences. For a species with 3
#samples (which most of them are), there should be
#6 sequences retrieved.


############################################
##                LOOP
############################################
#Here, the loop is created
#n is number of species, define start index
1 -> n

#define a list in which will be stored a list
#containing the objects with dna to each species
gene1_sp <- list()

while (n <= length(speciesall)) {
  
  #In the loop, substitute "Lupinus albicaulis" with an array variable
  #called "species" that is defined as the name in each row of speciesall
  
  sp <- NAP_SPECIES[grep(speciesall[n], NAP_SPECIES$Species), ]
  biocodes = sp$`Bioinf code`
  
  #Account for the fact that, for each sample, there
  #is now two copies d2r
  
  bio2codes <- c(rep(biocodes, each = 2))
  paste(bio2codes, rep(c("_1", "_2"), length(biocodes) * 1), sep = "") -> bio2codes
  
  #Extract the rows correspondent to samples
  #of the species of interest
  gene1_sp[[n]] <- pgene[bio2codes,]
  
  #Increment index
  
  n <- n + 1
  
}
  
#Make it so the name of the list is that of the
#species, for easier read

names(gene1_sp) <- speciesall

############################################
##                LOOP (NAA)
############################################

1 -> n
gene1_Asp <- list()

while (n <= length(Aspeciesall)) {
  
  sp <- NAA_SPECIES[grep(Aspeciesall[n], NAA_SPECIES$Species), ]
  biocodes = sp$`Bioinf code`
 
  bio2codes <- c(rep(biocodes, each = 2))
  paste(bio2codes, rep(c("_1", "_2"), length(biocodes) * 1), sep = "") -> bio2codes
  
  gene1_Asp[[n]] <- pgene[bio2codes,]
  
  n <- n + 1
  
}

names(gene1_Asp) <- Aspeciesall
############################################

gene_NAPsp <- gene1_sp
gene_NAAsp <- gene1_Asp

############################################
##                STATS
############################################

# calculate Watterson's theta for 1 species
thetaW_sp1 <- theta.s(gene1_sp1)

# calculate nucleotide diversity for 1 species
pi_sp1 <- nuc.div(gene1_sp1)

# calculate Tajima's D for 1 species
tajD_sp1<-tajima.test(gene1_sp1)

############################################
#Watterson's theta 
############################################

###
# NAP
###

# Initialize an empty vector to store theta values for each species
theta_NAPvalues <- numeric()

# Loop through each species in the list by names
for (species_name in names(gene_NAPsp)) {
  # Extract the DNAbin objects for the current species
  dna_seq <- gene_NAPsp[[species_name]]
  
  # Calculate theta for the current species using theta.s function
  # Here, we use the Watterson's estimator by default
  theta_NAPvalues[species_name] <- theta.s(dna_seq)
}

# Print or further analyze the theta values for each species
print(theta_NAPvalues)

###
# NAA
###

theta_NAAvalues <- numeric()

for (species_name in names(gene_NAAsp)) {

  dna_seq <- gene_NAAsp[[species_name]]
  
  theta_NAAvalues[species_name] <- theta.s(dna_seq)
}

print(theta_NAAvalues)

############################################
#Nucleotide Diversity
############################################

###
# NAP
###

ND_NAPvalues <- numeric()

for (species_name in names(gene_NAPsp)) {
  
  dna_seq <- gene_NAPsp[[species_name]]
  
  ND_NAPvalues[species_name] <- nuc.div(dna_seq)
}

print(ND_NAPvalues)

###
# NAA
###

ND_NAAvalues <- numeric()

for (species_name in names(gene_NAAsp)) {
  
  dna_seq <- gene_NAAsp[[species_name]]
  
  ND_NAAvalues[species_name] <- nuc.div(dna_seq)
}

print(ND_NAAvalues)

############################################
#Tajima's D [ERROR]
############################################

###
# NAP
###

TajD_NAPvalues <- matrix(NA, nrow = length(gene_NAPsp), ncol = 3,
                         dimnames = list(names(gene_NAPsp), c("D", "p-value", "p-value beta")))

for (species_name in names(gene_NAPsp)) {
  
  dna_seq <- gene_NAPsp[[species_name]]
  
  TajD_output <- tajima.test(dna_seq)
  
  TajD_NAPvalues[species_name, "D"] <- TajD_output$D
  TajD_NAPvalues[species_name, "p-value"] <- TajD_output$Pval.normal
  TajD_NAPvalues[species_name, "p-value beta"] <- TajD_output$Pval.beta
  
  }

print(TajD_NAPvalues)

###
# NAA
###

TajD_NAAvalues <- matrix(NA, nrow = length(gene_NAAsp), ncol = 3,
                         dimnames = list(names(gene_NAAsp), c("D", "p-value", "p-value beta")))

for (species_name in names(gene_NAAsp)) {
  
  dna_seq <- gene_NAAsp[[species_name]]
  
  TajD_output <- tajima.test(dna_seq)
  
  TajD_NAAvalues[species_name, "D"] <- TajD_output$D
  TajD_NAAvalues[species_name, "p-value"] <- TajD_output$Pval.normal
  TajD_NAAvalues[species_name, "p-value beta"] <- TajD_output$Pval.beta
  
}

print(TajD_NAAvalues)


