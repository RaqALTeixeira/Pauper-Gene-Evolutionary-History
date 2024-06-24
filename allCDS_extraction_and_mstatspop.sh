grep "ID=lmut" scaffold_9 | uniq | sort | awk '{print $9}' 
awk -F';' '{split($1,a,"="); print a[2]}' 


while read -r gene_id; do
    grep "Parent=${gene_id}" scaffold_9.gff | awk '$3 == "CDS"'
done < gene_ids.txt > cds_sequences.txt

#cleanfasta2 -infile scaffold_9.fas -outfile gene_Lmutsamples.fas -maxMissing 100 -format fasta -gff paupergeneannote_Lmut.gff -scaffold scaffold_9 -include 1

#For each line of cds_annots (each CDS), run the
#cleanfasta for all the samples in scaffold_9.fas
#the outfile should be composed of all the samples
#cut along the .gff intervals
#and each outfile should have the name of the CDS code
#and the corresponding gene code

#Save the names of all CDS codes
#awk '{print $9}' cds_annots > test3
#awk -F';' '{split($1,a,"="); split($2,b,"="); print a[2] "." b[2]}' test3 > CDScodes

#Save only the annotations which are unique
sort -u cds_annots > cds_annots.sort


speciesOfInterest.names.txt has a list of species from the NAP region with more than 3 samples
NAP_samplesTspecies contains a list of all the samples for the NAP region species 
	column 1 is sample names
	column 2 is genus name (Lupinus)
	column 3 is species names 
	separator " "
samplesOfInterest.names.txt contains a list of the samples (from NAP_samplesTspecies)

while IFS= read -r species; do
    # Filter lines based on the species and redirect the output to a new file
    awk -v species="$species" '$3 == species { print $1 }' NAP_samplesTspecies > "${species}.samples.txt"
done < speciesOfInterest.names.txt

directory samplesOfInterestPERspecies/ contains ${species}.samples.txt

|||||||||||||||||||||||||||||
### INPUT WE NEED

# list of sample names for eahc species, in file speciesOfInterest.names.txt
Lup2327
Lup2328
Lup2334


# gff with entries for scaffold_9
../scaffold_9

# list of gene names
allGenesCHR9.txt

directory gene_annots contains the .gff files for each gene

#while read species in speciesOfInterest.names.txt

mkdir fastas.scaffold9
mkdir fastas.cds

while read species; do
	currentsamples=samplesOfInterestPERspecies/$species.samples.txt

	#create a fasta file for the full scaffold_9 containing
	#for the 3 or 4 samples of this species

	cleanfasta2 -infile scaffold_9.fas -outfile fastas.scaffold9/scaffold_9.$species.fas -maxMissing 100 -format fasta -samples $currentsamples

	#for the genes in allGenesCHR9.txt

	while read gene; do

		#get the CDS 
		gff=gene_annots/$gene.t1.txt
		cleanfasta2 -infile fastas.scaffold9/scaffold_9.$species.fas -outfile fastas.cds/$species.$gene.cds.fas -maxMissing 100 -verbose 1 -format fasta -gff $gff -scaffold scaffold_9 -include 1 -feature CDS

	done < allGenesCHR9.txt

done < speciesOfInterest.names.txt


#Verify missingness

#see stats
for file in *.cds.fas ; do

	cleanfasta2 -infile "$file" -outfile "../gene.stats/$file.stats" -maxMissing 100 -format stats

done

#select for files with only missingness < 0.3 (30%)
#and 'ping' files that have more than one
#sample with > 0.3 missingness

nano script_missing3ness30
mkdir species.genes.fewgoodsamples #in allsequences
mkdir allGoodSamples.list #in allsequences
mkdir fastas.GoodCDS #in allsequences

script_missing3ness30
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#!/bin/bash
#runs from directory allsequences

echo "Samples Values Species.Gene" > "./poorsamples.txt"

for gene in gene.stats/*; do
        output=${gene%.cds.fas.stats}
        cds=${gene%.stats}

        good_samples_file="allGoodSamples.list/${output##*/}.goodsamples.txt"
        poor_samples_file="./poorsamples.txt"
        few_samples_link="species.genes.fewgoodsamples/${output##*/}.goodsamples.few"

        good_samples_count=0
        total_samples_count=0

        while read -r sample; do
                value=$(echo "$sample" | awk -F ' ' '{print $2}')
                label=$(echo "$sample" | awk -F ' ' '{print $1}')

                #check the value of mislsingness and if it is

                if (( $(echo "$value < 0.3" | bc -l) )); then
                        echo "$label $value" >> "$good_samples_file"
                        ((good_samples_count++))
                else
                        echo "$label $value ${output##*/}" >> "$poor_samples_file"
                fi

                ((total_samples_count++))

                #if more than one sample had > 30% missingness
                #by checking if nº samples in output =< nº samples in input - 2

        done < $gene

        if [[ -e "$good_samples_file" ]]; then

                line_count=$(wc -l < "$good_samples_file")

                if (( line_count < 3 )); then

                        mv "$good_samples_file" "$few_samples_link"

                fi
        else
                echo "$gene" >> "./species.genes.WithoutGoodSamples"
        fi

        #make a "ping" directory to store symbolic links to the species with more than one sample with more than 30% missingness
        #create cds fasta files containing only the good samples for each gene in each species
        #currentsamples=$(awk -F ' ' '{print $1}' $good_samples_file)
done

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

#for samples in allGoodSamples list
#extract the good cds into file fastas.GoodCDS

for gene in species.genes.fewgoodsamples/*; do
	sample=${gene%.goodsamples.few}
	list=$(awk '{print $1}' "$gene")

	echo "$list" > currentsamples.txt

	cds_file_name="fastas.cds/${sample##*/}.cds.fas"
	outfile="fastas.fewCDS/${sample##*/}.fewCDS.fas"

	cleanfasta2 -infile "$cds_file_name" -outfile "$outfile" -maxMissing 100 -samples "currentsamples.txt" 

done


#i need a list of genes like so:
lmut_g52410
lmut_g52411
lmut_g52412

(allGenesCHR9.txt)

#the loop is going to circle through the gene list
#using lmut_gXXXXX as a variable
#by doing *.$gene.* I should be able to locate
#every iteration of that gene in GoodCDS

#FOR EACH GENE
for gene in $(cat allGenesCHR9.txt); do 

	while read species; do #in speciesOfInterest.names.txt
		

		#instead of cleanfasta, I already have, for each gene
		#each species cds samples
		#so i just take the cds file corresponding to the
		#current gene iteration for the current species iteration

		currentfile="fastas.GoodCDS/$species.$gene.GoodCDS.fas"
		ref_fasta="fastas.ref/$species.$gene.fas.ref"
		ref_fai="ref.fai/$species.$gene.fas.ref.fai"

		if [[ -e "$currentfile" ]]; then #file exists

			#make a reference file

			head -n 2 $currentfile > $ref_fasta

			#create .fai for the reference file

			samtools faidx $ref_fasta

			#mstats -i infile (fasta)

			mstatspop -f fasta -i $currentfile -n $ref_fai -o 1 -N 1 3 -p 2 -u 1 -G 0 > mstats.$species.$genes 

		else

			echo "ERROR: $currentfile does not exist. This gene may not have sufficient samples of this species" 


	done < speciesOfInterest.names.txt

done

#NOTE: i forgot to make an output directory

mkdir results.mstats
mv mstats* results.mstats/


#create a full table (?)

cat results.mstats/mstats.* > mstats.GENERAL

#fix full table





















































	# Requires module load: bcftools/1.17 cleanfasta2/161123 mstatspop/20220915 samtools/1.17
	species_list=
	samples_list=
	input_fas= #fastas.GoodCDS

###SPECIESLIST (1st argument) is speciesNAP.list.txt
# angustiflorus
# arboreus
# etc (1 column)
###SAMPLELIST (2nd argument) is goodNAP.txt
# LupXXX species
# Lup XXY species
###INPUTFAS (3rd argument) is allExons.cnct.fas
# >LupXXX
# AGATCGATC (concatenated)

mkdir $OUTPUTDK

while read species
do
	grep "$species" $SAMPLELIST | cut -d ' ' -f 1 > currentSpecies.txt
	cleanfasta2 -infile $INPUTFAS -outfile $species.fas -maxMissing 100 -samples currentSpecies.txt
	head -n 2 $species.fas > $species.fas.ref

	#module load bcftools/1.17
	#module load samtools/1.17
	samtools faidx $species.fas.ref

	#module load mstatspop/20220915
	mstatspop -f fasta -i $species.fas -n $species.fas.ref.fai -o 1 -N 1 3 -p 2 -u 1 -G 0 mstats.$INPUTFAS.$SPECIESLIST

done < $SPECIESLIST





done
