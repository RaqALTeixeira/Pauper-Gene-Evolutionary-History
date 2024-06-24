#!/bin/bash

# This script takes the .phy and .tre files previously created by
# fasta2paml and runs them through the codeml comand from paml to obtain
# the parameters w p and the likelihood of the model

# It requires a template control file for codeml
# With all the conditions one wishes to apply WITHOUT the seqfile, outfile and treefile
# It also requires a list of the genes that will be cycled through
# this is a text file with the name of the genes

# Run from directory DNDS

TEMPLATERUN="../$1"
GENELIST="../$2"
# is given by xargs
# from Path file to each gene list in the parallel execution


#mkdir OUT_FULLcodemlAllGenes
mkdir -p oddGenes.codeml
RFOLDER="oddGenes.codeml/"
cd $RFOLDER

mkdir controlfiles4codeml
mkdir OUT_FULLcodemlAllGenes
#One of these directories is created within each
#temporary folder

errorfile="${GENELIST%.txt}.codeml2prmtrs.errorfile.txt"
echo "" > $errorfile #clear errorfile if there is any
#error file for each gene list will be stored in the temp folder

while read -r gene_name; do

	# create a control file for codeml for that gene by
	# adding the file names to the template

	if [[ -e "../forCODEML/$gene_name" ]]; then

		# cycles through the previously created directories for each gene, which contain the .phy and .tre

		control_file="../controlfiles4codeml/$gene_name.ctl"

		cp $TEMPLATERUN $control_file

		seqfile="../forCODEML/$gene_name/$gene_name.paml.phy"
		outfile="OUT_FULLcodemlAllGenes/$gene_name.FULLcodeml" # full out codeml is in temp folder
		treefile="../forCODEML/$gene_name/$gene_name.paml.tre"

		sed -i "1s|^|	seqfile = $seqfile\n|" $control_file
		sed -i "2s|^|	outfile = $outfile\n|" $control_file
		sed -i "3s|^|	treefile = $treefile\n|" $control_file #add the tree name files in a line
		sed -i "4s|^|\n|" $control_file
	else
		echo "Error: $gene_name is not a directory in forCODEML/" >> "$errorfile"

	fi

done < "$GENELIST"

echo "Finished preparing control files."
echo "Initiating CODEML"

while read -r gene_name; do

	# For the genes in the list, it looks for if the control file for codeml was created
	# proceeds to run codeml using the control file for that gene in the list

	control_file="../controlfiles4codeml/$gene_name.ctl"

	if [[ -e "$control_file" ]]; then

		codeml $control_file
	else
		echo "Error: couldn't find control file for the gene $gene_name" >> "$errorfile"
	fi

	outfile="OUT_FULLcodemlAllGenes/$gene_name.FULLcodeml" #SAME AS IN THE FIRST LOOP WHERE THE CONTROL FILE DEFINED THE OUTPUT FILE NAME
	file_name="../OUT_PrmtrsCodeml/$gene_name/$gene_name"

	mkdir ../OUT_PrmtrsCodeml/$gene_name # !!!! in DNDS directory

	if [[ -e "$outfile" ]]; then

		# Finds the section of the 'outfile' dedicated to each model
		# It looks for the common motifs "Model", which initiates model name
		# and "dN & dS for each branch" which all models display, after which theres a table of dN dS
		# In the intermediate space, the parameters p, w and likelihood can all be found

		awk -v prefix="$file_name" '
		BEGIN {
			RS="Model";
			ORS=""
		}

		NR > 1 {
			model = $1;
			model = gensub(/[^a-zA-Z0-9]/, "_", "g", model);
			if (match($0, /dN & dS for each branch/)) {
				print "Model" substr($0, 1, RSTART-1) > prefix "_" "model" model ".txt";
			}
		}
		' "$outfile"

		# Once intermediate files for each model have been created
		# this following part finds the common motifs of each parameters
		# and extracts them onto a new file prmters_

		for file in $file_name*; do

			prmtrs_file="../OUT_PrmtrsCodeml/$gene_name/prmtrs_$gene_name"

			grep "lnL" "$file" >> "$prmtrs_file"

			if [[ "$file" == *model0* ]]; then
				grep "omega" "$file" >> "$prmtrs_file"
			else
				awk '
				/for site classes/, /dN & dS for each branch/ {
					if (/dN & dS for each branch/) exit
					print
				}
				' "$file" >> "$prmtrs_file"
			fi

			rm "$file" #eliminates the intermediate 'section' file before the next loop
		done


	else
		echo "Error: couldn't find output codeml file for gene $gene_name" >> $errorfile
	fi

done < "$GENELIST"


echo "RUN COMPLETED"

#cp oputput file to outputfolder --> no need, things are already being redirected, hopefully

cd ..
