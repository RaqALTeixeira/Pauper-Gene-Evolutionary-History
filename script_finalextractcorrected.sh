#!/bin/bash

# Run from directory DNDS


for TEMPDIR in temp*/; do

	cd $TEMPDIR

	for file in OUT_FULLcodemlAllGenes/*; do

		gene=$(basename "$file")
		gene_name=${gene%.FULLcodeml} #the variable gene_name is now just lmut_XXXX

		FULLcodeml="$gene_name.FULLcodeml"
		file_name="../OUT_PrmtrsCodeml/$gene_name/$gene_name"

		#mkdir ../OUT_PrmtrsCodeml/$gene_name !!!! in DNDS directory

		if [[ -s "$FULLcodeml" ]]; then #FULLcodeml exists and is not empty

			# Finds the section of the 'FULLcodeml' dedicated to each model
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
			' "$FULLcodeml"

			# Once intermediate files for each model have been created
			# this following part finds the common motifs of each parameters
			# and extracts them onto a new file prmters_

			prmtrs_file="../OUT_PrmtrsCodeml/$gene_name/prmtrs_$gene_name" 

			rm -f $prmtrs_file #these already exist, and need to be replaced before initiating the corrected extraction 

			for file in "$file_name*"; do

				grep "Model" "$file" >> "$prmtrs_file"

				grep "lnL" "$file" >> "$prmtrs_file" ### GO BACK AND FIX

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

				rm "$file" # eliminates the intermediate 'section' file before the next loop
			done


		else
			echo "Error: $gene is an empty file"
			echo "$gene_name" >> "../2manyambiguity.codeml.txt"
		fi

	done

	echo "RUN COMPLETED for $TEMPDIR"

	#cp oputput file to outputfolder --> no need, things are already being redirected, hopefully

	cd ..
done