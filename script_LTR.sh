#!/bin/bash

# This code takes all prmtrs_files output from
# Codeml containing just the names, likelihoods and parameters
# of M1 and M2a

#Input prmtrs_files are in a directory with the name of the gene
#and look like thus (#1# is the line number, ignore for prmtrs_file disposition)

#1#Model 1: NearlyNeutral (2 categories)
#2#lnL(ntime: 49  np: 52):  -2933.673985      +0.000000
#3#dN/dS (w) for site classes (K=2)
#4#
#5#p:   0.75063  0.24937
#6#w:   0.00000  1.00000
#7#
#8#Model 2: PositiveSelection (3 categories)
#9#lnL(ntime: 49  np: 54):  -2812.710066      +0.000000
#10#dN/dS (w) for site classes (K=3)
#11#
#12#p:   0.51406  0.40113  0.08481
#13#w:   0.00000  1.00000 31.56577

TABLE="$1"
echo "Gene" "Model" "df" "LR" > $TABLE


mkdir -p M1vM2_Genes

for dir in OUT_PrmtrsCodeml/*; do

    sample_name=$(basename $dir)
    prmtrs_file="$dir/prmtrs_$sample_name"

    output="M1vM2_Genes/$sample_name.m1vm2.txt"

    lnL1=$(sed -n '2p' "$prmtrs_file" | awk '{print $5}')
    lnL2=$(sed -n '9p' "$prmtrs_file" | awk '{print $5}')

    lnL1=$(echo "$lnL1" | awk '{print ($1<0)?-$1:$1}')
    lnL2=$(echo "$lnL2" | awk '{print ($1<0)?-$1:$1}')

    LR=$(echo "2 * ($lnL1 - $lnL2)" | bc)
    critvalue=5.991 #df = 2 ; p-value = 0.05%

    if (( $(echo "$LR > $critvalue" | bc -l) )); then # Model 2 significantly fits the data better than model 1
        sed -n '8p' "$prmtrs_file" > "$output"
        sed -n '12p' "$prmtrs_file" >> "$output"
        sed -n '13p' "$prmtrs_file" >> "$output" # Append to output

        echo "$sample_name" "Model 2" "2" "$LR" >> $TABLE

    else
        sed -n '1p' "$prmtrs_file" > "$output"
        sed -n '5p' "$prmtrs_file" >> "$output"
        sed -n '6p' "$prmtrs_file" >> "$output" # Append to output

        echo "$sample_name" "Model 1" "2" "$LR" >> $TABLE

    fi

done
