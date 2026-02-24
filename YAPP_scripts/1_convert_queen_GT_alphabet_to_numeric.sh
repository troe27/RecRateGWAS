#!/bin/bash

# Check if the input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input.vcf>"
    exit 1
fi

input_file="$1"

# Process the VCF file
awk -v OFS="\t" '
BEGIN {
    # Initialize REF and ALT
    ref = ""
    alt = ""
}
{
    # Process header lines
    if ($1 ~ /^#/) {
        print $0
        next
    }

    # Store REF and ALT from the current data line
    ref = $4
    alt = $5

    # Convert genotypes starting from the queen column
    for (i = 1662; i <= NF; i += 2) {
        if (i + 1 <= NF) {
                # Determine the genotype based on REF and ALT
                if ($i == ref && $(i + 1) == ref) {
                    $i = "0|0";
                } else if ($i == ref && $(i + 1) == alt) {
                    $i = "0|1";
                } else if ($i == alt && $(i + 1) == ref) {
                    $i = "1|0";
                } else if ($i == alt && $(i + 1) == alt) {
                    $i = "1|1";
                } else if ($i == alt && $(i + 1) != ref && $(i + 1) != alt) {
                    $i = "1|.";
                } else if ($i == ref && $(i + 1) != ref && $(i + 1) != alt) {
                    $i = "0|.";
                } else if ($i != ref && $i != alt && $(i + 1) == alt) {
                    $i = ".|1";
                } else if ($i != ref && $i != alt && $(i + 1) == ref) {
                    $i = ".|0";
                } else {
                    $i = ".|.";
                }
                # Combine the two columns into one
                $(i + 1) = "";  # Clear the next column to avoid printing it
            }
        }

    # Print the modified line without empty columns
    output = "";
    for (j = 1; j <= NF; j++) {
        if ($(j) != "") {
            output = output (output == "" ? "" : OFS) $(j);
        }
    }
    print output;
}' "$input_file" 