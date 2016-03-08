#!/bin/bash

# THIS SCRIPT IS AN EXCERPT FROM AMPLICONDIVIDER_DRIVER.SH BY MATTHEW LAFAVE
# https://github.com/mlafave/ampliconDIVider

function parseBam {
    BAM=$1
    SAMPLE_NAME=$2
    LEFT=$3
    RIGHT=$4    

    samtools view $BAM \
    | awk '$6 ~ /[DI]/' \
    | perl/sam_all_div_linebreak.pl \
    | awk -F_ '{if($1 == "###"){print}else{print $0"\t"$2"\t"$2+($3-1)}}' \
    | awk -v LEFT="${LEFT}" -v RIGHT="${RIGHT}" '($3 >= LEFT && $3 <= RIGHT) || ($2 >= LEFT && $2 <= RIGHT) {print} $1 == "###" {print}' \
    | cut -f1 \
    | awk '{if($1 == "###"){if(out){print s"\t"out; out=""; delete a; s=0}}else{out=out$1"\t"; if(/^D/){split($0,a,"_"); s -= a[3]};if(/^I/){split($0,a,"_"); s += a[3]} }}END{if(out){print s"\t"out}}' \
    | awk '{ if($1%3){print "y\t"$0}else{print "n\t"$0} }' \
    > in-range-mutant-reads_${SAMPLE_NAME}

    cut -f1 in-range-mutant-reads_${SAMPLE_NAME} \
    | sort \
    | uniq -c \
    > fs_sums_${SAMPLE_NAME}
	
    # Identify the number of mutations with and without a frameshift	
    NOSHIFT=`head -1 fs_sums_${SAMPLE_NAME} | awk '{print $1}'`
    YESSHIFT=`tail -1 fs_sums_${SAMPLE_NAME} | awk '{print $1}'`
	# Count the number of total reads in the relevant BAM, and print a
	# summary line indicating total reads, reads without a variant near the
	# CRISPR site, inframe & frameshift variant reads, and calculations of
	# the overall mutation rate, the inframe rate, and the frameshift rate.
	
	samtools view ${BAM} \
	| wc -l \
	| awk -v target="${SAMPLE_NAME}" \
	-v noshift="$NOSHIFT" \
	-v yesshift="$YESSHIFT" \
	-v OFS="\t" \
	'{print target,$1,$1-(noshift+yesshift),noshift,yesshift,(noshift+yesshift)/$1,noshift/$1,yesshift/$1}' \
	>> frameshift_summary_${SAMPLE_NAME}

    rm in-range-mutant-reads_${SAMPLE_NAME} fs_sums_${SAMPLE_NAME}
}

