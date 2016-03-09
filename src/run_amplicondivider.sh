bam_template=".bam"
run_dir=../ampliconDIVider-master/
base_dir=../src
out_base="frameshift_summary_"
amplicon_div=../results/ampliconDIV_filtered/

source ${run_dir}ampliconDIV_minimal.sh

bam_dir=../bam/
#directories=( merged_split_pear
directories=( split_merged_pear
              merged_split_seqprep merged_split_l55_n90 
              merged_split_pear_tolerant merged_split_pear_strict)


while read name chr start end strand astart aend 
do
  for dr in ${directories[@]}
  do
    bam_fn=${name}${bam_template}
    cd ${run_dir}  
    samtools view -hb ${bam_dir}${dr}/${bam_fn} ${chr}":"${start}"-"${end} > ${bam_fn}
    parseBam ${bam_fn} ${name} $(($start - 5)) $(($end + 5))
    rm ${bam_fn}
    mv ${out_base}${name} ${amplicon_div}${dr}/${name}
    cd ${base_dir}
  done
done < ../annotation/Shah_guide_and_amplicon_locs.txt
