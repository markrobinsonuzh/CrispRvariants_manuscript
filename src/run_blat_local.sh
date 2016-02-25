primers=../annotation/shah_primers.bed
gen=../idx/danRer7.fa
amplicon_refs=../annotation/Shah_amplicon_refs/
psl2sam=psl2sam.pl
addseq=rackJ/scripts/samGetSEQfast.pl
header=test_header.sam
samtools view -H ../bam/original.bam > $header

while read chr start end name score strand; 
do 
  echo $chr $start $end
  amplicon=${amplicon_refs}${name}".fa"
  samtools faidx $gen $chr":"$start"-"$end | sed 's/:.*//g' > ${amplicon} 
  fa_in="../merged_split/"${name}"/SRR1769728_merged_pear.assembled.fastq.gz"
  base="../bam/blat_local/${name}"
  fa=${base}".fa"
  out="../bam/blat_local/"${name}"_temp.sam"
  outf="../bam/blat_local/"${name}".sam"
  cp ${header} ${out}

  # Convert fastq to fasta
  zcat ${fa_in} | awk '(NR % 4 == 1){print ">"substr($1, 2, length($1))}(NR % 4 == 2){print $1}' > ${fa}  
  
  # Map with blat 
  blat -extendThroughN -out=pslx $amplicon $fa ${base}.pslx
  pslReps -singleHit ${base}.pslx ${base}_f.pslx delme.psr
  
  # Convert to sam format, add the sequence
  ${psl2sam} ${base}_f.pslx >> ${out}
  ${addseq} ${out} ${fa} >> ${outf} && rm $out

  # Convert to bam
  samtools view -Sb $outf > ${base}"_temp.bam"
  rm $fa $outf ${base}.pslx ${base}_f.pslx
  (samtools sort ${base}"_temp.bam" ${base} && rm ${base}"_temp.bam") &

done < ${primers}
