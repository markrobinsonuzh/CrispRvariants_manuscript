gen=../idx/danRer7.fa
psl2sam=psl2sam.pl
header=test_header.sam
addseq=rackJ/scripts/samGetSEQfast.pl
outdir=../bam/blat_global/

for f in ../merged_split/*
do
  base=$(basename $f)
  base=${outdir}${base%%.*}  
  fq=$f"/SRR1769728_merged_pear.assembled.fastq.gz"
  fa=${base}".fa"
  out=${base}"_temp.sam"
  outf=${base}".sam"
  cp $header $out
  
  # Create fasta from fastq
  zcat $fq | awk '(NR % 4 == 1){print ">"substr($1, 2, length($1))}(NR % 4 == 2){print $1}' > $fa
  
  # run mapping
  blat -extendThroughN -out=pslx $gen $fa ${base}.pslx
  pslReps -singleHit ${base}.pslx ${base}_f.pslx delme_g.psr
  
  # Convert to sam format, add sequence, then convert to bam format
  ${psl2sam} ${base}_f.pslx >> ${out}
  ${addseq} ${out} ${fa} >> ${outf} && rm $out
  samtools view -Sb $outf > ${base}"_temp.bam"
  rm $fa $out ${base}.pslx ${base}_f.pslx
  (samtools sort ${base}"_temp.bam" ${base} && rm ${base}"_temp.bam") &
done
