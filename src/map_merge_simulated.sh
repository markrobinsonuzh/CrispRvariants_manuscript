index="../idx/danRer7.fa"
for r1 in ../simulation/*_sim1.fq
do
  r2="${r1/sim1/sim2}"
  out="${r1/_sim*/}"
  pear-0.9.4-64 -j 12 -f ${r1} -r ${r2} -o ${out}
  rm ../simulation/*unassembled* ../simulation/*discarded*
  bwa mem $index ${out}.assembled.fastq | samtools view -Sb - > ${out}_merged_temp.bam && rm ${out}.assembled.fastq
  samtools sort ${out}_merged_temp.bam ${out}_merged && rm ${out}_merged_temp.bam
  samtools index ${out}_merged.bam
done
