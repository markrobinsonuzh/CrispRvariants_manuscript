cd ../simulation
index="/home/helen/CrispRVariants_paper/idx/danRer7.fa"
for r1 in *_sim1.fq
do
  r2="${r1/sim1/sim2}"
  out="${r1/_sim*/}"
  bwa mem $index $r1 $r2 | samtools view -Sb - > ${out}_temp.bam
  samtools sort ${out}_temp.bam ${out} && rm ${out}_temp.bam
  samtools index ${out}.bam
done
cd ../src
