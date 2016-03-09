r1=../fastq/SRR1769728_1.fastq.gz
r2=../fastq/SRR1769728_2.fastq.gz

amplicons=../annotation/CRISPResso_pooled_layout.txt
idx=../idx/danRer7.fa
nm=crispresso_pooled_mixed

CRISPRessoPooled -r1 ${r1} -r2 ${r2} -f ${amplicons} -x ${idx} --name ${nm} -p 12 --window_around_sgrna 5

mv CRISPRessoPOOLED_on_crispresso_pooled_mixed ../results/CRISPResso_pooled_mixed 
