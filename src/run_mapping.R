map <- function(in_fq, out_prefix){
  index <- "../idx/danRer7.fa"
  cmd <- paste("bwa mem -t 12 %s %s | samtools view -Sb - > %s_temp.bam",
                "samtools sort %s_temp.bam %s && rm %s_temp.bam",
                "samtools index %s.bam", sep = "; ")
  cmd <- do.call(sprintf, c(cmd, as.list(c(index,in_fq,rep(out_prefix, 5)))))
  print(cmd)
  system(cmd)
}


# Map merged reads
map("../fastq/SRR1769728_merged_seqprep.fastq.gz", "../bam/merged_seqprep")


# Map merged split reads
guides <- list.files("../merged_split")
guides <- guides[! guides %in% c("tolerant_pear", "strict_pear")]

dummy <- lapply(guides, function(gd){
        pear_template <- "../merged_split/%s/SRR1769728_merged_pear.assembled.fastq.gz"
        seqprep_template <- "../merged_split/%s/SRR1769728_merged_seqprep.fastq.gz"
        seqprep55_template <- "../merged_split/%s/SRR1769728_merged_seqprep_l55.fastq.gz"
        pear_bam <- "../bam/merged_split_pear/%s"
        seqprep_bam <- "../bam/merged_split_seqprep/%s"
        seqprep55_bam <- "../bam/merged_split_l55_n90/%s"
        map(sprintf(pear_template, gd), sprintf(pear_bam, gd))
        map(sprintf(seqprep_template, gd), sprintf(seqprep_bam, gd))
        map(sprintf(seqprep55_template, gd), sprintf(seqprep55_bam, gd))
})

# Map split_merged (guides consistent across all conditions)
dummy <- lapply(guides, function(gd){
        pear_template <- "../split_merged/%s/SRR1769728_merged_pear.fastq.gz"
        pear_bam <- "../bam/split_merged_pear/%s"
        map(sprintf(pear_template, gd), sprintf(pear_bam, gd))
})

# Map primer comparison
dummy <- lapply(guides, function(gd){
   tolerant <- "../merged_split/tolerant_pear/%s.fastq"
   strict <- "../merged_split/strict_pear/%s.fastq"
   tolerant_out <- "../bam/merged_split_pear_tolerant/%s"
   strict_out <- "../bam/merged_split_pear_strict/%s"
   map(sprintf(tolerant, gd), sprintf(tolerant_out, gd))
   map(sprintf(strict, gd), sprintf(strict_out, gd))
})
