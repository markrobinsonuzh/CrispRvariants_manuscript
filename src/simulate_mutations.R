library("GenomicRanges")
library("Biostrings")


deleteNucleotides <- function(original, original_range, start, width){
    y <- IRanges(start, width = width)
    result <- original[setdiff(disjoin(c(y,original_range)), y)]
    result
}

insertNucleotides <- function(original, original_range, start, width){
  ins <- sample(c("A","C","T","G"), width, replace = TRUE)
  ins <- paste0(ins, collapse = "")
  rngs <- IRanges(c(1,start+1), end = c(start, width(original_range)))
  sqs <- as.character(Views(original, rngs))
  result <- DNAString(paste0(c(sqs[1], ins,sqs[2]), collapse = ""))
  result
}

mutateNucleotides <- function(original, fraction){
  nsample <- floor(fraction*length(original))
  mutate <- sample(c(1:length(original)), nsample)
  nucs <- strsplit(as.character(original),"")[[1]]
  replace_nucs <- list("A" = c("C","G","T"), "C" = c("A","G","T"),
                       "T" = c("A","C","G"), "G" = c("A","C","T"))
  replace_idx <- sample(c(1,2,3), nsample, replace = TRUE)
  replace_idx <- replace_idx + seq(0, 3*(nsample-1), by = 3)
  new_nucs <- unlist(replace_nucs[nucs[mutate]])[replace_idx]  
  nucs[mutate] <- new_nucs
  paste0(nucs, collapse = "")
}

set.seed(30)


freqs <- read.table("~/mutation_weights.txt", sep = "\t")
# Set a low value for large indels
freqs$x[freqs$Variant > 10] <- 1e-10
freqs$x <- freqs$x/sum(freqs$x) 

amplicons <- read.table("../annotation/Shah_cut_sites.txt", sep = "\t",
                stringsAsFactors = FALSE)
colnames(amplicons) <- c("name","original","target_loc")
amplicons[,"target_loc"] <- as.integer(amplicons[,"target_loc"])

# Remove amplicons where the cut site is too close to the read length (for 150bp)
# For read length 200, one read will still cover target
#keep <- ! (amplicons$target_loc >= 140 & amplicons$target_loc <= 160)
#amplicons <- amplicons[keep,]
amplicons <- amplicons[1:20,]

sample_seqs <- function(n_mut, n_original, n_offtarget, amplicons,
                        out_dir, sim_file, mut_frac = 0.3, 
                        log_file = NA, nvars = 10, crispresso_file =  NA,
                        read_len = 200){
                       
  dummy <- apply(amplicons, 1, function(a_rw){
    original <- DNAString(a_rw["original"])
    target_loc <- as.integer(a_rw["target_loc"])
    gd_name <- a_rw["name"]    
    guide <- substr(original, target_loc-17, target_loc + 5)
    print(unname(c(gd_name, as.character(guide))))
    original_range <- IRanges(1, length(original))
    amp_loc <- target_loc + freqs$Location - 1
    subsamp <- sample(1:nrow(freqs), nvars, prob = freqs$x)
    row_idxs <- sample(c(1:nrow(freqs))[subsamp], n_mut, 
                  prob = freqs$x[subsamp], replace = TRUE)
    
    subf <- freqs[subsamp,]
    var_labs <- paste(subf$Location, paste0(subf$Variant, subf$var_type), sep = ":")
    var_freqs <- table(factor(row_idxs, levels = subsamp))
    names(var_freqs) <- var_labs

    new_seqs <- lapply(row_idxs, function(i){
      rw <- freqs[i,]
      if (rw$var_type == "D") func <- deleteNucleotides
      if (rw$var_type == "I") func <- insertNucleotides
      new_seq <-  func(original, original_range, amp_loc[i], rw$Variant)
      new_seq                       
    })
    
    mut_seqs <- lapply(seq_len(n_offtarget), function(i) mutateNucleotides(original, mut_frac))

    out_fname <- file.path(out_dir, sprintf("%s_%smut_%swt_%sofftarget_%sreadlen.fa", 
                  gd_name, n_mut, n_original, n_offtarget, read_len))
    out <- file(out_fname)

    result <- c(new_seqs, mut_seqs, replicate(n_original, original))
    result_names <- c(rep(">var", n_mut), 
                      rep(">offtarget", n_offtarget),
                      rep(">original", n_original))


    print(sprintf("Length variant: %s Length unmutated: %s Length offtarget %s",
          length(new_seqs), n_original, length(mut_seqs)))
    stopifnot(length(result_names) == (length(new_seqs)+ length(mut_seqs) + n_original))

    result_names <- paste(result_names, c(seq_len(n_mut), 
                           seq_len(n_offtarget), seq_len(n_original)), sep = "_")

    new_seqs <- paste(result_names, sapply(result, as.character), sep = "\n")
    new_seqs <- paste0(new_seqs, collapse = "\n")
    
    write(new_seqs, out)
    close(out)
    
    # art_illumina command with seed = 30
    sim_template <- "art_illumina -amp -rs 30 -f 10 -l %s -p -ss MS -na -i %s -o %s\n"
    cat(sprintf(sim_template, read_len, out_fname, gsub(".fa", "_sim", out_fname)), 
        file = sim_file, append = TRUE)

    f1 <- gsub(".fa", "_sim1.fq", out_fname)
    f2 <- gsub(".fa", "_sim2.fq", out_fname)
    crispresso_dir <- "crispresso"
    crispresso_template <- "CRISPResso -r1 %s -r2 %s -a %s -g %s -o %s -w 5\n"
    cat(sprintf(crispresso_template, f1,f2, original, guide, crispresso_dir),
        file = crispresso_file, append = TRUE)
    
  })
}

sim_cmds <- "~/scratch/simulation_commands.sh"
crispresso_cmds <- "~/scratch/crispresso_commands.sh"

#for (read_len in c(150, 200))
for (nofftargets in c(100,33,0)){
    # 0% efficient  
    sample_seqs(0,300, nofftargets, amplicons, "~/scratch", 
               sim_file = sim_cmds, crispresso_file = crispresso_cmds)
    
    # 33% efficient
    sample_seqs(100,200, nofftargets, amplicons, "~/scratch", 
               sim_file = sim_cmds, crispresso_file = crispresso_cmds)
    
    # 66% efficient
    sample_seqs(200,100, nofftargets, amplicons, "~/scratch", 
               sim_file = sim_cmds, crispresso_file = crispresso_cmds)
               
    # 90% efficient
    sample_seqs(270,30, nofftargets, amplicons, "~/scratch", 
               sim_file = sim_cmds, crispresso_file = crispresso_cmds)
}

    




# cnsta
#original <- DNAString("CATGCACTGCTCTCAGAAAGCCTGAAGCACAATTACTGCTGAACTATGAGTCAGTTCAGTTAAAGGCCGGCTATTTCTTTACTGTTCCTGATTATTCTGTGTCTCCCAGAGGCAGAGCGGAGCGGCGTCCAGTGCCTGAGCGGATCTTGATGGATGACGGTCAGTGGCAGAGCGACGGCCCTGTAAGGAGAGGAGGAGCGACTGAGCCGGAGCTACAGAGCTGCTCTGAGTCCTTCCGCAGCCCTGATGACAACCACAACCGGCTGCT")
#original_range <- IRanges(1, length(original))
#target_loc <- 185 

#out <- file("~/Desktop/out_seqs.fa")

# generate - as for Sup 1, without min.freq = 2
# dat <- dat[,c(1:3)]
# x <- aggregate(dat$Frequency, dat[c(1,2)], mean)
# x$var_type  <- gsub("[0-9]", "", x$Variant)
# x$Variant <- gsub("I|D", "", x$Variant)
# write.table(x, file = "~/mutation_weights.txt", sep = "\t", quote = FALSE)


# amplicon sequencing simulation with paired-end reads
#  art_illumina -amp -f 100 -l 198 -p -ss MS -na -i amp_reference.fa -o amplicon_simulation



#bam <- "amplicon_simulation_bwa_s.bam"
#guides <- rtracklayer::import("~/crispRvariants_supplementary/annotation/guides.bed")
#guide <- guides[guides$name == "cnsta"]
#library(BSgenome.Drerio.UCSC.danRer7) 
#danRer7 <- BSgenome.Drerio.UCSC.danRer7
#guide <- guide + 5
#reference <- getSeq(danRer7, guide)
#library(CrispRVariants)
#cset <- readsToTarget(bam, target = guide, reference= reference, target.loc = 22, split.snv = FALSE)
