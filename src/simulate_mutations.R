library("GenomicRanges")
library("Biostrings")
library("rtracklayer")

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


freqs <- read.table("../annotation/Shah_mutation_weights.txt", sep = "\t")
# Set a low value for large indels
freqs$x[freqs$Variant > 10] <- 1e-10
freqs$x <- freqs$x/sum(freqs$x) 

amplicons <- read.table("../annotation/Shah_cut_sites.txt", sep = "\t",
                stringsAsFactors = FALSE)
colnames(amplicons) <- c("name","original","target_loc")
amplicons[,"target_loc"] <- as.integer(amplicons[,"target_loc"])
amplicons <- amplicons[1:20,]

# Get guides for ampliconDIVider
guides <- rtracklayer::import("../annotation/shah_guides.bed")
guides <- guides[match(amplicons$name, guides$name)]
guides <- guides + 5
adiv_out <- "amplicondivider_simulation_commands.sh"
cat('cd ../ampliconDIVider-master\nsource ampliconDIV_minimal.sh\n',
    file = adiv_out)


sample_seqs <- function(n_mut, n_original, n_offtarget, amplicons,
                        out_dir, sim_file, mut_frac = 0.3, 
                        log_file = NA, nvars = 10, crispresso_file =  NA,
                        read_len = 200){
                       
  dummy <- lapply(1:nrow(amplicons), function(i){
    a_rw <- amplicons[i, ]
    original <- DNAString(a_rw$original)
    target_loc <- as.integer(a_rw["target_loc"])
    gd_name <- a_rw["name"]    
    guide <- substr(original, target_loc-17, target_loc + 5)
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
    print(sprintf("Writing sequences to: %s", out_fname))
    cat(new_seqs, file = out_fname)
    
    # art_illumina command with seed = 30
    sim_template <- "art_illumina -amp -rs 30 -f 10 -l %s -p -ss MS -na -i %s -o %s\n"
    cat(sprintf(sim_template, read_len, out_fname, gsub(".fa", "_sim", out_fname)), 
        file = sim_file, append = TRUE)

    # Commands for CRISPResso
    f1 <- gsub(".fa", "_sim1.fq", out_fname)
    f2 <- gsub(".fa", "_sim2.fq", out_fname)
    crispresso_dir <- "crispresso"
    crispresso_template <- "CRISPResso -r1 %s -r2 %s -a %s -g %s -o %s -w 5\n"
    cat(sprintf(crispresso_template, f1,f2, original, guide, crispresso_dir),
        file = crispresso_file, append = TRUE)
        
    # Commands for ampliconDIVider
    adiv_dir <- "../simulation/amplicondivider"
    adiv_tmp1 <- 'samtools view -hb %s %s > temp.bam'
    adiv_tmp2 <- 'parseBam temp.bam %s %s %s; rm temp.bam'
    adiv_tmp3 <- 'mv frameshift_summary_%s %s\n\n' 
    bam <- gsub(".fa", "_merged.bam", out_fname)
    gd <- guides[i]
    gd_rng <- sprintf("%s:%s-%s", seqnames(gd), start(gd)-200, end(gd)+200)
    base <- gsub(".fa", "", basename(out_fname))
    a1 <- sprintf(adiv_tmp1, bam, gd_rng) 
    a2 <- sprintf(adiv_tmp2, base, start(gd) - 5, end(gd) + 5)
    a3 <- sprintf(adiv_tmp3, base, file.path(adiv_dir, paste0(base, ".txt")))
    cat(paste(a1,a2,a3, sep = "\n"), file = adiv_out, append = TRUE)
  })
}

sim_cmds <- "simulation_commands.sh"
crispresso_cmds <- "crispresso_simulation_commands.sh"

for (nofftargets in c(100,33,0)){
    # 0% efficient  
    sample_seqs(0,300, nofftargets, amplicons, "../simulation", 
               sim_file = sim_cmds, crispresso_file = crispresso_cmds)
    
    # 33% efficient
    sample_seqs(100,200, nofftargets, amplicons, "../simulation", 
               sim_file = sim_cmds, crispresso_file = crispresso_cmds)
    
    # 66% efficient
    sample_seqs(200,100, nofftargets, amplicons, "../simulation", 
               sim_file = sim_cmds, crispresso_file = crispresso_cmds)
               
    # 90% efficient
    sample_seqs(270,30, nofftargets, amplicons, "../simulation", 
               sim_file = sim_cmds, crispresso_file = crispresso_cmds)
}

    
# Write amplicons and guides into file for CRISPRessoPooled 
original <- amplicons$original
guide <- substr(original, amplicons$target_loc-17, amplicons$target_loc + 5)
amplicon_inf <- paste(amplicons$name, original, guide, sep = "\t", collapse = "\n")
cat(amplicon_inf, file = "../simulation/merged/crispresso_pooled_amplicons.txt")




cat("\n\nmv crispresso/* ../simulation/crispresso; rmdir crispresso\n", 
    file = crispresso_cmds, append = TRUE)