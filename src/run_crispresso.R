library(BSgenome.Drerio.UCSC.danRer7)
library(GenomicAlignments)
library(rtracklayer)

#____________________________________________________________
# Setup data

danRer7 <- BSgenome.Drerio.UCSC.danRer7
primers <- rtracklayer::import("../annotation/shah_primers.bed")
amplicons <- getSeq(danRer7, primers)
names(amplicons) <- primers$name

temp <- read.table("../annotation/shah_custom_amplicons.txt", sep = "\t")
custom <- temp$V2
names(custom) <- temp$V1

temp <- rtracklayer::import("../annotation/shah_guides.bed")
guides <- getSeq(danRer7, temp)
names(guides) <- temp$name 

merged_split <- list.files("../merged_split", full.names = TRUE)
out_dir <- "../results/test/CRISPResso_std_ref"
out_merged <- "../results/test/CRISPResso_std_ref_merged"
out_custom <- "../results/test/CRISPResso"

#____________________________________________________________

parse_results <- function(results_dir){
  results_f <- file.path(results_dir, "Quantification_of_editing_frequency.txt")
  system(paste0("echo '\n' >>", results_f))
  
  f <- file(results_f)
  lns <- readLines(f)
  close(f)
  
  nhej <- lns[grep("\tNHEJ:", lns)]
  total <- lns[grep("TOTAL:", lns)]
  counts <- as.numeric(gsub(".*:([0-9]+)\ .*", "\\1", c(nhej, total)))
  result <- counts[1]/counts[2]*100
  result
}

run_crispresso <- function(r1, name, amplicon, guide, results, r2 = NULL){
  # Run CRISPResso, remove unneeded files, rename folder
  
  single_cmd <- "CRISPResso -r1 %s -a %s -g %s -o %s"
  paired_cmd <- "CRISPResso -r1 %s -r2 %s -a %s -g %s -o %s --window_around_sgrna 5"
  if (is.null(r2)){
    x <- system(sprintf(single_cmd, r1, amplicon, guide, results),
           ignore.stdout = TRUE)
    out_dir <- paste("CRISPResso_on", basename(r1), sep = "_")
  } else {
    x <- system(sprintf(paired_cmd, r1, r2, amplicon, guide, results),
           ignore.stdout = TRUE)
    out_dir <- paste("CRISPResso_on", basename(r1), basename(r2), sep = "_") 
  }

  if (x == 2) return(NA)

  out_dir <- gsub(".fastq.gz", "", out_dir)
  out_dir <- file.path(results, out_dir)
  new_out_dir <- file.path(results, name)
  system(sprintf("mv %s %s", out_dir, new_out_dir))
    
  # cd to crispr dir, cd .., mv cripsr dir
  system(sprintf("rm %s/*vector*", new_out_dir))
  needle_out <- sprintf("%s/needle_output_SRR1769728_1_renamed_SRR1769728_2_renamed.txt.gz",
                     new_out_dir)
  system(sprintf("zcat %s | head -n 1000 > %s/needle_output.txt; rm %s",
   needle_out, new_out_dir, needle_out))
  result <- parse_results(new_out_dir)
  result
}


run_split_merged <- function(amplicons, guides, out){
  split_merged <- list.files("../split_merged", full.names = TRUE)
  result <- vector("numeric", length=length(split_merged))
  names(result) <- basename(split_merged)
  
  existing <- list.files(out)
  sm_base <- basename(split_merged)
  split_merged <- split_merged[! sm_base %in% existing]
  print(split_merged)
 
  for (dr in split_merged){
    r1 <- file.path(dr, "SRR1769728_1.fastq.gz")
    r2 <- file.path(dr, "SRR1769728_2.fastq.gz")
    nm <- basename(dr)
    amplicon <- amplicons[[nm]]
    guide <- as.character(guides[[nm]])
    result[[nm]] <- run_crispresso(r1 = r1,r2 = r2, name = nm,
                     amplicon = amplicon, guide = guide, results = out)
  }
  result 
}

run_merged_split <- function(amplicons, guides, out){
  merged_split <- list.files("../merged_split", full.names = TRUE)
  existing <- list.files(out)
  ms_base <- basename(merged_split)
  merged_split <- merged_split[! ms_base %in% existing]
  merged_split <- merged_split[!grepl("tolerant|strict", merged_split)]

  result <- vector("numeric", length=length(merged_split))
  names(result) <- basename(merged_split)
  for (dr in merged_split){
    r1 <- file.path(dr, "SRR1769728_merged_seqprep.fastq.gz")
    nm <- basename(dr)
    amplicon <- amplicons[[nm]]
    guide <- as.character(guides[[nm]])
    result[[nm]] <- run_crispresso(r1 = r1, name = nm,
                     amplicon = amplicon, guide = guide, results = out)
  }
  result
}

# Run CRISPResso with custom reference
custom <- run_merged_split(custom, guides, out = out_merged)

result <- run_split_merged(amplicons, guides, out = out_dir)
merged <- run_merged_split(amplicons, guides, out = out_merged)

# Run CRISPRessoPooled
system("run_crispresso_pooled.sh")
