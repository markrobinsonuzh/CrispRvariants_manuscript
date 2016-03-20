eibrary(httr)
library(rvest)


getNHEJ <- function(ref.seq, fq.file, ref = NA){
  # Submit the data
  dnaseq <- as.character(ref.seq)
  body <- list(uploaded_file = upload_file(fq.file),
               original = ref.seq)
  if (! is.na(ref)){
    body["target"] <- ref
  }               

  result <- POST("http://54.80.152.219/upload.php", body = body)
  
  # Get the results in xml format
  result_loc <- html_nodes(content(result), 'div a')[[4]] %>% html_attr("href")
  crisprga <- "http://54.80.152.219/"
  xml_result <- content(GET(paste0(crisprga, result_loc)))
  nhej <- as.numeric(html_nodes(xml_result, 'NHEJ') %>% html_text())
  print(nhej)
  return(nhej)
}


# Get the reference sequences
library(BSgenome.Drerio.UCSC.danRer7)
danRer7 <- BSgenome.Drerio.UCSC.danRer7

primers <- import("../annotation/shah_primers.bed")
start(primers) <- start(primers) + 1
sqs <- as.character(getSeq(danRer7, primers))
names(sqs) <- primers$name

# Run CRISPR-GA on the merging comparison data

directories <- c("../merged_split/tolerant_pear", "../merged_split/strict_pear")

for (dr in directories){
  # Get the matching fastqs
  fqs <- list.files(dr, full.names = TRUE)
  fqs <- fqs[match(primers$name, gsub(".fastq", "", basename(fqs)))]

  nhej <- vector("numeric", length=length(sqs))
  for (i in seq_along(sqs)){
    print(c(i, fqs[i]))
    eff <- getNHEJ(sqs[i], fqs[i])
    if (length(eff) == 0){
      nhej[i] <- NA
    } else {
      nhej[i] <- eff
    }
    print("sleeping")
    Sys.sleep(60)
  }
  results <- data.frame(name = primers$name, crisprga = nhej)
  write.table(results, sprintf("../results/CRISPRGA/%s_crisprga_refonly.txt", basename(dr)), sep = "\t", quote = FALSE)
}

# Run CRISPR-GA on the PEAR merged_split data

drs <- list.files("../merged_split", full.names = TRUE)
drs <- drs[!grepl("strict|tolerant", drs)]
fqs <- file.path(drs,"SRR1769728_merged_pear.assembled.fastq.gz")
fqs <- fqs[match(primers$name, basename(drs))]

nhej <- vector("numeric", length=length(sqs))
  for (i in seq_along(sqs)){
    print(c(i, fqs[i]))
    eff <- getNHEJ(sqs[i], fqs[i])
    if (length(eff) == 0){
      nhej[i] <- NA
    } else {
      nhej[i] <- eff
    }
    print("sleeping")
    Sys.sleep(60)
  }
results <- data.frame(name = primers$name, crisprga = nhej)
write.table(results, "../results/CRISPRGA/fwd_only_pear_crisprga_refonly.txt", sep = "\t", quote = FALSE)


# Analysis including guide sequence
# NOTE 7/3/16 - we were unable to repeat these results either through the 
# web form or by running through the website as the link to the xml file
# containing the results doesn't work

shah_results <- read.table("../annotation/Shah_metadata_edited.txt", sep = "\t", header = TRUE)
guide_nms <- gsub("\ ", "", shah_results$Gene)
guide_seqs <- gsub("\ ", "", shah_results$sgRNA)
refs <- guide_seqs[match(primers$name, guide_nms)]

nhej_w_ref <- vector("numeric", length=length(sqs))
for (i in seq_along(sqs)){
  print(c(i, fqs[i]))
  eff <- getNHEJ(sqs[i], fqs[i], refs[i])
  if (length(eff) == 0){
    nhej_w_ref[i] <- NA
  } else {
    nhej_w_ref[i] <- eff
  }
  print("sleeping")
  Sys.sleep(60)
}

write.table(nhej_w_ref, "../results/CRISPRGA/crisprga_w_guide_fwd_only_pear.txt", sep = "\t", quote = FALSE)

# Analysis of the seqprep originals

drs <- list.files("../merged_split", full.names = TRUE)
fqs <- file.path(drs,"SRR1769728_merged_seqprep.fastq.gz")
fqs <- fqs[match(primers$name, basename(drs))]

nhej <- vector("numeric", length=length(sqs))
  for (i in seq_along(sqs)){
    print(c(i, fqs[i]))
    eff <- getNHEJ(sqs[i], fqs[i])
    if (length(eff) == 0){
      nhej[i] <- NA
    } else {
      nhej[i] <- eff
    }
    print("sleeping")
    Sys.sleep(60)
  }
results <- data.frame(name = primers$name, crisprga = nhej)
write.table(results, "../results/CRISPRGA/seqprep_originals.txt", sep = "\t", quote = FALSE)




