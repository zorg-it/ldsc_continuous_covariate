#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(GenomicRanges))

args <- commandArgs(trailingOnly = TRUE)
bed_file   <- args[1]
bim_file   <- args[2]
annot_file <- args[3]
mode       <- ifelse(length(args)>=4, args[4], "full-annot")

# Leggi BED continua
annot_bed <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
colnames(annot_bed) <- c("chr", "start", "end", "value")
annot_bed$chr <- gsub("^chr", "", annot_bed$chr)
annot_bed$start <- annot_bed$start + 1

# Leggi BIM
bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

# Prepara GRanges
bim_gr <- GRanges(seqnames = bim$CHR, ranges = IRanges(start=bim$BP, end=bim$BP))
bed_gr <- GRanges(seqnames = annot_bed$chr, ranges = IRanges(start=annot_bed$start, end=annot_bed$end),
                  value = annot_bed$value)

# Sovrapposizione
ov <- findOverlaps(query = bim_gr, subject = bed_gr)
annot_values <- rep(0, nrow(bim))
annot_values[queryHits(ov)] <- annot_bed$value[subjectHits(ov)]

# Costruisci tabella
out_df <- data.frame(
  CHR  = bim$CHR,
  BP   = bim$BP,
  SNP  = bim$SNP,
  CM   = bim$CM,
  ANNOT = annot_values,
  stringsAsFactors = FALSE
)

# Scrittura file
write.table(out_df, gzfile(annot_file), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

message("✅ Continuous annotation saved: ", annot_file)
