
##SB sequencing 

library(tidyverse)
library(dada2)
library(here)
library(phyloseq)

##Create Folders
fastq = "fastq"        # raw fastq files
filtered = "filtered"  # dada2 trimmed fastq files
output = "output"      # output files
images = "images"       # output images


if(!dir.exists(filtered)) dir.create(filtered)
if(!dir.exists(output)) dir.create(output)
if(!dir.exists(images)) dir.create(images)

##Create List of File Names
fns = sort(list.files(fastq, full.names = TRUE))
fnFwds = fns[grep("R1.fastq", fns)]
fnRevs = fns[grep("R2.fastq", fns)]

sample_names = gsub("_R1.fastq", "", basename(fnFwds))


##Visualize Read Quality
ii = 1:length(sample_names)
pdf(paste0(images, "/plotQualityProfile.pdf"), width = 8, height = 8, pointsize = 12)
for(i in ii) {
  message(paste0("[", i ,"/", length(sample_names), "] ", sample_names[i]))
  print(plotQualityProfile(fnFwds[i]) + ggtitle("Fwd"))
  print(plotQualityProfile(fnRevs[i]) + ggtitle("Rev"))
}
invisible(dev.off())


##Filter and Trim
filtFwds = file.path(filtered, basename(fnFwds))
filtRevs = file.path(filtered, basename(fnRevs))

out = filterAndTrim(fnFwds, filtFwds, fnRevs, filtRevs,
                          truncLen = c(145,140), minLen = 100, maxN = 0, truncQ = 2, maxEE =c(2,2),
                          rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = FALSE)

out = as.data.frame(out)
rownames(out) = sample_names
head(out, 10)


##Learn the Error Rates
errFwd = learnErrors(filtFwds, multithread = FALSE)
errRev = learnErrors(filtRevs, multithread = FALSE)


pdf(paste0(images, "/plotErrors.pdf"), width = 10, height = 10, pointsize = 12)
plotErrors(errFwd, nominalQ = TRUE)
plotErrors(errRev, nominalQ = TRUE)
invisible(dev.off())


##Sample Inference
dadaFwds = dada(filtFwds, err = errFwd, pool = FALSE, multithread = FALSE)
dadaRevs = dada(filtRevs, err = errRev, pool = FALSE, multithread = FALSE)

##Merge Paired Reads
merged = mergePairs(dadaFwds, filtFwds, dadaRevs, filtRevs, verbose = TRUE)

##Construct the Sequence Table
seqtab <- makeSequenceTable(merged)

table(nchar(getSequences(seqtab)))

##Remove Chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=2, verbose=TRUE)
##Percent of sequences that are non-chimeric 
sum(seqtab_nochim)/sum(seqtab)

rownames(seqtab_nochim) <- sample_names

##Track Reads through the Pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFwds, getN), sapply(dadaRevs, getN), sapply(merged, getN), rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged","nonchim")
rownames(track) <- sample_names
head(track)


##Assign Taxonomy
seqs = getSequences(seqtab_nochim)
dbpath = "16S_DB/"
ref_db = paste0(dbpath, "silva_nr99_v138.1_train_set.fa.gz")

taxonomy_tab <- assignTaxonomy(seqs, refFasta = ref_db)

##Create Phyloseq Object
phyloseq_unfiltered <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE), 
                                     tax_table(taxonomy_tab), sample_data(metadata))


save(phyloseq_unfiltered, file="phyloseq_unfiltered.RData")



