library(GenomicFeatures)
txdb <- makeTranscriptDbFromUCSC(genome = "mm9", tablename = "knownGene")
ex0 <- exonsBy(txdb, "gene")
head(table(elementLengths(ex0)))
ex <- reduce(ex0)
BAMdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/BAM/"
outdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/edgeR/"
setwd(BAMdir)
## Read BAM files ##
fls <- list.files(file.path(BAMdir), "_sorted.bam$", full=TRUE)
names(fls) <-sub("_sorted.*", "", basename(fls))

for(i in 1:length(fls)){
	print(fls[i])
	aln <- readGappedAlignments(fls[i], param=ScanBamParam(what="mapq"))
	aln <- aln[values(aln)[["mapq"]] > 10]
	print(length(aln))
	print(table(width(aln)))
	hits <- countOverlaps(aln, ex)
	print(table(hits))
	
	alnGR <- granges(aln)
	print(length(alnGR))

	uniqueGR <- unique(alnGR)
	print(length(uniqueGR))
	strand(uniqueGR) <- "*"
	
	hits <- countOverlaps(uniqueGR, ex)
	print(table(hits))
#	cnt <- countOverlaps(range, aln[hits==1])
}

