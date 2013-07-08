### R code from vignette source 'RNASeq.Rnw'
library(Rsamtools)

#BAMdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/RawData/Time_course/Rep2/BAM"
#outdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/RawData/Time_course/Rep2/edgeR/"
#re <- "Rep2"

BAMdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/RawData/Time_course/Rep1/BAM"
outdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/RawData/Time_course/Rep1/edgeR/"
setwd(BAMdir)
re <- "Rep1"

## For Dosage and DN files ##
BAMdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/RawData/Suzanne_New_RNA/BAM"
outdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/RawData/Suzanne_New_RNA/edgeR/"

## Read BAM files ##
files <- tools::list_files_with_exts(BAMdir, "bam")
fls <- list.files(file.path(BAMdir), "_sorted.bam$", full=TRUE)
names(fls) <-sub("_sorted.*", "", basename(fls))
colNames <- unlist(strsplit(basename(fls), "\\_sorted.bam"))
bamFiles <- BamFileList(fls)



library(GenomicFeatures)
txdb <- makeTranscriptDbFromUCSC(genome = "mm9", tablename = "knownGene")
tx_by_gene=transcriptsBy(txdb,'gene')
ex_by_gene=exonsBy(txdb,'gene')


tx_by_gene_rg <- ranges(tx_by_gene)
ex_by_gene_rg <- ranges(ex_by_gene) ## u can pull out the exon size from here ##

temp_ex_size <- NULL
for(i in 1:length(ex_by_gene_rg)){
	temp_ex_size[i] <- sum(width(ex_by_gene_rg[i]))
}
#save(temp_ex_size, file=paste(outdir, "temp_ex_size_gene_level.rda", sep=""))
load(paste(outdir, "temp_ex_size_gene_level.rda", sep=""))

genome(ex_by_gene_rg) <- genome(ex_by_gene) <-  "mm9"

counter <- function(filePath, range){
	cat("File name is ",  filePath, "\n")
	aln <- readGappedAlignments(filePath)
	strand(aln) <- "*"
	cat("length of total reads is ", length(aln), "\n")
	genome(aln) <- "mm9"

	hits <- countOverlaps(aln, range)
#	temp <- countOverlaps(range, aln[hits==1])
	return(countOverlaps(range, aln[hits==1]))
	cat("length of unique hits is ", length(aln[hits==1]), "\n")
}

counts <- sapply(fls, counter, ex_by_gene)

## Counts is a matrix including the number of reads overlapped with exon regions ##
## The columns are data names and the row are exon IDs							 ##
#save(counts, file=paste(outdir, "total_data_exon_TimeCourse_", re, ".rda", sep=""))
save(counts, file=paste(outdir, "total_data_exon_Dosage_DN_", ".rda", sep=""))

#load(paste(outdir, "total_data_exon_TimeCourse_", re, ".rda", sep=""))	## call "counts" R object including the number of overlapped reads ##
load(paste(outdir, "temp_ex_size_gene_level.rda", sep=""))  ## call "temp_ex_size" object for exon-size for each gene ##
exon_size <- temp_ex_size
lib.size <- colSums(counts)

m <- 1e3*1e6 * t(t(counts/exon_size)/lib.size)
mm <- 1e7 * t(t(counts/exon_size)/lib.size)
#write.table(mm, file=paste(outdir, "total_data_exon_TimeCourse_", re, "_counts_ex_normalized_uniqueOverlap.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
#write.table(m, file=paste(outdir, "total_data_exon_TimeCourse_", re, "_counts_ex_RPKM_uniqueOverlap.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
write.table(mm, file=paste(outdir, "total_data_exon_Dosage_DN_", "_counts_ex_normalized_uniqueOverlap.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
write.table(m, file=paste(outdir, "total_data_exon_Dosage_DN_", "_counts_ex_RPKM_uniqueOverlap.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)

## Suzanne wants to look at log2-transformed RPKM data sets ##
#outdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/RawData/Time_course/Rep1/edgeR/"
#re <- "Rep1"

#outdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/RawData/Time_course/Rep2/edgeR/"
#re <- "Rep2"

m <- read.table(paste(outdir, "total_data_exon_Dosage_DN_", "_counts_ex_RPKM_uniqueOverlap.txt", sep=""), header=TRUE)
#m <- read.table(paste(outdir, "total_data_exon_TimeCourse_", re, "_counts_ex_RPKM_uniqueOverlap.txt", sep=""), header=TRUE)
temp <- log2(as.matrix(m))
temp1 <- ifelse(temp==-Inf, -99, temp)

#write.table(temp1, file=paste(outdir, "total_data_exon_TimeCourse_", re, "_counts_ex_Log2_RPKM_uniqueOverlap.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
write.table(temp1, file=paste(outdir, "total_data_exon_Dosage_DN", "_counts_ex_Log2_RPKM_uniqueOverlap.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)




















load(paste(outdir, "total_data_exon_TimeCourse_", re, ".rda", sep=""))
total_data_exon_IgG_DN2a <- total_data_exon
load(paste(outdir, "total_data_exon.rda", sep=""))
total_data_exon_D7 <- total_data_exon
total_data_exon  <- c(total_data_exon_D7, total_data_exon_IgG_DN2a)
save(total_data_exon, file=paste(outdir, "total_data_exon_merge.rda", sep=""))


load(paste(outdir, "total_data_trans_IgG_DN2a.rda", sep=""))
total_data_trans_IgG_DN2a <- total_data_trans
load(paste(outdir, "total_data_trans.rda", sep=""))
total_data_trans_D7 <- total_data_trans
total_data_trans  <- c(total_data_trans_D7, total_data_trans_IgG_DN2a)
save(total_data_trans, file=paste(outdir, "total_data_trans_merge.rda", sep=""))


total_data_exon_data <- data.frame(total_data_exon)
total_data_trans_data <- data.frame(total_data_trans)

save(total_data_exon_data, file=paste(outdir, "total_data_exon_data.rda", sep=""))
save(total_data_trans_data, file=paste(outdir, "total_data_trans_data.rda", sep=""))
#write.table(total_data_exon_data, file=paste(outdir, "total_data_exon_data.txt", sep=""), sep="\t")

load(paste(outdir, "total_data_exon_data.rda", sep="")) ## call "total_data_exon_data" object ##
lib.size <- colSums(total_data_exon_data)
m <- 1e7 * t(t(total_data_exon_data) / lib.size)
write.table(m, file=paste(outdir, "Jagged_lowDose_Dx_RNA_seq_counts_ex_normalized_03202013.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)

## Day 7 comparisons ##
coln <- c("D7_Dx_0.5_GGCTAC_L001", "D7_Dx_0.75_AGTCAA_L001_R2", "D7_Dx_1_AGTTCC_L002_R2")
coln <- c("D7_Dx_5_ATGTCA_L002_R2", "D7_Jx_20_CCGTCC_L003_R2")
temp <- total_data_exon_data[, coln]
grp <- coln
lib.size <- colSums(temp)
m <- 1e7 * t(t(temp) / lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)
temp <- temp[ridx,]
dif <- (log(temp[,grp[1]])-log(temp[,grp[2]]))/2
av <- (log(temp[,grp[1]])+log(temp[,grp[2]]))/2
pdf(file=paste(outdir,grp,"diff_ave_log.pdf", sep=""), width=12, height=12)
plot(dif~av, main=paste(grp[1], "vs", grp[2], sep=""), xlab="Ave", ylab="Diff")
dev.off()


library(biomaRt)
countdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/edgeR/"
load(paste(countdir, "gene_symbol.rda", sep=""))

###################################################
###################################################
###################################################

idx <- (rownames(total_data_exon_data) %in% gene_symbol$entrezgene)
total_data_exon_data_entrezegene <- total_data_exon_data[idx,]

#coln <- c("D7_Dx_5_ATGTCA_L002_R2", "D7_Jx_20_CCGTCC_L003_R2")
#coln <- c("D7_Dx_0.75_AGTCAA_L001_R2", "D7_Dx_5_ATGTCA_L002_R2")
#coln <- c("D7_Dx_0.75_AGTCAA_L001_R2", "D7_Jx_20_CCGTCC_L003_R2")

#coln <- c("D7_Dx_0.5_GGCTAC_L001", "D7_Dx_0.75_AGTCAA_L001_R2")
#coln <- c("D7_Dx_0.5_GGCTAC_L001", "D7_Dx_1_AGTTCC_L002_R2")

#coln <- c("D7_Dx_0.5_GGCTAC_L001", "D7_IgG_20_GTGAAA_L003")
#coln <- c("D7_Dx_0.75_AGTCAA_L001_R2", "D7_IgG_20_GTGAAA_L003")
#coln <- c("D7_Dx_5_ATGTCA_L002_R2", "D7_IgG_20_GTGAAA_L003")
coln <- c("D7_Jx_20_CCGTCC_L003_R2", "D7_IgG_20_GTGAAA_L003")




temp <- total_data_exon_data_entrezegene[, coln]

grp <- coln
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust(edgerRes_result$PValue, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)

print(table(idx))
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

#fileName <- paste("edgeR_exons_", coln[1], "_",coln[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", coln[1], "_",coln[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(outdir, fileName, sep=""))


edgerRes_result$keyGene <- ifelse(edgerRes_result$Gene %in% gene_list_RefSeq$Gene, 1, 0)
write.csv(edgerRes_result, file=paste(outdir, fileName, sep=""))
subset(edgerRes_result, keyGene==1)



#coln <- c("D7_Dx_0.5_GGCTAC_L001", "D7_IgG_20_GTGAAA_L003")
#coln <- c("D7_Dx_0.75_AGTCAA_L001_R2", "D7_IgG_20_GTGAAA_L003")


coln <- c("D7_Dx_0.75_AGTCAA_L001_R2", "D7_IgG_20_GTGAAA_L003")
fileName <- paste("edgeR_exons_", coln[1], "_",coln[2], ".csv", sep="")
edgerRes_result_075_IgG <- read.csv(paste(outdir, fileName, sep=""), header=TRUE)

coln <- c("D7_Jx_20_CCGTCC_L003_R2", "D7_IgG_20_GTGAAA_L003")
fileName <- paste("edgeR_exons_", coln[1], "_",coln[2], ".csv", sep="")
edgerRes_result_Jx20_IgG <- read.csv(paste(outdir, fileName, sep=""), header=TRUE)

coln <- c("D7_Dx_5_ATGTCA_L002_R2", "D7_IgG_20_GTGAAA_L003")
fileName <- paste("edgeR_exons_", coln[1], "_",coln[2], ".csv", sep="")
edgerRes_result_Dx5_IgG <- read.csv(paste(outdir, fileName, sep=""), header=TRUE)


edgerRes_result_075_IgG$logFC/edgerRes_result_Jx20_IgG$logFC
edgerRes_result_075_IgG$logFC/edgerRes_result_Dx5_IgG$logFC

edgerRes_result_Jx20_IgG$logFC/edgerRes_result_Dx5_IgG$logFC


fc <-seq(1.5, 10, by=0.5)
num_fc <- vector("list", length(fc))
names(num_fc) <- paste("FCcutoff", fc, sep="")

for(i in 1:length(fc)){
	edgerRes_result_075_IgG_1.5 <- edgerRes_result_075_IgG[edgerRes_result_075_IgG$logFC>fc[i],]
	edgerRes_result_Jx20_IgG_1.5 <- edgerRes_result_Jx20_IgG[edgerRes_result_Jx20_IgG$logFC>fc[i],]
	edgerRes_result_Dx5_IgG_1.5 <- edgerRes_result_Dx5_IgG[edgerRes_result_Dx5_IgG$logFC>fc[i],]

	num_fc[[i]] <- c(dim(edgerRes_result_Dx5_IgG_1.5)[1], dim(edgerRes_result_075_IgG_1.5)[1], dim(edgerRes_result_Jx20_IgG_1.5)[1])
}

tempp <- t(as.data.frame(num_fc))
colnames(tempp) <- c("Dx5", "Dx0.75", "Jx20")


pcut <-seq( 0.01, 10^(-5), by=-10^(-5)*5)

num_pcut <- vector("list", length(pcut))
names(num_pcut) <- paste("pcutoff", pcut, sep="")

for(i in 1:length(pcut)){
	edgerRes_result_075_IgG_1.5 <- edgerRes_result_075_IgG[edgerRes_result_075_IgG$"p.value.adj"<pcut[i],]
	edgerRes_result_Jx20_IgG_1.5 <- edgerRes_result_Jx20_IgG[edgerRes_result_Jx20_IgG$"p.value.adj"<pcut[i],]
	edgerRes_result_Dx5_IgG_1.5 <- edgerRes_result_Dx5_IgG[edgerRes_result_Dx5_IgG$"p.value.adj"<pcut[i],]
	
	num_pcut[[i]] <- c(dim(edgerRes_result_Dx5_IgG_1.5)[1], dim(edgerRes_result_075_IgG_1.5)[1], dim(edgerRes_result_Jx20_IgG_1.5)[1])
}


temppp <- t(as.data.frame(num_pcut))
colnames(temppp) <- c("Dx5", "Dx0.75", "Jx20")

edgerRes_result_075_IgG_1.5 <- edgerRes_result_075_IgG[edgerRes_result_075_IgG$logFC>fc,]
edgerRes_result_Jx20_IgG_1.5 <- edgerRes_result_Jx20_IgG[edgerRes_result_Jx20_IgG$logFC>fc,]
edgerRes_result_Dx5_IgG_1.5 <- edgerRes_result_Dx5_IgG[edgerRes_result_Dx5_IgG$logFC>fc,]




#coln <- c("D7_Dx_0.5_GGCTAC_L001", "D7_IgG_20_GTGAAA_L003")
#coln <- c("D7_Dx_0.75_AGTCAA_L001_R2", "D7_IgG_20_GTGAAA_L003")

coln <- c("D7_Dx_0.75_AGTCAA_L001_R2", "D7_Dx_5_ATGTCA_L002_R2")
#coln <- c("D7_Dx_0.75_AGTCAA_L001_R2", "D7_Jx_20_CCGTCC_L003_R2")
#coln <- c("D7_Jx_20_CCGTCC_L003_R2", "D7_Dx_5_ATGTCA_L002_R2")

temp <- total_data_exon_data_entrezegene[, coln]-total_data_exon_data_entrezegene[, "D7_IgG_20_GTGAAA_L003"]
temp[,1] <- ifelse(temp[,1]<0, 0, temp[,1])
temp[,2] <- ifelse(temp[,2]<0, 0, temp[,2])


coln <- c("D7_Dx_5_ATGTCA_L002_R2", "D7_Jx_20_CCGTCC_L003_R2")

fileName <- paste("edgeR_exons_", coln[1], "_",coln[2], ".csv", sep="")
edgerRes_result_Dx5_Jx20 <- read.csv(paste(outdir, fileName, sep=""), header=TRUE)
edgerRes_result_Dx5_Jx20 <- edgerRes_result_Dx5_Jx20[order(-edgerRes_result_Dx5_Jx20$logFC),]

coln <- c("D7_Dx_5_ATGTCA_L002_R2", "D7_Dx_0.75_AGTCAA_L001_R2")
fileName <- paste("edgeR_exons_", coln[1], "_",coln[2], ".csv", sep="")
edgerRes_result_Dx5_Dx075 <- read.csv(paste(outdir, fileName, sep=""), header=TRUE)
edgerRes_result_Dx5_Dx075 <- edgerRes_result_Dx5_Dx075[order(-edgerRes_result_Dx5_Dx075$logFC),]





grp <- coln
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust(edgerRes_result$PValue, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)

print(table(idx))
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

#fileName <- paste("edgeR_exons_", coln[1], "_",coln[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_IgG_adjustment", coln[1], "_",coln[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(outdir, fileName, sep=""))







reads=GRanges(seqnames=new_read_chr_names,ranges=IRanges(start=start(reads),end=end(reads)), 
strand=rep("*",length(reads)))
reads=GRanges(seqnames=new_read_chr_names,ranges=IRanges(start=start(reads),end=end(reads)), 
strand=strand(reads))
reads=GRanges(seqnames=rname(reads),ranges=IRanges(start=start(reads),end=end(reads)),
strand=rep("*",length(reads)))
   counts=countOverlaps(tx_by_gene,reads)
toc=data.frame(condition1_rep1=counts1.1,condition1_rep2=counts1.2,
condition2_rep1=counts2.1,condition2_rep2=counts2.2,stringsAsFactors=FALSE)
rownames(toc)=names(tx_by_gene)

library(edgeR)
norm_factors=calcNormFactors(as.matrix(toc))
   DGE=DGEList(toc,lib.size=norm_factors*colSums(toc),group=rep(c("Condition1","Condition2"),c(2,2)))
   disp=estimateCommonDisp(DGE)
   tested=exactTest(disp)



library(goseq)
genes = as.integer(p.adjust(tested$table$p.value, method = "BH") < 0.05)
names(genes) = row.names(tested$table)
pwf=nullp(genes,'mm9','ensGene')
   pwf=nullp(genes,bias.data=rowsum(counts[match(names(genes),rownames(counts))]))
   GO.pvals=goseq(pwf,'mm9','ensGene')

#The annotations have chromosomes called
names(seqlengths(tx_by_gene))
#The reads have chromosomes called
as.character(unique(rname(reads)))


#saveFeatures(txdb, paste(BAMdir, "my.dm3.ensGene.txdb.sqlite", sep=""))
#loadFeatures(paste(BAMdir, "my.dm3.ensGene.txdb.sqlite", sep=""))
head(getChromInfoFromUCSC("mm9"))

## Transcript regions from the first to last exon, contiguous ##
#library(GenomicFeatures)
#mm9_tx <- transcripts(mm9KG, columns="gene_id")
#gene_id <- values(mm9_tx)$gene_id
#all(elementLengths(gene_id) <= 1)	## no mapping or single mapping ##
#flat_gene_id <- character(length(mm9_tx))
#flat_gene_id[elementLengths(gene_id) == 1] <- unlist(gene_id)
#values(mm9_tx)$gene_id <- flat_gene_id

###################################################
## DF peaks b/w WT and N1n on each gene region ##
## promotor regions ##
#promoters <- flank(mm9_tx, 1000, both=TRUE)
#peakSummary$inPromoter <- peakSummary %in% promoters
#xtabs(~ inPromoter + change, peakSummary)

## Upstream or in a gene ##
#peakSummary$inUpstream <- peakSummary %in% flank(mm9_tx, 20000)
#peakSummary$inGene <- peakSummary %in% mm9_tx


## QC: quality control ##

aln <- readGappedAlignments(fls[1], param=ScanBamParam(what="mapq"))
aln <- aln[values(aln)[["mapq"]] > 10]
#strand(aln) <- "*"
alnGR <- granges(aln)
alnRD <- RangedData(IRanges(start=start(alnGR), end=end(alnGR)), space=seqnames(alnGR), strand=strand(alnGR))
alnData <- data.frame(space=space(alnRD), start=start(alnRD), end=end(alnRD), strand=strand(alnRD))
uniqueData <- unique(alnData)
dupData <- duplicated(alnData)
dupData1 <- duplicated(alnData[,-4])

hits <- countOverlaps(alnGR, ex)
table(hits)
cnt <- countOverlaps(range, aln[hits==1])

nrepeat <- as.vector(table(alnData[dupData=="TRUE",]))

nrepeats <- tabDuplReads(alnRD)

## Generate count matrix ##

###################################################
### code chunk number 29: bam-count-fun
###################################################
counter <- 
function(filePath, range)
{
#    aln <- readGappedAlignments(filePath)
	aln <- readGappedAlignments(filePath, param=ScanBamParam(what="mapq"))
	aln <- aln[values(aln)[["mapq"]] > 10]
	alnGR <- granges(aln)

#    strand(aln) <- "*"
#    hits <- countOverlaps(aln, ex)
#    cnt <- countOverlaps(range, aln[hits==1])

	uniqueGR <- unique(alnGR)

    strand(uniqueGR) <- "*"
    hits <- countOverlaps(uniqueGR, ex)
table(hits)
	cnt <- countOverlaps(range, uniqueGR[hits==1])

    names(cnt) <- names(range)
    cnt
}

for(i in 1:length(fls)){
	print(fls[i])
	p1 <- ScanBamParam(what=c("rname", "strand", "pos", "mapq", "qwidth"), flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE))
	rawData <- scanBam(fls[i], param=p1)[[1]]
	retained <-  rawData$mapq > 10  # filter out reads with bad quality scores

#	aln <- readGappedAlignments(fls[i], param=ScanBamParam(what="mapq"))
#	cat(length(aln))
#	aln1 <- aln[values(aln)[["mapq"]] > 10]
print(	length(rawData$mapq))
	cat("\n")
print(table(retained))
	addon <- rawData$qwidth - 1
	IP <- data.frame(space=rawData$rname[retained], start=(rawData$pos)[retained], end=(rawData$pos+addon)[retained], strand=rawData$strand[retained])
	IP$space <- as.character(IP$space) #change strand into characters
	IP <- unique(IP)
	print(dim(IP))
}
addon <- rawData$qwidth - 1
#addon[rawData$strand!="-"] <- 0  # to get end of "-" reads, seq length need to be added to their positions
retained <-  rawData$mapq > 10  # filter out reads with bad quality scores

## retrieve Known genes track of established exons, transcripts, and coding sequences of mouse ##
###################################################
### code chunk number 18: txdb
###################################################
#txdbFile <- list.files(bigdata(), "sqlite", full=TRUE)
#txdb <- loadFeatures(txdbFile)

## exons ##
exID0 <- exons(txdb, columns=c("exon_id","gene_id"))
exID <- reduce(exID0)
counts_exid <- sapply(fls, counter, exID)		## ex: gene exon ranges ##
save(counts_exid, file=paste(outdir, "counts_exons_mapq_10.rda", sep=""))


## Group exons into gene level ##
ex0 <- exonsBy(txdb, "gene")
head(table(elementLengths(ex0)))
ex <- reduce(ex0)
counts_ex <- sapply(fls, counter, ex)		## ex: gene exon ranges ##

save(counts_ex, file=paste(outdir, "counts_ex_mapq_10_unique.rda", sep=""))

## transcripts ##
tx0 <- transcriptsBy(txdb, "gene")
head(table(elementLengths(tx0)))
tx <- reduce(tx0)


###################################################
### code chunk number 30: bam-count-all
###################################################
counts_ex <- sapply(fls, counter, ex)		## ex: gene exon ranges ##

save(counts_ex, file=paste(outdir, "counts_ex_mapq_10.rda", sep=""))

outdir <-"/Users/sirius/Documents/POST_Doctoral_Fellow/Hamid_Notch/Suzanne_RNA_Seq/Counts_Data/"
write.table(counts_ex, file=paste(outdir, "counts_ex_mapq_10.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)

lib.size <- colSums(counts_ex)
m <- 1e7 * t(t(counts_ex) / lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)
temp <- m[ridx,]
write.table(m, file=paste(outdir, "counts_ex_mapq_10_normalized.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)
write.table(temp, file=paste(outdir, "counts_ex_mapq_10_normalized_filtered.txt", sep=""), sep="\t", quote=FALSE, row.names=TRUE)


#counts_tx <- sapply(fls, counter, tx)		## tx: gene transcripts ranges ##


library(biomaRt)

myENSEMBLids <- rownames(counts_ex)
myEnsembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
#listAttributes(myEnsembl)
gene_symbol <- getBM(attributes=c("entrezgene","mgi_symbol"), filters="entrezgene", values=myENSEMBLids, mart=myEnsembl, uniqueRows=TRUE)
save(gene_symbol, file=paste(countdir, "gene_symbol.rda"))
write.csv(gene_symbol, file=paste(countdir, "entrezgeneID_GeneSymbol.csv", sep=""), row.names=FALSE) 


###################################################
### code chunk number 1: setup
###################################################
BAMdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/BAM/"
countdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/edgeR/"
QCdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/QC/"
options(digits=2)
#library(SeattleIntro2011)
library(edgeR)
library(goseq)


###################################################
### code chunk number 2: counts
###################################################
load(paste(countdir,"counts_ex_mapq_10.rda", sep=""))
load(paste(countdir,"counts_ex_mapq_10.rda", sep=""))
#data(counts_ex)
dim(counts_ex)
#grp <- factor(sub("[1-4].*", "", colnames(counts_ex)), levels=c("untreated", "treated"))
colnames(counts_ex)

#tr_name <- rownames(counts_ex)
#tg_name <- gene_symbol[tr_name==gene_symbol$entrezgene, ]

## Day 7 comparisons ##
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d7", "RNA_Seq_N1n_Dx_d7")]
grp <- c("RNA_Seq_WT_Dx_d7", "RNA_Seq_N1n_Dx_d7")
lib.size <- colSums(temp)
m <- 1e7 * t(t(temp) / lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)
temp <- temp[ridx,]
dif <- (log(temp[,grp[1]])-log(temp[,grp[2]]))/2
av <- (log(temp[,grp[1]])+log(temp[,grp[2]]))/2
pdf(file=paste(QCdir,grp[1],"_",grp[2],"diff_ave_log.pdf", sep=""), width=12, height=12)
plot(dif~av, main=paste(grp[1], "vs", grp[2], sep=""), xlab="Ave", ylab="Diff")
dev.off()

temp <- counts_ex[, c("RNA_Seq_WT_Dx_d7", "RNA_Seq_WT_IgG_d7")]
grp <- c("RNA_Seq_WT_Dx_d7", "RNA_Seq_WT_IgG_d7")
lib.size <- colSums(temp)
m <- 1e7 * t(t(temp) / lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)
temp <- temp[ridx,]
dif <- (log(temp[,grp[1]])-log(temp[,grp[2]]))/2
av <- (log(temp[,grp[1]])+log(temp[,grp[2]]))/2
pdf(file=paste(QCdir,grp[1],"_",grp[2],"diff_ave_log.pdf", sep=""), width=12, height=12)
plot(dif~av, main=paste(grp[1], "vs", grp[2], sep=""), xlab="Ave", ylab="Diff")
dev.off()

temp <- counts_ex[, c("RNA_Seq_N1n_Dx_d7", "RNA_Seq_WT_IgG_d7")]
grp <- c("RNA_Seq_N1n_Dx_d7", "RNA_Seq_WT_IgG_d7")
lib.size <- colSums(temp)
m <- 1e7 * t(t(temp) / lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)
temp <- temp[ridx,]
dif <- (log(temp[,grp[1]])-log(temp[,grp[2]]))/2
av <- (log(temp[,grp[1]])+log(temp[,grp[2]]))/2
pdf(file=paste(QCdir,grp[1],"_",grp[2],"diff_ave_log.pdf", sep=""), width=12, height=12)
plot(dif~av, main=paste(grp[1], "vs", grp[2], sep=""), xlab="Ave", ylab="Diff")
dev.off()

## Day 14 comparisons ##
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d14", "RNA_Seq_N1n_Dx_d14")]
grp <- c("RNA_Seq_WT_Dx_d14", "RNA_Seq_N1n_Dx_d14")
lib.size <- colSums(temp)
m <- 1e7 * t(t(temp) / lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)
temp <- temp[ridx,]
dif <- (log(temp[,grp[1]])-log(temp[,grp[2]]))/2
av <- (log(temp[,grp[1]])+log(temp[,grp[2]]))/2
pdf(file=paste(QCdir,grp[1],"_",grp[2],"diff_ave_log.pdf", sep=""), width=12, height=12)
plot(dif~av, main=paste(grp[1], "vs", grp[2], sep=""), xlab="Ave", ylab="Diff")
dev.off()

temp <- counts_ex[, c("RNA_Seq_WT_Dx_d14", "RNA_Seq_WT_IgG_d14")]
grp <- c("RNA_Seq_WT_Dx_d14", "RNA_Seq_WT_IgG_d14")
lib.size <- colSums(temp)
m <- 1e7 * t(t(temp) / lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)
temp <- temp[ridx,]
dif <- (log(temp[,grp[1]])-log(temp[,grp[2]]))/2
av <- (log(temp[,grp[1]])+log(temp[,grp[2]]))/2
pdf(file=paste(QCdir,grp[1],"_",grp[2],"diff_ave_log.pdf", sep=""), width=12, height=12)
plot(dif~av, main=paste(grp[1], "vs", grp[2], sep=""), xlab="Ave", ylab="Diff")
dev.off()

temp <- counts_ex[, c("RNA_Seq_N1n_Dx_d14", "RNA_Seq_WT_IgG_d14")]
grp <- c("RNA_Seq_N1n_Dx_d14", "RNA_Seq_WT_IgG_d14")
lib.size <- colSums(temp)
m <- 1e7 * t(t(temp) / lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)
temp <- temp[ridx,]
dif <- (log(temp[,grp[1]])-log(temp[,grp[2]]))/2
av <- (log(temp[,grp[1]])+log(temp[,grp[2]]))/2
pdf(file=paste(QCdir,grp[1],"_",grp[2],"diff_ave_log.pdf", sep=""), width=12, height=12)
plot(dif~av, main=paste(grp[1], "vs", grp[2], sep=""), xlab="Ave", ylab="Diff")
dev.off()

## Day 21 comparisons ##
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d21", "RNA_Seq_N1n_Dx_d21")]
grp <- c("RNA_Seq_WT_Dx_d21", "RNA_Seq_N1n_Dx_d21")
lib.size <- colSums(temp)
m <- 1e7 * t(t(temp) / lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)
temp <- temp[ridx,]
dif <- (log(temp[,grp[1]])-log(temp[,grp[2]]))/2
av <- (log(temp[,grp[1]])+log(temp[,grp[2]]))/2
pdf(file=paste(QCdir,grp[1],"_",grp[2],"diff_ave_log.pdf", sep=""), width=12, height=12)
plot(dif~av, main=paste(grp[1], "vs", grp[2], sep=""), xlab="Ave", ylab="Diff")
dev.off()

temp <- counts_ex[, c("RNA_Seq_WT_Dx_d21", "RNA_Seq_WT_IgG_d21")]
grp <- c("RNA_Seq_WT_Dx_d21", "RNA_Seq_WT_IgG_d21")
lib.size <- colSums(temp)
m <- 1e7 * t(t(temp) / lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)
temp <- temp[ridx,]
dif <- (log(temp[,grp[1]])-log(temp[,grp[2]]))/2
av <- (log(temp[,grp[1]])+log(temp[,grp[2]]))/2
pdf(file=paste(QCdir,grp[1],"_",grp[2],"diff_ave_log.pdf", sep=""), width=12, height=12)
plot(dif~av, main=paste(grp[1], "vs", grp[2], sep=""), xlab="Ave", ylab="Diff")
dev.off()

temp <- counts_ex[, c("RNA_Seq_N1n_Dx_d21", "RNA_Seq_WT_IgG_d21")]
grp <- c("RNA_Seq_N1n_Dx_d21", "RNA_Seq_WT_IgG_d21")
lib.size <- colSums(temp)
m <- 1e7 * t(t(temp) / lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)
temp <- temp[ridx,]
dif <- (log(temp[,grp[1]])-log(temp[,grp[2]]))/2
av <- (log(temp[,grp[1]])+log(temp[,grp[2]]))/2
pdf(file=paste(QCdir,grp[1],"_",grp[2],"diff_ave_log.pdf", sep=""), width=12, height=12)
plot(dif~av, main=paste(grp[1], "vs", grp[2], sep=""), xlab="Ave", ylab="Diff")
dev.off()

###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
ridx <- rowSums(m > 1) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]


###################################################
### code chunk number 5: mds
###################################################
#plotMDS.DGEList(dgl)


countdir <- "/shared/silo_researcher/Gottardo_R/Sangsoon_working/Hamid_Notch/Suzanne_RNA_Seq/edgeR/"
load(paste(countdir,"counts_ex_mapq_10_unique.rda", sep=""))


###################################################
###################################################
###################################################
load(paste(countdir, "gene_symbol.rda", sep=""))
## Day 7 comparisons ##
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d7", "RNA_Seq_N1n_Dx_d7")]
grp <- c("RNA_Seq_WT_Dx_d7", "RNA_Seq_N1n_Dx_d7")

###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################
dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

#fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d7", "RNA_Seq_WT_IgG_d7")]
#cut_ind <- rowSums(temp1) > 10 
#temp <- temp1[cut_ind,]
grp <- c("RNA_Seq_WT_Dx_d7", "RNA_Seq_WT_IgG_d7")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
print(table(idx))
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

#fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_N1n_Dx_d7", "RNA_Seq_WT_IgG_d7")]
grp <- c("RNA_Seq_N1n_Dx_d7", "RNA_Seq_WT_IgG_d7")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)

edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

#fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))

###################################################
###################################################
###################################################

## Day 14 comparisons ##
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d14", "RNA_Seq_N1n_Dx_d14")]
grp <- c("RNA_Seq_WT_Dx_d14", "RNA_Seq_N1n_Dx_d14")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))




###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d14", "RNA_Seq_WT_IgG_d14")]
grp <- c("RNA_Seq_WT_Dx_d14", "RNA_Seq_WT_IgG_d14")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_N1n_Dx_d14", "RNA_Seq_WT_IgG_d14")]
grp <- c("RNA_Seq_N1n_Dx_d14", "RNA_Seq_WT_IgG_d14")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)

edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
## Day 21 comparisons ##
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d21", "RNA_Seq_N1n_Dx_d21")]
grp <- c("RNA_Seq_WT_Dx_d21", "RNA_Seq_N1n_Dx_d21")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d21", "RNA_Seq_WT_IgG_d21")]
grp <- c("RNA_Seq_WT_Dx_d21", "RNA_Seq_WT_IgG_d21")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))



###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_N1n_Dx_d21", "RNA_Seq_WT_IgG_d21")]
grp <- c("RNA_Seq_N1n_Dx_d21", "RNA_Seq_WT_IgG_d21")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################

## Day 7 vs Day 14 comparisons ##
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d7", "RNA_Seq_WT_Dx_d14")]
grp <- c("RNA_Seq_WT_Dx_d7", "RNA_Seq_WT_Dx_d14")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_WT_IgG_d7", "RNA_Seq_WT_IgG_d14")]
grp <- c("RNA_Seq_WT_IgG_d7", "RNA_Seq_WT_IgG_d14")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_N1n_Dx_d7", "RNA_Seq_N1n_Dx_d14")]
grp <- c("RNA_Seq_N1n_Dx_d7", "RNA_Seq_N1n_Dx_d14")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
## Day 14 vs Day 21 comparisons ##
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d14", "RNA_Seq_WT_Dx_d21")]
grp <- c("RNA_Seq_WT_Dx_d14", "RNA_Seq_WT_Dx_d21")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_WT_IgG_d14", "RNA_Seq_WT_IgG_d21")]
grp <- c("RNA_Seq_WT_IgG_d14", "RNA_Seq_WT_IgG_d21")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_N1n_Dx_d14", "RNA_Seq_N1n_Dx_d21")]
grp <- c("RNA_Seq_N1n_Dx_d14", "RNA_Seq_N1n_Dx_d21")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
## Day 7 vs Day 21 comparisons ##
temp <- counts_ex[, c("RNA_Seq_WT_Dx_d7", "RNA_Seq_WT_Dx_d21")]
grp <- c("RNA_Seq_WT_Dx_d7", "RNA_Seq_WT_Dx_d21")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))


###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_WT_IgG_d7", "RNA_Seq_WT_IgG_d21")]
grp <- c("RNA_Seq_WT_IgG_d7", "RNA_Seq_WT_IgG_d21")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))



###################################################
###################################################
###################################################
temp <- counts_ex[, c("RNA_Seq_N1n_Dx_d7", "RNA_Seq_N1n_Dx_d21")]
grp <- c("RNA_Seq_N1n_Dx_d7", "RNA_Seq_N1n_Dx_d21")
###################################################
### code chunk number 3: DGEList
###################################################
library(edgeR)

dgl <- DGEList(counts=temp, group=grp)
dgl$samples$lib.size <- colSums(dgl$counts)

dgl <- calcNormFactors(dgl, , method="TMM")

###################################################
### code chunk number 4: DEGList-filter
###################################################
m <- 1e7 * t(t(dgl$counts) / dgl$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
ridx <- rowSums(m) >= 2
table(ridx)   # number filtered / retained $
dgl <- dgl[ridx,]

###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
table(idx)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], ".csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))



###################################################
### code chunk number 6: design
###################################################
(design <- model.matrix( ~ grp ))


###################################################
### code chunk number 7: common.dispersion
###################################################

dgl <- estimateCommonDisp(dgl, design)
dgl$common.dispersion <- 0.1
sqrt(dgl$common.dispersion)

edgerRes_totalredas <- exactTest( dgl )
edgerRes_result <- edgerRes_totalredas$table
edgerRes_result$BH_Padj <- p.adjust( edgerRes_totalredas$table$p.value, method="BH" )
colnames(edgerRes_result) <- c("logConc", "logFC", "p.value", "p.value.adj")
#edgerRes_result$id <- RefSeqCounts_TSS$Gene[which_gene=="TRUE"]
edgerRes_result<- edgerRes_result[order(edgerRes_result$p.value.adj),]
idx <- (rownames(edgerRes_result) %in% gene_symbol$entrezgene)
edgerRes_result <- edgerRes_result[idx,]
for(i in 1:dim(edgerRes_result)[1]){
	edgerRes_result$Gene[i] <- gene_symbol[rownames(edgerRes_result[i,])==gene_symbol$entrezgene,2]
}

fileName <- paste("edgeR_exons_", grp[1], "_",grp[2], "_uniqueR.csv", sep="")
write.csv(edgerRes_result, file=paste(countdir, fileName, sep=""))




#> names(gene_symbol)
#[1] "entrezgene" "mgi_symbol"
#> gene_symbol[gene_symbol$mgi_symbol=="Hes1",]
#entrezgene mgi_symbol
#4742      15205       Hes1
#> gene_symbol[gene_symbol$mgi_symbol=="Il2ra",]
#entrezgene mgi_symbol
#11760      16184      Il2ra
#> gene_symbol[gene_symbol$mgi_symbol=="Gata3",]
#entrezgene mgi_symbol
#9393      14462      Gata3
#> gene_symbol[gene_symbol$mgi_symbol=="Notch1",]
#entrezgene mgi_symbol
#18896      18128     Notch1
#> gene_symbol[gene_symbol$mgi_symbol=="Nrarp",]
#entrezgene mgi_symbol
#10048      67122      Nrarp


###################################################
### code chunk number 8: glmFit
###################################################
fit <- glmFit(dgl, design, dispersion=dgl$common.dispersion)


###################################################
### code chunk number 9: lrt
###################################################
lrTest <- glmLRT(dgl, fit, coef=2)


###################################################
### code chunk number 10: topTags
###################################################
(tt <- topTags(lrTest))


###################################################
### code chunk number 11: sanity
###################################################
sapply(rownames(tt$table)[1:4], 
       function(x) tapply(counts_ex[x,], grp, mean))


