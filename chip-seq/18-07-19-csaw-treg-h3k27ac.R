library(csaw)
library(edgeR)
library(statmod)
bam.files <- c("/media/serlia/Storage\ 2/rora-treg-chipseq/mapping/srtd-fltrd/treg/INPUT-WT_S9.sorted.mrkdup.bam",
               "/media/serlia/Storage\ 2/rora-treg-chipseq/mapping/srtd-fltrd/treg/INPUT-KO_S10.sorted.mrkdup.bam",
               "/media/serlia/Storage\ 2/rora-treg-chipseq/mapping/srtd-fltrd/treg/WT-H-1_S1.sorted.mrkdup.bam",
               "/media/serlia/Storage\ 2/rora-treg-chipseq/mapping/srtd-fltrd/treg/WT-H-2_S2.sorted.mrkdup.bam",
               "/media/serlia/Storage\ 2/rora-treg-chipseq/mapping/srtd-fltrd/treg/WT-H-3_S3.sorted.mrkdup.bam",
               "/media/serlia/Storage\ 2/rora-treg-chipseq/mapping/srtd-fltrd/treg/KO-H-1_S4.sorted.mrkdup.bam",
               "/media/serlia/Storage\ 2/rora-treg-chipseq/mapping/srtd-fltrd/treg/KO-H-2_S5.sorted.mrkdup.bam",
               "/media/serlia/Storage\ 2/rora-treg-chipseq/mapping/srtd-fltrd/treg/KO-H-3_S6.sorted.mrkdup.bam")


# depends on the aligner? Check for bowtie 
# Also better to set dedup = TRUE for readParam
# at this stage

param <- readParam(minq = 30, dedup = TRUE)

# frag.length determined on 11-07-19 was with input samples,
# will recalculate without input, plot saved

bam.files.noinput <- bam.files[3:8]
bam.files.wt.input <- c(bam.files[1], bam.files[3:8])
# Was unable to subset Ranged
max.delay <- 500
x <- correlateReads(bam.files.noinput,
                    max.delay, param = param)
plot(0:max.delay, 
     x,
     type="l", ylab="CCF", xlab="Delay (bp)")
# Determine fragment length with highest  
max.CCF <- max(x)
frag.length <- which(x %in% max.CCF)
# Ah, there is a function for that
frag.length <- maximizeCcf(x)
frag.length <- 186
window.width <- 150

data.no.input <- windowCounts(bam.files.noinput,
                              ext = frag.length,
                              width = window.width,
                              param = param)

data.with.input <- windowCounts(bam.files.wt.input,
                                       ext=frag.length,
                                       width=window.width,
                                       param=param)
chip <- data.with.input[,2:7]
control <- data.with.input[,1]

# According to the author, input controls are irrelevant
# to use in DB analysis. But they might be used to filter
# windows with average abundance less than in the input control
# Binning with large windows is used to determine low abundance
# regions.

window.width.binning <- 10000
data.with.input.binned <- windowCounts(bam.files.wt.input,
                                       bin=TRUE,
                                       ext=frag.length,
                                       width=window.width.binning,
                                       param=param)
chip.binned <- data.with.input.binned[,2:7]

control.binned <- data.with.input.binned[,1]

scale.info <- scaleControlFilter(chip.binned, control.binned)
filter.stat <- filterWindows(chip, control, type = "control",
                                prior.count = 5, scale.info = scale.info)
keep <- filter.stat$filter > log2(15)

sum(keep)

data.no.input.binned <- windowCounts(bam.files.noinput, bin = TRUE,
                                     width = window.width.binning,
                                     param = param)

# Constructing filtered.data object
filtered.data <- chip[keep,]
filtered.data1 <- chip[keep1,]
# By now I do not understand how to deal with wt and ko inputs.
# I use wt input for filtering of all the data, I will see if 
# it makes sense later

# Now I proceed to normalization
#
# Normalization aims to generate normalization coefficient
# for each dataset. These might be generated based on either
# 1) Composition biases (uneven distribution of sequencing resources)
# or 2) Efficiency biases (Variable efficiencies of IPs)
# Using either 1) or 2) is mutually exclusive.
# Eliminating 1) is useful for designes where binary (1/0)-type differences
# are expected between experimental groups (WT-KO -- ChIP for mutated protein).
# Eliminating 2) is useful when differences expected are quite subtle
# as with histone marks. 
# So for H3K27ac of Treg I will use elimination of efficiency biases
#
# Eliminating efficiency biases
# filtered dataset is used (input control-based normalization performed above)

filtered.data <- normFactors(filtered.data, se.out = TRUE)
data.eff <- filtered.data$norm.factors
data.eff
# Check with MA plots
# Data is ploted on non-filtered reads with binned regions to include
# background in the plot

chip.comp <- normFactors(chip.binned, se.out = FALSE)
chip.comp
# If I use data.no.input.binned dataset
# > chip.comp
# [1] 1.0151403 0.9922660 0.9666523 1.0108136 0.9984452 1.0176073
# If I use chip.binned dataset
# > chip.comp
# [1] 1.0151403 0.9922660 0.9666523 1.0108136 0.9984452 1.0176073
# So it's okay -- it's the same
par(mfrow=c(2, 3))
bins <- data.no.input.binned
comp <- chip.comp
eff <- data.eff1
adjc <- cpm(asDGEList(bins), log = TRUE)
for (i in c(2:6)) {
  smoothScatter(x=rowMeans(adjc), y=adjc[,1]-adjc[,i],
                xlab="A", ylab="M", main=paste("1 vs", i))
  abline(h=log2(eff[1]/eff[i]), col="red")
  abline(h=log2(comp[1]/comp[i]), col="red", lty = 2)
}
# Differences for potentially DB cloud are obviously not huge
# between the samples. So I will probably proceed with
# efficiency bias elimination.
# I will also try what loess normalization gives
# Will use less stringent filtering (log2(4) intead of log2(15))
keep.loess <- filter.stat$filter > log2(4)
filtered.data.loess <- chip[keep.loess,]
filtered.data.loess <- normOffsets(filtered.data.loess,
                                   se.out = filtered.data.loess)
loess.offsets <- assay(filtered.data.loess,
                       "offset")
head(loess.offsets)
# Filtering log2(4)
# > head(loess.offsets)
# [,1]       [,2]        [,3]      [,4]      [,5]      [,6]
# [1,] -0.2830150 -0.4424788 -0.08704122 0.3257994 0.2648910 0.2218447
# [2,] -0.2849254 -0.4322047 -0.08078586 0.3157196 0.2697697 0.2124266
# [3,] -0.2860894 -0.4184608 -0.07510660 0.3059113 0.2719245 0.2018209
# [4,] -0.2868542 -0.4055504 -0.06953929 0.2994935 0.2704166 0.1920337
# [5,] -0.2860143 -0.4015795 -0.06840179 0.2971201 0.2701748 0.1887005
# [6,] -0.2871081 -0.4080238 -0.07053192 0.3008878 0.2707280 0.1940479
#
# Result is a bit alarming as it obviously corrects WT and KO
# samples differently in opposite directions.
# It is mentioned in the guide that loess correction may remove
# some truly DB regions, so I will first try efficiency TMM normalization
#
# Filtering log2(15)
# > head(loess.offsets)
# [,1]       [,2]        [,3]      [,4]      [,5]      [,6]
# [1,] -0.2680033 -0.4170878 -0.05581422 0.2776614 0.2989191 0.1643248
# [2,] -0.2661286 -0.4115538 -0.05607539 0.2760606 0.2972989 0.1603983
# [3,] -0.2658210 -0.4107074 -0.05611850 0.2757561 0.2970737 0.1598171
# [4,] -0.2678112 -0.4165738 -0.05584343 0.2775338 0.2987631 0.1639315
# [5,] -0.2665601 -0.4128829 -0.05603048 0.2765113 0.2976778 0.1612843
# [6,] -0.2646900 -0.4083786 -0.05641008 0.2748914 0.2965234 0.1580640
# 
# Check with MA plots

par(mfrow=c(5, 2))
data.y <- asDGEList(filtered.data.loess)
adjc <- cpm(data.y, log=TRUE)
abval <- aveLogCPM(data.y)library(statmod)
for (i in c(2:6)) {
  # Not normalized
  mval <- adjc[,1]-adjc[,i]
  fit <- loessFit(x=abval, y=mval)
  smoothScatter(abval, mval, ylab = "M", xlab = "Average logCPM",
                main="Raw", ylim = c(-2,2), xlim = c(0,7))
  o <- order(abval)
  lines(abval[o], fit$filtered[o], col="red")
  # Normalized
  re.adjc <- log2(assay(filtered.data.loess)+0.5) - loess.offsets/log(2)
  mval <- re.adjc[,1]-re.adjc[,i]
  fit <- loessFit(x=abval, y=mval)
  smoothScatter(abval, re.adjc[,1]-re.adjc[,i], ylab = "M", xlab = "Average logCPM",
                main="Normalized", ylim = c(-2,2), xlim = c(0,7))
  lines(abval[o], fit$filtered[o], col="red")
}
dev.off()

# Will try to exclude samples which seem outliers on MDS plot

filtered.data.excl <- subset(filtered.data, select = c(T,T,F,F,T,T))

# Testing for differential binding
y <- asDGEList(filtered.data.excl)
genotype <- c("WT", "WT",
              "KO", "KO")
design <- model.matrix(~factor(genotype))
colnames(design) <- c("intercept", "genotype")
design
# > design
# intercept genotype
# 1         1        1
# 2         1        1
# 3         1        1
# 4         1        0
# 5         1        0
# 6         1        0
# attr(,"assign")
# [1] 0 1
# attr(,"contrasts")
# attr(,"contrasts")$`factor(genotype)`
# [1] "contr.treatment"

y <- estimateDisp(y, design)
summary(y$trended.dispersion)
# > summary(y$trended.dispersion)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.007331 0.011404 0.014500 0.014951 0.018931 0.022936 
fit <- glmQLFit(y, design, robust = TRUE)
summary(fit$var.post)
# > summary(fit$var.post)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7044  0.9928  1.0044  1.0038  1.0210  1.2782
par(mfrow=c(1,2))
o <- order(y$AveLogCPM)
plot(y$AveLogCPM, sqrt(y$trended.dispersion[o]), type="l",
     lwd=2, ylim=c(0,1), xlab=expression("Ave."~Log[2]~"CPM"),
     ylab=("Biological Coefficient of variation"))
plotQLDisp(fit)
dev.off()

######################################################################
# Something is wrong here or job is too heavy, will try later
# relevant <- rowSums(assay(data.no.input)) >= 20
# yo <- asDGEList(data.no.input[relevant], norm.factors=data.eff1)
# yo <- estimateDisp(yo, design)
# oo <- order(yo$AveLogCPM)
# plot(yo$AveLogCPM[oo], sqrt(yo$trended.dispersion[oo]), type="l", lwd=2,
#     ylim=c(0, max(sqrt(yo$))))
# lines.....
# legend.....
# manual page 42-43
#######################################################################

summary(fit$df.prior)
# Here the degrees of freedom seem to be quite large, will need to shrink
# model more, I added prior.count = 170 to glmQLFit arguments
# according to documentation of glmFit. ??? I think it changes nothing
results <- glmQLFTest(fit, contrast=c(0, 1))
head(results$table)

# Replicate similarity
par(mfrow=c(2,2), mar=c(5,4,2,2))
adj.counts <- cpm(y, log=TRUE)
for (top in c(100, 500, 1000, 5000)) {
  out <- plotMDS(adj.counts, main=top, col=c("blue", "blue", "red", "red"),
                 labels=c("wt.1", "wt.2", "ko.1", "ko.2"),
                 top = top)
}

# Clustering windows into regions
# Seems that I should use ensembl genome
#
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# broads <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)


BiocManager::install("EnsDb.Mmusculus.v79")
library(EnsDb.Mmusculus.v79)
broads <- genes(EnsDb.Mmusculus.v79)
broads <- resize(broads, width(broads)+3000, fix="end")
head(broads)
suppressWarnings(keep <- overlapsAny(rowRanges(filtered.data.excl), broads))
sum(keep)
olap <- findOverlaps(broads, rowRanges(filtered.data.excl))
olap

# Now will assign p-values to regions instead of windows
tabbroad <- combineOverlaps(olap, results$table)
head(tabbroad[!is.na(tabbroad$PValue),])
# Find significantly DB regions
is.sig.region <- tabbroad$FDR <= 0.25
table(tabbroad$direction[is.sig.region])
# Tried to change filtering threshold to log2(5) and log2(30)
# But log2(15) gives more DB regions
# Quick and dirty clustering (alternative)
merged <- mergeWindows(rowRanges(filtered.data.excl), tol=1000L)
merged$region
tabcom <- combineTests(merged$id, results$table)
head(tabcom)

# Determine window size range
summary(width(merged$region))

# To eliminate large windows
merged.max <- mergeWindows(rowRanges(filtered.data.excl), tol=1000L,
                           max.width = 5000L)
summary(width(merged.max$region))

# DB
is.sig.region <- tabcom$FDR <= 0.05
table(tabcom$direction[is.sig.region])
# Gives a bit more results than clustering based on annotation
# > is.sig.region <- tabcom$FDR <= 0.25
# > table(tabcom$direction[is.sig.region])
# 
# down   up 
# 31   19 
# > # DB
#   > is.sig.region <- tabcom$FDR <= 0.05
# > table(tabcom$direction[is.sig.region])
# 
# down   up 
# 2    1 
#
# Another approach -- up/down based on the most significant window of the region
tab.best <- getBestTest(merged$id, results$table)
head(tab.best)
# Adding start location of the most significant window within the region is useful
tabcom$best.logFC <- tab.best$logFC
tabcom$best.start <- start(rowRanges(filtered.data.excl))[tab.best$best]
head(tabcom[,c("best.logFC", "best.start")])
#
# Similar thing may be used for overlaps of regions and genes
tab.best.broad <- getBestOverlaps(olap, results$table)
tabbroad$best.logFC <- tab.best.broad$logFC
tabbroad$best.start <- start(rowRanges(filtered.data.excl))[tab.best.broad$best]
head(tabbroad[!is.na(tabbroad$PValue), c("best.logFC", "best.start")])
#
# Changing window size in windowCounts may result in more DB events
# in csaw it is possible to combine results obtained using different window widths
# It is done using functions consolidateWindows, consolidateOverlaps
# Will try to combine my 150bp windows with 1000bp windows
#################################################################################
#######Sometthing is wrong here with filtering, will check later#################
#################################################################################
# data.with.input.1k <- windowCounts(bam.files.wt.input,
#                                 ext = frag.length,
#                                 width = 1000L,
#                                 spacing = 500L,
#                                 filter = 35,
#                                 param=param)
# chip.1k <- data.with.input.1k[,2:7]
# control.1k <- data.with.input.1k[,1]
# # Will use same binned windows and their scale.info
# filter.stat.1k <- filterWindows(chip.1k, control.1k, type = "control",
#                              prior.count = 5, scale.info = scale.info)
# keep.1k <- filter.stat.1k$filter > log2(15)
# 
# filtered.data.1k <- chip.1k[keep,]
# 
# filtered.data.1k <- normFactors(filtered.data.1k, se.out = TRUE)
# filtered.data.excl.1k <- subset(filtered.data.1k, select = c(T,T,F,F,T,T))
# 
# # Testing for differential binding
# y.1k <- asDGEList(filtered.data.excl.1k)
# 
# y.1k <- estimateDisp(y.1k, design)
# fit.1k <- glmQLFit(y.1k, design, robust = TRUE)
# 
# results.1k <- glmQLFTest(fit.1k, contrast=c(0, 1))
# cons.win <- consolidateWindows(data.list = list(filtered.data.excl,
#                                                 filtered.data.excl.1k),
#                                equiweight = TRUE,
#                                merge.args = list(tol = 1000))
# cons.res <- consolidateTests(id.list = cons.win$id,
#                              result.list = list(results, results.1k),
#                              weight.list = cons.win$weight)
# cons.olap <- consolidateWindows(data.list = list(filtered.data.excl,
#                                                  filtered.data.excl.1k),
#                                 equiweight = TRUE, region = broads)
# cons.broad <- consolidateOverlaps(olap.list = cons.olap$olap,
#                                   result.list = list(results, results.1k),
#                                   weight.list = cons.win$weight)
#############################################################################
#
# Post-hoc analysis only on significant windows. Good for diffuse broad marks
postclust <- clusterWindows(rowRanges(filtered.data.excl), results$table,
                            target = 0.25, tol = 100, max.width = 1000)
postclust$FDR
postclust$region

# Using empirical FDR. Useful for noisy data

empres <- empiricalFDR(merged$id, results$table)
head(empres)

mcols(broads) <- tabbroad
mcols(merged$region) <- tabcom
# Adding gene-based annotation
BiocManager::install("ensembldb")
library(ensembldb)
BiocManager::install("TxDb")
Tx <- transcripts(EnsDb.Mmusculus.v79)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
install.packages("RMariaDB")
library(RMariaDB)
makeTxDbFromEnsembl(organism="Mus musculus",
#                    release=NA,
#                    circ_seqs=DEFAULT_CIRC_SEQS,
                    server="ensembldb.ensembl.org",
                    username="anonymous", password=NULL, port=3306L,
#                    tx_attrib=NULL
)
makeTxDbFromEnsembl(organism="Mus musculus",server="ensembldb.ensembl.org",
                    username="anonymous", port=3337)
Tx <- toGRanges(EnsDb.Mmusculus.v79, feature ="transcript")
anno <- detailRanges(merged$region, txdb = Tx,
                     orgdb = org.Mm.eg.db, promoter = c(3000, 1000),
                     dist = 5000, name.field = )
# Will try another package
BiocManager::install("ChIPpeakAnno")
library(ChIPpeakAnno)
data(TSS.mouse.GRCm38)
minimal <- merged$region
elementMetadata(minimal) <- NULL
overlaps.anno <- annotatePeakInBatch(minimal, AnnotationData=TSS.mouse.GRCm38)
colnames(elementMetadata(overlaps.anno))

library(org.Mm.eg.db)
overlaps.anno <- addGeneIDs(overlaps.anno,
                            "org.Mm.eg.db",
                            IDs2Add = "entrez_id")
head(overlaps.anno)
write.csv(as.data.frame(unname(overlaps.anno)),
          "/media/serlia/Storage\ 2/rora-treg-chipseq/anno.csv")

anno.regions
anno.to.file <- as.data.frame(anno.regions@elementMetadata@listData)
write.csv(anno.to.file, file = "/media/serlia/Storage\ 2/rora-treg-chipseq/anno.csv")
pie1(table(overlaps.anno$insideFeature))
over <- getEnrichedGO(overlaps.anno, orgAnn="org.Mm.eg.db", 
                      maxP=.05, minGOterm=10, 
                      multiAdjMethod="BH", condense=TRUE)
# Try to save to BED file
is.sig <- merged$region$FDR <= 0.25
library(rtracklayer)
test <- merged$region[is.sig]
test$score <- -10*log10(merged$region$FDR[is.sig])
names(test) <- paste0("region", 1:sum(is.sig))
export(test, "clusters.bed")
merged$region$FDR[is.sig]
head(read.table("clusters.bed"))
test
is.sig
merged$region$FDR




# Try once again to create txdb object from ensembl
txdb <- makeTxDbFromUCSC(genome = "mm10", tablename = "ensGene")
gtfFile <- 
txdb2 <- makeTxDbFromGFF(file=gtfFile,
                         chrominfo=chrominfo,
                         dataSource="ensemblgenomes",
                         organism="Aedes aegypti",
                         metadata=metadata)

# visualize
cur.region <- GRanges("1", IRanges(36939550, 36939950))
extractReads(bam.files.wt.input[2], cur.region, param=param)

BiocManager::install("Gviz")
library(Gviz)
collected <- vector("list", 6)
for (i in c(1, 2, 3, 4, 5, 6)) {
  print(i)
  reads <- extractReads(bam.files.noinput[i], cur.region, param=param)
  adj.total <- chip$totals[i]/1e6
  pcov <- as(coverage(reads[strand(reads)=="+"])/adj.total, "GRanges")
  ncov <- as(coverage(reads[strand(reads)=="-"])/adj.total, "GRanges")
  ptrack <- DataTrack(pcov, type="histogram", lwd=0, fill=rgb(0,0,1,.4),
                      ylim=c(0,1.1), name=bam.files.noinput[i], col.axis="black",
                      col.title="black")
  ntrack <- DataTrack(pcov, type="histogram", lwd=0, fill=rgb(0,0,1,.4),
                      ylim=c(0,1.1))
  collected[[i]] <- OverlayTrack(trackList=list(ptrack, ntrack))
}
gax <- GenomeAxisTrack(col="black")
plotTracks(c(gax, collected), from=start(cur.region), to=end(cur.region))
dev.off()
