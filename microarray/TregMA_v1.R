setwd("/Volumes/BigHDD/RORaTregMA")
source("https://bioconductor.org/biocLite.R")
biocLite("oligo")
biocLite("mogene20sttranscriptcluster.db")
biocLite("pd.mogene.2.0.st")
biocLite("affycoretools")

library("oligo")
library("affycoretools")
library("pd.mogene.2.0.st")
library("genefilter")
library(RColorBrewer)
library(gplots)

# Reading Data, Correcting Background, annotating
celFiles <- list.celfiles(full.names=T)
rawData <- read.celfiles(celFiles)
eset <- rma(rawData, target="core")
eset <- getMainProbes(eset)
eset <- annotateEset(eset, "mogene20sttranscriptcluster.db")
eset <- eset[!is.na(fData(eset)$ENTREZID),]

# Expression table
exp <- exprs(eset)
exp <- as.data.frame(exp)
exp <- cbind.data.frame(fData(eset)$ENTREZID, fData(eset)$SYMBOL, exp)
colnames(exp) <- c("ENTREZ_ID", "SYMBOL", "invitr_IL2_WT","invitr_IL2_IL4_WT","invitr_IL2_KO",
                   "invitr_IL2_IL4_KO","invitr_IL2_WT","invitr_IL2_IL4_WT",
                   "invitr_IL2_WT","invitr_IL2_IL4_WT","invitr_IL2_KO",
                   "invitr_IL2_IL4_KO","invitr_IL2_KO","invitr_IL2_IL4_KO",
                   "invitr_IL2_WT","invitr_IL2_IL4_WT","invitr_IL2_WT",
                   "invitr_IL2_IL4_WT","invitr_IL2_KO","invitr_IL2_IL4_KO",
                   "invitr_IL2_KO","invitr_IL2_IL4_KO","eWATHFD_WT","eWATHFD_WT",
                   "eWATHFD_WT","eWATHFD_KO","eWATHFD_KO","eWATHFD_KO","SPLCD_WT",
                   "SPLCD_WT","SPLCD_KO","SPLCD_KO","SPLCD_WT","SPLCD_WT","SPLCD_KO",
                   "SPLCD_KO","LUNG_WT","LUNG_WT","LUNG_WT","LUNG_WT","LUNG_KO",
                   "LUNG_KO","LUNG_KO","LUNG_KO","MLN_WT","MLN_WT","MLN_WT",
                   "MLN_WT","MLN_KO","MLN_KO","MLN_KO","MLN_KO")
geno <- c("WT","WT","KO", "KO", "WT","WT","WT","WT","KO","KO","KO",
	"KO","WT","WT","WT","WT","KO","KO","KO","KO","WT","WT","WT","KO","KO",
	"KO","WT","WT","KO","KO","WT","WT","KO","KO","WT","WT","WT","WT","KO",
	"KO","KO","KO","WT","WT","WT","WT","KO","KO","KO","KO")
treatment <- c("IL2","IL2_IL4","IL2","IL2_IL4","IL2","IL2_IL4","IL2",
	"IL2_IL4","IL2","IL2","IL2","IL2_IL4","IL2","IL2_IL4","IL2","IL2_IL4",
	"IL2","IL2_IL4","IL2","IL2_IL4","HFD","HFD","HFD","HFD","HFD","HFD","CD",
	"CD","CD","CD","CD","CD","CD","CD","HDM","HDM","HDM","HDM","HDM","HDM",
	"HDM","HDM","HDM","HDM","HDM","HDM","HDM","HDM","HDM","HDM")
tissue <- c("invitr","invitr","invitr","invitr","invitr","invitr","invitr",
	"invitr","invitr","invitr","invitr","invitr","invitr","invitr","invitr",
	"invitr","invitr","invitr","invitr","invitr","eWAT","eWAT","eWAT","eWAT",
	"eWAT","eWAT","SPL","SPL","SPL","SPL","SPL","SPL","SPL","SPL","LUNG","LUNG",
	"LUNG","LUNG","LUNG","LUNG","LUNG","LUNG","MLN","MLN","MLN","MLN","MLN",
	"MLN","MLN","MLN")
exp <- rbind.data.frame(c("","",geno), c("","",tissue), c("","",treatment), exp)
write.table(x=exp, file="RORaTregMAexp.txt", sep = "\t", row.names = F, col.names = T)

# Add variables to eSet
eset@phenoData@data$genotype <-
	factor(c("WT","WT","KO", "KO", "WT","WT","WT","WT","KO","KO","KO",
	"KO","WT","WT","WT","WT","KO","KO","KO","KO","WT","WT","WT","KO","KO",
	"KO","WT","WT","KO","KO","WT","WT","KO","KO","WT","WT","WT","WT","KO",
	"KO","KO","KO","WT","WT","WT","WT","KO","KO","KO","KO"))
eset@phenoData@data$treatment <-
	factor(c("IL2","IL2_IL4","IL2","IL2_IL4","IL2","IL2_IL4","IL2",
	"IL2_IL4","IL2","IL2","IL2","IL2_IL4","IL2","IL2_IL4","IL2","IL2_IL4",
	"IL2","IL2_IL4","IL2","IL2_IL4","HFD","HFD","HFD","HFD","HFD","HFD","CD",
	"CD","CD","CD","CD","CD","CD","CD","HDM","HDM","HDM","HDM","HDM","HDM",
	"HDM","HDM","HDM","HDM","HDM","HDM","HDM","HDM","HDM","HDM"))
eset@phenoData@data$tissue <-
	factor(c("invitr","invitr","invitr","invitr","invitr","invitr","invitr",
	"invitr","invitr","invitr","invitr","invitr","invitr","invitr","invitr",
	"invitr","invitr","invitr","invitr","invitr","eWAT","eWAT","eWAT","eWAT",
	"eWAT","eWAT","SPL","SPL","SPL","SPL","SPL","SPL","SPL","SPL","LUNG","LUNG",
	"LUNG","LUNG","LUNG","LUNG","LUNG","LUNG","MLN","MLN","MLN","MLN","MLN",
	"MLN","MLN","MLN"))
eset@phenoData@varMetadata <-
	rbind(eset@phenoData@varMetadata, data.frame(row.names = c("genotype", 
	"treatment","tissue"),labelDescription = c("genotype", "treatment", 
	"tissue"),channel = c("_ALL_","_ALL","_ALL_")))

# PCA
color <- c("black","blue","grey",
                   "cyan","black","blue",
                   "black","blue","grey",
                   "cyan","grey","cyan",
                   "black","blue","black",
                   "blue","grey","cyan",
                   "grey","cyan","yellow","yellow",
                   "yellow","orange","orange","orange","red",
                   "red","brown","brown","red","red","brown",
                   "brown","green","green","green","green","darkgreen",
                   "darkgreen","darkgreen","darkgreen","deeppink","deeppink","deeppink",
                   "deeppink","deeppink4","deeppink4","deeppink4","deeppink4")
data.PC <- prcomp(t(exprs(eset)), scale.=T)
plot(data.PC$x, col=color)
leg <- paste(pData(eset)$tissue, " ", pData(eset)$genotype, " ", pData(eset)$treatment)
leg <- leg[!duplicated(leg)]
col_leg <- color[!duplicated(color)]
legend("topright", col = col_leg, legend = leg, pch = 20, bty = "n", cex=.75) 
                   
# Heatmap.2 for WT spleen vs WT lung
# Correlation distance and clustering functions
# https://bioramble.wordpress.com/2015/08/03/heatmaps-part-3-how-to-create-a-microarray-heatmap-with-r/
dist_cor <- function(x) {
  as.dist(1 - cor(t(x), method = "pearson"))
}
clus_wd2 <- function(x) {
  hclust(x, method = "ward.D2")
}
# colour scheme for the heatmap
redblackgreen <- colorRampPalette(c("green", "black", "red"))(n = 100)

# subset eSet WT spl-lung. top 200 variable genes
WT <- eset@phenoData@data$genotype %in% c("WT")
spl_lung <- eset@phenoData@data$tissue %in% c("SPL", "LUNG")

esetWTspl_lung <- eset[, WT & spl_lung]
esetWTspl_lung$genotype <- droplevels(esetWTspl_lung$genotype)
esetWTspl_lung$treatment <- droplevels(esetWTspl_lung$treatment)
esetWTspl_lung$tissue <- droplevels(esetWTspl_lung$tissue)
WTspl_lung_sd <- rowSds(exprs(esetWTspl_lung))
top200 <- names(sort(WTspl_lung_sd, decreasing = TRUE))[1:200]
WTspl_lung_var <- esetWTspl_lung[top200, ]
WTspl_lung_var@featureData@data$SYMBOL <-
  droplevels(WTspl_lung_var@featureData@data$SYMBOL)
class_labels <- ifelse(WTspl_lung_var$tissue == "SPL", "blue", "orange")
geneSym <- fData(WTspl_lung_var)[3]
geneSym <- unlist(geneSym)
heatmap.2(exprs(WTspl_lung_var),
          # clustering
          distfun = dist_cor, 
          hclust = clus_wd2,
          # scaling (genes are in rows)
          scale = "row",
          # color
          col = redblackgreen, 
          # labels
          labRow = geneSym, 
          #labCol = " ",
          ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none",
          key = 1)
          
# subset eSet KO spl-lung
KO <- eset@phenoData@data$genotype %in% c("KO")
spl_lung <- eset@phenoData@data$tissue %in% c("SPL", "LUNG")
#NArm <- !is.na(eset@featureData@data$SYMBOL)

esetKOspl_lung <- eset[, KO & spl_lung]
esetKOspl_lung$genotype <- droplevels(esetKOspl_lung$genotype)
esetKOspl_lung$treatment <- droplevels(esetKOspl_lung$treatment)
esetKOspl_lung$tissue <- droplevels(esetKOspl_lung$tissue)

#esetWTspl_lung1 <- nsFilter(esetWTspl_lung, require.entrez=TRUE,remove.dupEntrez=FALSE,
#                            var.filter=FALSE,feature.exclude="^AFFX")
KOspl_lung_sd <- rowSds(exprs(esetKOspl_lung))
top200 <- names(sort(KOspl_lung_sd, decreasing = TRUE))[1:200]
KOspl_lung_var <- esetKOspl_lung[top200, ]
KOspl_lung_var@featureData@data$SYMBOL <-
  droplevels(KOspl_lung_var@featureData@data$SYMBOL)
class_labels <- ifelse(KOspl_lung_var$tissue == "SPL", "blue", "orange")
geneSym <- fData(KOspl_lung_var)[3]
geneSym <- unlist(geneSym)
heatmap.2(exprs(KOspl_lung_var),
          # clustering
          distfun = dist_cor, 
          hclust = clus_wd2,
          # scaling (genes are in rows)
          scale = "row",
          # color
          col = redblackgreen, 
          # labels
          labRow = geneSym, 
          #labCol = " ",
          ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none",
          key = 1)



library(limma)
Conditions <- factor(c("invitr_IL2_WT","invitr_IL2_IL4_WT","invitr_IL2_KO",
                       "invitr_IL2_IL4_KO","invitr_IL2_WT","invitr_IL2_IL4_WT",
                       "invitr_IL2_WT","invitr_IL2_IL4_WT","invitr_IL2_KO",
                       "invitr_IL2_IL4_KO","invitr_IL2_KO","invitr_IL2_IL4_KO",
                       "invitr_IL2_WT","invitr_IL2_IL4_WT","invitr_IL2_WT",
                       "invitr_IL2_IL4_WT","invitr_IL2_KO","invitr_IL2_IL4_KO",
                       "invitr_IL2_KO","invitr_IL2_IL4_KO","eWATHFD_WT","eWATHFD_WT",
                       "eWATHFD_WT","eWATHFD_KO","eWATHFD_KO","eWATHFD_KO","SPLCD_WT",
                       "SPLCD_WT","SPLCD_KO","SPLCD_KO","SPLCD_WT","SPLCD_WT","SPLCD_KO",
                       "SPLCD_KO","LUNG_WT","LUNG_WT","LUNG_WT","LUNG_WT","LUNG_KO",
                       "LUNG_KO","LUNG_KO","LUNG_KO","MLN_WT","MLN_WT","MLN_WT",
                       "MLN_WT","MLN_KO","MLN_KO","MLN_KO","MLN_KO"))

design <- model.matrix(~0+Conditions)
colnames(design) <- levels(Conditions)
fit <- lmFit(eset, design)

# With WT spleen as control
contrast.matrix <- makeContrasts(Lung_WTvsSpl_WT = LUNG_WT - SPLCD_WT,
								Lung_KOvsSpl_WT = LUNG_KO - SPLCD_WT,
								Diff = (LUNG_WT - SPLCD_WT) - (LUNG_KO - SPLCD_WT),
                                 levels = design)
# Just KO-WT in lungs
contrast.matrix <- makeContrasts(Lung_KOvsLung_WT = LUNG_KO - LUNG_WT, levels = design)
                                                                 
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
table <- topTable(fit2, coef=1, adjust="BH", number = 10000) #BH -- Benjamini-Hochberg adjusted P value
write.table(x = table, file = "WTvsKO_Lung_limma.txt", sep =",", col.names = NA, row.names=T)
results <- decideTests(fit2)
volcanoplot(fit2, coef=1, highlight=20, names = fit2$genes$SYMBOL)

vennDiagram(results)
top.iv <- as.numeric(rownames(topTable(fit2,n=50)))

# Genes upregulated in lung Treg vs spleen WT
genesUP_SPLvsLUNGwt <- c("Vps37b","Dusp4","Gadd45b","Ell2","Hilpda","Dgat1","Mxd1","Epas1","Ccr5","Chil3","Ccl4","Gas2l3","Ccl3","Car4","Ctla2b","Cd63","Plaur","Xist","Gm29719","Ramp3","Ier3","Lmna","Plet1","Clec4d","Emp1","Thbs1","Ier5l","Mxi1","Gpnmb","Tubb6","Nr4a3","Il10","Areg","Dusp5","Ets2","Dusp10","Nr4a1","Got1","Nfil3","Rgs2","Egln3","Sik1","Fosl2","Nr4a2","Cxcl2","Gem","Wdyhv1","Errfi1","Atf3","Cd14","Klk1","Gm19723","Fgl2","Cdkn1a","Cpd","Slc15a3","Ctla2a","Kctd12","Il1rl1","Mt2","Ccl6","Ccr8","Pparg","Furin","AA467197","Mt1","Ccr1","Cxcl3","Gzmb","Coq10b","Ltb4r1","Calca","Klf4","Il17a","Il13","Ccl17","Wfdc17","Fabp4","Il1b","Spp1","Acod1","Gprc5a","Dyrk3","Retnlg","Cxcr6","Clec7a","Retnla","Scgb1a1","Fpr2","Ear1","Ear2","Ccl8","Snora34","Asb2","Stfa2l1","Mir1946b","Snora62","Cmc2","Ifnk","Galr3","Hist1h4h","Fcor","Rnase2a","Igkv14-126","Cox7a1","Lrr1","Bambi","Bex1","Igh-V3609N")

genesDOWN_SPLvsLUNGwt <- c("Gm3264","Gm5796","Gm3373","Igk","Tespa1","Mir103-2","Gm5796
Clec2i","Dennd2d","1810006J02Rik","4930431P03Rik","Gm10408","D830030K20Rik","Gm8271","Gm8271","Naip5","Igkv1-135","Igk-V28","Igh-VX24","Ddx43","Ighm","Ighm","LOC382693","Ighm","Igh-V7183","Nup37","Mterf1b","Kdm5d","Uty","Gm4767","Eif2s3y","Igh-VJ558","Treml2","Ddx3y","Sostdc1","Igh-VJ558","Gm4956","Igk","Gm13707","Gucy1a1","Ighv1-19","Gm17764","Igh-VJ558","Iglv1","Gm4951","Trim34a","Naip6","Igkv4-72","Cr2","LOC544905","Igh-VJ558","Ighv10-3","Ighg3","Spic","Igkv14-111","Ighv1-62","Igk","Slc40a1","Igkv4-55","Igk-V8","Art2b","Igkv10-96","Igkv1-117","Igh-VJ558","Mir1955","Akr1c13","Mir20a","Stfa3","Igh-VX24","Mela","Igkv3-4","Snora69","Vmn1r127","Tnpo3","Tnpo3
Asb2","Snord110","Ighg","Igkv5-43","Igkv15-103","Igk","Igk-V1","Mir669k","Mir669k","Mir669k","Ttc36","Mir669i","Ighm","Igkv4-91")

genesGSEA <- c("IL13","CCR1","CXCL3","CCL8","IL17A","CCL17","CCR5","CALCA","IL10","CLEC7A","THBS1","LTB4R","IFNK","MT2A","PPARG","CLEC4D","GZMB","CXCR6","FPR2","ASB2","IL1RL1","NFIL3","KLF4","EPAS1","FURIN","BAMBI","SCGB1A1","PLAUR","HIST1H4H","KCTD12","IER3","TUBB6","AREG","GPNMB","FABP4","GPRC5A","IER5L","LRR1","CA4","COX7A1","GAS2L3","WDYHV1","SNORA34","SNORA62")

genesUPinKOlung <- c("Mela","Rnase2a","Ccl24","Mmp12","Mir302b","LOC102638993","Gzmb","Ccr3","Depdc1a","Ermn","Cops7b","Reg3g","Magohb","9230102O04Rik","Serpinb2","Clca1","Rnf113a1","Trim12a","Cd36","Med21","Gm4767","Mir302a","Cd80","Gm4737","Ctsk","Bhlhe41","Gpr25","Dusp16","Gm7160","Calca","Melk","Prr11","Ptgir","LOC105244102","Akip1","Gm17745","Itgax","Alg8","Arg1","2510002D24Rik","Hsh2d","Tmem106a","Gnb4","Cacybp","4930433I11Rik","Commd9","Chrna1os","Wdr76","Tpx2","Rhoc","Katnbl1","Zfp438","Parp11","Mut","Ctbs","Ssbp1","Ap3s2","Uba3","Fundc1","Pde3b","Fra10ac1","Larp4b","Cdc25b","Dram2","Mmp19","Sec22a","Nkrf","Dlgap5","Plk4","Plet1","Nfyb","Vcpkmt","Gm4890","Lztfl1","E2f7","Tmem9","Mcm7","Gm14005","Rad17","AI662270","Mirt2","Cdca2","Arel1","Ahsa1","Gatad1","Cdc23","Slc25a19","Cyp2ab1","Pold3","Lman1","Rbp4","Msl1","Socs6","Fgfr1op2")

genesDOWNinKOlung <- c("Gjb1","Gm10666","Gm10666","Polr3b","Gm10665","Gm10665","Myo5c","Npdc1","Tmem125","Bmf","Gm14199","Muc19","Peg13","Nop14","Trim65","Mlh1","Sh2b3","Endov","Gm3948","Ankrd6","Myoc","Whrn","Gm1987","Cnst","Nebl","Plxna1","Zfp956","Nrip1","Urb1","Cryab","Sdc4","Rab3ip","Hlf","Coq10a","Pcsk1","Hdhd5","Gm1987","Emx1","Zfp446","St6galnac3","Arntl","Gm11938","Fam124b","March9","Snord14a","Ubtd2","Gm26567","Gm24463","Gm24463","Gm6410","Gm24463","Gm24463","Rasl11a","Nr1d1","Ccl7","Il23r","Snord83b","Rnu73b","Fam160a1","Jchain","Mir16-2","Snora7a","Rrp7a","Ighg","Igh-V3609N")




eset1 <- eset[(casefold(fData(eset)$SYMBOL) %in% casefold(genesGSEA)),]

mln_lung <- eset@phenoData@data$tissue %in% c("MLN", "LUNG")
eset1 <- eset1[, mln_lung]
eset1$genotype <- droplevels(eset1$genotype)
eset1$treatment <- droplevels(eset1$treatment)
eset1$tissue <- droplevels(eset1$tissue)
eset1@featureData@data$SYMBOL <-
  droplevels(eset1@featureData@data$SYMBOL)
class_labels <- ifelse(pData(eset1)$genotype == "WT", "grey40", "grey80")
geneSym <- fData(eset1)[3]
geneSym <- unlist(geneSym)
heatmap.2(exprs(eset1),
          # clustering
          #distfun = dist_cor, 
          #hclust = clus_wd2,
          # scaling (genes are in rows)
          scale = "row",
          # color
          col = redblackgreen, 
          # labels
          labRow = geneSym, 
          #labCol = " ",
          ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none",
          key = 1,
          margins = c(9,9))



# subset eSet WT spl-eWAT. top 200 variable genes

WT <- eset@phenoData@data$genotype %in% c("WT")
spl_ewat <- eset@phenoData@data$tissue %in% c("SPL", "eWAT")

esetWTspl_ewat <- eset[, WT & spl_ewat]
esetWTspl_ewat$genotype <- droplevels(esetWTspl_ewat$genotype)
esetWTspl_ewat$treatment <- droplevels(esetWTspl_ewat$treatment)
esetWTspl_ewat$tissue <- droplevels(esetWTspl_ewat$tissue)
WTspl_ewat_sd <- rowSds(exprs(esetWTspl_ewat))
top200 <- names(sort(WTspl_ewat_sd, decreasing = TRUE))[1:200]
WTspl_ewat_var <- esetWTspl_ewat[top200, ]
WTspl_ewat_var@featureData@data$SYMBOL <-
  droplevels(WTspl_ewat_var@featureData@data$SYMBOL)
class_labels <- ifelse(WTspl_ewat_var$tissue == "SPL", "blue", "orange")
geneSym <- fData(WTspl_ewat_var)[3]
geneSym <- unlist(geneSym)
heatmap.2(exprs(WTspl_ewat_var),
          # clustering
          distfun = dist_cor, 
          hclust = clus_wd2,
          # scaling (genes are in rows)
          scale = "row",
          # color
          col = redblackgreen, 
          # labels
          labRow = geneSym, 
          #labCol = " ",
          ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none",
          key = 1)
          

# subset eSet KO spl-eWAT

KO <- eset@phenoData@data$genotype %in% c("KO")
spl_ewat <- eset@phenoData@data$tissue %in% c("SPL", "eWAT")

esetKOspl_ewat <- eset[, KO & spl_ewat]
esetKOspl_ewat$genotype <- droplevels(esetKOspl_ewat$genotype)
esetKOspl_ewat$treatment <- droplevels(esetKOspl_ewat$treatment)
esetKOspl_ewat$tissue <- droplevels(esetKOspl_ewat$tissue)
KOspl_ewat_sd <- rowSds(exprs(esetKOspl_ewat))
top200 <- names(sort(KOspl_ewat_sd, decreasing = TRUE))[1:200]
KOspl_ewat_var <- esetKOspl_ewat[top200, ]
KOspl_ewat_var@featureData@data$SYMBOL <-
  droplevels(KOspl_ewat_var@featureData@data$SYMBOL)
class_labels <- ifelse(KOspl_ewat_var$tissue == "SPL", "blue", "orange")
geneSym <- fData(KOspl_ewat_var)[3]
geneSym <- unlist(geneSym)
heatmap.2(exprs(KOspl_ewat_var),
          # clustering
          distfun = dist_cor, 
          hclust = clus_wd2,
          # scaling (genes are in rows)
          scale = "row",
          # color
          col = redblackgreen, 
          # labels
          labRow = geneSym, 
          #labCol = " ",
          ColSideColors = class_labels, 
          # tweaking
          trace = "none",
          density.info = "none",
          key = 1)


library(limma)
Conditions <- factor(c("invitr_IL2_WT","invitr_IL2_IL4_WT","invitr_IL2_KO",
                       "invitr_IL2_IL4_KO","invitr_IL2_WT","invitr_IL2_IL4_WT",
                       "invitr_IL2_WT","invitr_IL2_IL4_WT","invitr_IL2_KO",
                       "invitr_IL2_IL4_KO","invitr_IL2_KO","invitr_IL2_IL4_KO",
                       "invitr_IL2_WT","invitr_IL2_IL4_WT","invitr_IL2_WT",
                       "invitr_IL2_IL4_WT","invitr_IL2_KO","invitr_IL2_IL4_KO",
                       "invitr_IL2_KO","invitr_IL2_IL4_KO","eWATHFD_WT","eWATHFD_WT",
                       "eWATHFD_WT","eWATHFD_KO","eWATHFD_KO","eWATHFD_KO","SPLCD_WT",
                       "SPLCD_WT","SPLCD_KO","SPLCD_KO","SPLCD_WT","SPLCD_WT","SPLCD_KO",
                       "SPLCD_KO","LUNG_WT","LUNG_WT","LUNG_WT","LUNG_WT","LUNG_KO",
                       "LUNG_KO","LUNG_KO","LUNG_KO","MLN_WT","MLN_WT","MLN_WT",
                       "MLN_WT","MLN_KO","MLN_KO","MLN_KO","MLN_KO"))

design <- model.matrix(~0+Conditions)
colnames(design) <- levels(Conditions)
fit <- lmFit(eset, design)

# With WT spleen as control
contrast.matrix <- makeContrasts(Lung_WTvsSpl_WT = LUNG_WT - SPLCD_WT,
								Lung_KOvsSpl_WT = LUNG_KO - SPLCD_WT,
								Diff = (LUNG_WT - SPLCD_WT) - (LUNG_KO - SPLCD_WT),
                                 levels = design)
# Just KO-WT in eWAT
contrast.matrix <- makeContrasts(eWAT_KOvseWAT_WT = eWATHFD_KO - eWATHFD_WT, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
table <- topTable(fit2, coef=1, adjust="BH", number = 10000) #BH -- Benjamini-Hochberg adjusted P value
write.table(x = table, file = "WTvsKO_eWAT_limma.txt", sep =",", col.names = NA, row.names=T)

genesUPewatWTvsSPL <- c("Vps37b","Dgat1","Ell2","Hilpda","Gadd45b","Furin","AA467197","Lilr4b","Ltb4r1","Syne1","Syne1","Syne1","","Dpt","Ccr1","Fn1","Pcsk1","Syne1","Slc24a1","Cyp11a1","Sytl3","Mfap5","Dab2","LOC102642252","Syne1","","Fbn1","Atf3","Tnfaip6","Igkv4-68","Adam8","Rab4a","Syne1","Bgn","Syne1","Sparc","Il5","Syne1","Syne1","Cxcl1","Mmp25","Syne1","Efemp1","Syne1","Ccl22","Cyr61","Syne1","Areg","Dusp5","Irs2","Nr4a1","Fosl2","Nr4a3","Nr4a2","Il10","Egln3","Got1","Klf4","Cxcl2","Syne1","Gem","Igkv4-91","Syne1","Gas2l3","Syne1","Mgp","Ccl2","Syne1","Lyz1","Syne1","Acod1","Syne1","Syne1","Ccl9","Igh-VX24","Syne1","Ighg","Calca","Ighm","Errfi1","Il13","Igkv1-117","Tnfrsf10b","Igkv19-93","LOC544905","Lsmem1","Igkv2-112","Igk","","Ier3","Fam183b","Ccl11","Dennd4a","Cxcl3","Cpm","Cxcr6","Dcn","Igh-V11","Lmna","Dyrk3","Il1rl1","Igkv16-104","Myadm","Gm29721","Mt2","Anxa1","Cdkn1a","Retnla","Plin2","Zdhhc23","Pcyt1a","Ndrg1","Kcna4","Igk","Gm29719","Ctla2a","Frmd5","Kctd12","Clec9a","Ptafr","Serping1???,???Igh-V7183","Mt1","Ctla2b","C3","Fabp4","Raph1","Ccl6","Pparg","Slc15a3","Fgl2","Ifi205","Cd36","Igkv14-126","Klrg1","Iglv2","6330407A03Rik","Igkv4-53","Igh-V7183","Igkv4-70","Alox5ap","Lyn","Cd14","Igk-V1","Ighm","Il1b","Igh-V3660","Jag1","Olfr1258","Syne1","Igk","Igkv6-14","Ighv9-4","4930412O13Rik","Syne1","Syne1","Rnase2a","Syne1","Ighv7-3???,???A630038E17Rik???,","???Plac9b","Plac9b","Fam71a","Klrb1b","Syne1","LOC105245453")

genesUPKOvsWT_eWAT <- table[table$logFC >= 0.5 & table$adj.P.Val <= 0.05, ]$SYMBOL
geneOverlap <- genesUPewatWTvsSPL[genesUPewatWTvsSPL %in% genesUPKOvsWT_eWAT]
write.table(x=genesUPKOvsWT_eWAT, file="genesUPKOvsWT_eWAT_limma.txt")

genesDOWNKOvsWT_eWAT <- table[table$logFC <= -0.5 & table$adj.P.Val <= 0.05, ]$SYMBOL
geneOverlap <- genesUPewatWTvsSPL[genesUPewatWTvsSPL %in% genesDOWNKOvsWT_eWAT]
write.table(x=genesDOWNKOvsWT_eWAT, file="genesDOWNKOvsWT_eWAT_limma.txt")
